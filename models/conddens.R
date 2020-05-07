
estimateJointDensities <- function(
  data = data.frame(), densFunSuffix, ignoreGeneration = TRUE,
  featuresPdfPmf = c(), featuresCdf = c())
{
  if (!is.data.frame(data) || ncol(data) < 1) {
    stop("The data given is not a data.frame or too sparse.")
  }
  
  if (length(featuresPdfPmf) > 0 && length(featuresCdf) > 0) {
    warning("Estimating PDF and CDF, be careful using them together.")
  }
  
  estFuncs <- list()
  cn <- colnames(data)
  fn <- c(featuresPdfPmf, featuresCdf)
  
  if (ignoreGeneration) {
    cn <- gsub("_t_\\d+$", "", cn)
    fn <- gsub("_t_\\d+$", "", fn)
    colnames(data) <- cn
  }
  
  for (feat in fn) {
    if (!(feat %in% cn)) {
      stop(paste("The feature", feat, "is not in the data.frame."))
    }
    
    densFunName <- paste(feat, densFunSuffix, sep = "__@dens@__")
    isDiscrete <- !is.numeric(data[[feat]])
    doEcdf <- feat %in% featuresCdf
    
    if (isDiscrete) {
      estFuncs[[densFunName]] <- (function() {
        temp <- as.character(data[[feat]])
        return(function(x) mmb::getProbForDiscrete(temp, x))
      })()
    } else if (doEcdf) {
      estFuncs[[densFunName]] <- (function() {
        tryCatch({
          temp <- data[[feat]]
          ecdf(temp)
        }, error=function(cond) function(x) 0)
      })()
    } else {
      estFuncs[[densFunName]] <- (function() {
        pdf <- mmb::estimatePdf(data[[feat]])
        return(function(x) pdf$fun(x))
      })()
    }
  }
  
  return(estFuncs)
}



condDens_phi1j_default <- function(O_1, initProbs, states, stateColumn, densities) {
  b_j_O_1 <- computeDensitiesSum(
    O_t = O_1, states = states, densities = densities)
  
  predTemp <- c()
  for (state in states) {
    predTemp[state] <- initProbs[state] * b_j_O_1[state]
  }
  O_1[[stateColumn]] <- names(which.max(predTemp))
  return(O_1[[stateColumn]])
}


#' 1st-order conditional density model (model A).
condDens_A_1stOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE)
{
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  
  # initialize the two lists of densities, one for d_j, one for b_j:
  # d_j is for observations that had j as successor, and it builds
  # densities over the variables ending in _t_1.
  # b_j is like in the depmix-models; i.e., it builds conditional
  # densities over observations that had state j.
  # Also, in this model, d_j, b_j (and v_j in 2nd-order) all use the
  # same segmented dataset; the difference is over which generations'
  # variables they estimate their densities.
  
  densities_d <- list()
  for (state in states) {
    # Data at t-0 that had j(=state), which means that any data
    # at t-1 had this j as successor:
    condData <- mmb::conditionalDataMin(
      df = df,
      features = mmb::createFeatureForBayes(name = stateColumn, value = state),
      selectedFeatureNames = c(stateColumn))
    
    # Now we want conditional densities over the t_1-variables:
    featNames <- colnames(condData)[grepl("_t_1$", colnames(condData))]
    densities_d <- append(
      densities_d,
      estimateJointDensities(
        data = condData[, featNames],
        densFunSuffix = state,
        ignoreGeneration = TRUE,
        featuresPdfPmf = if (doEcdf) c() else featNames,
        featuresCdf = if (doEcdf) featNames else c()
      ))
  }
  
  # The densities for b_j we can get from the depmix-estimator.
  # However, we want only the t_0 features:
  featNames <- colnames(df)[grepl("_t_0$", colnames(df))]
  densities_b <- estimateDepmixDensities(
    states = states, stateColumn = stateColumn, data = df,
    featuresPdfPmf = if (doEcdf) c() else featNames,
    featuresCdf = if (doEcdf) featNames else c())
  
  
  # Estimate O_1 using \phi_1(j) and \pi_j:
  O_1 <- observations[1, ]
  initProbs <- getInitialStateProbs(
    data = df, stateColumn = stateColumn)
  
  predTemp <- condDens_phi1j_default(
    O_1 = O_1, initProbs = initProbs, states = states,
    stateColumn = stateColumn, densities = densities_b)
  O_1[[stateColumn]] <- predTemp
  pred <- c(predTemp)
  
  if (nrow(observations) > 1) {
    # Now assign the label of the greatest likelihood to each
    # observation:
    O_t_1 <- O_1
    for (j in 2:nrow(observations)) { # we start at 2, as O_1 is done already
      # Attention: We use O_t_1 for d_j!
      d_j <- computeDensitiesSum(
        O_t = O_t_1, states = states,
        densities = densities_d, ignoreGeneration = TRUE) # IMPORTANT!
      
      O_t <- observations[j, ]
      b_j <- computeDensitiesSum(
        O_t = O_t, states = states, densities = densities_b)
      
      predTemp <- c()
      for (state in states) {
        predTemp[state] <- d_j[state] * b_j[state]
      }
      
      O_t[stateColumn] <- names(which.max(predTemp))
      pred <- c(pred, O_t[[stateColumn]])
      
      # Also, move O_t_1:
      O_t_1 <- O_t
    }
  }
  
  mmb::setWarnings(w)
  return(pred)
}


#' 2nd-order conditional density model (model A). The underlying
#' model is $$\phi_t^2(j) = d_j^2(O_{t-2}) * d_j^1(O_{t-1}) * b_j(O_t) $$.
condDens_A_2ndOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE)
{
  # BEWARE: In this function, we use the NEW notation already!
  # d_j, v_j become d_j^1(O_{t-1}) and d_j^2(O_{t-2}).
  
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  
  # initialize the three lists of densities, one for d_j, one for b_j:
  # d_j is for observations that had j as successor, and it builds
  # densities over the variables ending in _t_1.
  # b_j is like in the depmix-models; i.e., it builds conditional
  # densities over observations that had state j.
  # Also, in this model, d_j, b_j (and v_j in 2nd-order) all use the
  # same segmented dataset; the difference is over which generations'
  # variables they estimate their densities.
  
  densities_d_2 <- list()
  densities_d_1 <- list()
  for (state in states) {
    condData <- mmb::conditionalDataMin(
      df = df, selectedFeatureNames = c(stateColumn),
      features = mmb::createFeatureForBayes(name = stateColumn, value = state))
    
    cn <- colnames(condData)
    featNames_2 <- cn[grepl("_t_2$", cn)]
    featNames_1 <- cn[grepl("_t_1$", cn)]
    
    densities_d_2 <- append(densities_d_2, estimateJointDensities(
      data = condData[, featNames_2],
      densFunSuffix = state,
      ignoreGeneration = TRUE,
      featuresPdfPmf = if (doEcdf) c() else featNames_2,
      featuresCdf = if (doEcdf) featNames_2 else c()
    ))
    
    densities_d_1 <- append(densities_d_1, estimateJointDensities(
      data = condData[, featNames_1],
      densFunSuffix = state,
      ignoreGeneration = TRUE,
      featuresPdfPmf = if (doEcdf) c() else featNames_2,
      featuresCdf = if (doEcdf) featNames_2 else c()
    ))
  }
  
  # The densities for b_j we can get from the depmix-estimator.
  # However, we want only the t_0 features:
  featNames <- colnames(df)[grepl("_t_0$", colnames(df))]
  densities_b <- estimateDepmixDensities(
    states = states, stateColumn = stateColumn, data = df,
    featuresPdfPmf = if (doEcdf) c() else featNames,
    featuresCdf = if (doEcdf) featNames else c())
  
  
  # Estimate O_1 using \phi_1(j) and \pi_j:
  O_1 <- observations[1, ]
  initProbs <- getInitialStateProbs(
    data = df, stateColumn = stateColumn)
  
  predTemp <- condDens_phi1j_default(
    O_1 = O_1, initProbs = initProbs, states = states,
    stateColumn = stateColumn, densities = densities_b)
  O_1[[stateColumn]] <- predTemp
  pred <- c(predTemp)
  
  if (nrow(observations) > 1) {
    O_t_2 <- 0 # will only be available after the next round
    O_t_1 <- O_1
    
    for (j in 2:nrow(observations)) {
      use_d_2 <- j > 2
      d_2 <- if (!use_d_2) NULL else computeDensitiesSum(
        O_t = O_t_2, states = states,
        densities = densities_d_2, ignoreGeneration = TRUE)
      
      d_1 <- computeDensitiesSum(
        O_t = O_t_1, states = states,
        densities = densities_d_1, ignoreGeneration = TRUE)
      
      O_t <- observations[j, ]
      b_j <- computeDensitiesSum(
        O_t = O_t, states = states, densities = densities_b)
      
      predTemp <- c()
      for (state in states) {
        if (use_d_2) {
          predTemp[state] <- d_2[state] * d_1[state] * b_j[state]
        } else {
          predTemp[state] <- d_1[state] * b_j[state]
        }
      }
      
      O_t[[stateColumn]] <- names(which.max(predTemp))
      pred <- c(pred, O_t[[stateColumn]])
      
      # Also move O_t_2 and O_t_1:
      O_t_2 <- O_t_1
      O_t_1 <- O_t
    }
  }
  
  
  mmb::setWarnings(w)
  return(pred)
}






#' 1st-order conditional density model (model B). The underlying
#' model is $$\frac_{\sum_{i=1}^{N} c_i^1(O_t) + b_j(O_t) }{N *  b_j(O_t)}$$
condDens_B_1stOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE)
{
  # Here, we use the c_i functions, which fix the data at t-1=j,
  # and then segment over the t_0 variables.
  # This model is special, because it does not need to use an
  # extra \phi_1(j) function.
  
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  # We need to fix this column in each iteration:
  fixColName <- gsub("_t_0$", "_t_1", stateColumn)
  densities_c_1 <- list()
  for (state in states) {
    condData <- mmb::conditionalDataMin(
      df = df, selectedFeatureNames = c(fixColName),
      features = mmb::createFeatureForBayes(name = fixColName, value = state))
    
    cn <- colnames(condData)
    featNames <- cn[grepl("_t_0$", cn)]
    
    densities_c_1 <- append(densities_c_1, estimateJointDensities(
      data = condData[, featNames],
      densFunSuffix = state,
      ignoreGeneration = TRUE,
      featuresPdfPmf = if (doEcdf) c() else featNames,
      featuresCdf = if (doEcdf) featNames else c()
    ))
  }
  
  
  # The densities for b_j we can get from the depmix-estimator.
  # However, we want only the t_0 features:
  featNames <- colnames(df)[grepl("_t_0$", colnames(df))]
  densities_b <- estimateDepmixDensities(
    states = states, stateColumn = stateColumn, data = df,
    featuresPdfPmf = if (doEcdf) c() else featNames,
    featuresCdf = if (doEcdf) featNames else c())
  
  
  pred <- c()
  for (j in 1:nrow(observations)) {
    O_t <- observations[j, ]
    b_j_O_t <- computeDensitiesSum(
      O_t = O_t, states = states, densities = densities_b)
    
    predTemp <- c()
    
    for (state_j in states) {
      c_i_O_t <- computeDensitiesSum(
        O_t = O_t, states = states,
        densities = densities_c_1, ignoreGeneration = TRUE)
      
      nominator <- 0
      denominator <- 0
      
      for (state_i in states) {
        nominator <- nominator + c_i_O_t[[state_i]] + b_j_O_t[[state_j]]
        denominator <- denominator + b_j_O_t[[state_j]]
      }
      
      predTemp[state_j] <- nominator / denominator
    }
    
    pred <- c(pred, names(which.max(predTemp)))
  }
  
  mmb::setWarnings(w)
  return(pred)
}
