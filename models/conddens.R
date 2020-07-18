
estimateJointDensities <- function(
  data = data.frame(), densFunSuffix, ignoreGeneration = TRUE,
  featuresPdfPmf = c(), featuresCdf = c(), normalizePdf = FALSE,
  ecdfMinusOne = FALSE)
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
          cdf <- stats::ecdf(data[[feat]])
          return(function(x) if (ecdfMinusOne) 1 - cdf(x) else cdf(x))
        }, error=function(cond) function(x) 0)
      })()
    } else {
      estFuncs[[densFunName]] <- (function() {
        pdf <- mmb::estimatePdf(data[[feat]], densFun = stats::density)
        normFac <- if (normalizePdf) max(pdf$y) else 1
        return(function(x) pdf$fun(x) / normFac)
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
#' model is $$\frac_{\sum_{i=1}^{N} c_i^1(O_t) * b_j(O_t) }{N * b_j(O_t)}$$
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
        nominator <- nominator + c_i_O_t[[state_i]] * b_j_O_t[[state_j]]
        denominator <- denominator + b_j_O_t[[state_j]]
      }
      
      predTemp[state_j] <- nominator / denominator
    }
    
    pred <- c(pred, names(which.max(predTemp)))
  }
  
  mmb::setWarnings(w)
  return(pred)
}


#' 2nd-order conditional density model (model B). The underlying
#' model is $$\frac_{\sum_{i=1}^{N} c_i^2(O_t) * c_i^1(O_t) * b_j(O_t) }{N * b_j(O_t)}$$
condDens_B_2ndOrder <- function(
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
  fixColName_c_1 <- gsub("_t_0$", "_t_1", stateColumn)
  densities_c_1 <- list()
  fixColName_c_2 <- gsub("_t_0$", "_t_2", stateColumn)
  densities_c_2 <- list()
  
  
  for (state in states) {
    condData <- mmb::conditionalDataMin(
      df = df, selectedFeatureNames = c(fixColName_c_1),
      features = mmb::createFeatureForBayes(name = fixColName_c_1, value = state))
    
    cn <- colnames(condData)
    featNames <- cn[grepl("_t_0$", cn)]
    
    densities_c_1 <- append(densities_c_1, estimateJointDensities(
      data = condData[, featNames],
      densFunSuffix = state,
      ignoreGeneration = TRUE,
      featuresPdfPmf = if (doEcdf) c() else featNames,
      featuresCdf = if (doEcdf) featNames else c()
    ))
    
    
    condData <- mmb::conditionalDataMin(
      df = df, selectedFeatureNames = c(fixColName_c_2),
      features = mmb::createFeatureForBayes(name = fixColName_c_2, value = state))
    
    cn <- colnames(condData)
    featNames <- cn[grepl("_t_0$", cn)]
    
    densities_c_2 <- append(densities_c_2, estimateJointDensities(
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
      c_i_2_O_t <- computeDensitiesSum(
        O_t = O_t, states = states,
        densities = densities_c_2, ignoreGeneration = TRUE)
      
      c_i_1_O_t <- computeDensitiesSum(
        O_t = O_t, states = states,
        densities = densities_c_1, ignoreGeneration = TRUE)
      
      nominator <- 0
      denominator <- 0
      
      for (state_i in states) {
        nominator <- nominator + c_i_2_O_t[[state_i]] * c_i_1_O_t[[state_i]] * b_j_O_t[[state_j]]
        denominator <- denominator + b_j_O_t[[state_j]]
      }
      
      predTemp[state_j] <- nominator / denominator
    }
    
    pred <- c(pred, names(which.max(predTemp)))
  }
  
  mmb::setWarnings(w)
  return(pred)
}





#' 3rd-order conditional density model (model B). The underlying
#' model is $$\frac_{\sum_{i=1}^{N} c_i^3(O_t) * c_i^2(O_t) * c_i^1(O_t) * b_j(O_t) }{N * b_j(O_t)}$$
condDens_B_3rdOrder <- function(
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
  fixColName_c_1 <- gsub("_t_0$", "_t_1", stateColumn)
  densities_c_1 <- list()
  fixColName_c_2 <- gsub("_t_0$", "_t_2", stateColumn)
  densities_c_2 <- list()
  fixColName_c_3 <- gsub("_t_0$", "_t_3", stateColumn)
  densities_c_3 <- list()
  
  
  for (state in states) {
    condData <- mmb::conditionalDataMin(
      df = df, selectedFeatureNames = c(fixColName_c_1),
      features = mmb::createFeatureForBayes(name = fixColName_c_1, value = state))
    
    cn <- colnames(condData)
    featNames <- cn[grepl("_t_0$", cn)]
    
    densities_c_1 <- append(densities_c_1, estimateJointDensities(
      data = condData[, featNames],
      densFunSuffix = state,
      ignoreGeneration = TRUE,
      featuresPdfPmf = if (doEcdf) c() else featNames,
      featuresCdf = if (doEcdf) featNames else c()
    ))
    
    
    condData <- mmb::conditionalDataMin(
      df = df, selectedFeatureNames = c(fixColName_c_2),
      features = mmb::createFeatureForBayes(name = fixColName_c_2, value = state))
    
    cn <- colnames(condData)
    featNames <- cn[grepl("_t_0$", cn)]
    
    densities_c_2 <- append(densities_c_2, estimateJointDensities(
      data = condData[, featNames],
      densFunSuffix = state,
      ignoreGeneration = TRUE,
      featuresPdfPmf = if (doEcdf) c() else featNames,
      featuresCdf = if (doEcdf) featNames else c()
    ))
    
    
    condData <- mmb::conditionalDataMin(
      df = df, selectedFeatureNames = c(fixColName_c_3),
      features = mmb::createFeatureForBayes(name = fixColName_c_3, value = state))
    
    cn <- colnames(condData)
    featNames <- cn[grepl("_t_0$", cn)]
    
    densities_c_3 <- append(densities_c_3, estimateJointDensities(
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
      c_i_3_O_t <- computeDensitiesSum(
        O_t = O_t, states = states,
        densities = densities_c_3, ignoreGeneration = TRUE)
      
      c_i_2_O_t <- computeDensitiesSum(
        O_t = O_t, states = states,
        densities = densities_c_2, ignoreGeneration = TRUE)
      
      c_i_1_O_t <- computeDensitiesSum(
        O_t = O_t, states = states,
        densities = densities_c_1, ignoreGeneration = TRUE)
      
      nominator <- 0
      denominator <- 0
      
      for (state_i in states) {
        nominator <- nominator + c_i_3_O_t[[state_i]] * c_i_2_O_t[[state_i]] * c_i_1_O_t[[state_i]] * b_j_O_t[[state_j]]
        denominator <- denominator + b_j_O_t[[state_j]]
      }
      
      predTemp[state_j] <- nominator / denominator
    }
    
    pred <- c(pred, names(which.max(predTemp)))
  }
  
  mmb::setWarnings(w)
  return(pred)
}




balanceDatasetSmote <- function(data, stateColumn) {
  lvls <- if (is.factor(data[[stateColumn]])) levels(data[[stateColumn]]) else NULL
  d <- table(data[[stateColumn]])
  m <- names(which.max(d))
  # We'll sample all other classes until we reach this for each:
  targetAmount <- d[[m]]
  
  # Get the other classes:
  otherClasses <- names(d)[!(names(d) %in% m)]
  
  # Add the over-represented class already to the final data:
  dataLargestClass <- data[data[[stateColumn]] == m, ]
  dataFinal <- dataLargestClass[, ]
  dataFinal[[stateColumn]] <- as.character(dataFinal[[stateColumn]])
  
  # Now, for each class, over-sample it and add to final frame:
  for (oc in otherClasses) {
    dataOtherClass <- data[data[[stateColumn]] == oc, ]
    temp <- rbind(dataLargestClass, dataOtherClass)
    
    # SMOTE requires factor-labels:
    temp[[stateColumn]] <- factor(temp[[stateColumn]])
    
    overSampled <- DMwR::SMOTE(
      form = formula(paste0(stateColumn, "~.")),
      data = temp,
      perc.over = 100 * ceiling(nrow(dataLargestClass) / nrow(dataOtherClass)),
      perc.under = 100
    )
    
    # Since we rounded up, let's only sample what we need:
    overSampled <- overSampled[overSampled[[stateColumn]] == oc, ]
    overSampled <- overSampled[sample(
      x = rownames(overSampled), size = min(nrow(overSampled), nrow(dataLargestClass))), ]
    
    # .. change to character again:
    overSampled[[stateColumn]] <- as.character(overSampled[[stateColumn]])
    dataFinal <- rbind(dataFinal, overSampled)
  }
  
  if (is.character(lvls)) {
    dataFinal[[stateColumn]] <- factor(dataFinal[[stateColumn]], levels = lvls)
  }
  
  return(dataFinal)
}




data2Densities_1stOrder <- function(
  states, data, stateColumn, split = 0.85, doEcdf = FALSE,
  normalizePdfs = FALSE, ecdfMinusOne = FALSE,
  returnDataOnly = FALSE
) {
  # First, we create a joined label for the data.
  stateColumn_t_1 <- gsub("_t_0$", "_t_1", stateColumn)
  stateColumn_joined <- paste0(stateColumn, "_d2d1")
  
  data[[stateColumn_joined]] <- ""
  for (rn in rownames(data)) {
    data[rn, stateColumn_joined] <- paste0(
      data[rn, stateColumn], data[rn, stateColumn_t_1])
  }
  
  
  # Second, make a split across the joined labels:
  temp <- caret::createDataPartition(
    y = data[[stateColumn_joined]], p = split, list = FALSE)
  # .. remove the joined label:
  data[[stateColumn_joined]] <- NULL
  
  # Temp will have no rows of split = 0:
  useDatasets <- c("train")
  if (nrow(temp) == 0) {
    train <- data[, ]
    valid <- data[c(), ]
  } else {
    train <- data[temp, ]
    valid <- data[-temp, ]
    useDatasets <- c(useDatasets, "valid")
    
    # As a help for inferncing, we will add an ID to every
    # observation in the valid-dataframe:
    valid$obsId <- 1:nrow(valid)
  }
  
  # Also add an observation-ID to train as helper:
  train$obsId <- 1:nrow(train)
  
  # Sometimes, we want the split functionality but
  # skip the densities, so just return the data:
  if (returnDataOnly) {
    return(list(
      train_data = train,
      valid_data = valid
    ))
  }
  
  
  # Third, estimate densities and PMFs over the train-data only:
  densities_ij <- list()
  
  for (state_i in states) {
    for (state_j in states) {
      condData <- mmb::conditionalDataMin(
        df = train, selectedFeatureNames = c(stateColumn, stateColumn_t_1),
        features = rbind(
          mmb::createFeatureForBayes(name = stateColumn, value = state_j),
          mmb::createFeatureForBayes(name = stateColumn_t_1, value = state_i)
        ))
      
      cn <- colnames(condData)
      
      densities_ij[[paste0(state_i, state_j)]] <- estimateJointDensities(
        data = condData,
        densFunSuffix = paste0(state_i, state_j),
        ignoreGeneration = FALSE,
        normalizePdf = normalizePdfs,
        ecdfMinusOne = ecdfMinusOne,
        featuresPdfPmf = if (doEcdf) c() else cn,
        featuresCdf = if (doEcdf) cn else c()
      )
    }
  }
  
  
  # Fourth, pass all data through the estimated functions. We need to
  # do this for the training data and the validation data as well.
  # We will create simple binary labels for when a sample matches
  # the current estimator's conditional data.
  
  # We gonna store the whole joined datasets in these:
  train_final <- NULL
  valid_final <- NULL
  
  for (state_i in states) {
    for (state_j in states) {
      for (useDataset in useDatasets) {
        # We are going to pass them through these:
        densFuns <- densities_ij[[paste0(state_i, state_j)]]
        # .. using this data:
        df <- if (useDataset == "train") train else valid
        
        # In 1st-order models like here, we can separate all data into
        # 4 classes: those where t-1,t-0 are correct, then those where
        # these are incorrect, and those where just either one is correct.
        # In the following, we'll call these data_11, data_00, data_01 etc.
        
        
        org_11 <- df[df[[stateColumn_t_1]] == state_i & df[[stateColumn]] == state_j, ]
        org_00 <- df[df[[stateColumn_t_1]] != state_i & df[[stateColumn]] != state_j, ]
        org_01 <- df[df[[stateColumn_t_1]] != state_i & df[[stateColumn]] == state_j, ]
        org_10 <- df[df[[stateColumn_t_1]] == state_i & df[[stateColumn]] != state_j, ]
        
        # Now pass all these through the density-estimation:
        temp <- colnames(org_11)[!(colnames(org_11) %in% "obsId")]
        data_11 <- computeLikelihoods(observations = org_11[, temp], densities = densFuns)
        data_00 <- computeLikelihoods(observations = org_00[, temp], densities = densFuns)
        data_01 <- computeLikelihoods(observations = org_01[, temp], densities = densFuns)
        data_10 <- computeLikelihoods(observations = org_10[, temp], densities = densFuns)
        
        # Rename all columns, by prefixing them with "d_" (for density):
        colnames(data_11) <- paste0("d_", colnames(data_11))
        colnames(data_00) <- paste0("d_", colnames(data_00))
        colnames(data_01) <- paste0("d_", colnames(data_01))
        colnames(data_10) <- paste0("d_", colnames(data_10))
        
        # Store which estimator we used:
        data_11$est_t1 <- which(states == state_i)
        data_00$est_t1 <- which(states == state_i)
        data_01$est_t1 <- which(states == state_i)
        data_10$est_t1 <- which(states == state_i)
        data_11$est_t0 <- which(states == state_j)
        data_00$est_t0 <- which(states == state_j)
        data_01$est_t0 <- which(states == state_j)
        data_10$est_t0 <- which(states == state_j)
        # .. and initialize the binary Y-values:
        data_11$y_t1 <- 1
        data_11$y_t0 <- 1
        data_00$y_t1 <- 0
        data_00$y_t0 <- 0
        data_01$y_t1 <- 0
        data_01$y_t0 <- 1
        data_10$y_t1 <- 1
        data_10$y_t0 <- 0
        
        # .. we change this to use true one-hot activation, with the following
        # mapping describing which 'bit' is set to 1 (1 through 4):
        # 11: 1, 00: 2, 01: 3, 10: 4
        #data_11$y_oh1 <- 1
        #data_11$y_oh2 <- 0
        #data_11$y_oh3 <- 0
        #data_11$y_oh4 <- 0
        # data-00:
        #data_00$y_oh1 <- 0
        #data_00$y_oh2 <- 1
        #data_00$y_oh3 <- 0
        #data_00$y_oh4 <- 0
        # data-01:
        #data_01$y_oh1 <- 0
        #data_01$y_oh2 <- 0
        #data_01$y_oh3 <- 1
        #data_01$y_oh4 <- 0
        # data-10:
        #data_10$y_oh1 <- 0
        #data_10$y_oh2 <- 0
        #data_10$y_oh3 <- 0
        #data_10$y_oh4 <- 1
        
        
        temp <- rbind(
          cbind(org_11, data_11),
          cbind(org_00, data_00),
          cbind(org_01, data_01),
          cbind(org_10, data_10))
        
        # .. and append to correct data-frame:
        if (useDataset == "train") {
          train_final <- rbind(train_final, temp)
        } else {
          valid_final <- rbind(valid_final, temp)
        }
      }
    }
  }
  
  return(list(
    densities = densities_ij,
    train = train_final,
    valid = valid_final,
    train_data = train,
    valid_data = valid
  ))
}




condDens_D_1stOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, normalizePdfs = TRUE,
  logTol = 1e-3
) {
  w <- mmb::getWarnings()
  
  lvls <- if (is.factor(data[[stateColumn]])) levels(data[[stateColumn]]) else NULL
  
  # Creates aa, ac, .., pp etc.
  states_1stOrder <- sort(apply(expand.grid(list(
    a = states,
    b = states
  )), 1, function(x) paste0(x[1], x[2])))
  
  d2d <- data2Densities_1stOrder(
    states = states, data = data, stateColumn = stateColumn,
    split = 0, doEcdf = doEcdf, normalizePdfs = normalizePdfs)
  
  pred <- c()
  for (rno in rownames(observations)) {
    est <- c()
    for (s in states_1stOrder) {
      est[[s]] <- computeDensitiesSum(
        O_t = observations[rno, ],
        states = states_1stOrder,
        densities = d2d$densities[[s]])[[s]]
    }
    
    ######## IGNORE EVERYTHING, JUST PICK THE MAX!
    #pred <- c(pred, stringi::stri_sub(names(which.max(est)), from = -1))
    #next
    
    
    
    
    
    # Now combine all states that end in a specific t-0, so that
    # we can do a group-voting. We will do a data-frame, in which
    # we store (per column) state at t-0, with the first row
    # holding the sum of all voters in the group, the second row
    # holding the standard-deviation of votes within the group.
    # In the 3rd row we store the log10 of the sum of votes.
    voteMat <- matrix(nrow = 3, ncol = length(states))
    for (s in states) {
      grpData <- est[grepAnyOrAll(paste0(s, "$"), names(est))]
      grpIdx <- which(states == s)
      
      voteMat[1, grpIdx] <- sum(grpData)
      voteMat[2, grpIdx] <- sd(grpData)
      voteMat[3, grpIdx] <- log10(voteMat[1, grpIdx])
    }
    
    # Now that we have the vote-matrix, we regard all votes as
    # equal that have their log10 closer than some threshold
    # to the log10 of the highest group's sum. If there is more
    # than one such group, we the group with the lowest std-dev,
    # indicating the most consistent votes within the group.
    maxLog <- voteMat[3, which.max(voteMat[3, ])]
    equalLogs <- which(voteMat[3, ] >= (maxLog - logTol))
    minSd <- which.min(voteMat[2, equalLogs])
    
    pred <- c(pred, states[equalLogs[minSd]])
  }
  
  if (!is.null(lvls)) {
    pred <- factor(pred, levels = lvls)
  }
  
  
  mmb::setWarnings(w)
  return(pred)
}




#' 1st-order conditional density model (model C). The underlying
#' model is $$\sum_{i=1}^{N} e_{ij}(O_t)$$
condDens_C_1stOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE
) {
  # Here, we use the e_i functions, which fix the data at multiple
  # h[,i[,j]] points, but always segment over the t-0 variables.
  
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  
  densities_e_ij <- list()
  
  for (state_i in states) {
    for (state_j in states) {
      
      stateColumn_e_1 <- gsub("_t_0$", "_t_1", stateColumn)
      
      condData <- mmb::conditionalDataMin(
        df = df, selectedFeatureNames = c(stateColumn, stateColumn_e_1),
        features = rbind(
          mmb::createFeatureForBayes(name = stateColumn, value = state_j),
          mmb::createFeatureForBayes(name = stateColumn_e_1, value = state_i)
        ))
      
      cn <- colnames(condData)
      featNames <- cn[grepl("_t_0$", cn)]
      
      densities_e_ij <- append(densities_e_ij, estimateJointDensities(
        data = condData[, featNames],
        densFunSuffix = paste0(state_i, state_j),
        ignoreGeneration = TRUE,
        featuresPdfPmf = if (doEcdf) c() else featNames,
        featuresCdf = if (doEcdf) featNames else c()
      ))
    }
  }
  
  
  pred <- c()
  for (j in 1:nrow(observations)) {
    
    O_t <- observations[j, ]
    predTemp <- c()
    
    for (state_j in states) {
      temp <- 0
      for (state_i in states) {
        s <- paste0(state_i, state_j)
        temp <- temp + computeDensitiesSum(
          O_t = O_t, densities = densities_e_ij,
          states = c(s), ignoreGeneration = TRUE)[s]
      }
      
      predTemp[state_j] <- temp
    }
    
    pred <- c(pred, names(which.max(predTemp)))
  }
  
  mmb::setWarnings(w)
  return(pred)
}



#' 1st-order conditional density model (model C). The underlying
#' model is $$\sum_{i=1}^{N} e_{ij}(O_t)$$ and we are using its
#' output in a 2nd stage by attaching weights.
condDens_C_1stOrder_weighted <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE
) {
  # Here, we use the e_i functions, which fix the data at multiple
  # h[,i[,j]] points, but always segment over the t-0 variables.
  
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  # We store the data as likelihoods in this one:
  dataLLH <- data.frame()
          # In this one, we store the likelihoods of the data we are
          # about to predict. The prediction result will be a one-hot
          # encoded matrix which we can then decode into labels.
          #predLLH <- data.frame()
  
  densities_e_ij <- list()
  
  for (state_i in states) {
    for (state_j in states) {
      
      stateColumn_e_1 <- gsub("_t_0$", "_t_1", stateColumn)
      
      condData <- mmb::conditionalDataMin(
        df = df, selectedFeatureNames = c(stateColumn, stateColumn_e_1),
        features = rbind(
          mmb::createFeatureForBayes(name = stateColumn, value = state_j),
          mmb::createFeatureForBayes(name = stateColumn_e_1, value = state_i)
        ))
      
      cn <- colnames(condData)
      featNames <- cn[grepl("_t_0$", cn)]
      
      tempDens <- estimateJointDensities(
        data = condData[, featNames],
        densFunSuffix = paste0(state_i, state_j),
        ignoreGeneration = TRUE,
        featuresPdfPmf = if (doEcdf) c() else featNames,
        featuresCdf = if (doEcdf) featNames else c(),
        normalizePdf = TRUE # important
      )
      
      densities_e_ij <- append(densities_e_ij, tempDens)
      
      # Now let's compute the relative likelihoods:
      llhTemp <- computeLikelihoods(
        observations = condData[, featNames],
        densities = tempDens,
        ignoreGenerations = TRUE)
      
      # Also, we need to add the one-hot coding of state_j:
      llhTemp <- cbind(llhTemp, nominalToOneHot(
        labels = condData[[stateColumn]], states = states))
      
      dataLLH <- rbind(dataLLH, llhTemp)
    }
  }
  
  
  # Let's train a neural network over our data:
  dataLLH <- dataLLH[, !(colnames(dataLLH) %in% caret::nzv(dataLLH, names = TRUE))]
  nnet <- neuralnet::neuralnet(
    formula = paste0(paste(paste0(paste0("`__@onehot@__", states), "`"), collapse = " + "), " ~."),
    data = dataLLH,
    #hidden = floor((length(colnames(dataLLH)) - length(states)) / 2),
    hidden = floor(log2(ncol(dataLLH) - length(states))),
    threshold = 0.01,
    stepmax = 5e5,
    # approx-reLu; < -10 not working; https://stackoverflow.com/questions/34532878/package-neuralnet-in-r-rectified-linear-unit-relu-activation-function
    #act.fct = function(x) x / (1 + exp(-10 * x)),
    #act.fct = function(x) (1 / 1 + exp(-1 * x)),
    # Currently, this does not work. I am not sure why, maybe for 'ce' the
    # labels must be categorical. The default, SSE ('sse'), will work, too.
    #err.fct = 'ce',
    #err.fct = function(t, s) {
    #  # https://gombru.github.io/2018/05/23/cross_entropy_loss/
    #  return(-1 * (t[1] * log(s[1]) + t[2] * log(s[2]) + t[3] * log(s[3])))
    #},
    lifesign = "full"
  )
  
  plot(nnet)
  
  colnames(observations) <- gsub("_t_\\d+", "", colnames(observations))
  pred <- matrix(nrow = 0, ncol = 5)
  colnames(pred) <- c("i", "j", states)
  for (rn in rownames(observations)) {
    obs <- observations[rn, colnames(observations) %in% colnames(dataLLH)]
    for (state_i in states) {
      for (state_j in states) {
        densFuns <- densities_e_ij[grepl(
          paste0(state_i, state_j, "$"), names(densities_e_ij))]
        
        asLLH <- computeLikelihoods(
          observations = obs,
          densities = densFuns,
          ignoreGenerations = TRUE)
        
        pred <- rbind(pred, c(
          which(states == state_i),
          which(states == state_j),
          predict(nnet, newdata = asLLH)[1, ]
        ))
      }
    }
  }
  
  mmb::setWarnings(w)
  return(pred)
}




#' 2nd-order conditional density model (model C). The underlying
#' model is $$\sum_{i=1}^{N} e_{hij}(O_t)$$
condDens_C_2ndOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE
) {
  # Here, we use the e_i functions, which fix the data at multiple
  # h[,i[,j]] points, but always segment over the t-0 variables.
  
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  
  densities_e_hij <- list()
  
  for (state_h in states) {
    for (state_i in states) {
      for (state_j in states) {
        
        stateColumn_e_1 <- gsub("_t_0$", "_t_1", stateColumn)
        stateColumn_e_2 <- gsub("_t_0$", "_t_2", stateColumn)
        
        condData <- mmb::conditionalDataMin(
          df = df, selectedFeatureNames = c(stateColumn, stateColumn_e_1, stateColumn_e_2),
          features = rbind(
            mmb::createFeatureForBayes(name = stateColumn, value = state_j),
            mmb::createFeatureForBayes(name = stateColumn_e_1, value = state_i),
            mmb::createFeatureForBayes(name = stateColumn_e_2, value = state_h)
          ))
        
        cn <- colnames(condData)
        featNames <- cn[grepl("_t_0$", cn)]
        
        densities_e_hij <- append(densities_e_hij, estimateJointDensities(
          data = condData[, featNames],
          densFunSuffix = paste0(state_h, state_i, state_j),
          ignoreGeneration = TRUE,
          featuresPdfPmf = if (doEcdf) c() else featNames,
          featuresCdf = if (doEcdf) featNames else c()
        ))
      }
    }
  }
  
  
  pred <- c()
  for (j in 1:nrow(observations)) {
    
    O_t <- observations[j, ]
    predTemp <- c()
    
    for (state_j in states) {
      temp <- 0
      for (state_i in states) {
        for (state_h in states) {
          s <- paste0(state_h, state_i, state_j)
          temp <- temp + computeDensitiesSum(
            O_t = O_t, densities = densities_e_hij,
            states = c(s), ignoreGeneration = TRUE)[s]
        }
      }
      
      predTemp[state_j] <- temp
    }
    
    pred <- c(pred, names(which.max(predTemp)))
  }
  
  mmb::setWarnings(w)
  return(pred)
}




