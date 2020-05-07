
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




condDens_A_1stOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE)
{
  if (nrow(observations) < 2) {
    stop("Need at least 2 observations to predict.")
  }
  
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
  b_j_O_1 <- computeDensitiesSum(
    O_t = O_1, states = states, densities = densities_b)
  
  predTemp <- c()
  for (state in states) {
    predTemp[state] <- initProbs[state] * b_j_O_1[state]
  }
  O_1[[stateColumn]] <- names(which.max(predTemp))
  pred <- c(O_1[[stateColumn]])
  
  # Now assign the label of the greatest likelihood to each
  # observation:
  O_t_1 <- O_1
  for (j in 2:nrow(observations)) { # we start at 2, as O_1 is done already
    
    # Now we have a bit of trouble. d_j and b_j were using the same
    # segmented dataset, but b_j was estimated over all t_0 variables,
    # and d_j over t_1 variables, meaning it stores the conditional
    # densities of such observations, that had t_0=j as parent. Now
    # when we feed the preceding observation to d_j, it would look at
    # those t_1 variables. However, this is wrong because the actual
    # data resides in O_t_1's t_0 variables! We need to delete those
    # t_1 variables and rename t_0 to t_1, as computeDensitiesSum
    # expects t_1-named variables.
    #O_temp <- O_t_1[1, !grepl("t_1$", colnames(O_t_1))]
    #colnames(O_temp) <- gsub("t_0$", "t_1", colnames(O_temp))
    #O_t_1 <- O_temp
    
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
  
  mmb::setWarnings(w)
  return(pred)
}

condDens_A_2ndOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE)
{
}