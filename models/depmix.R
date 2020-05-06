
depmixMsg <- TRUE
depmixSetMessages <- function(enable = TRUE) depmixMsg <-- enable
depmixGetMessages <- function() depmixMsg


estimateDensities <- function(
  states, stateColumn, data = data.frame(), featuresPdfPmf = c(), featuresCdf = c()
) {
  if (length(states) == 0) {
    stop("No states were given.")
  }
  if (!is.character(stateColumn) || length(stateColumn) == 0) {
    stop("The state's column-name is not correct.")
  }
  if (!is.data.frame(data) || ncol(data) < 1) {
    stop("The data given is not a data.frame or too sparse.")
  }
  
  if (length(featuresPdfPmf) > 0 && length(featuresCdf) > 0) {
    warning("Estimating PDF and CDF, be careful using them together.")
  }
  
  estFuncs <- list()
  cn <- colnames(data)
  
  for (state in states) {
    # We always segment all densities by the current label (state) of State t0:
    condData <- mmb::conditionalDataMin(
      df = data,
      features = mmb::createFeatureForBayes(name = stateColumn, value = state),
      selectedFeatureNames = c(stateColumn)
    )
    
    for (feat in c(featuresPdfPmf, featuresCdf)) {
      if (!(feat %in% cn)) {
        stop(paste("The feature", feat, "is not in the data.frame."))
      }
      
      densFunName <- paste(feat, state, sep = "_@dens@_")
      isDiscrete <- !is.numeric(data[[feat]])
      doEcdf <- feat %in% featuresCdf
      
      if (isDiscrete) {
        estFuncs[[feat]] <- (function() {
          temp <- as.character(condData[[feat]])
          return(function(x) mmb::getProbForDiscrete(temp, x))
        })()
      } else if (doEcdf) {
        estFuncs[[feat]] <- stats::ecdf(condData[[feat]])
      } else {
        estFuncs[[feat]] <- (function() {
          pdf <- mmb::estimatePdf(condData[[feat]])
          return(function(x) pdf$fun(x))
        })()
      }
    }
  }
  
  return(estFuncs)
}

getInitialStateProbs <- function(data, stateColumn) {
  
  probs <- c()
  labelData <- table(as.character(data[[stateColumn]]))
  labels <- names(labelData)
  
  for (label in labels) {
    probs[[label]] <- labelData[[label]] / sum(labelData)
  }
  
  return(probs)
}


getTransprobs_1stOrder <- function(states, data, stateColumn) {
  numStates <- length(states)
  transprobs_1stOrder <- matrix(data = 0, nrow = numStates, ncol = numStates)
  colnames(transprobs_1stOrder) <- states
  rownames(transprobs_1stOrder) <- states
  
  label_t_0 <- stateColumn
  label_t_1 <- paste(stateColumn, "t_1", sep = "_")
  
  df <- data[, ]
  df[label_t_0] <- as.character(df[[label_t_0]])
  df[label_t_1] <- as.character(df[[label_t_1]])
  
  for (t_1 in states) { # column-wise
    for (t_0 in states) {
      # Sum how often we went from t_1 to t_0
      transprobs_1stOrder[t_0, t_1] <- transprobs_1stOrder[t_0, t_1] +
        sum(df[label_t_1] == t_1 & df[label_t_0] == t_0)
    }
    
    # Normalize all options for ending up in t_0 coming from t_1:
    transprobs_1stOrder[, t_1] <- transprobs_1stOrder[, t_1] /
      sum(transprobs_1stOrder[, t_1])
  }
  
  return(transprobs_1stOrder)
}



#' Computes the sum of densities, by taking a vector of density
#' functions and computes the likelihood using an observation's
#' random variables' values (similar to a dot-product).
computeDensitiesSum <- function(O_t, states, densities, smoothing = 0.1) {
  # We have to compute the likelihood for all possible states
  # in order to find the maximum later.
  allJs <- c()
  for (state in states) {
    densFactors <- c()
    for (densFunName in names(densities)) {
      featName <- strsplit(densFunName, split = "_@dens@_")[[1]][1]
      densFun <- densities[[densFunName]]
      densFactors[[featName]] <- smoothing + densFun(O_t[[featName]])
    }
    
    allJs[[state]] <- sum(densFactors)
  }
  
  return(allJs)
}



#' 1st-order Forward-algorithm for maximum-likelihood estimation of
#' hidden states. Builds initial state- and transition-probabilities
#' from the given data.
#' 
#' @param states character vector with all possible states (labels).
#' @param data data.frame with all known observations/states. This is
#' labeled data, and the label is to be found in the column designated
#' by the argument \code{stateColumn}. The data should ideally be ar-
#' bitrary many arbitrary long chains of consecutively observed data-
#' points. Also, it is expected the data has a second column, called
#' the same as the \code{stateColumn}, followed by "_t_1", to identify
#' the state/label of the previous observation. This is needed to build
#' transition probabilities.
#' @param observations data.frame with observations, where the oldest
#' observation is the first row, and the most recent is the last row.
#' Uses the forward-algorithm to assign labels. If this data.frame
#' consists of exactly two observations and the first has a label,
#' then the the label for the second observation is estimated using
#' the short form: $\phi_t^1(i,j) = A_{ij} * b_j(O_t)$.
#' @param stateColumn character, the name of the column holding the
#' states. These need to be discrete.
#' @param smoothing double, defaults to 0.1. This is Laplacian/Lidstone
#' smoothing and a constant factor added to each density.
#' @param doEcdf default FALSE, whether to use empirical CDF instead of
#' empirical PDF for continuous features.
#' @param returnLogLikelihood default FALSE, if true, returns the log-
#' likelihood of the entire sequence.
#' @return character vector with most likely labels for the given ob-
#' servations, in the same order.
depmixForward_1stOrder <- function(
  states, data, observations, stateColumn, smoothing = 0.1,
  doEcdf = FALSE, returnLogLikelihood = FALSE)
{
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  initProbs <- getInitialStateProbs(
    data = df, stateColumn = stateColumn)
  transProbs <- getTransprobs_1stOrder(
    states = states, data = df, stateColumn = stateColumn)
  
  if (depmixGetMessages()) {
    print(initProbs)
    print(transProbs)
  }
  
  densityCols <- colnames(df)
  densityCols <- densityCols[!(densityCols %in% stateColumn)]
  densities <- estimateDensities(
    states = states,
    stateColumn = stateColumn,
    data = df,
    featuresPdfPmf = if (doEcdf) c() else densityCols,
    featuresCdf = if (doEcdf) densityCols else c()
  )
  
  # Before we go to the forward-algorithm, let's check the special
  # case: 2 observations, and t-1 is labeled:
  if (nrow(observations) == 2 &&
      (utils::head(observations, 1)[[stateColumn]] %in% states))
  {
    warning("Using short-form for special case.")
    O_t_1 <- utils::head(observations, 1)
    O_t <- utils::tail(observations, 1)
    # This is a named vector, we need to multiply each entry still
    # with the transition probability.
    b_j_O_t <- computeDensitiesSum(
      O_t = O_t, states = states, densities = densities, smoothing = smoothing)
    state_t_1 <- O_t_1[[stateColumn]]
    
    for (state in states) {
      # Let's reuse this vector..
      b_j_O_t[[state]] <- b_j_O_t[[state]] * transProbs[state, state_t_1]
    }
    
    if (depmixGetMessages()) {
      print(b_j_O_t)
    }
    
    if (returnLogLikelihood) {
      return(log(sum(b_j_O_t)))
    }
    
    # Return the name of the state with highest likelihood:
    return(names(which.max(b_j_O_t)))
  }
  
  

  # Now we got everything ready, we can go ahead and start the
  # forward algorithm. For it, we need a matrix to store the
  # likelihoods of t-1 states.
  prevLh <- matrix(nrow = length(states), ncol = nrow(observations))
  rownames(prevLh) <- states
  
  for (i in 1:nrow(observations)) {
    # Those two have always to be done:
    O_t <- observations[i, ]
    b_j_O_t <- computeDensitiesSum(
      O_t = O_t, states = states, densities = densities, smoothing = smoothing)
    
    # These are constant for every possible i (if i > 1):
    sumPrevLh <- 0
    if (i > 1) {
      sumPrevLh <- sum(prevLh[, i - 1])
    }
    
    for (state in states) {
      if (i == 1) {
        # Compute \phi_1(j)
        prevLh[state, i] <- initProbs[state] * b_j_O_t[[state]]
      } else {
        # Compute \phi_t(j)
        temp <- 0
        for (state_t_1 in states) {
          temp <- temp +
            (prevLh[state_t_1, i - 1] * transProbs[state, state_t_1] * b_j_O_t[[state]])
        }
        prevLh[state, i] <- temp / sumPrevLh
      }
    }
  }
  
  mmb::setWarnings(w)
  
  if (depmixGetMessages()) {
    print(prevLh)
  }
  
  if (returnLogLikelihood) {
    # The log-likelihood of the entire forecast is the sum over the
    # logs of sums of \phi_t(j).
    return(sum(apply(prevLh, 2, function(col) log(sum(col)))))
  }
  
  return(sapply(unname(apply(prevLh, 2, function(col) which.max(col))),
                function(i) rownames(prevLh)[i]))
}




#' 2nd-order Forward-algorithm for maximum-likelihood estimation of
#' hidden states. Builds initial state- and transition-probabilities
#' from the given data.
depmixForward_2ndOrder <- function(
  states, data, observations, stateColumn, smoothing = 0.1, doEcdf = FALSE)
{
  
}



























