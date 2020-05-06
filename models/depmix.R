#' Creates a closure over a variable and returns its getter and setter.
#' 
#' @author Sebastian Hönel <sebastian.honel@lnu.se>
#' @param initVarVal the initial value of the closed variable.
#' @keywords internal
#' @return list with entries 'get' and 'set' which are getter/setter for the
#' variable that a closure was made over.
make.varClosure <- function(initVarVal = NULL) {
  varClosed <- initVarVal
  
  return(list(
    get = function() varClosed,
    set = function(val) varClosed <<- val
  ))
}


varMsg <- make.varClosure(TRUE)
depmixSetMessages <- function(enable = TRUE) varMsg$set(!!enable)
depmixGetMessages <- function() varMsg$get()


#' Creates conditional density functions for each feature of the data,
#' segmented for by each state of the data. For continuous features,
#' this estimates a PDF or eCDF, and for discrete features, it estimates
#' a PMF.
#' 
#' @author Sebastian Hönel <sebastian.honel@lnu.se>
#' @param states character vector of available states
#' @param data data.frame with observations
#' @param featuresPdfPmf character vector with names of features to build
#' a PDF or PMF for.
#' @param featuresCdf character vector with names of features to build an
#' eCDF for (only continuous features allowed).
#' @return list of density functions, that follow the name-scheme
#' \code{paste(featureName, state, sep = "_@dens@_")}.
estimateDepmixDensities <- function(
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
        estFuncs[[feat]] <- (function() {
          tryCatch({
            temp <- condData[[feat]]
            ecdf(temp)
          }, error=function(cond) function(x) 0)
        })()
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


#' Compute the initial state probabilities.
#' 
#' @author Sebastian Hönel <sebastian.honel@lnu.se>
#' @param data data.frame with observations and states assigned (i.e.,
#' labeled data).
#' @param stateColumn character, the name of the column holding the
#' state/label.
#' @return named vector with probabilities (sums to 1).
getInitialStateProbs <- function(data, stateColumn) {
  
  probs <- c()
  labelData <- table(as.character(data[[stateColumn]]))
  labels <- names(labelData)
  
  for (label in labels) {
    probs[[label]] <- labelData[[label]] / sum(labelData)
  }
  
  return(probs)
}



#' Creates transition probabilities for a 2nd-order dependency model.
#' 
#' @author Sebastian Hönel <sebastian.honel@lnu.se>
#' @param states character vector of available states
#' @param data data.frame with all known observations/states. This is
#' labeled data, and the label is to be found in the column designated
#' by the argument \code{stateColumn}. The data should ideally be ar-
#' bitrary many arbitrary long chains of consecutively observed data-
#' points. Also, it is expected the data has a second column, called
#' the same as the \code{stateColumn}, followed by "_t_1", to identify
#' the state/label of the previous observation. This is needed to build
#' transition probabilities.
#' @param stateColumn the name of the column in the data that holds
#' the state (label).
#' @return matrix with column- and row-names according to the states.
#' The first index is t-0, the 2nd is t-1.
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



#' Creates transition probabilities for a 2nd-order dependency model.
#' 
#' @author Sebastian Hönel <sebastian.honel@lnu.se>
#' @param states character vector of available states
#' @param data data.frame with observations. Like for @seealso
#' \code{getTransprobs_1stOrder()}, the data is expected to have additonal
#' columns for t-1 and t-2 state-labels (see description there).
#' @param stateColumn the name of the column in the data that holds
#' the state (label).
#' @param returnTransprobsOnly default FALSE, returns a list with the
#' two entries 'mapping', 'transprobs' by default. The mapping is created
#' by indexing the states from \code{1:length(states)} and is useful to
#' access the tensor by characters for indices. If TRUE, then this function
#' only returns transition-tensor.
#' @return list or transition-tensor (\code{tensorr::dtensor}). Depends on
#' how \code{returnTransprobsOnly} was set. The 1st index of the tensor is
#' t-0, the 2nd is t-1, and the 3rd then is t-2. E.g., to go from a over b
#' to c, the index would be A[c,b,a].
getTransprobs_2ndOrder <- function(states, data, stateColumn, returnTransprobsOnly = FALSE) {
  install.packagesCond("tensorr")
  library("tensorr")
  
  L <- length(states)
  # Create a dense LxLxL tensor
  transprobs_2ndOrder <- dtensor(array(data = 0, dim = c(L,L,L)))
  dimnames(transprobs_2ndOrder) <- list(states, states, states)
  
  label_t_0 <- stateColumn
  label_t_1 <- paste(stateColumn, "t_1", sep = "_")
  label_t_2 <- paste(stateColumn, "t_2", sep = "_")
  
  df <- data[, ]
  df[label_t_0] <- as.character(df[[label_t_0]])
  df[label_t_1] <- as.character(df[[label_t_1]])
  df[label_t_2] <- as.character(df[[label_t_2]])
  
  # Also, let's create a mapping between states and their numeric index:
  m <- c()
  for (i in 1:L) {
    m[[states[i]]] <- i
  }
  
  for (t_2 in states) {
    i2 <- m[t_2]
    for (t_1 in states) {
      i1 <- m[t_1]
      for (t_0 in states) {
        i0 <- m[t_0]
        
        # Sum how often we went from t_2, over t_1, to t_0
        transprobs_2ndOrder[i0, i1, i2] <- transprobs_2ndOrder[i0, i1, i2] +
          sum(df[label_t_2] == t_2 &
              df[label_t_1] == t_1 &
              df[label_t_0] == t_0)
      }
    }
    
    # Normalize each 3x3x1 tensor:
    n <- transprobs_2ndOrder[,, i2] / sum(transprobs_2ndOrder[,, i2])
    transprobs_2ndOrder[,, i2] <- array(n, dim = dim(n))
  }
  
  if (returnTransprobsOnly) {
    return(transprobs_2ndOrder)
  }
  
  return(list(
    mapping = m,
    transprobs = transprobs_2ndOrder
  ))
}


#' Computes the sum of densities, by taking a vector of density
#' functions and computes the likelihood using an observation's
#' random variables' values (similar to a dot-product).
#' 
#' @author Sebastian Hönel <sebastian.honel@lnu.se>
#' @param O_t data.frame with one row, holding a value for each feature.
#' @param states character vector of available states
#' @param densities list of density functions, as obtained by @seealso
#' \code{estimateDepmixDensities()}.
#' @return named vector with computed and summed up densities
#' for each state.
computeDensitiesSum <- function(O_t, states, densities) {
  # We have to compute the likelihood for all possible states
  # in order to find the maximum later.
  allJs <- c()
  for (state in states) {
    densFactors <- c()
    for (densFunName in names(densities)) {
      featName <- strsplit(densFunName, split = "_@dens@_")[[1]][1]
      densFun <- densities[[densFunName]]
      densFactors[[featName]] <- densFun(O_t[[featName]])
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
#' @param doEcdf default FALSE, whether to use empirical CDF instead of
#' empirical PDF for continuous features.
#' @param returnLogLikelihood default FALSE, if true, returns the log-
#' likelihood of the entire sequence.
#' @return character vector with most likely labels for the given ob-
#' servations, in the same order.
depmixForward_1stOrder <- function(
  states, data, observations, stateColumn,
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
  densities <- estimateDepmixDensities(
    states = states,
    stateColumn = stateColumn,
    data = df,
    featuresPdfPmf = if (doEcdf) c() else densityCols,
    featuresCdf = if (doEcdf) densityCols else c()
  )
  
  # Before we go to the forward-algorithm, let's check the special
  # case: 2 observations, and t-1 is labeled:
  if (nrow(observations) == 2 &&
      (observations[1, stateColumn] %in% states))
  {
    warning("Using short-form for special case.")
    O_t_1 <- observations[1, ]
    O_t <- observations[2, ]
    # This is a named vector, we need to multiply each entry still
    # with the transition probability.
    b_j_O_t <- computeDensitiesSum(
      O_t = O_t, states = states, densities = densities)
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
      O_t = O_t, states = states, densities = densities)
    
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
#' 
#' @param states character vector with all possible states (labels).
#' @param data data.frame with all known observations/states. This is
#' labeled data, and the label is to be found in the column designated
#' by the argument \code{stateColumn}. The data should ideally be ar-
#' bitrary many arbitrary long chains of consecutively observed data-
#' points. Also, it is expected the data has a second column, called
#' the same as the \code{stateColumn}, followed by "_t_1", to identify
#' the state/label of the previous observation. This is needed to build
#' transition probabilities. Also, since this is a 2nd-order model, we
#' need to have additional columns ending in "_t_2".
#' @param observations data.frame with observations, where the oldest
#' observation is the first row, and the most recent is the last row.
#' Uses the forward-algorithm to assign labels. If this data.frame
#' consists of exactly two observations and the first has a label,
#' then the the label for the second observation is estimated using
#' the short form: $\phi_t^1(i,j) = A_{ij} * b_j(O_t)$.
#' @param stateColumn character, the name of the column holding the
#' states. These need to be discrete.
#' @param doEcdf default FALSE, whether to use empirical CDF instead of
#' empirical PDF for continuous features.
#' @param returnLogLikelihood default FALSE, if true, returns the log-
#' likelihood of the entire sequence.
#' @return character vector with most likely labels for the given ob-
#' servations, in the same order.
depmixForward_2ndOrder <- function(
  states, data, observations, stateColumn,
  doEcdf = FALSE, returnLogLikelihood = FALSE)
{
  w <- mmb::getWarnings()
  # Because otherwise, mmb will likely warn a lot about scarce data.
  mmb::setWarnings(FALSE)
  
  df <- data[, ]
  df[stateColumn] <- as.character(df[[stateColumn]])
  
  initProbs <- getInitialStateProbs(
    data = df, stateColumn = stateColumn)
  temp <- getTransprobs_2ndOrder(
    states = states, data = df, stateColumn = stateColumn)
  transProbs <- temp[["transprobs"]]
  transProbs_mapping <- temp[["mapping"]]
  
  if (depmixGetMessages()) {
    print(initProbs)
    print(transProbs)
  }
  
  densityCols <- colnames(df)
  densityCols <- densityCols[!(densityCols %in% stateColumn)]
  densities <- estimateDepmixDensities(
    states = states,
    stateColumn = stateColumn,
    data = df,
    featuresPdfPmf = if (doEcdf) c() else densityCols,
    featuresCdf = if (doEcdf) densityCols else c()
  )
  
  # Before we go to the forward-algorithm, let's check the special
  # case: 3 observations, and t-2,t-1 are labeled:
  if (nrow(observations) == 3 &&
      (observations[1:2, stateColumn] %in% states))
  {
    warning("Using short-form for special case.")
    O_t_2 <- observations[1, ]
    O_t_1 <- observations[2, ]
    O_t <- observations[3, ]
    # This is a named vector, we need to multiply each entry still
    # with the transition probability tensor later.
    b_j_O_t <- computeDensitiesSum(
      O_t = O_t, states = states, densities = densities)
    
    m <- transProbs_mapping
    state_t_1 <- m[O_t_1[[stateColumn]]]
    state_t_2 <- m[O_t_2[[stateColumn]]]
    
    for (state in states) {
      # Let's reuse this vector..
      b_j_O_t[[state]] <- b_j_O_t[[state]] *
        transProbs[m[state], state_t_1, state_t_2]
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
  
  # Also, for a 2nd-order model, we need to be able to compute
  # \phi_2(j), so we need an extra 2D transition matrix for that.
  transprobs_1stOrder <-getTransprobs_1stOrder(
    states = states, data = df, stateColumn = stateColumn)
  m <- transProbs_mapping
  
  for (i in 1:nrow(observations)) {
    # Those two have always to be done:
    O_t <- observations[i, ]
    b_j_O_t <- computeDensitiesSum(
      O_t = O_t, states = states, densities = densities)
    
    # These are constant for every possible i > 1:
    sumPrevLh <- 0 # inference first observation
    if (i == 2) {  # one prev. observation known/labeled 
      sumPrevLh <- sum(prevLh[, i - 1])
    } else if (i > 2) { # two or more previous observations known
      for (state_t_2 in states) {
        for (state_t_1 in states) {
          sumPrevLh <- sumPrevLh +
            (prevLh[state_t_2, i - 2] * prevLh[state_t_1, i - 1])
        }
      }
    }
    
    for (state in states) {
      if (i == 1) {
        # Compute \phi_1(j)
        prevLh[state, i] <- initProbs[state] * b_j_O_t[[state]]
      } else if (i == 2) {
        # Compute \phi_2(j)
        temp <- 0
        for (state_t_1 in states) {
          temp <- temp +
            (prevLh[state_t_1, i - 1] *
               transprobs_1stOrder[state, state_t_1] *
               b_j_O_t[[state]])
        }
        prevLh[state, i] <- temp / sumPrevLh
      } else {
        # Compute \phi_t(j)
        temp <- 0
        for (state_t_2 in states) {
          st_2 <- m[state_t_2]
          for (state_t_1 in states) {
            st_1 <- m[state_t_1]
            temp <- temp + (b_j_O_t[[state]] *
                            prevLh[state_t_1, i - 1] *
                            prevLh[state_t_2, i - 2] *
                            transProbs[m[[state]], st_1, st_2])
          }
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

