# As returned by a tryCatch-construct on error
is.error <- function(e) is.list(e) && all(c("message", "call") %in% names(e))

install.packagesCond <- function(pkgs) {
  inst <- rownames(installed.packages())
  for (pkg in c(pkgs)) {
    if (!(pkg %in% inst)) {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
  }
}

library.silent <- function(libName) {
  suppressWarnings(suppressMessages(library(package=libName, character.only = TRUE)))
}



if (interactive()) {
  if (length(grep("anacon", file.path(R.home("bin"), "R"), ignore.case = TRUE)) > 0) {
    if (rstudioapi::showQuestion("", "Switch library path exclusively to Anaconda Environment?", "Yes", "No")) {
      #
      for (p in .libPaths()) {
        if (length(grep("anacon", p, ignore.case = TRUE)) > 0) {
          rstudioapi::showDialog("", paste("Switching to", p))
          .libPaths(p)
          break
        }
      }
    }
  }
}


prepareDataForDepmix <- function(
  data, orderedFactor2Numeric = TRUE, sparseInteger2Nominal = TRUE
) {
  if (!is.data.frame(data)) {
    stop("data is not a data.frame")
  }
  
  l <- nrow(data)
  cn <- colnames(data)
  for (c in cn) {
    d <- data[[c]]
    if (is.numeric(d)) {
      # Check if is integers:
      ni <- sum(d %% 1)
      if (ni == 0) { # all are integers
        if (sparseInteger2Nominal && length(unique(d)) <= base::ceiling(log2(l))) {
          d <- as.character(d)
          data[[c]] <- d
        }
      }
    } else if (is.logical(d)) {
      data[[c]] <- as.character(d)
    } else if (is.character(d)) {
      next
    } else if (is.factor(d)) {
      if (orderedFactor2Numeric && is.ordered(d)) {
        data[[c]] <- as.numeric(d)
      } else {
        data[[c]] <- as.character(d)
      }
    } else {
      data[[c]] <- as.character(d)
    }
  }
  
  return(data)
}


getExperimentConn <- function() {
  cnfFile <- normalizePath("../my.cnf")
  if (!file.exists(cnfFile)) {
    stop(paste("Create a my.cnf in", cnfFile, "that specifies a connection a MySQL/MariaDB server to connect to."))
  }
  
  return(dbConnect(
    RMariaDB::MariaDB(),
    default.file=cnfFile,
    group="experiments")
  )
}

get1stOrderCommits <- function() {
  conn <- getExperimentConn()
  result <- dbSendQuery(conn, "SELECT t0.*, t1.* FROM `rq5_treex` as t0 inner join `rq5_treex` as t1 on t0.ParentCommitSHA1s = t1.commitId;")
  ds <- dbFetch(result)
  
  # We get duplicate names that end in ..\d+,
  # and we replace those names by "_t_1" (t-1)
  colnames(ds) <- gsub("\\.\\.\\d+$", "_t_1", colnames(ds))
  
  # Also, rename the other columns to have
  # the suffix "_t_0":
  temp <- colnames(ds)
  m <- !grepl("_t_1", temp)
  for (i in 1:length(temp)) {
    if (m[i]) {
      temp[i] <- paste(temp[i], "_t_0", sep = "")
    }
  }
  colnames(ds) <- temp
  
  dbClearResult(result)
  dbDisconnect(conn)
  return(ds)
}

get2ndOrderCommits <- function() {
  conn <- getExperimentConn()
  result <- dbSendQuery(conn, "SELECT t0.*, t1.* ,t2.* FROM `rq5_treex` as t0 inner join `rq5_treex` as t1 on t0.ParentCommitSHA1s = t1.commitId inner join `rq5_treex` as t2 on t1.ParentCommitSHA1s = t2.commitId;")
  ds <- dbFetch(result)
  
  # Now we get order+1 times the columns, we split the colnames
  # by order+1 and then rename using the schema from 1st-order.
  cn <- colnames(ds)
  l <- length(cn) / 3
  names1st <- gsub("\\.\\.\\d+$", "_t_1", cn[(l+1):(l*2)])
  names2nd <- gsub("\\.\\.\\d+$", "_t_2", cn[(l*2+1):(l*3)])
  colnames(ds) <- c(cn[1:l], names1st, names2nd)
  
  # Also, rename the other columns to have
  # the suffix "_t_0":
  temp <- colnames(ds)
  m <- !grepl("_t_\\d", temp)
  for (i in 1:length(temp)) {
    if (m[i]) {
      temp[i] <- paste(temp[i], "_t_0", sep = "")
    }
  }
  colnames(ds) <- temp
  
  dbClearResult(result)
  dbDisconnect(conn)
  return(ds)
}

get3rdOrderCommits <- function() {
  conn <- getExperimentConn()
  result <- dbSendQuery(conn, "SELECT t0.*, t1.*, t2.*, t3.* FROM `rq5_treex` as t0 inner join `rq5_treex` as t1 on t0.ParentCommitSHA1s = t1.commitId inner join `rq5_treex` as t2 on t1.ParentCommitSHA1s = t2.commitId inner join `rq5_treex` as t3 on t2.ParentCommitSHA1s = t3.commitId;")
  ds <- dbFetch(result)
  
  # Now we get order+1 times the columns, we split the colnames
  # by order+1 and then rename using the schema from 1st-order.
  cn <- colnames(ds)
  l <- length(cn) / 4
  names1st <- gsub("\\.\\.\\d+$", "_t_1", cn[(l+1):(l*2)])
  names2nd <- gsub("\\.\\.\\d+$", "_t_2", cn[(l*2+1):(l*3)])
  names3rd <- gsub("\\.\\.\\d+$", "_t_3", cn[(l*3+1):(l*4)])
  colnames(ds) <- c(cn[1:l], names1st, names2nd, names3rd)
  
  # Also, rename the other columns to have
  # the suffix "_t_0":
  temp <- colnames(ds)
  m <- !grepl("_t_\\d", temp)
  for (i in 1:length(temp)) {
    if (m[i]) {
      temp[i] <- paste(temp[i], "_t_0", sep = "")
    }
  }
  colnames(ds) <- temp
  
  dbClearResult(result)
  dbDisconnect(conn)
  return(ds)
}

getDataset <- function(dsName, removeUnwantedColums = TRUE) {
  conn <- getExperimentConn()
  result <- dbSendQuery(conn, paste("SELECT * FROM ", dsName))
  ds <- dbFetch(result)
  
  if (removeUnwantedColums) {
    removeNames <- c("SHA1",
                     #"RepoPathOrUrl",
                     "AuthorName", "CommitterName", "AuthorTime",
                     "CommitterTime", "MinutesSincePreviousCommit", "Message",
                     "AuthorEmail", "CommitterEmail",
                     "AuthorNominalLabel", "CommitterNominalLabel",
                     "ParentCommitSHA1s")
    
    ds <- ds[, !names(ds) %in% removeNames]
  }
  
  dbClearResult(result)
  dbDisconnect(conn)
  return(ds)
}


#' "data" needs to be data.frame with commits (i.e., not a list)
#' @param data data.frame with observations. Each observation must have
#' its own ID and the ID of a parent. For each set of observations in
#' the data that can be concatenated by these two IDs, a data.frame of
#' consecutive commits is produced then.
#' @param returnOnlyInitial default FALSE, if TRUE, returns a data.frame
#' with only those observations that do not have any parent (i.e.,
#' observations that are initial in an ordered set of consecutive ob-
#' servations).
#' @return list where the key is an observation ID, and the value is a
#' data.frame with consecutive observations. The first entry/row in this
#' data.frame is the observation with that ID, and all following obser-
#' vations are consecutive children. In other words, the first observation
#' in that data.frame is the one that has no further parents, and the last
#' observation is the "youngest", meaning there are no further children.
buildObservationChains <- function(data, idCol, parentIdCol, returnOnlyInitial = F) {
  observIds <- unique(data[[idCol]])
  parentIds <- unique(data[[parentIdCol]])
  chains <- list()
  
  initialIds <- observIds[!(data[[parentIdCol]] %in% observIds)]
  initialObs <- data[data[[idCol]] %in% initialIds, ]
  
  if (returnOnlyInitial) {
    return(initialObs)
  }
  
  dataByParentId <- list()
  for (i in 1:nrow(data)) {
    obs <- data[i, ]
    dataByParentId[[obs[[parentIdCol]]]] <- obs
  }
  
  for (i in 1:nrow(initialObs)) {
    parent <- initialObs[i, ]
    chain <- data.frame(parent)
    
    parentId <- parent[[idCol]]
    while (TRUE) {
      child <- data[data[[parentIdCol]] == parentId, ]
      
      if (nrow(child) != 1) {
        break
      }
      
      chain <- rbind(chain, child)
      parentId <- child[[idCol]]
    }
    
    chains[[parent[[idCol]]]] <- chain
  }
  
  return(chains)
}

grepAnyOrAll <- function(patterns, subjects, allMustMatch = FALSE) {
  return(sapply(subjects, function(subject) {
    res <- sapply(patterns, function(pattern) {
      return(grepl(pattern = pattern, x = subject))
    })
    
    if (allMustMatch) {
      return(all(res))
    }
    return(sum(res) > 0) # return if at least one matches
  }))
}


#' Returns a list of seeds used in parallel training with caret. For
#' repeatability, we need deterministic seeds. The amount depends on
#' the amounts of hyperparamenters, and number of folds/repeats.
#' @param nh integer, the number of hyperparameters
#' @param amount integer, the number of seeds, usually this is number
#' of folds * number of repeats.
#' @param seed integer used in \code{set.seed()}. Given an identical
#' seed, this function produces the same seeds (idempotent).
#' @return list with seeds that can be used in caret's trainControl
get_seeds <- function(nh, amount, seed = 42) {
  set.seed(seed)
  
  seeds <- vector(mode = "list", length = amount + 1)
  for(i in 1:amount) seeds[[i]] <- sample.int(.Machine$integer.max, nh)
  # For the last model:
  seeds[[amount + 1]] <- sample.int(.Machine$integer.max, 1)
  return(seeds)
}


doRfe <- function(data, y, yColName = "label_t", ignoreCols, number = 10, repeats = 10, maxSize = 50) {
  set.seed(1337)
  
  cn <- colnames(data)
  data <- data[, !grepAnyOrAll(c(ignoreCols, yColName), cn)]
  
  # We are training models now using RandomForests:
  control <- caret::rfeControl(
    functions = caret::rfFuncs, method = "cv", number = number, repeats = number)
  
  resultsRfe <- caret::rfe(
    x = data, y = y,
    # attempt to try sets of attributes between 1 and all attributes of size
    sizes = 1:min(maxSize, length(colnames(data))), rfeControl = control)
  
  return(resultsRfe)
}


evaluateStatelessModel <- function(
  ds, modelName, yColName, trP = 0.8, rcvNumber = 5, rcvRepeats = 3
) {
  isMMB <- modelName == "bayesCaret"
  
  y <- ds[, yColName]
  x <- ds[, !(colnames(ds) %in% yColName)]
  
  trMethod <- "repeatedcv"
  tuneLength <- ifelse(trMethod == "none", 1, 3)
  numHypers <- nrow(if (isMMB) mmb::bayesCaret$grid(x = x, y = y) else caret::getModelInfo(modelName)[[1]]$grid(len = tuneLength, x = x, y = y))
  
  set.seed(23)
  trControl <- caret::trainControl(
    method = trMethod, number = rcvNumber, repeats = rcvRepeats, p = trP,
    seeds = get_seeds(numHypers, rcvNumber * rcvRepeats, seed = 23))
  
  set.seed(371)
  modelObj <- if (isMMB) mmb::bayesCaret else modelName
  model <- caret::train(x = x, y = y, method = modelObj, trControl = trControl)
  
  return(model)
}


doWithParallelCluster <- function(expr, errorValue = NULL, numCores = parallel::detectCores()) {
  cl <- parallel::makePSOCKcluster(numCores)
  doSNOW::registerDoSNOW(cl)
  
  result <- tryCatch(expr, error=function(cond) {
    if (!is.null(errorValue)) {
      return(errorValue)
    }
    return(cond)
  }, finally = {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    cl <- NULL
    gc()
  })
  return(result)
}


extractChain <- function(data, idColData, idColChain, t, chains, chainId) {
  chainIds <- chains[[chainId]][[idColChain]]
  if (t > 0) {
    chainIds <- utils::tail(chainIds, -t)
  }
  
  obs <- data.frame()
  for (cId in chainIds) {
    dRow <- data[data[[idColData]] == cId, ]
    if (nrow(dRow) != 1) {
      stop(paste("Cannot extract row with ID", cId, nrow(dRow), nrow(data)))
    }
    obs <- rbind(obs, dRow)
  }
  
  return(list(
    data = data[!(data[[idColData]] %in% chainIds), ],
    obs = obs
  ))
}


detectCommitKeywords <- function(messages, use.binary = TRUE) {
  # (1) add (2) allow (3) bug (4) chang (5) error (6) fail (7) fix (8) implement
  # (9) improv (10) issu (11) method (12) new (13) npe
  # (14) refactor (15) remov (16) report (17) set (18) support
  # (19) test (20) use
  words <- c("add", "allow", "bug", "chang", "error", "fail", "fix", "implement",
             "improv", "issu", "method", "new", "npe",
             "refactor", "remov", "report", "set", "support",
             "test", "use")
  
  trueVal <- if (use.binary) 1 else TRUE
  falseVal <- if (use.binary) 0 else FALSE
  
  kw <- matrix(nrow = length(messages), ncol = length(words))
  for (i in 1:length(messages)) {
    mesVec <- strsplit(as.character(messages[i]), "\\s+")[[1]]
    kw[i, ] <- sapply(words, function(word) {
      if (any(grepl(x = mesVec, pattern = word, ignore.case = TRUE))) trueVal else falseVal
    })
  }
  
  kw <- data.frame(kw)
  colnames(kw) <- paste0("kw_", words)
  return(kw)
}


confmat2Dataframe <- function(cm, prefix = NULL) {
  clazzes <- colnames(cm$table)
  nClass <- length(clazzes)
  # We store from overall:
  # (4) Accuracy, Kappa, AccuracyNull (ZeroR), McnemarPValue,
  # then, for each class we store:
  # (N*5) Precision, Recall, F1, Prevalence, Balanced Accuracy
  df <- data.frame(matrix(nrow = 1, ncol = 4 + nClass * 5))
  
  colsO <- c("Accuracy", "Kappa", "AccuracyNull", "McnemarPValue")
  colsC <- c("Precision", "Recall", "F1", "Prevalence", "BalancedAccuracy")
  colsCExpanded <- c(sapply(clazzes, function(clz) {
    paste0(colsC, paste0("_", clz))
  }))
  
  colnames(df) <- c(colsO, colsCExpanded)

  o <- cm$overall
  for (c in colsO) {
    df[1, c] <- o[[c]]
  }
  
  for (clz in clazzes) {
    bc <- as.data.frame(cm$byClass)[paste("Class:", clz), ]
    colnames(bc) <- gsub("\\s+", "", colnames(bc))
    
    for (c in colsC) {
      df[1, paste0(c, "_", clz)] <- bc[[c]]
    }
  }
  
  # Let's conditionally prepend the prefix:
  if (!is.null(prefix)) {
    colnames(df) <- paste0(paste0(prefix, "_"), colnames(df))
  }
  
  return(df)
}


firstOrderPredToFlatConfmats <- function(groundTruth, pred, states) {
  # Creates 'aa', 'ac', .., 'pp'
  statesCombined <- apply(expand.grid(
    a = states,
    b = states
  ), 1, function(r) paste0(r[1], r[2]))
  
  gt_L <- sapply(groundTruth, function(gt) stringi::stri_sub(gt, from = 1, to = 1))
  gt_R <- sapply(groundTruth, function(gt) stringi::stri_sub(gt, from = 2, to = 2))
  pr_L <- sapply(pred, function(pr) stringi::stri_sub(pr, from = 1, to = 1))
  pr_R <- sapply(pred, function(pr) stringi::stri_sub(pr, from = 2, to = 2))
  
  cm <- caret::confusionMatrix(
    data = factor(pred, levels = statesCombined),
    reference = factor(groundTruth, levels = statesCombined),
    mode = "everything")
  res_both <- confmat2Dataframe(cm = cm, prefix = "both")
  
  cm <- caret::confusionMatrix(
    data = factor(pr_L, levels = states),
    reference = factor(gt_L, levels = states),
    mode = "everything")
  res_left <- confmat2Dataframe(cm = cm, prefix = "left")
  
  cm <- caret::confusionMatrix(
    data = factor(pr_R, levels = states),
    reference = factor(gt_R, levels = states),
    mode = "everything")
  res_right <- confmat2Dataframe(cm = cm, prefix = "right")
  
  return(list(
    both = res_both,
    left = res_left,
    right = res_right
  ))
}


generateDataForHpOrderedLoss <- function(hp, data, dataKw, states, stateColumn) {
  
  # First, we create a joined label for the data.
  stateColumn_t_1 <- gsub("_t_0$", "_t_1", stateColumn)
  
  # Let's check if we need to append the keywords to the data first:
  if (hp$useKeywords) {
    data <- cbind(data, dataKw)
  }
  
  ignoreCols <- c(
    "project", # because of scarce data, we train and evaluate only cross-project
    "commitId",
    "labelledBefore",
    "reasonOrRule",
    "RepoPathOrUrl",
    "CommitterName",
    "CommitterTime",
    "CommitterEmail",
    "CommitterNominalLabel",
    "ParentCommitSHA1s",
    "Message",
    "branch",
    "SHA1",
    "AuthorName",
    "AuthorTime",
    "AuthorEmail",
    "AuthorNominalLabel"
  )
  
  # Remove features we ALWAYS ignore:
  data <- data[, !grepAnyOrAll(ignoreCols, colnames(data))]
  
  # If this is true, then no densities are used. We can still
  # use data2densities, but let's tell it to skip the densities.
  useDens <- hp$useDatatype %in% c("density", "both", "density_scaled", "both_scaled")
  scaleDens <- hp$useDatatype %in% c("density_scaled", "both_scaled")
  useData <- hp$useKeywords || hp$useDatatype %in% c("data", "both", "both_scaled")
  useOnlyData <- hp$useDatatype == "data"
  
  set.seed(hp$resampleSeed)
  trainValid <- data2Densities_1stOrder(
    states = states, data = data, stateColumn = stateColumn,
    split = 0.8, returnDataOnly = useOnlyData,
    doEcdf = hp$dens.fct == "ecdf",
    normalizePdfs = hp$dens.fct == "epdf_scaled",
    ecdfMinusOne = hp$dens.fct == "1-ecdf")
  
  train <- if (useOnlyData) trainValid$train_data else trainValid$train[, ]
  train_Y_obsId <- train$obsId
  train$obsId <- NULL
  
  yCols <- if (useOnlyData) c(stateColumn, stateColumn_t_1) else c("y_t0", "y_t1")
  yColsInt <- !useOnlyData # y_t1, y_t0 are integers
  
  # Sometimes we need to oversample the training data to account
  # for the non-even distribution of class labels.
  if (hp$oversample) {
    # What we do here is to create a combined label and then balance on that.
    # After balancing, we split the combined label again and reassign the
    # correct classes for each state.
    train$y__ <- paste0(train[[yCols[1]]], "__|__", train[[yCols[2]]])
    
    train <- balanceDatasetSmote(train, "y__")
    train[[yCols[1]]] <- sapply(train$y__, function(l) {
      temp <- strsplit(l, "__|__")[[1]][1]
      return(if (yColsInt) as.integer(temp) else temp)
    })
    train[[yCols[2]]] <- sapply(train$y__, function(l) {
      temp <- strsplit(l, "__|__")[[1]][3]
      return(if (yColsInt) as.integer(temp) else temp)
    })
    train$y__ <- NULL
    
    if (!useOnlyData) {
      train$est_t1 <- as.integer(round(train$est_t1))
      train$est_t0 <- as.integer(round(train$est_t0))
    }
  }
  
  train_Y <- as.matrix(train[, sort(colnames(train)[grepAnyOrAll(paste0("^", yCols, "$"), colnames(train))])])
  train_est <- if (useOnlyData) NULL else train[, colnames(train) %in% c("est_t0", "est_t1")]
  
  # Let's remove columns we ALWAYS ignore:
  train <- train[, !grepAnyOrAll(ignoreCols, colnames(train))]
  
  # Let's remove columns with zero variance:
  nzv <- caret::nzv(train, names = TRUE, saveMetrics = TRUE)
  train <- train[, !nzv$zeroVar]
  # Then remove labels or what kind of estimator was used:
  train <- train[, !grepAnyOrAll(c("label", "est_t", "y_t"), colnames(train))]
  
  # Now we need to make some final selection, and choose between
  # data, density(data), or both:
  if (useOnlyData) {
    # Remove density data:
    train <- train[, !grepAnyOrAll("^d_", colnames(train))]
  } else if (hp$useDatatype %in% c("density", "density_scaled")) {
    # Remove the ordinary data (but NOT the keywords):
    train <- train[, grepAnyOrAll(c("^d_", "^kw_"), colnames(train))]
  } # else keep both, no sub-selection required.
  
  cols_data <- !grepAnyOrAll("^d_", colnames(train))
  cols_dens <- !cols_data
  
  scaler_train_data_X <- if (useData) {
    caret::preProcess(train[, cols_data], method = c("center", "scale")) } else { NULL }
  scaler_train_dens_X <- if (useDens && scaleDens) {
    caret::preProcess(train[, cols_dens], method = c("center", "scale")) } else { NULL }
  
  # The validation data
  valid <- if (useOnlyData) trainValid$valid_data else trainValid$valid[, ]
  valid_est <- if (useOnlyData) NULL else valid[, colnames(valid) %in% c("est_t0", "est_t1")]
  valid_Y_obsId <- valid$obsId
  valid_Y <- as.matrix(valid[, colnames(train_Y)])
  valid <- valid[, colnames(valid) %in% colnames(train)]
  
  if (useData) {
    if (useDens) {
      if (scaleDens) {
        train_X <- cbind(
          as.matrix(stats::predict(scaler_train_data_X, train[, cols_data])),
          as.matrix(stats::predict(scaler_train_dens_X, train[, cols_dens]))
        )
        valid_X <- cbind(
          as.matrix(stats::predict(scaler_train_data_X, valid[, cols_data])),
          as.matrix(stats::predict(scaler_train_dens_X, valid[, cols_dens]))
        )
      } else {
        train_X <- cbind(
          as.matrix(stats::predict(scaler_train_data_X, train[, cols_data])),
          as.matrix(train[, cols_dens])
        )
        valid_X <- cbind(
          as.matrix(stats::predict(scaler_train_data_X, valid[, cols_data])),
          as.matrix(valid[, cols_dens])
        )
        
      }
    } else {
      # Only data, no density
      train_X <- as.matrix(stats::predict(scaler_train_data_X, train[, cols_data]))
      valid_X <- as.matrix(stats::predict(scaler_train_data_X, valid[, cols_data]))
    }
  } else {
    # Only density
    if (scaleDens) {
      train_X <- as.matrix(stats::predict(scaler_train_dens_X, train[, cols_dens]))
      valid_X <- as.matrix(stats::predict(scaler_train_dens_X, valid[, cols_dens]))
    } else {
      train_X <- train[, cols_dens]
      valid_X <- valid[, cols_dens]
    }
  }
  
  if (!useOnlyData) {
    # These belong to the input data and designate from which density
    # estimator the data is coming.
    train_X <- cbind(train_X, train_est)
    valid_X <- cbind(valid_X, valid_est)
  }
  
  return(list(
    train_X = train_X,
    valid_X = valid_X,
    train_Y = train_Y,
    valid_Y = valid_Y,
    train_Y_obsId = train_Y_obsId,
    valid_Y_obsId = valid_Y_obsId,
    scaler_train_data_X = scaler_train_data_X,
    scaler_train_dens_X = scaler_train_dens_X
  ))
}


#' In this function, we perform a voting for each scheme, depending
#' on the data (not all schemes can always be performed). The result
#' is a list with all applicable results. If the voting for one scheme
#' failed, the result list contains an error message for that scheme.
#' If everything worked, the result contains a flattened confusion
#' matrix. If not applicable, the entry will be NULL.
evaluateHpOrderedLoss <- function(hp, data, dataKw, states, stateColumn) {
  if (!(hp$useDatatype %in% c("density", "density_scaled", "both", "both_scaled"))) {
    stop(paste0(hp$useDatatype, " not supported - need density."))
  }
  
  data <- generateDataForHpOrderedLoss(
    hp = hp, data = data, dataKw = dataKw,
    states = states, stateColumn = stateColumn)
  
  isOversampled <- hp$oversample
  unqObsIds <- sort(unique(data$valid_Y_obsId))
  
  # In this case, we cannot use schemes one and two, as they require
  # the density data to be strictly positive.
  isDensScaled <- hp$useDatatype %in% c("density_scaled", "both_scaled")
  
  # For the following predictive scenarios and voting schemes it is
  # important to understand what we are interested in and what data
  # we have. If we use the density of the observations in any way,
  # then we always try to predict the correct sequence of density
  # estimators, which results in attempting to predict an n-times one-
  # hot encoding of correctly chosen estimators (where n is the order
  # of the model). Each sample exists n^2 times and the overall goal
  # in this scenario is to find the best-matching succession of density
  # estimators for each sample.
  
  # Columns end in t_1, t_0 where t_0 is the most important and must
  # come first. Using sort() ensures this.
  yCols <- sort(colnames(data$valid_Y))
  validJoined <- data.frame(cbind(
    data$valid_Y, data$valid_X$est_t1, data$valid_X$est_t0, data$valid_Y_obsId))
  colnames(validJoined) <- c(yCols, "est_t1", "est_t0", "obsId")
  
  # This is the ground-truth for our 2-times one-hot encoding; out of all density-
  # activations it determines which was the correct one. Every observation has an
  # ID and there are n^2 activations per ID. Select the estimators where both both,
  # y_t1 and y_t0 are correct.
  groundTruth <- sapply(unqObsIds, function(oId) {
    temp <- validJoined[validJoined$obsId == oId & validJoined$y_t1 == 1 & validJoined$y_t0 == 1, ]
    if (nrow(temp) > 1) stop("This must not happen.")
    paste0(states[temp$est_t1], states[temp$est_t0])
  })
  
  # Left is t-1, right is t0
  groundTruth_L <- sapply(groundTruth, function(gt) stringi::stri_sub(gt, from = 1, to = 1))
  groundTruth_R <- sapply(groundTruth, function(gt) stringi::stri_sub(gt, from = 2, to = 2))
  
  # Some schemes produce a result for each of these:
  logTols <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8)
  
  # Creates 'aa', 'ac', .., 'pp'
  statesCombined <- apply(expand.grid(
    a = states,
    b = states
  ), 1, function(r) paste0(r[1], r[2]))
  
  results <- list(
    One = NULL,
    Two = NULL,
    Three = NULL,
    Four = NULL,
    # Include some metadata before we return:
    meta = list()
  )

  
  
  # Let's do the voting, scheme by scheme:
  
  #########################################################################
  ############  Scheme One
  #########################################################################
  
  if (!isDensScaled) {
    
    # We can do the schemes 'One' and 'Two', they require density. For scheme
    # One, we produce a matrix where each column corresponds to one specific
    # estimator. The estimator the data was generated with is stored in the
    # data in the columns 'est_t1', 'est_t0' as integers and refers to the
    # index in the 'states'-vector. A fully-qualified estimator thus is the
    # concat of both states.
    
    # Techniques a/b/c are the first three columns, and the next 2*n columns
    # are technique d with a different logTol every time (prod/sum). The rows
    # are the predictions. (a/b/c) are: max(prod), max(sum), min(sd).
    pred_1 <- matrix(ncol = 3 + 2 * length(logTols), nrow = length(unqObsIds))
    colnames(pred_1) <- c(paste0("simple_", c("a", "b", "c")),
                          paste0("d_prod_", formatC(logTols, format = "e", digits = 0)),
                          paste0("d_sum_", formatC(logTols, format = "e", digits = 0)))
    rownames(pred_1) <- unqObsIds
    
    for (obsId in unqObsIds) {
      
      # Cols is estimators, rows is sum, prod, sd, log10(prod), log10(sum)
      votemat <- matrix(nrow = 5, ncol = length(states)**2)
      colnames(votemat) <- statesCombined
      validSample <- data$valid_X[data$valid_Y_obsId == obsId, ]
      
      for (state_i in states) {
        for (state_j in states) {
          estimatorName <- paste0(state_i, state_j)
          estimatorData <- validSample[validSample$est_t1 == which(states == state_i) & validSample$est_t0 == which(states == state_j), grepAnyOrAll("^d_", colnames(validSample))]
          estimatorData <- estimatorData
          
          votemat[1, estimatorName] <- prod(estimatorData + 1/3)
          votemat[2, estimatorName] <- sum(estimatorData)
          votemat[3, estimatorName] <- sd(estimatorData)
          # log10's of prod and sum
          votemat[4, estimatorName] <- log10(votemat[1, estimatorName])
          votemat[5, estimatorName] <- log10(votemat[2, estimatorName])
        }
      }
      
      # The vote-matrix is complete now, let's pick the vote for the current observation.
      # a) (max(prod))
      pred_1[obsId, "simple_a"] <- names(which.max(votemat[1, ]))
      # b) (max(sum))
      pred_1[obsId, "simple_b"] <- names(which.max(votemat[2, ]))
      # c) (min(sd))
      pred_1[obsId, "simple_c"] <- names(which.min(votemat[3, ]))
      
      # d-1) (prod) Let's pick the votes using each log-tolerance:
      # d-2) (sum)
      for (lt in logTols) {
        ltStr <- formatC(lt, format = "e", digits = 0)
        
        maxLogProd <- votemat[4, which.max(votemat[4, ])]
        maxLogSum <- votemat[5, which.max(votemat[5, ])]
        
        equalLogsProd <- which(votemat[4, ] >= (maxLogProd - lt))
        equalLogsSum <- which(votemat[5, ] >= (maxLogSum - lt))
        
        minSdProd <- which.min(votemat[3, equalLogsProd])
        minSdSum <- which.min(votemat[3, equalLogsSum])
        
        pred_1[obsId, paste0("d_prod_", ltStr)] <- names(minSdProd)
        pred_1[obsId, paste0("d_sum_", ltStr)] <- names(minSdSum)
      }
    }
    
    resultsOne <- list()
    for (resType in colnames(pred_1)) {
      resultsOne[[resType]] <- firstOrderPredToFlatConfmats(
        groundTruth = groundTruth, pred = pred_1[, resType], states = states)
    }
    
    results$One <- resultsOne
  }

  
      
  #########################################################################
  ############  Scheme Two
  #########################################################################
  
  if (!isDensScaled) {
    # For scheme Two, we combine the estimators by grouping them. Each group
    # represents all estimators that end in a state t_0 (group estimators
    # with same state for t_0). Also, make a votematrix grouping the estimators
    # that start in a specific state t-1, so that we can test whether this also
    # works backward.
    # Again, we have a couple of simple techniques a/b/c and a few that depend
    # on a log-tolerance and standard deviation for votes within a group.
    
    
    # This is exactly like pred_1, except that we do everything twice: first
    # forward, then backward.
    pred_2 <- matrix(ncol = 2 * (3 + 2 * length(logTols)), nrow = length(unqObsIds))
    cn_pred_2 <- c(paste0("simple_", c("a", "b", "c")),
                   paste0("d_prod_", formatC(logTols, format = "e", digits = 0)),
                   paste0("d_sum_", formatC(logTols, format = "e", digits = 0)))
    colnames(pred_2) <- c(
      paste0("forward_", cn_pred_2),
      paste0("backward_", cn_pred_2))
  
    rownames(pred_2) <- unqObsIds
    
    for (obsId in unqObsIds) {
      
      validSample <- data$valid_X[data$valid_Y_obsId == obsId, ]
      
      # Cols is estimators, rows is sum, prod, sd, log10(sum), log10(prod)
      # 'to_t0' is forward, and 'to_t1' is backward.
      votemat_to_t0 <- matrix(nrow = 5, ncol = length(states))
      colnames(votemat_to_t0) <- states
      votemat_to_t1 <- matrix(nrow = 5, ncol = length(states))
      colnames(votemat_to_t1) <- states
      
      for (state in states) {
        estimatorData_to_t0 <- validSample[validSample$est_t0 == which(states == state), grepAnyOrAll("^d_", colnames(validSample))]
        estimatorData_to_t1 <- validSample[validSample$est_t1 == which(states == state), grepAnyOrAll("^d_", colnames(validSample))]
        
        votemat_to_t0[1, state] <- prod(apply(estimatorData_to_t0, 1, sum))
        votemat_to_t0[2, state] <- sum(estimatorData_to_t0)
        votemat_to_t0[3, state] <- sd(apply(estimatorData_to_t0, 1, sum))
        votemat_to_t0[4, state] <- log10(votemat_to_t0[1, state])
        votemat_to_t0[5, state] <- log10(votemat_to_t0[2, state])
        
        votemat_to_t1[1, state] <- prod(apply(estimatorData_to_t1, 1, sum))
        votemat_to_t1[2, state] <- sum(estimatorData_to_t1)
        votemat_to_t1[3, state] <- sd(apply(estimatorData_to_t1, 1, sum))
        votemat_to_t1[4, state] <- log10(votemat_to_t1[1, state])
        votemat_to_t1[5, state] <- log10(votemat_to_t1[2, state])
      }
      
      
      # The vote-matrix is complete now, let's pick the vote for the current observation.
      # a) (max(prod))
      pred_2[obsId, "forward_simple_a"] <- names(which.max(votemat_to_t0[1, ]))
      pred_2[obsId, "backward_simple_a"] <- names(which.max(votemat_to_t1[1, ]))
      # b) (max(sum))
      pred_2[obsId, "forward_simple_b"] <- names(which.max(votemat_to_t0[2, ]))
      pred_2[obsId, "backward_simple_b"] <- names(which.max(votemat_to_t1[2, ]))
      # c) (min(sd))
      pred_2[obsId, "forward_simple_c"] <- names(which.min(votemat_to_t0[3, ]))
      pred_2[obsId, "backward_simple_c"] <- names(which.min(votemat_to_t1[3, ]))
      
      # d-1) (prod) Let's pick the votes using each log-tolerance:
      # d-2) (sum)
      for (lt in logTols) {
        ltStr <- formatC(lt, format = "e", digits = 0)
        
        maxLogProd_fw <- votemat_to_t0[4, which.max(votemat_to_t0[4, ])]
        maxLogProd_bw <- votemat_to_t1[4, which.max(votemat_to_t1[4, ])]
        maxLogSum_fw <- votemat_to_t0[5, which.max(votemat_to_t0[5, ])]
        maxLogSum_bw <- votemat_to_t1[5, which.max(votemat_to_t1[5, ])]
        
        equalLogsProd_fw <- which(votemat_to_t0[4, ] >= (maxLogProd_fw - lt))
        equalLogsProd_bw <- which(votemat_to_t1[4, ] >= (maxLogProd_bw - lt))
        equalLogsSum_fw <- which(votemat_to_t0[5, ] >= (maxLogSum_fw - lt))
        equalLogsSum_bw <- which(votemat_to_t1[5, ] >= (maxLogSum_bw - lt))
        
        minSdProd_fw <- which.min(votemat_to_t0[3, equalLogsProd_fw])
        minSdProd_bw <- which.min(votemat_to_t1[3, equalLogsProd_bw])
        minSdSum_fw <- which.min(votemat_to_t0[3, equalLogsSum_fw])
        minSdSum_bw <- which.min(votemat_to_t1[3, equalLogsSum_bw])
        
        pred_2[obsId, paste0("forward_d_prod_", ltStr)] <- names(minSdProd_fw)
        pred_2[obsId, paste0("backward_d_prod_", ltStr)] <- names(minSdProd_bw)
        pred_2[obsId, paste0("forward_d_sum_", ltStr)] <- names(minSdSum_fw)
        pred_2[obsId, paste0("backward_d_sum_", ltStr)] <- names(minSdSum_bw)
      }
    }
  
    
    resultsTwo <- list()
    for (resType in colnames(pred_2)) {
      isForward <- length(grep("^forward_", resType)) > 0
      useTruth <- if (isForward) groundTruth_R else groundTruth_L
      resultsTwo[[resType]] <- confmat2Dataframe(
        cm = caret::confusionMatrix(
          data = factor(pred_2[, resType], levels = states),
          reference = factor(useTruth, levels = states),
          mode = "everything"))
    }
    
    results$Two <- resultsTwo
  }

    
  
  #########################################################################
  ############  Scheme Three
  #########################################################################
  
  # We have covered the first two schemes, now it's time for the main scheme:
  # Training a customized neural network. This network supports any kind of
  # data, and will even train with scaled density.
  set.seed(hp$resampleSeed)
  
  lHl <- hp$neurons
  weights_hidden <- glorot_weights(ncol(data$train_X),lHl)
  biases_hidden <- glorot_weights(1, lHl)
  weights_output <- glorot_weights(lHl, 2)
  biases_output <- glorot_weights(1, ncol(data$train_Y))
  
  
  act.fct.derive <- relu_d1
  if (hp$act.fct == "Lrelu") {
    act.fct <- lrelu
    act.fct.derive <- lrelu_d1
  } else if (hp$act.fct == "swish") {
    act.fct <- swish
    act.fct.derive <- swish_d1
  } else if (hp$act.fct == "gelu") {
    act.fct <- gelu
    act.fct.derive <- gelu_d1
  } else if (hp$act.fct == "relu") {
    act.fct <- relu
    act.fct.derive <- relu_d1
  } else {
    stop(hp$act.fct)
  }
  
  
  err.fct.derive <- RSS_d1
  if (hp$err.fct == "wRSS") {
    err.fct <- wRSS
    err.fct.derive <- wRSS_d1
  } else if (hp$err.fct == "omeRSS") {
    err.fct <- omeRSS
    err.fct.derive <- omeRSS_d1
  } else if (hp$err.fct == "ome2RSS") {
    err.fct <- ome2RSS
    err.fct.derive <- ome2RSS_d1
  } else if (hp$err.fct == "RSS") {
    err.fct <- RSS
    err.fct.derive <- RSS_d1
  } else {
    stop(hp$err.fct)
  }
  
  start_three <- as.numeric(Sys.time())
  m1Fit <- gradient_descent_m1(
    X = as.matrix(data$train_X), Y = as.matrix(data$train_Y),
    w_h_0 = weights_hidden, b_h_0 = biases_hidden,
    w_o_0 = weights_output, b_o_0 = biases_output,
    epochs = hp$epochs, learning_rate = hp$learningRate,
    precision = 1e-3, batch_size = hp$batchSize,
    act.fn = act.fct, act.fn.derive = act.fct.derive,
    err.fn = err.fct, err.fn.derive = err.fct.derive
  )
  duration_three <- as.numeric(Sys.time()) - start_three
  
  # Store some meta for the result obtained:
  tempLm <- lm(
    formula = y~x,
    data = data.frame(x = 1:length(m1Fit$hist_loss), y = m1Fit$hist_loss))
  results$meta$hist_length <- length(m1Fit$hist_loss) - 1
  results$meta$epochs_used_perc <- (length(m1Fit$hist_loss) - 1) / hp$epochs
  results$meta$hist_intercept <- tempLm$coefficients[[1]]
  results$meta$hist_slope <- tempLm$coefficients[[2]]
  results$meta$duration_three <- duration_three
  
  # Again, we have 3 simple a/b/c schemes (pick max prod/sum; min sd) and
  # then we have the minimum sd according to the log-tolerance, both for
  # products and sums.
  pred_3 <- matrix(ncol = 3 + 2 * length(logTols), nrow = length(unqObsIds))
  colnames(pred_3) <- c(paste0("simple_", c("a", "b", "c")),
                        paste0("d_prod_", formatC(logTols, format = "e", digits = 0)),
                        paste0("d_sum_", formatC(logTols, format = "e", digits = 0)))
  rownames(pred_3) <- unqObsIds
  
  for (obsId in unqObsIds) {
    
    # Rows is estimators, cols are pred_t0, pred_t1, prod, sum, sd, log10(prod), log10(sum)
    votemat <- matrix(nrow = length(states)**2, ncol = 7)
    rownames(votemat) <- statesCombined
    
    validSample <- data$valid_X[data$valid_Y_obsId == obsId, ]
    
    for (state_i in states) {
      for (state_j in states) {
        # Now for each estimator, we inference the sample and produce estimates
        # for state t-1 and t0.
        est_t1 <- which(states == state_i)
        est_t0 <- which(states == state_j)
        
        # Returns a matrix with one row:
        predRaw <- m1(x_i = as.matrix(validSample[validSample$est_t1 == est_t1 & validSample$est_t0 == est_t0, ]),
                      w_h = m1Fit$w_h, b_h = m1Fit$b_h, w_o = m1Fit$w_o, b_o = m1Fit$b_o,
                      act.fn = act.fct, act.fn.derive = act.fct.derive)
        
        # Note the order here; j corresponds to t0, and i to t-1.
        # This is how it was passed to the network as Y.
        estimatorName <- paste0(state_i, state_j)
        votemat[estimatorName, 1:2] <- predRaw[1, ]
        votemat[estimatorName, 3:5] <- c(prod(predRaw[1, ]), sum(predRaw[1, ]), sd(predRaw[1, ]))
        votemat[estimatorName, 6:7] <- log10(votemat[estimatorName, 3:4])
      }
    }
    
    pred_3[obsId, "simple_a"] <- names(which.max(votemat[, 3]))
    pred_3[obsId, "simple_b"] <- names(which.max(votemat[, 4]))
    pred_3[obsId, "simple_c"] <- names(which.min(votemat[, 5]))
    
    for (lt in logTols) {
      ltStr <- formatC(lt, format = "e", digits = 0)
      
      maxLogProd <- votemat[which.max(votemat[, 6]), 6]
      maxLogSum <- votemat[which.max(votemat[, 7]), 7]
      
      equalLogsProd <- which(votemat[, 6] >= (maxLogProd - lt))
      equalLogsSum <- which(votemat[, 7] >= (maxLogSum - lt))
      
      minSdProd <- which.min(votemat[equalLogsProd, 5])
      minSdSum <- which.min(votemat[equalLogsSum, 5])
      
      pred_3[obsId, paste0("d_prod_", ltStr)] <- names(minSdProd)
      pred_3[obsId, paste0("d_sum_", ltStr)] <- names(minSdSum)
    }
  }
  
  resultsThree <- list()
  for (resType in colnames(pred_3)) {
    resultsThree[[resType]] <- firstOrderPredToFlatConfmats(
      groundTruth = groundTruth, pred = pred_3[, resType], states = states)
  }
  
  results$Three <- resultsThree

  
  
  #########################################################################
  ############  Scheme Four
  #########################################################################
  
  # Scheme 4 is conditional on the data NOT be oversampled. That is
  # because we flatten an observation and thus require a fixed length.
  # With oversampling, we have too many or few representations of the same
  # sample.
  # Scheme 4 reuses the model we fit in scheme 3, and roughly does this:
  # - Pass the training data through the network (not the validation data)
  #   to obtain the network's activations
  # - Flatten each of the activation's for both states and create a true
  #   one-hot encoding of which estimators were picked correctly.
  if (!hp$oversample) {
    unqTrainIds <- sort(unique(data$train_Y_obsId))

    train_X_four <- matrix(nrow = length(unqTrainIds), ncol = 2 + 2*length(states)**2)
    train_Y_four <- matrix(nrow = length(unqTrainIds), ncol = length(states)**2)
    rownames(train_X_four) <- unqTrainIds
    
    valid_X_four <- matrix(nrow = length(unqObsIds), ncol = 2 + 2*length(states)**2)
    valid_Y_four <- matrix(nrow = length(unqObsIds), ncol = length(states)**2)
    rownames(valid_X_four) <- unqObsIds
    
    # Attention: We take the training IDs:
    for (obsId in unqTrainIds) {
      trainSample_X <- as.matrix(data$train_X[data$train_Y_obsId == obsId, ])
      trainSample_Y <- as.matrix(data$train_Y[data$train_Y_obsId == obsId, ])
      idx11 <- which(trainSample_Y[, yCols[1]] == 1 & trainSample_Y[, yCols[2]] == 1)
      vec11 <- rep(0, length(states)**2)
      vec11[idx11] <- 1
      
      # Let's predict all of the densitiy-versions:
      predRaw <- matrix(nrow = nrow(trainSample_X), ncol = ncol(trainSample_Y))
      for (i in 1:nrow(trainSample_X)) {
        predRaw[i, ] <- m1(
          x_i = as.matrix(trainSample_X[i, ]),
          w_h = m1Fit$w_h, b_h = m1Fit$b_h, w_o = m1Fit$w_o, b_o = m1Fit$b_o,
          act.fn = act.fct, act.fn.derive = act.fct.derive)
      }
      
      train_X_four[obsId, ] <- c(predRaw[, 1], sd(predRaw[, 1]), predRaw[, 2], sd(predRaw[, 2]))
      train_Y_four[obsId, ] <- vec11
    }
    
    # Let's do the same for the validation data:
    for (obsId in unqObsIds) {
      validSample_X <- as.matrix(data$valid_X[data$valid_Y_obsId == obsId, ])
      validSample_Y <- as.matrix(data$valid_Y[data$valid_Y_obsId == obsId, ])
      idx11 <- which(validSample_Y[, "y_t1"] == 1 & validSample_Y[, "y_t0"] == 1)
      vec11 <- rep(0, length(states)**2)
      vec11[idx11] <- 1
      
      predRaw <- matrix(nrow = nrow(validSample_Y), ncol = ncol(validSample_Y))
      for (i in 1:nrow(validSample_X)) {
        predRaw[i, ] <- m1(
          x_i = as.matrix(validSample_X[i, ]),
          w_h = m1Fit$w_h, b_h = m1Fit$b_h, w_o = m1Fit$w_o, b_o = m1Fit$b_o,
          act.fn = act.fct, act.fn.derive = act.fct.derive)
      }
      
      valid_X_four[obsId, ] <- c(predRaw[, 1], sd(predRaw[, 1]), predRaw[, 2], sd(predRaw[, 2]))
      valid_Y_four[obsId, ] <- vec11
    }
    
    # Don't forget to scale the data:•
    train_X_four <- as.data.frame(train_X_four)
    valid_X_four <- as.data.frame(valid_X_four)
    
    scaler_X <- caret::preProcess(train_X_four, method = c("center", "scale"))
    train_X_four <- as.matrix(predict(scaler_X, newdata = train_X_four))
    valid_X_four <- as.matrix(predict(scaler_X, newdata = valid_X_four))
    
    last.warning <- NULL
    # Now we can train the 2nd-stage network and make predictions.
    set.seed(hp$resampleSeed)
    nnet <- tryCatch({
      neuralnet::neuralnet(
        formula = formula(paste0(paste(paste0("V", (ncol(train_X_four) + 1):(ncol(train_X_four) + ncol(train_Y_four))), collapse = " + "), "~.")),
        data = as.data.frame(cbind(train_X_four, train_Y_four)), # so we get col-names
        hidden = 7,
        stepmax = 15e4,
        startweights = glorot_weights(numIn = length(states)**2 + 2, numOut = 8),
        lifesign = "none")
    }, error=function(cond) {
      cond
    })
    
    if (("weights" %in% names(nnet)) && length(grep("did not converge", names(last.warning))) == 0) {
      # There is no strategy other than picking the maximum activation out
      # of the true one-hot vector.
      pred_4 <- matrix(nrow = length(unqObsIds), ncol = 1)
      rownames(pred_4) <- unqObsIds
      
      for (obsId in unqObsIds) {
        validSample_X <- as.matrix(data$valid_X[data$valid_Y_obsId == obsId, ])
        
        predOneHotIdx <- which.max(predict(nnet, newdata = t(as.data.frame(valid_X_four[obsId, ]))))
        
        estimatorPred <- paste0(
          states[validSample_X[predOneHotIdx, "est_t1"]], states[validSample_X[predOneHotIdx, "est_t0"]])
        
        pred_4[obsId, ] <- estimatorPred
      }
      
      results$Four <- firstOrderPredToFlatConfmats(
        groundTruth = groundTruth, pred = pred_4[, 1], states = states)
    } else {
      results$Four <- NULL
      # do not set an extra error string, we know when we're
      # supposed to have the fourth result.
      #results$Four <- "Error: did not converge"
    }
  }
  
  return(results)
}


writeOrAppend.csv <- (function() {
  padLockFile <- tempfile()

  return(function(csvFile, data) {
    existed <- file.exists(csvFile)
    if (!existed) {
      file.create(csvFile)
    }
    
    padLock <- filelock::lock(padLockFile, exclusive = TRUE, timeout = Inf)
    suppressWarnings(write.table(x = data, file = csvFile, append = TRUE, sep = ";", row.names = FALSE, col.names = !existed))
    filelock::unlock(padLock)
  })
})()



