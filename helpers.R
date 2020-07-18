install.packagesCond <- function(pkg) {
  if (pkg %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
}

library.silent <- function(libName) {
  suppressWarnings(suppressMessages(library(package=libName, character.only = TRUE)))
}

install.packagesCond("rstudioapi")
library.silent("rstudioapi")
install.packagesCond("RMariaDB")
library.silent("RMariaDB")


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


doWithParallelCluster <- function(expr, errorValue = NULL) {
  cl <- parallel::makePSOCKcluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  
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


confmat2Dataframe <- function(cm, prefix) {
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
    bc <- as.data.frame(temp$byClass)[paste("Class:", clz), ]
    colnames(bc) <- gsub("\\s+", "", colnames(bc))
    
    for (c in colsC) {
      df[1, paste0(c, "_", clz)] <- bc[[c]]
    }
  }
  
  # Let's prepend the prefix:
  colnames(df) <- paste0(paste0(prefix, "_"), colnames(df))
  
  return(df)
}


generateDataForHpOrderedLoss <- function(hp, data, dataKw, states, stateColumn) {
  
  # First, we create a joined label for the data.
  stateColumn_t_1 <- gsub("_t_0$", "_t_1", stateColumn)
  
  # Let's check if we need to append the keywords to the data first:
  if (hp$useKeywords) {
    data <- cbind(data, dataKw)
  }
  
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
  
  yCols <- if (useOnlyData) c(stateColumn_t_1, stateColumn) else c("y_t1", "y_t0")
  
  # Sometimes we need to oversample the training data to account
  # for the non-even distribution of class labels.
  if (hp$oversample) {
    # What we do here is to create a combined label and then balance on that.
    # After balancing, we split the combined label again and reassign the
    # correct classes for each state.
    train$y__ <- paste0(train[[yCols[1]]], "__|__", train[[yCols[2]]])
    
    train <- balanceDatasetSmote(train, "y__")
    train[[yCols[1]]] <- sapply(train$y__, function(l) as.integer(strsplit(l, "__|__")[[1]][1]))
    train[[yCols[2]]] <- sapply(train$y__, function(l) as.integer(strsplit(l, "__|__")[[1]][3]))
    train$y__ <- NULL
    
    if (!useOnlyData) {
      train$est_t1 <- as.integer(round(train$est_t1))
      train$est_t0 <- as.integer(round(train$est_t0))
    }
  }
  
  train_Y <- as.matrix(train[, grepAnyOrAll(paste0("^", yCols, "$"), colnames(train))])
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


evaluateHpOrderedLoss <- function(hp, data, dataKw, states, stateColumn) {
  data <- generateDataForHpOrderedLoss(hp = hp, data = data, dataKw = dataKw, states = states, stateColumn = stateColumn)
  
}





