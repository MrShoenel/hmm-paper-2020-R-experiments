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
  cl <- makePSOCKcluster(detectCores())
  registerDoParallel(cl)
  
  result <- tryCatch(expr, error=function(cond) {
    if (!is.null(errorValue)) {
      return(errorValue)
    }
    return(cond)
  }, finally = {
    stopCluster(cl)
    registerDoSEQ()
    cl <- NULL
    gc()
  })
  return(result)
}




