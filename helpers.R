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

