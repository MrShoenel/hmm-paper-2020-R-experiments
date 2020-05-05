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
  # and we replace those names by "_t1" (t-1)
  colnames(ds) <- gsub("\\.\\.\\d+$", "_t_1", colnames(ds))
  
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
