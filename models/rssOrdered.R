RSS <- function(e, j, y_true, y_pred) {
  return((y_true - y_pred)^2)
}

RSS_d1 <- function(e, j, y_true, m, m_d1) {
  return(2 * m_d1 * (m - y_true))
}

wRSS <- function(e, j, y_true, y_pred) {
  return(e[j] * ((y_true - y_pred)^2))
}

wRSS_d1 <- function(e, j, y_true, m, m_d1) {
  return(2 * e[j] * m_d1 * (m - y_true))
}

omeRSS <- function(e, j, y_true, y_pred) {
  res <- y_true - y_pred
  exponent <- (e[j] - 2) + 2 * e1071::sigmoid(res^2)
  return(2^exponent)
}

omeRSS_d1 <- function(e, j, y_true, m, m_d1) {
  res <- y_true - m
  return(log(2) * res * m_d1 * sigmoid_d1(res^2) *
           (-2^(e[j] + 2 * e1071::sigmoid(res^2))))
}

ome2RSS <- function(e, j, y_true, y_pred) {
  res <- y_true - y_pred
  exponent <- e[j] * 2 * e1071::sigmoid(res^2)
  return(2^exponent)
}

ome2RSS_d1 <- function(e, j, y_true, m, m_d1) {
  res <- y_true - m
  return(e[j] * log(2) * res * m_d1 * sigmoid_d1(res^2) *
           (-2^(2 + 2 * e[j] * e1071::sigmoid(res^2))))
}


rssOrdered <- function(yTrue = matrix(), yPred = matrix()) {
  yTrue <- as.matrix(yTrue)
  yPred <- as.matrix(yPred)
  
  N <- nrow(yTrue)
  lOl <- ncol(yTrue)
  
  err <- 0
  
  for (i in 1:N) {
    for (j in 1:lOl) {
      residual <- (yPred[i, j] - yTrue[i, j])**2
      if (residual == 0) {
        next
      }
      err <- err + 2**(lOl - j + 2 * e1071::sigmoid(residual))
    }
  }
  
  return(err)
}


rssOrdered2 <- function(yTrue = matrix(), yPred = matrix()) {
  yTrue <- as.matrix(yTrue)
  yPred <- as.matrix(yPred)
  
  N <- nrow(yTrue)
  lOl <- ncol(yTrue)
  
  err <- 0
  
  for (i in 1:N) {
    for (j in 1:lOl) {
      residual <- (yPred[i, j] - yTrue[i, j])**2
      if (residual == 0) {
        next
      }
      
      err <- err + 2**((lOl - j + 1) * (2 * e1071::sigmoid(residual)))
    }
  }
  
  return(err)
}
