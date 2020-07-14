relu <- function(x) {
  return(max(0.1 * x, x))
}

relu_d1 <- function(x) {
  if (x <= 0) return(0.1)
  return(1)
}

swish <- function(x) {
  return(x * sigmoid(x))
}

swish_d1 <- function(x) {
  ex <- exp(x)
  return(ex * (x + ex + 1) / (ex + 1)**2)
}

sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

sigmoid_d1 <- function(x) {
  return(sigmoid(x) * sigmoid(-x))
}

glorot_weights <- function(numIn, numOut) {
  n <- numIn * numOut
  l <- sqrt(6 / (numIn + numOut))
  return(runif(n, min = -l, max = l))
}


m1 <- function(x_i, w_h, b_h, w_o, b_o, act.fn = relu, act.fn.derive = relu_d1) {
  num_inputs <- length(x_i)
  num_hidden <- length(w_h) / num_inputs
  w_per_hidden <- length(w_h) / num_hidden
  num_output <- length(b_o)
  w_per_output <- length(w_o) / num_output

  if (num_hidden != length(b_h)) {
    stop("Number of hidden units not equal to number of hidden biases.")
  }
  if (num_output != length(b_o)) {
    stop("Number of output units not equal to number of output biases.")
  }
  
  # Compute dot-product and activation for each hidden unit:
  H <- matrix(nrow = 1, ncol = num_hidden)
  for (i in 1:num_hidden) {
    off <- (i-1) * w_per_hidden + 1
    w <- w_h[off:(off + w_per_hidden - 1)]
    H[1, i] <- swish(w %*% x_i + b_h[i])
  }
  
  O <- matrix(nrow = 1, ncol = num_output)
  for (i in 1:num_output) {
    off <- (i-1) * w_per_output + 1
    w <- w_o[off:(off + w_per_output - 1)]
    O[1, i] <- sigmoid(w %*% H[1, ] + b_o[i])
  }
  
  return(O[1, ])
}


loss_wRSS_m1 <- function(X, Y, w_h, b_h, w_o, b_o) {
  N <- nrow(X)
  lOl <- ncol(Y)
  e <- lOl:1
  
  l <- 0
  for (i in 1:N) {
    for (j in 1:lOl) {
      temp <- m1(X[i, ], w_h = w_h, b_h = b_h, w_o = w_o, b_o = b_o)
      l <- l +
        (e[j] * Y[i, j]**2) -
        (2 * e[j] * Y[i, j] * temp[j]) +
        e[j] * temp[j]**2
    }
  }
  
  return(l)
}


compute_grad_wRSS_m1 <- function(X, Y, w_h, b_h, w_o, b_o) {
  num_inputs <- ncol(X)
  num_hidden <- length(w_h) / num_inputs
  w_per_hidden <- length(w_h) / num_hidden
  num_output <- length(b_o)
  w_per_output <- length(w_o) / num_output
  
  if (num_hidden != length(b_h)) {
    stop("Number of hidden units not equal to number of hidden biases.")
  }
  if (num_output != length(b_o)) {
    stop("Number of output units not equal to number of output biases.")
  }
  
  
  N <- nrow(X)
  lOl <- ncol(Y)
  e <- lOl:1
  
  l_w_o <- length(w_o)
  l_b_o <- length(b_o)
  l_w_h <- length(w_h)
  l_b_h <- length(b_h)
  
  grad_w_o <- rep(0, l_w_o)
  grad_b_o <- rep(0, l_b_o)
  grad_w_h <- rep(0, l_w_h)
  grad_b_h <- rep(0, l_b_h)
  
  modp <- function(x, m) {
    temp <- x %% m
    if (temp == 0) return(m)
    return(temp)
  }
  
  for (i in 1:N) {
    for (j in 1:lOl) {
      # In each gradient, we need to compute the hidden layer
      # and the following first derivate of sigmoid:
      H <- matrix(nrow = 1, ncol = num_hidden)
      for (n in 1:num_hidden) {
        off <- (n-1) * w_per_hidden + 1
        w <- w_h[off:(off + w_per_hidden - 1)]
        H[1, n] <- swish(w %*% X[i, ] + b_h[n])
      }
      
      w <- w_o[((j-1) * num_hidden + 1):(j * num_hidden)]
      grad_all <- sigmoid_d1(w %*% H[1, ] + b_o[j])
      m1_all <- m1(X[i, ], w_h = w_h, b_h = b_h, w_o = w_o, b_o = b_o)
      
      
      # Now, for each of the four types of parameters, we'll
      # have a bit of a different gradient. We do them type
      # by type. We start with each k-th output-weight, then
      # continue with each k-th output bias, k-th hidden weight,
      # and finally each k-th hidden bias.
      
      # 1)
      for (k in 1:l_w_o) {
        w <- w_h[((modp(k, num_hidden) - 1) * num_inputs + 1):(modp(k, num_hidden) * num_inputs)]
        m1_prime <- grad_all * swish(w %*% X[i, ] + b_h[modp(k, num_hidden)])
        
        grad_w_o[k] <- grad_w_o[k] +
          2 * e[j] * (m1_all[j] * m1_prime - Y[i, j] * m1_prime)
      }
      
      # 2)
      # The partial derivative for each output-bias has no 'k' -
      # it depends (and corresponds) only on j.
      m1_prime <- grad_all # They are equivalent in this case, as
      # we'd only multiply by 1.
      grad_b_o[j] <- grad_b_o[j] +
        2 * e[j] * (m1_all[j] * m1_prime - Y[i, j] * m1_prime)
      
      # 3)
      for (k in 1:l_w_h) {
        single_w <- w_o[ceiling(k / num_inputs) + num_hidden * (j-1)]
        single_x_i <- X[i, modp(k, num_inputs)]
        
        w <- w_h[((ceiling(k / num_inputs) - 1) * num_inputs + 1):(ceiling(k / num_inputs) * num_inputs)]
        b <- b_h[ceiling(k / num_inputs)]
        m1_prime <- grad_all * swish_d1(w %*% X[i, ] + b) * single_x_i * single_w
        
        grad_w_h[k] <- grad_w_h[k] +
          2 * e[j] * (m1_all[j] * m1_prime - Y[i, j] * m1_prime)
      }
      
      # 4)
      for (k in 1:l_b_h) {
        single_w <- w_o[(j-1) * num_hidden + k]
        
        w <- w_h[((k-1) * num_inputs + 1):(k * num_inputs)]
        b <- b_h[k]
        m1_prime <- grad_all * swish_d1(w %*% X[i, ] + b) * single_w
        
        grad_b_h[k] <- grad_b_h[k] +
          2 * e[j] * (m1_all[j] * m1_prime - Y[i, j] * m1_prime)
      }
    }
  }
  
  return(list(
    grad_w_o = grad_w_o,
    grad_b_o = grad_b_o,
    grad_w_h = grad_w_h,
    grad_b_h = grad_b_h
  ))
}



gradient_descent_m1 <- function(X, Y, w_h_0, b_h_0, w_o_0, b_o_0, epochs = 1e2, learning_rate = 1e-3, precision = 1e-4, batch_size = -1) {
  hist_loss <- c(loss_wRSS_m1(X, Y, w_h = w_h_0, b_h = b_h_0, w_o = w_o_0, b_o = b_o_0))
  
  w_o <- w_o_0
  b_o <- b_o_0
  w_h <- w_h_0
  b_h <- b_h_0
  
  for (i in 1:epochs) {
    use_X <- X
    use_Y <- Y
    if (batch_size > 0) {
      rns <- sample(1:nrow(X), size = batch_size)
      use_X <- matrix(data = X[rns, ], nrow = length(rns))
      use_Y <- matrix(data = Y[rns, ], nrow = length(rns))
    }
    
    grads <- compute_grad_wRSS_m1(
      X = use_X, Y = use_Y,
      w_h = w_h, b_h = b_h, w_o = w_o, b_o = b_o)
    
    step_w_o <- grads$grad_w_o * learning_rate
    step_b_o <- grads$grad_b_o * learning_rate
    step_w_h <- grads$grad_w_h * learning_rate
    step_b_h <- grads$grad_b_h * learning_rate
    
    w_o <- w_o - step_w_o
    b_o <- b_o - step_b_o
    w_h <- w_h - step_w_h
    b_h <- b_h - step_b_h
    
    hist_loss <- c(
      hist_loss, loss_wRSS_m1(X, Y, w_h = w_h, b_h = b_h, w_o = w_o, b_o = b_o))
    
    if (all(c(step_w_o, step_b_o, step_w_h, step_b_h) < precision)) {
      print("Stopping early, as the improvement is less than the threshold.")
      break
    }
  }
  
  return(list(
    hist_loss = hist_loss,
    w_o = w_o,
    b_o = b_o,
    w_h = w_h,
    b_h = b_h
  ))
}












