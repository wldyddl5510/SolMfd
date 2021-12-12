# create quadratic function from given phi and Lambda
pd_function = function(phi, Lambda) {
  quadratic_f = function(x) {
    return(crossprod(phi(x), Lambda %*% phi(x)))
  }
  return(quadratic_f)
}

# calculate the gradient of quadratic function
grad_of_quadratic_f = function(phi, Lambda) {
  grad_f = function(x) {
    return(crossprod(2 * phi(x), Lambda %*% rootSolve::gradient(phi, x)))
  }
  return(grad_f)
}

# calculate the gradient of general function.
grad_of_f = function(f) {
  grad_f = function(x) {
    return(rootSolve::gradient(f, x))
  }
  return(grad_f)
}

# measure distance between rows of two matrices.
between_row_dist = function(X, M) {
  mat_dist = outer(rowSums(X), rowSums(M^2), '+') - tcrossprod(X, 2 * M)
  return(mat_dist)
}

# sampling from given prior. Also conduct compatibility check.
sampling_from_dist = function(d, prior, ...) {
  params = list(...)
  if(prior == "gaussian") {
    # set mu
    if(exists('mu', params)) {
      mu = params$mu
      # compatibility check
      if(length(mu) != d) {
        stop("dimension of mu mismatches with d.")
      }
    } else{ # default
      print("No mu declaration. Set to default 0.")
      mu = rep(0, d)
    }
    # set sigma
    if(exists('sigma', params)) {
      sigma = params$sigma
      # compatibility check
      if((nrow(sigma) != d) || (ncol(sigma) != d)) {
        stop("dimension of sigma mismatches with d.")
      }
      if(!matrixcalc::is.positive.definite(sigma)) {
        stop("sigma must be positive definite.")
      }
    } else{ # default
      print("No sigma declaration. Set to default identity matrix.")
      sigma = diag(d)
    }
    # sampling function
    sampling_function = function(N) {
      (mvtnorm::rmvnorm(N, mean = mu, sigma = sigma))
    }
  } else if (prior == "uniform") { # set parameter for uniform distribution
      # set lower
      if(exists('min', params)) {
        lower = params$min
        if(!is.numeric(lower)) {
          stop("min must be numeric.")
        }
      } else {
        print("No min declaration. Set to default 0.")
        lower = 0
      }
      # set upper
      if(exists('max', params)) {
        upper = params$upper
        if(!is.numeric(upper)) {
          stop("max must be numeric.")
        }
      } else {
        print("No max declaration. Set to default 1")
        upper = 1
      }
      # compatibility check
      if(lower >= upper) {
        stop("min must be smaller than max.")
      }
      # sampling function
      sampling_function = function(N) {
      return(matrix(stats::runif(N * d, min = lower, max = upper), N, d))
      }
    } else{
      stop("Not implemented yet. Please use gaussian or uniform as a prior.")
    }
  return(sampling_function)
}

# density from given prior. Also conduct compatibility check.
density_from_dist = function(d, prior, ...) {
  params = list(...)
  if(prior == "gaussian") {
    # set mu
    if(exists('mu', params)) {
      mu = params$mu
      # compatibility check
      if(length(mu) != d) {
        stop("dimension of mu mismatches with d.")
      }
    } else{ # default
      print("No mu declaration. Set to default 0.")
      mu = rep(0, d)
    }
    # set sigma
    if(exists('sigma', params)) {
      sigma = params$sigma
      # compatibility check
      if((nrow(sigma) != d) || (ncol(sigma) != d)) {
        stop("dimension of sigma mismatches with d.")
      }
      if(!matrixcalc::is.positive.definite(sigma)) {
        stop("sigma must be positive definite.")
      }
    } else{ # default
      print("No sigma declaration. Set to default identity matrix.")
      sigma = diag(d)
    }
    # sampling function
    density_function = function(x) {
      (mvtnorm::dmvnorm(x, mean = mu, sigma = sigma))
    }
  } else if (prior == "uniform") { # set parameter for uniform distribution
    # set lower
    if(exists('min', params)) {
      lower = params$min
      if(!is.numeric(lower)) {
        stop("min must be numeric.")
      }
    } else {
      print("No min declaration. Set to default 0.")
      lower = 0
    }
    # set upper
    if(exists('max', params)) {
      upper = params$upper
      if(!is.numeric(upper)) {
        stop("max must be numeric.")
      }
    } else {
      print("No max declaration. Set to default 1")
      upper = 1
    }
    # compatibility check
    if(lower >= upper) {
      stop("min must be smaller than max.")
    }
    # sampling function
    density_function = function(x) {
      return(stats::dunif(x, min = lower, max = upper))
    }
  } else{
    stop("Not implemented yet. Please use gaussian or uniform as a prior.")
  }
  return(density_function)
}
