require(matrixcalc)
require(mvtnorm)

# for given phi and lambda, construct positive definite function
# @params: phi: target function, d: input dim, s: otput dim, Lambda: pd matrix
pd_function = function(phi, d, s, Lambda = NULL) {
  # compatibility check
  if(is.null(Lambda)) {
    Lambda = diag(s)
  } else {
    if(nrow(Lambda) != s || ncol(Lambda) != s) {
      stop("Incorrect dimension for Lambda. Lambda must be s x s dim matrix")
    }
    if(!is.positive.definite(Lambda)) {
      stop("Lambda must be positive definite matrix")
    }
  }
  f = function(x) {
    return(crossprod(phi(x), Lambda %*% phi(x)))
  }
  return(f)
}

# Sampling from various distributions
# @params: N: num of samples, d; dim of space, distribution: sampling dist
# @return: X:matrix[N, d]
dist_sampling = function(N, d, distribution = "uniform", ...) {
  x = switch(distribution,
         uniform = {
           x = runif(N * d, ...)
           matrix(x, nrow = N, ncol = d)
         },
         gaussian = mvtnorm::rmvnorm(N, ...),
         stop("Not implemented yet.")
      )
  return(x)
}
