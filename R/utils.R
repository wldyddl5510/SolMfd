pd_function = function(phi, Lambda) {
  quadratic_f = function(x) {
    return(crossprod(phi(x), Lambda %*% phi(x)))
  }
  return(quadratic_f)
}

grad_of_quadratic_f = function(phi, Lambda) {
  grad_f = function(x) {
    return(crossprod(2 * phi(x), Lambda %*% rootSolve::gradient(phi, x)))
  }
  return(grad_f)
}

grad_of_f = function(f) {
  grad_f = function(x) {
    return(rootSolve::gradient(f, x))
  }
  return(grad_f)
}

between_row_dist = function(X, M) {
  mat_dist = outer(rowSums(X), rowSums(M^2), '+') - tcrossprod(X, 2 * M)
  return(mat_dist)
}
