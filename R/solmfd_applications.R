#' solving constraint likelihood function for the solution manifold.
#' @param X mat: data matrix
#' @param L function: likelihood
#' @param C function: constraint function
#' @param Theta
#' @param iter int: number of maximum iteration
#' @param tol: double: convergence threshold
#'
#' @return theta: list of parameters of interests that maximize the likelihood under constraint
#' @export
#' @examples
#'
const_likelihood = function(X, L, C, d, s, Theta, alpha, iter = 100, tol = 1e-07) {
  # TODO: solve constraint likelihood on solution manifold
  for(i in 1:iter) {
    # likelihood update
    theta = theta + alpha * grad(L)
    # descent to manifold
    f = pd_function(C, d, s)
    while(f(theta) > tol) {
      theta = theta - gamma * grad(f)
    }
  }
  return(theta)
}
