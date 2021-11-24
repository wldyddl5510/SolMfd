#' solving constraint likelihood function for the solution manifold.
#' @param X mat: data matrix
#' @param L function: log-likelihood
#' @param C function: constraint function
#' @param Theta
#' @param alpha: grad descent step for likelihood
#' @param gamma: grad descent step for sol_mfd
#' @param iter int: number of maximum iteration
#' @param tol: double: convergence threshold
#'
#' @return theta: list of parameters of interests that maximize the likelihood under constraint
#' @export
#' @examples
#'
const_likelihood = function(X, L, C, d, s, Theta, alpha = 0.1, gamma = 0.1, iter = 100, tol = 1e-07) {
  # TODO: solve constraint likelihood on solution manifold
  for(i in 1:iter) {
    # likelihood update
    theta = grad_descent(-L, theta, alpha, tol)
    # descent to manifold
    f = pd_function(C, d, s)
    theta = grad_descent(f, theta, gamma, tol)
  }
  return(theta)
}
