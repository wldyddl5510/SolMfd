#' solving constraint likelihood function using the solution manifold.
#'
#' @param nll function: negative log-likelihood given data X
#' @param C function: constraint function
#' @param theta vector: initial parameter
#' @param prior prior distribution for solution manifold algorithm. Currently "gaussian" or "uniform".
#' @param alpha gradient descent step for nll update
#' @param gamma gradient descent step for solution manifold algorithm
#' @param Lambda positive definite matrix for solution manifold algorithm. Default is identity
#' @param tol tolerance.
#' @param num_iter number of iteration.
#'
#' @return theta: a vector of parameters that maximize the likelihood under constraint. Since the algorithm is stochastic, the result will depends on init value.
#' @export
#'
#' @examples
#' n = 100
#' X = rnorm(n, mean = 1.5, sd = 3)
#' nll = function(theta) {return(dnorm(X, theta[1], theta[2]))}
#' C = function(mu, sd) {return(pnorm(2, mu, sd) - pnorm(-5, mu, sd) - 0.5)}
#' theta = rnorm(1)
#' theta_updated = constraint_likelihood(nll, C, theta)
#' C(theta)
#' C(theta_updated)
constraint_likelihood = function(nll, C, theta, prior = "gaussian", alpha = 0.1, gamma = 0.1, Lambda = NULL, tol = 1e-07, num_iter = 1000) {
  # compatibility check
  if(is.null(Lambda)) {
    Lambda = diag(s)
  } else {
    if(nrow(Lambda) != s || ncol(Lambda) != s) {
      stop("Incorrect dimension for Lambda. Lambda must be s x s dim matrix")
    }
    if(!matrixCalc:: is.positive.definite(Lambda)) {
      stop("Lambda must be positive definite matrix")
    }
  }
  grad_nll = grad_of_f(nll)
  quadratic_f = pd_function(C, Lambda)
  grad_f = grad_of_quadratic_f(C, Lambda)
  for(i in 1:num_iter) {
    # likelihood update
    # gradient descent to negative log likelihood (i.e. gradient ascenting to log-likelihood.)
    theta = grad_descent_cpp(theta, grad_nll, alpha, tol, num_iter)
    # descent to manifold
    while(quadratic_f(theta) > tol) {
      theta = grad_descent_cpp(theta, grad_f, gamma, tol, num_iter)
    }
    # stopping criterion: grad_nll(theta) \in span(row(grad_C))
    grad_C_mat = gradient(C, theta)
    if(rankMatrix(grad_C_mat)[1] == rankMatrix(rbind(grad_C_mat, grad_nll(theta)))[1]) {
      break
    }
  }
  return(theta)
}
