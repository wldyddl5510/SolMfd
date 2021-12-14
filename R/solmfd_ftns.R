#' solution manifold points sampling
#'
#' @param N int: number of output data
#' @param phi target function: R^d -> R^s
#' @param d int: input dimension
#' @param prior str: type of prior distribution. "Gaussian", "t" are currently provided
#' @param s int: output dimension
#' @param gamma double: parameter of gradient descent step size
#' @param Lambda mat (s * s): p.d. function. Identity by default
#' @param tol1 double: convergence threshold for manifold convergence
#' @param tol2 double: convergence threhold for gradient descent algorithm. Should be smaller than tol1
#' @param num_iter int: if algorithm does not converge, it iterates for num_iter times.
#' @param ... additional parameter for prior distribution if needed. If not provided, each distribution has its default
#'
#' @return final_points: mat (N * d): N number of points in R^d, which are points in solution manifold.
#' @export
#'
#' @examples
#' set.seed(10)
#' N = 10
#' phi = function(x) {return(pnorm(2, x[[1]], x[[2]]) - pnorm(-5, x[[1]], x[[2]]) - 0.5)}
#' d = 2
#' s = 1
#' res_points = sol_mfd_points(N, phi, d, s, prior = "uniform")
#' head(res_points)
#' phi(res_points[1, ]) # are they on solution manifold?
#' plot(res_points, xlab = "mean", ylab = "sigma", type = 'o') # how they are distributed
sol_mfd_points = function(N, phi, d, s, prior = "gaussian", gamma = 0.005, Lambda = NULL, tol1 = 1e-07, tol2 = 1e-15, num_iter = 100000, ...) {
  # construct positive definite function
  # compatibility check
  if(is.null(Lambda)) {
    Lambda = diag(s)
  } else {
    if(nrow(Lambda) != s || ncol(Lambda) != s) {
      stop("Incorrect dimension for Lambda. Lambda must be s x s dim matrix")
    }
    if(!matrixcalc::is.positive.definite(Lambda)) {
      stop("Lambda must be positive definite matrix")
    }
  }
  # result saving matrix
  final_points = matrix(0, nrow = N, ncol = d)
  iter = 0
  # quadratic form and grad of this to run and evaluate algorithm
  quadratic_f = pd_function(phi, Lambda)
  grad_f = grad_of_quadratic_f(phi, Lambda)
  # designate sampling function
  sampling_function = sampling_from_dist(d, prior, ...)
  final_points = sol_mfd_grad_descent(N, d, quadratic_f, grad_f, sampling_function, gamma, tol1, tol2, num_iter)
  return(final_points)
}



#' solving constraint likelihood function using the solution manifold.
#'
#' @param nll function: negative log-likelihood given data X
#' @param C function: constraint function
#' @param theta vector: initial parameter
#' @param s int: output dim of function C
#' @param alpha gradient descent step for nll update
#' @param gamma gradient descent step for solution manifold algorithm
#' @param Lambda positive definite matrix for solution manifold algorithm. Default is identity
#' @param tol1 double: convergence threshold for manifold convergence
#' @param tol2 double: convergence threhold for gradient descent algorithm. Should be smaller than tol1
#' @param num_iter maximum number of iterations for gradient descent.
#' @param num_iter2 number of iteration for all processes.
#'
#' @return theta_traj: matrix (num_iter2 * length(theta)) containing trajactory of theta updates. Last row is a final result.
#' @export
#'
#' @examples
#' # init value
#' set.seed(10)
#' # num of samples
#' n = 100
#' # data distribution
#' X = rnorm(n, mean = 1.5, sd = 3)
#' # negative log likelihood
#' nll = function(theta) {return(-sum(dnorm(X, theta[[1]], theta[[2]], log = TRUE)))}
#' # constraint
#' C = function(x) {return(pnorm(2, x[[1]], x[[2]]) - pnorm(-5, x[[1]], x[[2]]) - 0.5)}
#' theta = runif(2, 1, 3)
#' theta_updated = constraint_likelihood(nll, C, theta, 1)
#' # plot the convergences
#' const_val = apply(theta_updated, 1, C)
#' plot(x = seq(1, nrow(theta_updated)), const_val, xlab = "step", ylab = "constraint", type = 'o')
#' nll_val = apply(theta_updated, 1, nll)
#' plot(x = seq(1, nrow(theta_updated)), nll_val, xlab = "step", ylab = "nll", type = 'o')
#' plot(theta_updated, xlab = "mean", ylab = "sigma", type = 'o')
#' lines(theta_updated[seq(3, nrow(theta_updated), by = 2), ],  col = 'red')
constraint_likelihood = function(nll, C, theta, s, alpha = 0.005, gamma = 0.005, Lambda = NULL, tol1 = 1e-07, tol2 = 1e-15, num_iter = 100000, num_iter2 = 20) {
  d = length(theta)
  # compatibility check
  if(is.null(Lambda)) {
    Lambda = diag(s)
  } else {
    if(nrow(Lambda) != s || ncol(Lambda) != s) {
      stop("Incorrect dimension for Lambda. Lambda must be s x s dim matrix")
    }
    if(!matrixcalc::is.positive.definite(Lambda)) {
      stop("Lambda must be positive definite matrix")
    }
  }
  grad_nll = grad_of_f(nll)
  quadratic_f = pd_function(C, Lambda)
  grad_f = grad_of_quadratic_f(C, Lambda)
  grad_c = grad_of_f(C)
  theta_traj = matrix(NA, 2 * num_iter2 + 1, d)
  theta_traj[1, ] = theta
  # likelihood update
  for(i in 1:num_iter2) {
    # gradient descent to negative log likelihood (i.e. gradient ascending to log-likelihood.)
    # only one step
    theta = theta - alpha * grad_nll(theta)
    theta_traj[2 * i, ] = theta
    # descent to manifold, until it reaches manifold.
    theta = grad_descent(theta, grad_f, d, gamma, tol2, num_iter)
    if(quadratic_f(theta) > tol1) {
      stop("this theta does not converge to manifold. pick different theta.")
    }
    theta_traj[2 * i + 1, ] = theta
    # stopping criterion: grad_nll(theta) \in span(row(grad_C))
    # grad_nll \in row_grad_c
    grad_c_mat = grad_c(theta)
    # check grad_nll can be approximated by linear projection of rows of grad_c
    linear_projection = crossprod(grad_c_mat, solve(tcrossprod(grad_c_mat, grad_c_mat)) %*% grad_c_mat)
    target = grad_nll(theta)
    row_space_error = sum((target %*% linear_projection - target)^2)
    # error < tol1
    if (row_space_error < tol2) {
      break
    }
  }
  # return with removing NA's
  return(theta_traj[!rowSums(!is.finite(theta_traj)),])
}

#' Posterior density of solution manifold points
#'
#' @param X matrix: (n * m). input data. n is sample number and m is dimension of data.
#' @param prob_density function: P(X, Z). It is in fact P(X|Z)
#' @param points matrix (N * d). point clouds on manifold.
#' @param k function: kernel function
#' @param h normalizing factor. Default = 1
#' @param prior str: type of prior distribution. Extra argument ... will be passed to a distribution parameter
#' @param ... parameters of prior
#'
#' @return omega_i vec (N): density of each row in matrix: points's row.
#' @export
#'
#' @examples
#' set.seed(10)
#' k = get("dnorm", mode = 'function')
#' prob_density = function(x, theta) {return(dnorm(x, mean = theta[[1]], sd = theta[[2]]))}
#' n = 100
#' X = rnorm(n, 1.5, 3)
#' points = sol_mfd_points(N, phi, d, s, prior = "uniform")
#' res_with_density = post_density_solmfd(X, prob_density, points, k)
#' head(res_with_density)
post_density_solmfd = function(X, prob_density, points, k, h = 1, prior = "gaussian", ...) {
  if(!is.matrix(points)) {
    stop("data must be matrix.")
  }
  # input dim
  d = ncol(points)
  # kernel
  dist_between_row = between_row_dist(points, points)
  normalized = dist_between_row / h
  kerneled = k(normalized)
  # mean by N
  rho_i = rowMeans(kerneled)
  # calculate the \hat(pi_i) = \pi(Z_i) \prod(p(X_j | Z+i))
  # calculate pi_i
  # prior setting
  # designate sampling function
  density_function = density_from_dist(d, prior, ...)
  pi_z = density_function(points)
  # conditioning X ftn for vectoization
  z_conditioning_x = function(z) { return(prob_density(X, z)) }
  pi_i = exp(log(pi_z) + colSums(log(apply(points, 1, z_conditioning_x))))
  omega_i = pi_i / rho_i
  return(omega_i)
}

