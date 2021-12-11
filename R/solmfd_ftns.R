#' solution manifold points sampling
#'
#' @param N int: number of output data
#' @param phi target function: R^d -> R^s
#' @param d int: input dimension
#' @param prior str: type of prior distribution. "Gaussian", "t" are currently provided
#' @param s int: output dimension
#' @param gamma double: parameter of gradient descent step size
#' @param Lambda mat[s, s]: p.d. function. Identity by default
#' @param tol double: convergence threshold
#' @param num_iter int: if algorithm does not converge, it iterates for num_iter times.
#' @param ... additional parameter for prior distribution if needed. If not provided, each distribution has its default
#'
#' @return final_points: mat[N, d]: N number of points in R^d, which are points in solution manifold.
#' @export
#'
#' @examples
#' N = 10
#' phi = function(mu, sd) {return(pnorm(2, mu, sd) - pnorm(-5, mu, sd) - 0.5)}
#' d = 2
#' s = 1
#' res_points = sol_mfd_points(N, phi, d, s, mu = c(2, 2.5), sigma = 0.2 * diag(2))
#' head(res_points)
#' phi(res_points[1, ]) # are they on solution manifold?
sol_mfd_points = function(N, phi, d, s, prior = "gaussian", gamma = 0.1, Lambda = NULL, tol = 1e-07, num_iter = 1000, ...) {
  # construct positive definite function
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
  final_points = matrix(0, nrow = N, ncol = d)
  iter = 0
  params = list(...)
  quadratic_f = pd_function(phi, Lambda)
  grad_f = grad_of_quadratic_f(phi, Lambda)
  if(prior == "gaussian") {
    # set mu
    if(exists('mu', params)) {
      mu = params$mu
      # compatibility check
      if(length(mu) != d) {
        stop("dimension of mu mismatches with d.")
      }
    } else{
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
      if(!matrixCalc:: is.positive.definite(sigma)) {
        stop("sigma must be positive definite.")
      }
    } else{
      print("No sigma declaration. Set to default identity matrix.")
      sigma = diag(d)
    }
    final_points = sol_mfd_grad_descent_mvn_cpp(N, d, quadratic_f, grad_f, gamma, tol, num_iter, mu, sigma)
  } else if (prior == "uniform") {
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
    final_points = sol_mfd_grad_descent_munif_cpp(N, d, quadratic_f, grad_f, gamma, tol, num_iter, lower, upper)
  } else{
    stop("Not implemented yet. Please use gaussian or uniform as a prior.")
  }
  return(final_points)
}

#' solution manifold points sampling with Posterior density
#'
#' @param X matrix[n, m]. input data. n is sample number and m is dimension of data.
#' @param prob_density function: P(X, Z). It is in fact P(X|Z)
#' @param N int: number of data to ptoduce in solution manifold
#' @param phi function from R^d -> R^s
#' @param d int: input dimension
#' @param s int: output dimension
#' @param k function: kernel function
#' @param prior str: type of prior distribution. Extra argument ... will be passed to a distribution parameter
#' @param gamma double: parameter of gradient descent step size
#' @param Lambda positive definite function. Default identity
#' @param tol tolerance
#' @param num_iter number of maximum iteration
#' @param ... parameters of prior
#'
#' @return final_points: mat[N, d + 1]: N number of points in R^d+1, which are (R^d dim vector = points in solution manifold, density)
#' @export
#' @examples
#' k = function(x) {return(dnorm(x / 2))}
#' N = 10
#' phi = function(mu, sd) {return(pnorm(2, mu, sd) - pnorm(-5, mu, sd) - 0.5)}
#' d = 2
#' s = 1
#' prob_density = function(x, theta) {return(dnorm(x, mean = theta[1], sd = theta[2]))}
#' n = 100
#' X = rnorm(n, 1.5, 3)
#' post_density_solmfd(X, prob_density, N, phi, d, s, k, mu = c(2, 2.5), sigma = 0.2 * diag(2))
post_density_solmfd = function(X, prob_density, N, phi, d, s, k, h = 1, prior = "gaussian", gamma = 0.1, Lambda = NULL, tol = 1e-07, num_iter = 1000, ...) {
  # obtain points
  points = sol_mfd_points(N, phi, d, prior, s, gamma, Lambda, tol, num_iter, ...)
  # kernel
  dist_between_row = between_row_dist(points, points)
  normalized = dist_between_row / h
  kerneled = k(normalized)
  # mean by N
  rho_i = rowMeans(kerneled)
  # calculate the \hat(pi_i) = \pi(Z_i) \prod(p(X_j | Z+i))
  params = list(...)
  if(prior == "gaussian") {
    if(exists('mu', params)) {
      mu = params$mu
    } else {
      mu = rep(0, d)
    }
    if(exists('sigma', params)) {
      sigma = params$sigma
    } else {
      sigma = diag(d)
    }
    pi_z = mvtnorm::dmvnorm(points, mu, sigma)
  } else if(prior == "uniform") {
    if(exists('min', params)) {
      lower = params$min
    } else {
      lower = 0
    }
    if(exists("max", params)) {
      upper = params$max
    } else {
      upper = 1
    }
    pi_z = dunif(points, lower, upper)
  }
  # conditioning X ftn for vectoization
  z_conditioning_x = function(z) { return(prob_density(X, z)) }
  apply(points, 1, z_conditioning_x)
  pi_i = exp(log(pi_z) + colSum(log(apply(points, 1, z_conditioning_x))))
  omega_i = pi_i / rho_i
  return(cbind(points, omega_i))
}

