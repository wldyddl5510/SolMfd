#' solution manifold points sampling
#' @param N int: number of data
#' @param phi function from R^d -> R^s
#' @param d int: input dimension
#' @param prior str: type of prior distribution. Extra argument ... will be passed to a distribution parameter
#' @param s int: output dimension
#' @param gamma double: parameter of gradient descent step size
#' @param Lambda: mat[s, s]: positive definite matrix. Default is identity function
#' @param tol: double: convergence threshold
#'
#' @return final_points: mat[N, d]: N number of points in R^d, which are points in solution manifold.
#' @export
#' @examples
#' N = 10
#' phi = function(mu, sd) {return(pnorm(2, mu, sd) - pnorm(-5, mu, sd))}
#' d = 2
#' s = 1
#' sol_mfd_points(N, phi, d, s)
sol_mfd_points = function(N, phi, d, prior = "uniform", s, gamma = 0.1, Lambda = NULL, tol = 1e-07, num_iter = 1000, ...) {
  # TODO: Sampling point clouds in the solution manifold
  # construct positive definite function
  f = pd_function(phi, d, s, Lambda)
  num_points = 0
  final_points = matrix(0, nrow = N, ncol = d)
  iter = 0
  # iterate until obtain desired num of points, or limit iteration reached
  # TODO: Fix this part with properly designed cpp function
  while(num_points < N || iter < num_iter) {
    points = dist_sampling((N - num_points), d, distribution = prior)
    # gradient descent until converge to solution manifold
    # TODO: implement gradient descent (cpp? R?)
    points = grad_descent(f, points, gamma, tol)
    # keep or discard the point
    if(f(points) < tol) {
      num_points = num_points + nrow(points)
      final_points = points
    }
  }
  return(final_points)
}

#' solution manifold points sampling with Posterior density
#' @param h function: normalizer
#' @param N int: number of data
#' @param phi function from R^d -> R^s
#' @param d int: input dimension
#' @param prior str: type of prior distribution. Extra argument ... will be passed to a distribution parameter
#' @param s int: output dimension
#' @param gamma double: parameter of gradient descent step size
#' @param Lambda: mat[s, s]: positive definite matrix. Default is identity function
#' @param tol: double: convergence threshold
#'
#' @return final_points: mat[N, d + 1]: N number of points in R^d+1, which are c(R^d dim vector, density)
#' @export
#' @examples
#' h = function(x) {return(pnorm(x / 2))}
#' N = 10
#' phi = function(mu, sd) {return(pnorm(2, mu, sd) - pnorm(-5, mu, sd))}
#' d = 2
#' s = 1
#' post_density_solmfd(h, N, phi, d, s)
post_density_solmfd = function(h, N, phi, d, prior, s, gamma, Lambda, tol, num_iter, ...) {
  # obtain points
  points = sol_mfd_points(N, phi, d, prior, s, gamma, Lambda, tol, num_iter, ...)
  # call cpp function that performs density calculation
  score_vec = density_cpp(points, h)
  # TODO: Implement computing posterior
}

