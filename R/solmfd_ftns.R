source("utils.R")

# point clouds sampling in the solution manifold
# @params: N: data input, phi: target function, d:input dim, support: input space, s: output dim, gamma: step size, Lambda: p.d. matrix
# @returns: N data size
sol_mfd_points = function(N, phi, d, prior = "uniform", s, gamma = 0.1, Lambda = NULL, tol = 1e-07, num_iter = 1000, ...) {
  # TODO: Sampling point clouds in the solution manifold
  # construct positive definite function
  f = pd_function(phi, d, s, Lambda)
  num_points = 0
  final_points = matrix(0, nrow = N, ncol = d)
  iter = 0
  # iterate until obtain desired num of points, or limit iteration reached
  while(num_points < N || iter < num_iter) {
    points = dist_sampling((N - num_points), d, distribution = prior)
    # gradient descent until converge to solution manifold
    while(f(points) > tol) {
      # TODO: implement gradient descent
      points = points - gamma * grad(f)
    }
    # keep or discard the point
    if(f(points) < tol) {
      num_points = num_points + nrow(points)
      final_points = points
    }
  }
  return(final_points)
}

# prior distribution sampling in the solution manifold
# @params: h: normalizer function, rest: same with sol_mfd points algorithm
# @returns: N tuples of (point, density)
density_on_solmfd = function(h, N, phi, d, prior, s, gamma, Lambda, tol, num_iter, ...) {
  # obtain points
  points = sol_mfd_points(N, phi, d, prior, s, gamma, Lambda, tol, num_iter, ...)

}

