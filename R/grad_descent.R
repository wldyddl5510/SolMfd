# gradient descent in solution manifold.
sol_mfd_grad_descent = function(N, d, quadratic_f, grad_f, sampling_function, gamma = 0.1, tol1 = 1e-07, tol2 = 1e-15, num_iter = 1000) {
  final_mat = matrix(0, N, d)
  i = 1
  iter = 0
  while(i <= N) {
    # sampling required points
    curr_points = sampling_function(N - i + 1)
    curr_points = grad_descent(curr_points, grad_f, d, gamma, tol2, num_iter)
    # whether obtained points are in sol_mfd
    eval_points = apply(curr_points, 1, quadratic_f)
    # zero elements
    points_obtained = curr_points[which(eval_points < tol1), ]
    # points in sol_mfd this time
    num_of_points = nrow(points_obtained)
    final_mat[i: (i + num_of_points - 1), ] = points_obtained
    # updated number of samples this time
    i = i + num_of_points
    # reached maximum number of iteration
    if (num_iter < iter) {
      stop("Reached maximum iteration before getting N points. Try larger num_iter or different hyperparams.")
    }
    iter = iter + 1
  }
  return(final_mat)
}

# gradient descending in general case.
grad_descent = function(curr_points, grad_f, d, gamma, tol, num_iter) {
  # for each row in matrix
  if(!is.matrix(curr_points)) {
    curr_points = matrix(curr_points, ncol = d)
  }
  for(i in 1:nrow(curr_points)) {
    curr_point = curr_points[i, ]
    # error: to evaluate the convergence of grad_descent
    error = sum((curr_point)^2)
    new_point = curr_point
    iter = 0
    # perform grad_descent
    while(error > tol) {
      # update
      new_point = curr_point - gamma * grad_f(curr_point)
      # checking convergence
      error = sum((new_point - curr_point)^2)
      # update
      curr_point = new_point
      # check maximum iter
      if(iter > num_iter) {
        stop("Reached maximum iteration before grad_descent converge. Try larger num_iter or different hyperparams.")
      }
      iter = iter + 1
    }
    # update to matrix
    curr_points[i, ] = curr_point
  }
  return(curr_points)
}
