# gradient descent in solution manifold.
sol_mfd_grad_descent = function(N, d, quadratic_f, grad_f, sampling_function, gamma = 0.1, tol1 = 1e-07, tol2 = 1e-15, num_iter = 1000) {
  final_mat = matrix(0, N, d)
  i = 1
  iter = 0
  while(i <= N) {
    # sampling required points
    curr_point = sampling_function(1)
    curr_point = grad_descent(curr_point, grad_f, d, gamma, tol2, num_iter)
    # whether obtained points are in sol_mfd
    eval_point = as.numeric(quadratic_f(curr_point))
    # no convergence to manifold.
    if(eval_point > tol1) {
      next
    }
    # conv to manifold
    final_mat[i, ] = curr_point
    # updated number of samples this time
    i = i + 1
    # reached maximum number of iteration
    if (num_iter < iter) {
      stop("Reached maximum iteration before getting N points. Try larger num_iter or different hyperparams.")
    }
    iter = iter + 1
  }
  return(final_mat)
}

# gradient descending in general case.
grad_descent = function(curr_point, grad_f, d, gamma, tol, num_iter) {
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
  return(curr_point)
}
