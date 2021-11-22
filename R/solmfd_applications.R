source("utils.r")

# solving constraint likelihood function for the solution manifold.
# @params: X: data, L: likelihood, C: constraint, Theta: parameter space, iter: num of iteration
# @params: alpha: update rate of likelihood, gamma: update of manifold
# @params: d, s: input-output
# @return: theta: maximal likelihood parameter under constraint
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
}
