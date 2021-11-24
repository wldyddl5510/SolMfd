#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

/*
 * @params: x - init vector, f - quadratic form
 */
// [[Rcpp::export]]
arma::colvec solmfd_descent(arma::colvec& x, const arma::colvec& (*f)(arma::colvec&), double gamma, double tol) {
  // iterate until converge
  while(f(x) > tol) {
    // TODO: Implement gradient descent, especially consider gradient calculation.
    x += -(gamma * grad(f)(x));
  }
}

/*
 * @params: x - init vector, f - quadratic form, h: normalizer
 */
// [[Rcpp::export]]
arma::colvec density_cpp(const arma::mat& points, const arma::colvec& (*h)(double)) {
  // iterate
  N = n_row(points);
  d = n_col(points);
  for(int i = 0; i < N ; i++) {
    for(int j = 0 ; j < d ; j++) {
      // TODO: Implement this inner part
    }
  }
}
