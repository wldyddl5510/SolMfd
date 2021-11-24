#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

/*
 * @params
 */
// [[Rcpp::export]]
arma::colvec grad_descent_cpp(arma::colvec& (*f)(arma::colvec&), arma::colvec& init, double gamma, double tol) {
  arma::colev x_old(arma::size(x), fill::arma.zeors);
  arma::colvec x = init;
  double diff = arma::norm(x - x_old);
  // iterate until converge
  while(diff > tol) {
    x_old = x;
    // TODO: calculate grad of f: named grad_f
    grad_f = f(x);
    x += -gramma * grad_f;
    diff = arma::norm(x - x_old);
  }
  return(x)
}
