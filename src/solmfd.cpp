#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

/*
 * @params: x - init vector, f - quadratic form
 */
// [[Rcpp::export]]
solmfd_descent(arma::colvec& x, const arma::colvec& (*f)(arma::colvec&), double gamma, double tol) {
  // iterate until converge
  while(f(x) > tol) {
    // TODO: Implement gradient descent, especially consider gradient calculation.
    x += -(gamma * grad(f)(x))
  }
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# timesTwo(42)
*/
