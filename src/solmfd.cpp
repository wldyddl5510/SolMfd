#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;


arma::mat sol_mfd_grad_descent_mvn_cpp(const int N, const int d, Rcpp::Function quadratic_f, Rcpp::Function grad_f, double gamma, double tol, int num_iter, arma::vec& mu, arma::vec& sigma) {
  arma::mat final_mat(N, d);
  int main_iter = 0;
  int i = 0;
  // until reaches desired number of points or maximum iteration.
  while(i < N) {
    arma::rowvec curr_point = mv_norm_cpp(d, mu, sigma);
    curr_point = grad_descent_cpp(curr_point, d, grad_f, gamma, tol, num_iter)
      /*
       arma::rowvec new_point(d);
       double error = arma::accu(square(curr_point));
       int inner_iter = 0;
       // gradient descent
       while((error > tol) && (inner_iter++ < num_iter)) {
       new_point -= gamma % grad_f(curr_point);
       error = arma::accu(square(new_point - curr_point));
       curr_point = new_point;
       }
       // if maximum iteration
       if (num_iter < main_iter++) {
       Rcpp::stop("Reached maximum iteration before getting N points. Try larger num_iter or different prior.")
       }
       */
    // accept or reject the sample
    if(quadratic_f(curr_point) > 0) { //reject
      continue;
    }
    else { //accept. move to next row
      final_mat.row(i) = curr_point;
      i++;
    }
  }
  return final_mat;
}

arma::rowvec mv_norm_cpp(int d, arma::vec mu, arma::mat sigma) {
  arma::mat Y = arma::randn<vec>(d);
  return mu + Y * arma::chol(sigma);
}

arma::mat sol_mfd_grad_descent_munif_cpp(const int N, const int d, Rcpp::Function quadratic_f, Rcpp::Function grad_f, double gamma, double tol, int num_iter, const double lower, const double upper) {
  arma::mat final_mat(N, d);
  int main_iter = 0;
  int i = 0;
  // until reaches desired number of points or maximum iteration.
  while(i < N) {
    arma::rowvec curr_point = multi_uniform_cpp(d, lower, upper);
    curr_point = grad_descent_cpp(curr_point, d, grad_f, gamma, tol, num_iter)
    /*
    arma::rowvec new_point(d);
    double error = arma::accu(square(curr_point));
    int inner_iter = 0;
    // gradient descent
    while((error > tol) && (inner_iter++ < num_iter)) {
      new_point -= gamma % grad_f(curr_point);
      error = arma::accu(square(new_point - curr_point));
      curr_point = new_point;
    }
    // if maximum iteration
    if (num_iter < main_iter++) {
      Rcpp::stop("Reached maximum iteration before getting N points. Try larger num_iter or different prior.")
    }
    */
    // accept or reject the sample
    if(quadratic_f(curr_point) > 0) { //reject
      continue;
    }
    else { // accept. move to next row
      final_mat.row(i) = curr_point;
      i++;
    }
  }
  return final_mat;
}

arma::vec multi_uniform_cpp(int d, double lower = 0, double upper = 1) {
  return (lower + ((upper - lower) % arma::randu<vec>(d)));
}

// perform gradient descent for given gradient function.
arma::vec grad_descent_cpp(arma::vec& point, Rcpp::Function grad_f, double gamma, double tol, int num_iter) {
  arma::vec curr_point = point;
  int d = curr_point.n_elem;
  arma::vec new_point(d);
  int iter = 0;
  double error = arma::accu(square(curr_point));
  while((error > tol) && (iter++ < num_iter)) {
    new_point -= gamma % grad_f(curr_point);
    error = error = arma::accu(square(new_point - curr_point));
    curr_point = new_point;
  }
  // if maximum iteration
  if (num_iter < main_iter++) {
    Rcpp::stop("Reached maximum iteration before getting N points. Try larger num_iter or different prior.")
  }
  return curr_point;
}
