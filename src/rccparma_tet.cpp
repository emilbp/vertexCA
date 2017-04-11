#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat matrix_multiplication(int n, arma::mat A) {
  return A * n;
}
