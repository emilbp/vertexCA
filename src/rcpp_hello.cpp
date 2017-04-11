#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector convolveCpp(NumericVector a, NumericVector b) {
  int na = a.size(), nb = b.size();
  int nab = na + nb - 1;
  NumericVector xab(nab);
  for (int i = 0; i < na; i++)
    for (int j = 0; j < nb; j++)
      xab[i + j] += a[i] * b[j];
  return xab;
}

// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
//
//
// // This is a simple function using Rcpp that creates an R list
// // containing a character vector and a numeric vector.
// //
// // Learn more about how to use Rcpp at:
// //
// //   http://www.rcpp.org/
// //   http://adv-r.had.co.nz/Rcpp.html
// //
// // and browse examples of code using Rcpp at:
// //
// //   http://gallery.rcpp.org/
// //
//
// // [[Rcpp::export]]
// List rcpp_hello() {
//   //CharacterVector x = CharacterVector::create("foo", "bar");
//   //NumericVector y   = NumericVector::create(0.0, 1.0);
//   //List z            = List::create(x, y);
//   arma::mat A = arma::randn(5, 4);
//   List z = List::create(A, A * 10)
//   return z;
// }
