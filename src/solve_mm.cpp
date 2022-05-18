// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Matrix Inversion in RcppArmadillo.
//'
//' @param Am numeric matrix
//' @return Inversion of Am
//' @import RcppArmadillo
// [[Rcpp::export(solve_mm)]]
Rcpp::List flowCalcCpp(const arma::mat &Am) {
  arma::mat B = inv(Am);
  return Rcpp::List::create( Rcpp::Named("Imp") = B);
}
