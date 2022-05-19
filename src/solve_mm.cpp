// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

//' Matrix Inversion in RcppArmadillo.
//'
//' @param Am numeric matrix
//' @return Inversion of Am
// [[Rcpp::export(solve_mm)]]
Rcpp::List flowCalcCpp(const arma::mat &Am) {
  arma::mat B = inv(Am);
  return Rcpp::List::create( Rcpp::Named("Imp") = B);
}
