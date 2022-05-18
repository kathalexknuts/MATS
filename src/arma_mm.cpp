// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

//' Matrix Product in RcppArmadillo.
//'
//' @param m numeric matrix
//' @param m2 numeric matrix
//' @return matrix product of m and m2
// [[Rcpp::export(arma_mm)]]
arma::mat arma_mm(const arma::mat& m, const arma::mat& m2) {
  return m * m2;
};
