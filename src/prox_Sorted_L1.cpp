#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
#include "proxSortedL1.h"
}

// [[Rcpp::export]]
NumericVector prox_sorted_L1_C(NumericVector y, NumericVector lambda, IntegerVector order) {
  size_t n = y.size();
  NumericVector x(n);
  evaluateProx(y.begin(), lambda.begin(), x.begin(), n, order.begin());
  return x;
}
