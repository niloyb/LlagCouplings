#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
double logistic_logtarget_c(NumericVector chain_state, const NumericVector & Y, const NumericMatrix & X, double sigma2){
  int p = X.cols();
  int n = X.rows();
  // double logsigma2 = chain_state(p+1);
  double logsigma2 = log(sigma2);
  NumericVector xbeta(n);
  double eval = 0;
  NumericVector logvalues(n);
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < p; j ++){
      xbeta(i) += X(i,j) * chain_state(j);
    }
    eval -= xbeta(i) * (1. - Y(i));
    logvalues(i) = log(1 + exp(- xbeta(i)));
  }
  double maxlogvalues = max(logvalues);
  eval -= (sum(logvalues - maxlogvalues) + n * maxlogvalues);
  eval += - 0.5 * (sum(chain_state*chain_state)) / sigma2;
  return eval;
}


// [[Rcpp::export]]
NumericVector logistic_gradlogtarget_c(NumericVector chain_state, const NumericVector & Y, const NumericMatrix & X, double sigma2){
  int p = X.cols();
  int n = X.rows();
  double logsigma2 = log(sigma2);
  NumericVector xbeta(n);
  // NumericVector precomputed_vec(n);
  double precomp = 0.;
  NumericVector gradient(p);
  // std::fill(gradient.begin(), gradient.end(), 0);
  // gradient(0) = - alpha / sigma2;
  gradient = -chain_state / sigma2;
  // for (int j = 0; j < p; j++){
  //   gradient(j) = - chain_state(j) / sigma2;
  // }
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < p; j ++){
      xbeta(i) += X(i,j) * chain_state(j);
    }
    // precomputed_vec(i) = (1. + exp(xbeta(i)));
    precomp = (1. + exp(xbeta(i)));
    for (int j = 0; j < p; j++){
      gradient(j) += (Y(i)-1) * X(i,j) + X(i,j) / precomp;
    }
  }
  return gradient;
}
