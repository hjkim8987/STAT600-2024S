#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Create a design matrix including intercept term
// 
// @param x numeric vector
// @param n integer
// @return design matrix n by 2 numeric matrix, (1, x)
// @export
arma::mat get_degisn_matrix(arma::vec x, int n){
  arma::mat X(n, 2, arma::fill::ones);

  for(int i = 0; i < n; i++){
    X(i, 1) = x(i);
  }

  return X;
}

// Calculate estimated regression coefficients
// 
// @param X n by 2 numeric matrix
// @param y numeric vector
// @return coefficient estimates 2 by 1 numeric matrix
arma::mat get_coefficients(arma::mat X, arma::vec y){
  return (X.t() * X).i() * X.t() * y;
}

// Calculate fitted values
// 
// @param X n by 2 numeric matrix
// @param beta 2 by 1 numeric matrix
// @return fitted values n by 1 numeric vector
arma::vec get_predicted_values(arma::mat X, arma::mat beta){
  return X * beta;
}

// Calculate residuals
// 
// @param y n by 1 numeric vector
// @param yhat n by 1 numeric vector
// @return residuals n by 1 numeric vector
arma::vec get_residuals(arma::vec y, arma::vec yhat){
  return y - yhat;
}

// Calculate Sum of Squares of Residuals
// 
// @param resids n by 1 numeric vector
// @return sse sum of squared errors, double
double get_sse(arma::vec resids){
  return arma::as_scalar(resids.t() * resids);
}

// Calculate \hat{Sigma}^2
// an estimated variance matrix of regression coefficients
// 
// @param X n by 2 numeric matrix
// @param sse double
// @param n integer, number of observations
// @return estimated variance matrix 2 by 2 numeric matrix
arma::mat get_variance_estimate(arma::mat X, double sse, int n){
  return (X.t() * X).i() * sse / (double) (n - X.n_cols);
}

// Calculate (1 - \alpha) * 100% confidence interval of an estimate
// 
// @param estimate double
// @param sig_level double, significance level (\alpha)
// @param variance double, estimated variance of the estimate
// @return confidence interval lower and upper bounds, 2 by 1 numeric vector
arma::vec get_confidence_interval(double estimate,
                                  double sig_level,
                                  double variance){
  arma::vec ci(2, arma::fill::zeros);
  double q = R::qnorm(1.0 - sig_level / 2.0, 0.0, 1.0, true, false);

  ci(0) = estimate - q * sqrt(variance);
  ci(1) = estimate + q * sqrt(variance);

  return ci;
}

//' Cpp function to obtain
//'   estimated regression coefficients
//'   standard errors of the estimates
//'   95% confidence interval of each coefficient
//'   predicted values
//'   residuals
//' 
//' @param x predictor, n by 1 numeric vector
//' @param y response, n by 1 numeric vector
//' @param n number of observations, integer
//' @param sig_level significance level of the confidence interval, double
//' @return result of fit List
//' @export
// [[Rcpp::export]]
List SimpLinCpp(arma::vec x, arma::vec y, int n, double sig_level){
  arma::mat X = get_degisn_matrix(x, n);
  arma::mat beta = get_coefficients(X, y);
  arma::vec yhat = get_predicted_values(X, beta);
  arma::vec residuals = get_residuals(y, yhat);
  arma::mat sigma2hat = get_variance_estimate(X, get_sse(residuals), n);
  arma::mat ci_mat(beta.n_rows, 2, arma::fill::zeros);
  for(int i = 0; i < beta.n_rows; i++){
    ci_mat.col(i) =
      get_confidence_interval(beta(i), sig_level, sigma2hat(i, i));
  }

  return List::create(
    Named("coefficients") = beta,
    Named("residuals") = residuals,
    Named("fitted.values") = yhat,
    Named("vcov") = sigma2hat,
    Named("confidence.intervals") = ci_mat.t()
  );
}

