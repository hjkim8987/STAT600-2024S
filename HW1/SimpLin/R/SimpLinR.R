#' Fit a simple linear regression model
#' 
#' This function fits a simple linear regression model using an n by 1
#' predictor variable (x) and an n by 1 response variable (y).
#' It throws an error if
#'    i) length of x and y is different
#'    ii) x is not numeric
#'    iii) y is not numeric
#' It returns a fitted result, which includes
#'    estimated regression coefficients
#'    standard errors of the estimates
#'    95% confidence interval of each coefficient
#'    predicted values
#'    residuals
#' 
#' @param x A predictor variable, n by 1 numeric vector
#' @param y A response variable, n by 1 numeric vector
#' @param significance_level Significance level for confidence interval
#' @return A list of the fit
#' @export
#' 
#' @examples
#' n <- 100
#' x <- rnorm(n)
#' beta0 <- 1
#' beta1 <- -1
#' y <- beta0 + beta1 * x + rnorm(n)
#' SimpLinR(x, y)
SimpLinR <- function(x, y, significance_level = 0.05){
  n <- length(y)

  if(n != length(x)){
    stop("Length!")
  }

  if(!is.numeric(x)){
    stop("Type of x")
  }

  if(!is.numeric(y)){
    stop("Type of y")
  }

  return(SimpLinCpp(x, y, n, significance_level))
}
