# a)
library(Rcpp)
sourceCpp("utils.cpp")

data_ <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 1.09,
           1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 6.83)

x_grid <- seq(-20, 50, length.out = 1000)
y_value <- numeric(1000)

for(i in 1:1000){
  y_value[i] <- get_log_likelihood(x_grid[i], data_, get_ell_prime)
}

plot(x_grid, y_value, type = "l",
     main = "The First Derivative of Log-Likelihood",
     xlab = expression(theta),
     ylab = expression(paste("l'(", theta, ")")))
abline(h = 0, lty = 2, col = "grey")


# b)
# i) Bisection
res_bisection <- do_bisection(get_ell_prime, data_, 0.0, 5.0, 0.01)

# ii) Newton-Raphson
res_newton <- do_newton_raphson(get_ell_prime,
                                get_ell_double_prime,
                                data_, 1.0, 0.01)

# iii) Fisher Scoring
res_fisher <- do_fisher_scoring(get_ell_prime,
                                get_ell_double_prime,
                                data_, 1.0, 0.01)

# iv) Secant
res_secant <- do_secant(get_ell_prime, data_,
                        0.0, 2.0, 0.01)


# c)
library(dplyr)

knitr::kable(
  rbind(data.frame(res_bisection),
        data.frame(res_newton),
        data.frame(res_fisher),
        data.frame(res_secant)) %>% 
    `rownames<-`(c("Bisection", "Newton-Raphson",
                   "Fisher Scoring", "Secant")) %>% 
    `colnames<-`(c("$\\hat{\\theta}$", "#(Iterations)", "Convergence"))
  )


# d)
# I used the relative convergence criterion to stop the iteration.
# \begin{equation*}
# \frac{\mid x^{\left(t+1\right)}-x^{\left(t\right)}\mid}{\mid x^{\left(t\right)}\mid}<\epsilon=0.01
# \end{equation*}
# This criterion is chosen because the root of the first derivative of the log-likelihood function is small, so I did not want to be worried about the unit.

# e) 
# Using the fact that $\text{Var}\left(\theta\right)=-\mathbb{E}\left[\ell''\left(\theta\right)\right]^{-1}$,
# we can obtain the variance of our estimate by calculating
# \begin{equation*}
# \text{Var}\left(\hat{\theta}\right)=-\left[\ell''\left(\hat{\theta}\right)\right]^{-1}
# \end{equation*}
# In my opinion, the estimate calculated by Newton-Raphson algorithm is the best even though its standard error is slightly bigger than that of secant method.
# The reason is that Newton-Raphson algorithm gives a more stable result.
# So, the standard error of my estimate is
# \begin{equation*}
# \text{Var}\left(\hat{\theta}\right)=-\left[\ell''\left(\hat{\theta}\right)\right]^{-1}
# \end{equation*}
sqrt(-1 / get_log_likelihood(res_newton$root, data_, get_ell_double_prime))

# f)
# To initialize, I first checked the graph and set the values near the suspected root as the inital values.
# More specifically, if I choose initial values of $10$ and $20$, then both of the evaluated values are negative, so bisection method cannot be used.
# In addition to that, secant method cannot be used because it will diverge.
# And suppose that my initial value is $10$, then both Newton-Raphson and Fisher Scoring cannot be used as the root diverges.

do_bisection(get_ell_prime, data_, 10.0, 20.0, 0.01)
do_secant(get_ell_prime, data_, 10.0, 20.0, 0.01)
do_newton_raphson(get_ell_prime, get_ell_double_prime, data_, 10.0, 0.01)
do_fisher_scoring(get_ell_prime, get_ell_double_prime, data_, 10.0, 0.01)

# g) 

data_ <- c(data_,
           c(-8.34, -1.73, -0.40, -0.24, 0.60, 0.94, 1.05, 1.06, 1.45, 1.50,
             1.54, 1.72, 1.74, 1.88, 2.04, 2.16, 2.39, 3.01, 3.01, 3.08,
             4.66, 4.99, 6.01, 7.06, 25.45))

res_bisection <- do_bisection(get_ell_prime, data_, 0.0, 5.0, 0.01)

res_newton <- do_newton_raphson(get_ell_prime,
                                get_ell_double_prime,
                                data_, 1.0, 0.01)

res_fisher <- do_fisher_scoring(get_ell_prime,
                                get_ell_double_prime,
                                data_, 1.0, 0.01)

res_secant <- do_secant(get_ell_prime, data_,
                        0.0, 2.0, 0.01)

knitr::kable(
  rbind(data.frame(res_bisection),
        data.frame(res_newton),
        data.frame(res_fisher),
        data.frame(res_secant)) %>% 
    mutate(se = 
             c(sqrt(-1 / get_log_likelihood(res_bisection$root,
                                            data_, get_ell_double_prime)),
               sqrt(-1 / get_log_likelihood(res_newton$root,
                                            data_, get_ell_double_prime)),
               sqrt(-1 / get_log_likelihood(res_fisher$root,
                                            data_, get_ell_double_prime)),
               sqrt(-1 / get_log_likelihood(res_secant$root,
                                            data_, get_ell_double_prime)))) %>% 
    `rownames<-`(c("Bisection", "Newton-Raphson",
                   "Fisher Scoring", "Secant")) %>% 
    `colnames<-`(c("$\\hat{\\theta}$", "#(Iterations)",
                   "Convergence", "Standard Error"))
)

# 3.

# male = 1; female = 0;
df <- data.frame(gender = c(1, 1, 1, 1, 0, 0, 0, 0),
                 n = c(41, 213, 127, 142, 67, 211, 133, 76),
                 r = c(9, 94, 53, 60, 11, 59, 53, 28),
                 x1 = c(0, 2, 4, 5, 0, 2, 4, 5))

df <- rbind(
  as.data.frame(lapply(df, rep, df$r)) %>% 
    mutate(y = 1),
  as.data.frame(lapply(df, rep, df$n - df$r)) %>% 
    mutate(y = 0)) %>% 
  select(y, x1, gender)

fit <- logit(df$y, as.matrix(cbind(1, df[, 2:3])),
             c(-1.0, 0.1, 0.3), 0.01)

knitr::kable(
  data.frame(coef = fit$coefs,
             se = sqrt(diag(solve(fit$hessian)))) %>% 
    mutate(z = coef / se,
           p = round(pnorm(abs(coef / se),
                           lower.tail = FALSE), 3)) %>% 
    `rownames<-`(c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$")) %>% 
    `colnames<-`(c("Estimate", "Std. Error",
                   "z-value", "Pr(>\\|z\\|)"))
)


X <- cbind(1, rep(seq(0, 5, length.out = 100), 2),
           c(rep(1, 100), rep(0, 100)))
Xbeta <- X %*% fit$coefs
prob <- exp(Xbeta) / (1 + exp(Xbeta))
plot(X[X[, 3] == 1, 2], prob[X[, 3] == 1],
     ylim = c(0, 1),
     type = "l", col = "red", lwd = 2,
     xlab = "Amount of Coffee Consumption",
     ylab = "P(Cancer)",
     main = "Probability of getting cancer")
lines(X[X[, 3] == 0, 2], prob[X[, 3] == 0],
      col = "blue", lwd = 2, lty = 2)
points(df[df$gender == 1, 2] + runif(sum(df$gender == 1), 0, 0.1),
       df[df$gender == 1, 1], col = "red")
points(df[df$gender == 0, 2] + runif(sum(df$gender == 0), 0, 0.1),
       df[df$gender == 0, 1] + 0.01, 
       col = "blue")
legend(x = "topright", inset = 0.08, legend = c("Male", "Female"),
       lty = c(1, 2), col = c("red", "blue"), lwd = 2)

sourceCpp("utils.cpp")
runif()
