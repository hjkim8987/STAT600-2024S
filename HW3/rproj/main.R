library(Rcpp)
sourceCpp("EM.cpp")

##
generate_data <- function(n, p, lambda, mu){
  y <- numeric(n)
  for(i in 1:n){
    if(runif(1) < p){
      y[i] <- rexp(1, rate = lambda)
    }else{
      y[i] <- rexp(1, rate = mu)
    }
  }
  return(y)
}

set.seed(2024)

n <- 100
p <- 0.25
lambda <- 1
mu <- 2

nsim <- 100
y <- list()
for(i in 1:nsim){
  y[[i]] <- generate_data(n, p, lambda, mu)
}

##
res <- data.frame(matrix(nrow = nsim, ncol = 3))
colnames(res) <- c("p", "lambda", "mu")
for(i in 1:nsim){
  res[i, ] <- exponential_mixture_EM(y[[i]],
                                     c(0.2, 0.5, 1.5))
}

res <- rbind(round(colMeans(res), 4),
             round(apply(res, 2, sd), 4))
rownames(res) <- c("Mean", "StDev")
res

##
y100 <- y[[100]]
res_y100 <- exponential_mixture_EM(y100, c(0.2, 0.5, 1.5))

res_boot <- data.frame(matrix(nrow = 10000, ncol = 3))
colnames(res_boot) <- c("p", "lambda", "mu")
for(i in 1:10000){
  y_bootstrap <- sample(y100, n, replace = TRUE)
  res_boot[i, ] <- exponential_mixture_EM(y_bootstrap, res_y100)
}
res_boot <- rbind(round(apply(res_boot, 2, mean), 4),
                  round(sqrt(apply(res_boot, 2, var) * 9999 / 10000), 4))
rownames(res_boot) <- c("EM Estimates", "SE (Bootstrap)")
res_boot

##
estimates <- data.frame(matrix(nrow = nsim, ncol = 3))
ses <- data.frame(matrix(nrow = nsim, ncol = 3))
res_boot <- data.frame(matrix(nrow = 10000, ncol = 3))
coverage <- data.frame(matrix(nrow = nsim, ncol = 3))
for(i in 1:nsim){
  estimates_each <- exponential_mixture_EM(y[[i]], c(0.2, 0.5, 1.5))
  for(b in 1:10000){
    y_bootstrap <- sample(y[[i]], n, replace = TRUE)
    res_boot[b, ] <- exponential_mixture_EM(y_bootstrap, estimates_each)
  }
  estimates[i, ] <- estimates_each
  ses[i, ] <- sqrt(apply(res_boot, 2, var) * 9999 / 10000)
  for(l in 1:3){
    ci <- quantile(res_boot[, l], probs = c(0.025, 0.975))
    coverage[i, l] <- (ci[1] <= p) & (ci[2] >= p)
  }
}
bias <- estimates - matrix(rep(c(p, lambda, mu), nsim),
                           nrow = nsim, byrow = TRUE)

res6 <- rbind(apply(estimates, 2, mean),
              apply(bias, 2, mean),
              apply(ses, 2, mean),
              apply(coverage, 2, mean))
rownames(res6) <- c("Estimate", "Bias", "SE", "Coverage")
colnames(res6) <- c("$$p$$", "$$\\lambda$$", "$$\\mu$$")
knitr::kable(round(res6, 4))



