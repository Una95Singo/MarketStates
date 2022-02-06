#-----------------------------------------------------------------
#---------------------- Geometric brownian motion---------------------
#----------------------------------------------------------------
gbm_loop <- function(nsim = 100, t = 25, mu = 0, sigma = 0.1, S0 = 100, dt = 1./365) {
  gbm <- matrix(ncol = nsim, nrow = t)
  for (simu in 1:nsim) {
    gbm[1, simu] <- S0
    for (day in 2:t) {
      epsilon <- rnorm(1)
      dt = 1 / 365
      gbm[day, simu] <- gbm[(day-1), simu] * exp((mu - sigma * sigma / 2) * dt + sigma * epsilon * sqrt(dt))
    }
  }
  return(gbm)
}

library(tidyverse)
nsim <- 100
t <- 100
mu <- 0
sigma <- 0.1
S0 <- 100
gbm <- gbm_loop(nsim, t, mu, sigma, S0)

plot(y=gbm[,4], x=1:100, type = 'l')

Ret = diff(log(gbm))

C1 = t(Ret)%*%(Ret) # integrated - realised volatility
C2 = cov(Ret)
solve(C1)
is.positive.semi.definite(C1)
1/prod(eigen(C1)$values)


Decomp = chol(cov(Ret)) %*% t(chol(cov(Ret)))

det(LoGo.solve(cov(Ret)))

library(matrixcalc)
is.positive.definite(cov(GRet))
det(cov(GRet))

shrunk = 0.1*cov(GRet)+0.*diag(cov(GRet))

det(shrunk)
is.positive.semi.definite()
Ret


D = 100
N = 1000
Ret = matrix(NA, nrow = N, ncol = D)
for(i in 1:D){
  Ret[, i] = rnorm(N, mean = 0, 0.014)
}
var = cov(Ret)
var.inv = solve(var)
dim(var)

plot(y = (100*cumsum(rowMeans(Ret))+100), x = 1:N, type='l')


is.positive.semi.definite(var)
is.positive.definite(var)

det(var)
log(det(var.inv))
log((det(chol(var.inv)))^2)
sum(2*log(diag(chol(var.inv))))


alpha = 0.5

shrunk = alpha*var+(1-alpha)*diag(1, nrow = 50)
shrunk = var

det(var*10^3)
det(var*10^3)
det(solve(var*10^3))

det(solve((diag(diag(var, names = T), nrow = 100))))

is.singular.matrix(var*10^5)
is.singular.matrix(cov(GRet))
diag()

