

mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50")

C = 2
s = 365
N = C*s
spinLabels = c(replicate(365,1), replicate(365,2))
D = 400 # 400 stocks in the market

eta = matrix(c(rnorm(D, mean = 0.5, sd = 1), rnorm(D,mean = -0.5, sd = 1)), nrow = C, ncol = D, byrow = T) 
eta

hist(eta[1,],main = "Histogram of state effects", col = 'lightblue',freq = F)
mycol <- t_col("pink", perc = 50, name = "lt.pink")
hist(eta[2,], add= T, col = mycol,freq = F)

rowMeans(eta) # state means
apply(eta, 1, FUN = sd) # state standard deviation



eps = matrix(c(rnorm(N, mean = 0, sd = 1)), nrow = N, ncol = D, byrow = T)
eps

rowMeans(eps) # state means
apply(eps, 1, FUN = sd) # state standard deviation

gs = c(0.5,0.2)

xi = matrix(NA, N,D)


# double for loop for simualation
for (n in 1:N){
  for (day in 1:D){
    cluster = spinLabels[n]
    xi[n, day] = (sqrt(gs[cluster]) * eta[cluster, day] + eps[n, day])/(sqrt(1+gs[cluster]))
  }
}


plot(xi[1,], type='l', main = 'Single day', xlab = 'Stocks', ylab='returns')

plot(t(xi)[1,], type='l', main = 'Single stock', xlab = 'Stocks', ylab='returns')
plot(t(xi)[20,], type ='l')
plot(cumsum(t(xi)[20,]), type='l')
hist(t(xi)[50,])

DailyReturns = t(xi)

DailyAvgReturns = colMeans(DailyReturns)


plot(DailyAvgReturns, type ='l')
plot(cumsum(DailyAvgReturns), type ='l')

# timeseries plots. 
xi.ts = ts(t((DailyReturns)), start=c(2009, 1), end=c(2014, 12), frequency = 12)
ts.plot(xi.ts, plot.type ='single', col = spinLabels)

