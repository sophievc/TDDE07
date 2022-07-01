library('rstan')

########## A ##########
mu <- 10
sigma_sq <- 2
T <- 200 
x <- mu

par(mfrow=c(3,1))

AR <- function(mu, prev_x, phi, sigma_sq) {
  error <- rnorm(1, 0, sqrt(sigma_sq))
  return (mu + phi * (prev_x - mu) + error)
}

phis <- c(-1,0.1,0.95)

for(phi in phis) {
  result = c(x)
  for (i in 2:T) {
    x <- AR(mu, x, phi, sigma_sq)
    result <- c(result, x)
  }
  plot(result, type = 'line', main = bquote(phi ~ ' =' ~ .(phi)))
}


########## B ##########

mu <- 10
sigma_sq <- 2
T <- 200 
x <- mu

phis <- c(0.3, 0.95)

X = c()
for(phi in phis) {
  result <- c(x)
  for (i in 2:T) {
    x <- AR(mu, x, phi, sigma_sq)
    result <- c(result, x)
  }
  X <- rbind(X, result)
}

model <- stan_model('StanNormalModel.stan')

warmup=1000
iter=2000
niter = 4*(iter-warmup)

fitX <- sampling(model, data = list(T=200, x=X[1,]), iter = iter, warmup = warmup)
fitY <- sampling(model, data = list(T=200, x=X[2,]), iter = iter, warmup = warmup)

# Print the fitted model
print(fitX,digits_summary=3) # Extract posterior samples
print(fitY,digits_summary=3) # Extract posterior samples

postDrawsX <- extract(fitX)
postDrawsY <- extract(fitX)

par(mfrow=c(1,2))
plot(postDrawsX$mu, postDrawsX$phi, main = expression(phi ~ '= 0.3'), ylab=expression(phi), xlab=expression(mu))
lines(mean(postDrawsX$mu), mean(postDrawsX$phi), type='p', col="red")
# legend('topright', legend=expression(mean(postDrawsX$mu) ~ ',' ~ mean(postDrawsX$phi)))
plot(postDrawsY$mu, postDrawsY$phi, main = expression(phi ~ '= 0.95'), ylab=expression(phi), xlab=expression(mu))
lines(mean(postDrawsY$mu), mean(postDrawsY$phi), type='p', col="red")
# legend('topright', legend=expression(().mean(postDrawsY$mu) ~ ',' ~ ().mean(postDrawsY$phi)))

# plot convergence
traceplot(model)

########## C ##########
library(ggplot2)
data <- read.table('campy.txt', header=TRUE)

N <- length(data$c) 

model <- stan_model('StanPoissonModel.stan')
fit <- sampling(model, data = list(N=N, c=data$c), iter = 2000, warmup = 1000)

# Print the fitted model
print(fit,digits_summary=3) # Extract posterior samples
postDraws <- extract(fit)

x <- postDraws$x

xMean <- c()
thetaUp <- c()
thetaLow <- c()

for (i in 1:length(x[1,])) {
  xMean <- c(xMean, exp(mean(x[,i])))
  thetaUp <- c(thetaUp, exp(quantile(x[,i], probs = 0.975)))
  thetaLow <- c(thetaLow, exp(quantile(x[,i], probs = 0.025)))
}

# draws <- c()
# for(i in 1:140) {
#   draws <- c(draws, rpois(1, exp(xMean[i])))
# }

draws <- rpois(140, exp(xMean))

plot(data$c, pch=16, col='gray77', ylab='', main='Data plot')
lines(xMean, col='lightcyan3', lwd=2)
lines(thetaUp, col='lightcyan3', lwd=2)
lines(thetaLow, col='coral3', lwd=2)
legend('topleft',legend = c('Data', 'Mean','Interval'), 
       col = c('gray77', 'coral3', 'lightcyan3'), lwd=2)
dev.off()
# geom_ribbon(mapping = aes(seq(1,140,1), ymin = thetaLow, ymax = thetaUp), alpha = 0.25)

# Do traceplots of the first chain
par(mfrow = c(1,1))
# plot(postDraws$mu[1:(1000)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
# traceplot(fit)
# Bivariate posterior plots
pairs(fit)
# plot(fit)

########## D ##########

data <- read.table('campy.txt', header=TRUE)

N <- length(data$c) 

model <- stan_model('StanPoissonModelD.stan')
fit <- sampling(model, data = list(N=N, c=data$c), iter = 2000, warmup = 1000)

# Print the fitted model
print(fit,digits_summary=3) # Extract posterior samples
postDraws <- extract(fit)

x <- postDraws$x

xMean <- c()
thetaUp <- c()
thetaLow <- c()

for (i in 1:length(x[1,])) {
  xMean <- c(xMean, mean(exp(x[,i])))
  thetaUp <- c(thetaUp, quantile(exp(x[,i]), probs = 0.975))
  thetaLow <- c(thetaLow, quantile(exp(x[,i]), probs = 0.025))
}

plot(data$c, pch=16, col='gray77', ylab='', main='Data plot')
lines(xMean, col='coral3', lwd=2)
lines(thetaUp, col='lightcyan3', lwd=2)
lines(thetaLow, col='lightcyan3', lwd=2)
legend('topleft',legend = c('Data', 'Mean','Interval'), 
       col = c('gray77', 'coral3', 'lightcyan3'), lwd=2)


dev.off()
