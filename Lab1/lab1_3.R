
# A: Plot the posterior distribution of k
# posterior(k | y, mu) = likelihood * prior(k) = prod_prob * prior
# Plot of the normalized posterior with AUC = 1

y <- c(-2.44, 2.14, 2.54, 1.83, 2.01, 2.33, -2.79, 2.23, 2.07, 2.02)
gridWidth <- 0.01
kSeq <- seq(0,10, by=gridWidth)
mu <- 2.39
lambda <- 1
# kappa <- dexp(kSeq, lambda)

mises <- function(k, y, mu) {
  I <- besselI(k,0)
  return ((exp(k * cos(y - mu))) / (2 * pi * I))
}

kPos <- function(k, mu,y) {
  #prod since independent
  return ( prod( mises(k, y, mu) ) * dexp(k))
}

posterior = c()
for (k in kSeq){
  posterior = c(posterior, c(kPos(k, mu, y)))
}

# Plot the normalized posterior with AUC = 1
normConst <- sum(posterior)
plot(kSeq, posterior/(normConst*gridWidth),type='l', ylab = 'Normalized Density')
legend('topright', legend='Mode', fill='red')

# B: Compute the posterior mode of k

kPosMode <- kSeq[which.max(posterior)]
abline(v=kPosMode, col='red', lwd=1)
