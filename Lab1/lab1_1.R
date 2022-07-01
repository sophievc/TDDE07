alpha <- 2
beta <- 2
s <- 14
f <- 6
n <- 20
f <- n-s
nDraws <- c(10, 100, 1000)

# 1a - Random draws from posterior of θ | y ~ Beta(α + s, β + f)

xGrid <- seq(0.001, 0.999, by=0.001)
posterior = dbeta(xGrid, alpha + s, beta + (n-s))
trueMean <- (alpha + s)/(alpha + beta + s + f)
trueStdv <- sqrt(((alpha + s) * ( beta + f)) / ((alpha + beta + n)^2 * (alpha + beta + n + 1)))
par(mfrow=c(1,3))
for (draws in nDraws) {
  random = rbeta(draws, alpha+s, beta+(n-s))
  hist(random, xlim = c(0,1), freq = FALSE, breaks=10)
  lines(xGrid, posterior, type='l', col ='red')
}
dev.off() 

par(mfrow=c(1,2))
# Mean and standard deviation
mean <- c()
stdv <- c()
for (i in 1:1000) {
  random = rbeta(i, alpha+s, beta+(n-s))
  mean <- c(mean, mean(random))
  stdv <- c(stdv, sd(random))
}

# Plot mean to see convergence
plot(mean, type='l', col='green', xlab = 'n draws')
abline(h = trueMean)
legend('topright', legend = 'True mean', col = 'black', lwd = 1)

# Plot standard diviation to see convergence
plot(stdv, type='l', col='black', xlab='n draws')
abline(h = trueStdv, col='red')
legend('topright', legend = 'True stdv', col = 'red', lwd = 1)


# 1b - Simulation of Pr(θ < 0.4|y), compared with the real value using pbeta
nDraws <- 10000
random = rbeta(nDraws, alpha+s, beta+(n-s))
print(sum(random < 0.4)/length(random))
real <- pbeta(.4, alpha + s, beta + (n-s))


# 1c - Posterior distribution of the log-odds by simulation
nDraws <- 10000
random = rbeta(nDraws, alpha+s, beta+(n-s))
phi <- log(random/(1 - random))
plot(density(phi), lwd=1, type='l', col='red', main = expression('Density of' ~ phi))
