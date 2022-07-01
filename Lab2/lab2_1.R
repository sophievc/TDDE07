library(mvtnorm)
library(invgamma)
data = read.table("TempLinkoping.txt", header=TRUE)

######## A ########
I <- diag(3)
mu <- c(-10, 100, -90)
omega <- 10 * I
v <- 4
s2 <- 8
x <- seq(0, 10, by=0.01) 


# Inverse chisquare from lab1  
draw_sigma <- function(v, s2) {
  X_draw <- rchisq(1, v)
  return (v * s2 / X_draw)
}

# Given temperature regression function
temp <- function(beta, time) {
  return (beta[1] + beta[2] * time + beta[3] * (time ^ 2)) + rnorm(0,1)
}

plot(data$time, data$temp, pch=16)
for (i in 1:8) {
  sigma_sq = draw_sigma(v, s2)
  beta = rmvnorm(1, mu, sigma_sq*solve(omega))
  y = temp(beta[1,], data$time)
  lines(data$time, y, col='red')
}

dev.off()


######## B ########
#Plot the marginal posterior also for the sigma2 parameter
n_draws <- 1000
ones <- c(rep(1, length(data$time)))
x2 <- data$time^2
X <- cbind(ones, data$time, x2)
beta_hat <- solve(t(X)%*%X)%*%t(X)%*%data$temp
omega_n <- t(X) %*% X + omega
mu_n <- solve(omega_n) %*% (t(X) %*% X %*% beta_hat + omega %*% mu)

v_n <- v + length(data$time)
s2_n <- (v*s2 + (t(data$temp)%*%data$temp + t(mu)%*%omega%*%mu - t(mu_n)%*%omega_n%*%mu_n))/v_n

sigma_post = draw_sigma(v_n, s2_n)
beta_post = rmvnorm(n_draws, mu_n, solve(omega_n)*sigma_post[1])

Y = matrix(nrow=n_draws, ncol=366)
for (i in 1:n_draws) {
  Y[i,] = temp(beta_post[i,], data$time)
}

y_med = c()
y_up = c()
y_low = c()

for (i in 1:366) {
  y_med = c(y_med, median(Y[,i]))
  y_low = c(y_low, quantile(Y[,i], 0.025))
  y_up = c(y_up, quantile(Y[,i], 0.975))
}

par(mfrow=c(1,3))
hist(beta_post[,1], freq = FALSE, xlab = "Beta 0", main='')
hist(beta_post[,2], freq = FALSE, xlab = "Beta 1", main='')
hist(beta_post[,3], freq = FALSE, xlab = "Beta 2", main='')
dev.off()

sigma_draws <- c()
for (i in 1:n_draws) {
  sigma_draws <- c(draw_sigma(v_n, s2_n), sigma_draws)
}
hist(sigma_draws, breaks = 50)

plot(data$temp, type='p', col='lightgray', main = 'Data Plot', xlab = "Time", ylab = "Temperature")
lines(y_med, type='l')
lines(y_low, type='l', col='red')
lines(y_up, type='l', col='blue')

legend("bottomright", 
       legend = c("Mean","Lower", "Upper"),
       fill = c("black", "blue", "red"))

######## C ########

x_max = - 366 * beta_post[,2]/(2 * beta_post[,3])
hist(x_max, main='', xlab='Time')
dev.off()

######## D ########

# Prior
# Mu_0 = [...mu, 0, 0, 0, 0] due to suspicion that the introduced variables might not be needed
# Omega_0, want low variance (high bias to avoid overfitting) => High diagonal values in omega_0for the new variables, with the rest set to zero since we can't say anything about the covariance between the new variables


#Good! You should use mean 0 and high omega diagonal values 
#for shrinkage of the higher order polynomial terms. 
#However, using your previously estimated posterior parameters 
#for beta0-beta2 is a bit weird, unless you have some new data. 
#Otherwise I would say: use the same prior for beta0-beta2 as you did in 1b).

new_mu = c(mu_0, rep(0, 5))
new_omega <- 100 * diag(8)
new_omega[1:3,1:3] = omega_n




