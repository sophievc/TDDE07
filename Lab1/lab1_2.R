library(invgamma)

# -----------A-------------
income <- c(14,25,45,25,30,33,19,50,34,67)
#income <- c(1,1,1,1,1,1,1,1,1,1)
m <- 3.5
n <- length(income)
n_draws = 10000

t <- 0

for (v in income) {
  t <- t + (log(v) - m) ^ 2
}

tSquared <- t / n

# simulate postarior draws

X_draw <- rchisq(n_draws, n)
sigma_sq = n*tSquared / X_draw

#theoretical

theoretical <- function(theta, v, s) {
  return ( ((v/2) ^ (v/2)) / gamma(v/2)  * (s ^ v) * (theta ^ (-(v/2+1))) * exp((-v * s ^ 2) / (2 * theta) ))
}

range <- seq(0,10,by=0.01)
mx <- 0 * range

for (i in 1:length(range)) {
  mx[i] <- theoretical(range[i], length(income), sqrt(tSquared))
}

# plot
hist(sigma_sq, 100, freq=FALSE) 
lines(range, mx, lwd=3, type='l', col='red', xlab = 'sigma')

# -----------B-------------
z <- sqrt(sigma_sq/2)
G <- 2*pnorm(z)-1

hist(G, 100, freq = FALSE, xlab="posterior of Gini coefficient", ylab="") # Plot where y = values and x = index of the value in the vector

# -----------C-------------  

cred_int <- quantile(G, probs = c(0.025, 0.975))
G_dens = density(G)
y_ordered = G_dens$y[order(-G_dens$y)]
x_ordered = G_dens$x[order(-G_dens$y)]
dens_mass = sum(G_dens$y)
sum <- 0
current_mass <- 0

for (i in 1:length(y_ordered)) {
  current_mass <- y_ordered[i] + sum
  if ((current_mass/dens_mass) > 0.95) {
    break
  } else {
    sum <- current_mass
  }
}

a <- min(x_ordered[1:i])
b <- max(x_ordered[1:i])

plot(density(G), 
     col='blue', 
     xlim=c(0.1,0.8), 
     ylim=c(0,9),
     main='')
legend('topright',
       legend = c('Highest Posterior Density', 'Equal Tail Interval'),
       fill = c('green', 'red'))
abline(v=cred_int[1], col='red')
abline(v=cred_int[2], col='red')
abline(v=a, col='green')
abline(v=b, col='green')