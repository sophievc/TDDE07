
data = read.table("rainfall.txt", col.names = "x")

########### A ###########
# parameter init values
mu0 <- 30
tau0 <- 1
v0 <- 0
sigma0 <- 10

n <- nrow(data)
dataMean <- mean(data$x)
vn <- v0 + n

iter <- 1000

# draw mu 
drawMu <- function(prevMu, prevSigma) {
  tauSq <- 1/( (n/prevSigma) + (1/tau0^2) ) 
  w <- (n/prevSigma)/((n/prevSigma) + (1/tau0^2))
  mu <- w*dataMean + (1-w)*mu0
  draw <- rnorm(1, mu, sqrt(tauSq))
  return (draw)
}
#inv chi square
invChiSquare <- function(v, s) {
  return(v*s / rchisq(1,v))
}

# draw sigma 
drawSigma <- function(mu) {
  sum <- 0 
  for (i in 1:n) {
    sum <- sum + (data[i,1] - mu)^2
  }
  s <- (v0*sigma0 + sum)/(n+v0)
  return(invChiSquare(vn, s))
}

mu <- c()
sigma2 <- c()

currMu <- 32 
currSigma <- sigma0
for (i in 1:iter) {
  if(i %% 2 == 0) {
    currMu <- drawMu(currMu, currSigma)
  } else {
    currSigma <- drawSigma(currMu) 
  }
  mu <- c(mu, currMu)
  sigma2 <- c(sigma2, currSigma)
}

## plot trajectories of sampled mu and sigma
plot(mu, sqrt(sigma2),type='l')

# Also consider plotting the trajectories (the sampled values of mu and sigma2) over the iterations.
plot(mu, type = 'l')
plot(sqrt(sigma2), type = 'l')

########## C #########
# Your plot is very small and hard to see. 
# The green density looks incorrect, maybe just using the wrong grid.
# Your mixture density is plotted using the default output from Mattias' 
# code which computes the mean of the mixture densities computed at each Gibbs iteration. 
# This is ok, but not exactly the same as the density requested in the problem, 
# which is based on the posterior mean of the parameters.

densityData = density(data$x)

xGrid = seq(min(densityData$x),max(densityData$x),length = length(densityData$x))
ndens = dnorm(xGrid, mean(mu), mean(sqrt(sigma2)))

hist(data$x, 50, freq = FALSE)

lines(xGrid,
      ndens, 
      col='blue')

lines(mixDens, 
      col = 'green')

legend("topright", 
       box.lty = 1, 
       legend = c("Data",'Normal model', 'Mixed models'), 
       col = c("black",'blue', 'green'), 
       lwd = 2)

