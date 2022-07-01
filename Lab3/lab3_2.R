library(mvtnorm)

data <- read.table("eBayNumberOfBidderData.txt", header = TRUE)


########## A ##########

#The significance of a parameter is not governed by its estimated 
#absolute value, but by e.g. the p-value. Use e.g. summary(fit).

fit <- glm(nBids ~ 0 + ., data, family = poisson)

#Obtain p-values
p_values <- coef(summary(fit))[,4]
# Firstly, the p-value given for the Z-statistic would have to 
# be interpreted as how likely it is that a result as extreme or 
# more extreme than that observed would have occured under the 
# null hypothesis. I.e. 0.96 would in principle mean that the 
# data are providing very little evidence that the variable is 
# needed (while small values such as, say, 𝑝≤0.05
# would provide evidence for the likely relevance of the variable, 
# as pointed out by others already).

coeff <- t(fit$coefficients)
plot(abs(coeff), type='h', 
     lwd=2, 
     xlab = "coefficient index", 
     main='Significance of covariates',
     ylab='absolute value of coefficient')

X <- as.matrix(data[,2:10])
## The most significant covariate is minBidShare. 

########## B ##########

# Your estimated beta values should be close to the ones from glm() if 
# you did it right. Double check this line:

sigmaPrior <- 100 * solve(t(X)%*%X)

logPois <- function(beta, y, X, ...) {
  # log likelihood of poisson model
  if(!is.null(dim(beta))) {
    beta <- beta[1,]
  }
  
  lambda <- exp(X%*%beta)
  logLik <- sum(dpois(y, lambda, log=TRUE))
  
  if (logLik == -Inf) {
    logLik <- -2000
  }
  # logLik <- logLik + y[i] * t(beta)%*%x[i,] - exp( t(beta)%*%x[i,] - log(factorial(y[i])))
  # log of prior
  logPrior <- dmvnorm(beta, mean = rep(0, 9), sigma = sigmaPrior, log=TRUE)
  # add 
  return(logLik + logPrior)
}

OptimResults<-optim(coeff,logPois,gr=NULL, y = data$nBids,X = X, method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

postCov <- -solve(OptimResults$hessian)
st_div <- sqrt(diag(postCov))
betaMode <- OptimResults$par

########## C ##########
# Also plot the MCMC trajectories for the different betas. 
# Your results will be different when you fix 2b) also.

gaussianSample <- function(theta, sigma, c) {
  val <- rmvnorm(1, theta, c*sigma)
  return (val)
}

RWMSampler <- function(c, it, initBeta, fn, ...) {
  accRate <- 0
  sample <- c()
  prev <- gaussianSample(initBeta, postCov, c)
  for (i in 1:it) {
    candidate <- gaussianSample(prev, postCov, c)
    
    alpha <- min(1, exp( fn(candidate, ...) - fn(prev, ...) ))
    u <- runif(1, 0, 1)
    if (alpha >= u) {
      # accept candidate
      prev <- candidate
      accRate <- accRate + 1
      # as matrix
    }
    sample <- rbind(sample, prev)
  }
  return (sample)
}

# sample with random walk metropolis
burnin <- 500
it <- 10000
nIter <- burnin + it
sample <- RWMSampler(1,nIter, rep(0, 9), logPois, as.vector(data$nBids), X)
hist(sample[,9])
plot(sample[,1],
     sample[,2],
     type='l',
     xlab = expression(beta[1]),
     ylab = expression(beta[2]),
     main = expression("Samples of" ~ beta[1] ~ "and" ~ beta[2])
     )

# plot trajectoroies for betas 
for (i in 1:dim(sample)[2]) {
  file <- paste("lab3_2c_",i,".png", sep="")
  png(file, width = 1000, height = 650)
  plot(sample[,i], type = 'l', main=paste("Beta", i))
  abline(h=mean(sample[burnin+1:it,i]), col="red")
  dev.off()
}


