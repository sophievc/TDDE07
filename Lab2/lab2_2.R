womenData = read.table("WomenWork.txt", header=TRUE)
library(mvtnorm)

######## A ########
# Added a zero in the model formula so that R doesn’t add an extra intercept
# A . in the model formula means to add all other variables in the dataset as features
glmModel <- glm(Work ~ 0 + ., data = womenData, family = binomial)


######## B ########
#You have some error in the way you compute the credible interval. 
#PostDraws is never computed in the code you submitted. 
#Also, note that you can compute the credible interval based 
#on your normal approximation of the posterior analytically.

# Param setup
tau <- 10
nParams <- 8
mu <- as.vector(rep(0, nParams))
X <- as.matrix(womenData[,2:9])
y <- as.matrix(womenData[,1])
initVal <- as.vector(rep(0,nParams))
sigma <-  tau^2 * diag(nParams)
initBetas <- rep(0,nParams)


# Calculate the Log Posterior Logistic
logPost <- function(beta,y,X,mu,Sigma) {
  
  # Calculate predictions by multiplying data with betas
  pred <- X%*%beta
  
  # Log likelihood function, derived by taking the log of the likelihood function prod(e^pred*y/(1+e^pred))
  logLik <- sum(y * pred - log(1 + exp(pred)))
  
  # calculate log likelihood prior
  logPrior <- dmvnorm(beta, matrix(0,nParams,1), Sigma, log=TRUE);
  
  return (logPrior + logLik)
}

OptimResults<-optim(initBetas,logPost,gr=NULL,y,X,mu,sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

#compute hessian etc as in lecture 6 slide 4
postCov <- -solve(OptimResults$hessian)
st_div <- sqrt(diag(postCov))
betaMode <- OptimResults$par
nSmallChild <- betaMode[7]
postDist <- dnorm(x = seq(0,10,0.01), mean = betaMode[7], sd = st_div[7])

#get the posterior distribution of NSmallChild
vec <- seq(0,10, 0.01)
postDraws = rmvnorm(1000, mean = betaMode, postCov)
int = qnorm(c(0.025, 0.975), betaMode[7], st_div[7])
#compute a grid of beta_values around the posterior mode
beta_grid = seq(betaMode[7] - 3*st_div[7], betaMode[7] + 3 * st_div[7], length=1000)
#compute the density for each beta-value
dn = dnorm(x=beta_grid, mean=betaMode[7], st_div[7])
#plot the parameter for NSmall child
plot(beta_grid, dn, 
     type='l', 
     col='darkgrey', 
     lwd='2', 
     xlab=expression(beta[7]), 
     ylab='Density')
abline(v=int[1], col='coral')
abline(v=int[2], col='coral')

########## C ###########

# Good! Note that you were also asked to plot the 
# predictive distribution, but your computed probability
# p is a sufficient summary of the predictive distribution.

y <- c(constant = 1, 
       husbandInc = 10,
       educYears = 8, 
       expYears = 10,
       expYears2 = 1,
       age = 40, 
       nSmallChild = 1,
       nBigChild = 1)

n_draws <- 1000

predDist <- function(n, beta, sigma, y) {
  y_draws = c()
  for (i in 1:n) {
    #beta draw
    betaDraw = as.vector(rmvnorm(1, mean = beta, sigma = sigma))
    #probability
    p <- exp(y%*%betaDraw)/(1+exp(y%*%betaDraw))
    #prediction of y
    y_draw = rbinom(1, 1, p)
    y_draws = c(y_draws, y_draw)
  }
  return (y_draws)
}

draws = predDist(n_draws, betaMode, postCov, y)
workProb <- sum(draws == 1)
p <- workProb/n_draws

barplot(height = c(p, 1-p), 
        col = c('coral', 'lightcyan'),
        legend =c('working', 'not working'), 
        main = 'Probability distribution of predicted y')

