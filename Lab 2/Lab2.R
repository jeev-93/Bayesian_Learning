#Question 1a
library(mvtnorm)
library(readxl)


tempdf <- read_xlsx("Linkoping2022.xlsx") 
intercept <- rep(1,365)
data1 <- cbind(tempdf, "intercept"=intercept)
time <- c(1:365)/365
data1$time <- time
time2 <- time^2
time2 <- data1$time^2
data1 <- cbind(data1, "time2"=time2)
time_mtx <- as.matrix(data.frame(const,time,time2))

#Assign hyperparameters for the prior
mu0=matrix(c(0,100,-100))
omega0=diag(x=0.01, nrow=3, ncol=3)
nu0=1
sigmasq0=1

#Joint prior function
priorfunc = function(mu0,omega0,nu0,sigmasq0){
  set.seed(12345)
  for(i in 1:20){
    #using chi_sq to sample sigma^2
    chi_sample = rchisq(n=1, df=nu0)
    sigma2 = nu0*sigmasq0/chi_sample
    
    #using mvtnorm sample beta
    beta = rmvnorm(n=1, mean=mu0, sigma=sigmasq0*solve(omega0))
    
    #quadratic regression  
    quad_regre= beta[1]+beta[2]*data1$time+beta[3]*(data1$time^2)+rnorm(1,mean=0, sd=sqrt(sigma2))
    lines(x=data1$time, y=quad_regre,col="red",lwd=2)
  }
}

### Check the given hyperpara
plot(data1$time,data1$temp, col="blue", 
     main="Predicted Temperature with given hyperparameters", ylab="Temperature", 
     xlab="Time", type="l")
priorfunc( mu0, omega0, nu0, sigmasq0)

### change the hyperpara nu
plot(data1$time,data1$temp, col="blue",
     main="Predicted Temperature with given hyperparameters", ylab="Temperature", 
     xlab="Time", type="l")
priorfunc( mu0, omega0, nu0, sigmasq0=0.03)


# Change the hyperpara sigma
plot(data1$time,data1$temp, col="blue",
     main="Predicted Temperature with changed hyperparameters", ylab="Temperature", 
     xlab="Time", type="l")
priorfunc(mu0=matrix(c(-10,110,-105)), omega0, nu0, sigmasq0=0.03)


#Question 1b
y <- as.matrix(tempdf$temp)
n <- 365
beta_hat <- solve(t(time_mtx) %*% time_mtx) %*% t(time_mtx) %*% y

mu_n <- solve(t(time_mtx) %*% time_mtx + omega0) %*% (t(time_mtx) %*% time_mtx %*% beta_hat + omega0 %*% mu0)
omega_n <- t(time_mtx) %*% time_mtx + omega0
nu_n <- nu0 + n
sigmasqn <- (nu0 * sigmasq0 + (t(y) %*% y + t(mu0) %*% omega0 %*% mu0 - t(mu_n) %*% omega_n %*% mu_n))/nu_n
sigmasqn <- as.numeric(sigmasqn)

#simulate
beta <- c()
sigma_square <- c()

for(i in 1:100){
  X <- rchisq(1,nu_n)
  sigmasqtemp <- nu_n * sigmasqn / X
  beta_temp <- rmvnorm(1,mu_n,sigmasqtemp * solve(omega_n))
  beta <- rbind(beta,beta_temp)
  sigma_square <- c(sigma_square,sigmasqtemp)
}

#1.2.1
# plot the histogram for the first beta
hist(beta[,1],main = "Histogram of beta0",breaks = 20)
# plot the histogram for the second beta
hist(beta[,2],main = "Histogram of beta1",breaks = 20)
# plot the histogram for the third beta
hist(beta[,3],main = "Histogram of beta2",breaks = 20)
#plot the histogram of sigma2
hist(sigma_square, main = "Histogram of sigma", breaks = 20)

#1.2.2
beta <- as.matrix(beta)
vector_epsilon <- unlist(lapply(sqrt(sigma_square),rnorm,n=1,mean=0))
temp <- time_mtx %*% t(beta)
temp <- + t(vector_epsilon + t(temp))# This is used for adding epsilon.
temp_median <- apply(temp,1,median)
temp_lower <- apply(temp,1,quantile,prob = 0.05)
temp_upper <- apply(temp,1,quantile,prob = 0.95)
plot(time,tempdf$temp,col="blue",xlab = "time",ylab = "temp",main = "Posterior median of Regression function.")
lines(time,temp_upper,type = "l",col="orange",lwd=2)
lines(time,temp_median,type = "l",col="red",lwd=2)
lines(time,temp_lower,type = "l",col="green",lwd=2)

#Question 1c
x_t <- -beta[,2] /(2 * beta[,3])
x_t



#Assignment 2
#2.a

library("mvtnorm") # reads the mvtnorm package into R's memory. We can now use the necessary function dmvnorm.

### Prior and data inputs ###
Covs <- c(2:8) # Select which covariates/features to include
standardize <- TRUE # If TRUE, covariates/features are standardized to mean 0 and variance 1
tau <- 2# scaling factor for the prior of beta 

# Loading the women dataset
WomenData <- data.frame(read.csv("WomenAtWork.dat",sep = " ")) # read data from file
Nobs <- dim(WomenData)[1] # number of observations
y <- as.matrix(WomenData$Work) # y=1 if the quality of wine is above 5, otherwise y=0.

X <- as.matrix(WomenData[,Covs]);
Xnames <- colnames(X)
#if (standardize){
#  Index <- 2:(length(Covs)-1)
#  X[,Index] <- scale(X[,Index])
#}
Npar <- dim(X)[2]

# Setting up the prior
mu <- matrix(0,Npar,1) # Prior mean vector
Sigma <- tau^2*diag(Npar) # Prior covariance matrix

# Functions that returns the log posterior for the logistic regression.
# First input argument of this function must be the parameters we optimize on, 
# i.e. the regression coefficients beta.

LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}


# Select the initial values for beta
initVal <- matrix(0,Npar,1)

# The argument control is a list of options to the optimizer optim, where fnscale=-1 means that we minimize 
# the negative log posterior. Hence, we maximize the log posterior.  
OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Printing the results to the screen
PostCov <- solve(-OptimRes$hessian)
colnames(PostCov) <-Xnames
row.names(PostCov) <-Xnames
PostMode <- OptimRes$par[1:7]
names(PostMode) <- Xnames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(-solve(OptimRes$hessian))) # Computing approximate standard deviations.
names(approxPostStd) <- Xnames # Naming the coefficient by covariates
print('The posterior mode is:\n')
print(PostMode)
print("The covariance matrix is :\n")
print(PostCov)
print('The approximate posterior standard deviation is:\n')
print(approxPostStd)

#glm model evaluation for comparison purpose
glmModel <- glm(Work ~  0+., data = WomenData, family = binomial)

summary(glmModel)

#Calculating the credible interval of equal tail  95% for NSmallChild
womenDataPost <- rmvnorm(10000, mean = PostMode, sigma = PostCov)
NSmallChildPost <- womenDataPost[,6]
NSmallChildCI <- quantile(NSmallChildPost,c(0.025,0.975))
print("The credible interval of equal tail 95% for NsmallChild is :\n")
print(NSmallChildCI)



#2.b
postPredDistZero  <- function(sampleSize, mode, cov, xvalues) {
  sampleData <- rmvnorm(sampleSize, mean = mode, sigma = cov)
  yOfSample <- xvalues%*%t(sampleData)
  yProbZero <- 1-(exp(yOfSample)/(1+exp(yOfSample)))
  
  return(yProbZero)
  
}

xvals <- c(1,18,11,7,40,1,1)

probability <- postPredDistZero(10000, PostMode,PostCov, xvals)
plot(density(probability), type = 'l', lwd = 3, col = "blue",main = "Density distribution of Pr(Y=0|x)")


#2.c
#Additional function to find mode
FindMode <- function(listvalues) {
  uniqueValues <- unique(listvalues)
  modeValue <- uniqueValues[which.max(tabulate(match(listvalues, uniqueValues)))]
  return(modeValue)
}

#Modified function
postPredDistZeroBinomial  <- function(sampleSize, mode, cov, xvalues) {
  sampleData <- rmvnorm(sampleSize, mean = mode, sigma = cov)
  yOfSample <- xvalues%*%t(sampleData)
  yProbZero <- 1-(exp(yOfSample)/(1+exp(yOfSample)))
  yProbZero <- yProbZero[1,]
  #Finding Mode
  Mode <- FindMode(yProbZero)
  returnList <- list("Mode" = Mode ,"Probability" =yProbZero)
  return(returnList)
  
}

xvals <- c(1,18,11,7,40,1,1)

set.seed(1245)
ModeAndProbability <- postPredDistZeroBinomial(10000, PostMode,PostCov, xvals)
#setting the mode value of the Pr(Y=0|x) to be the probability of women do not work
set.seed(456)
BinomialDistribution <- rbinom(100000,13, ModeAndProbability$Mode)
hist(BinomialDistribution, col = "red",breaks = 13, main = "Density distribution of women who are not working")
