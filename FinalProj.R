
pf <- function(x, beta){
  out <- 1/(1 + exp(-t(x)%*%beta))
  return(out)
}


#need a matrix with first column as predictors, and second column as response as input
beta_ls <- function(x, y, beta){
  sum <- 0
  for(i in 1:nrow(x)){
    xi <- as.matrix(unlist(x[i,]), ncol = 1) #format it so dimensions work out later
    yi <- y[i]
    pi <- pf(xi,as.matrix(beta))
    
    res <- (-1*yi)*(log10(pi)) - (1 - yi)*log10(1-pi)
    sum <- sum + res
  }
  return(sum)

}





##generate random data 
#set.seed(123)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
X <- cbind(1,x1, x2, x3) #1 for the intercept term
beta0 <- -1
beta1 <- 2
beta2 <- -0.5
beta3 <- 1.5
true_beta <- matrix(c(beta0, beta1, beta2, beta3), ncol = 1)
p <- 1/(1 + exp(-(X%*%true_beta)))
y <- rbinom(100,1,p)
data <- data.frame(X, y)

B_initial <- solve(t(X)%*%X)%*%t(X)%*%y

test <- beta_ls(X,y,B_initial)
result <- optim(par = B_initial, fn = beta_ls, x = X, y = y)

#provides the coefficient estimates
result$par[1,] #betazero
result$par[2,] #beta1
result$par[3,] #beta2
result$par[4,] #beta3



## OUTLINE for bootstrap confidence intervals -- will need to be updated once logistic regression is complete



bootstrap_conf_intervals <- function(X, y, alpha = 0.05, n_bootstraps = 20) {
  n <- nrow(X)
  beta_bootstraps <- matrix(NA, ncol = n_bootstraps, nrow = length(name_log_regressionfxn(X, y)))
  
  for (i in 1:n_bootstraps) {
    # Sample with replacement
    indices <- sample(1:n, replace = TRUE)
    X_bootstrap <- X[indices, ]
    y_bootstrap <- y[indices]
    
    # Run logistic regression on the bootstrap sample
    beta_bootstraps[, i] <- namelogregressionfxn(X_bootstrap, y_bootstrap)
  }
  
  # Calculate confidence intervals
  lower <- apply(beta_bootstraps, 1, function(col) quantile(col, alpha / 2))
  upper <- apply(beta_bootstraps, 1, function(col) quantile(col, 1 - alpha / 2))
  
  intervals <- data.frame(lower = lower, upper = upper)
  return(intervals)
}
