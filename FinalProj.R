
pf <- function(x, beta){
  out <- 1/(1 + exp(-t(x)%*%beta))
  return(out)
}


#need a matrix with first column as predictors, and second column as response as input
beta_ls <- function(beta, x, y){
  sum <- 0
  for(i in 1:nrow(x)){
    xi <- as.matrix(unlist(x[i,]), ncol = 1) #format it so dimensions work out later
    yi <- y[i]
    pi <- pf(xi,as.matrix(beta))
   
    
    res <- (-1*yi)*(log(pi)) - (1 - yi)*log(1-pi)
    sum <- sum + res
  }
  return(sum)

}

log_reg <- function(X, y){
  B_initial <- solve(t(X)%*%X)%*%t(X)%*%y
  result <- optim(par = B_initial, fn = beta_ls, x = X, y = y)
  out <- list("int" = result$par[1,], "b1" = result$par[2,], "b2" = result$par[3,], "b3" = result$par[4,])
  class(out) <- "my_b"
  return(result$par)
}

plot.my_b <- function(obj){
  plot(obj$X ~ obj$y)
}


##generate random data 
set.seed(123)
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

B_initial <- solve(t(X)%*%X)%*%t(X)%*%y

test <- beta_ls(B_initial,X,y)
result <- optim(par = B_initial, fn = beta_ls, x = X, y = y)

log_reg(X,y)



## OUTLINE for bootstrap confidence intervals -- will need to be updated once logistic regression is complete


bootstrap_conf_intervals <- function(X, y, alpha = 0.05, n_bootstraps = 20) {
  n <- nrow(X)
  beta_bootstraps <- matrix(NA, ncol = n_bootstraps, nrow = length(log_reg(X, y)))
  
  for (i in 1:n_bootstraps) {
    # Sample with replacement
    indices <- sample(1:n, replace = TRUE)
    X_bootstrap <- X[indices, ]
    y_bootstrap <- y[indices]
    
    # Run logistic regression on the bootstrap sample
    beta_bootstraps[, i] <- log_reg(X_bootstrap, y_bootstrap)
  }
  
  # Calculate confidence intervals
  lower <- apply(beta_bootstraps, 1, function(row) quantile(row, alpha / 2))
  upper <- apply(beta_bootstraps, 1, function(row) quantile(row, 1 - alpha / 2))
  
  intervals <- data.frame(lower = lower, upper = upper)
  return(intervals)
}

# test bootstrap
set.seed(123)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
X <- cbind(1, x1, x2, x3) 
beta0 <- -1
beta1 <- 2
beta2 <- -0.5
beta3 <- 1.5
true_beta <- matrix(c(beta0, beta1, beta2, beta3), ncol = 1)
p <- 1 / (1 + exp(-(X %*% true_beta)))
y <- rbinom(100, 1, p)

# Bootstrap confidence intervals
intervals <- bootstrap_conf_intervals(X, y)

# Print the resulting intervals
print(intervals)


##logistic regression plot
##plot dots from before transformation of XB (or y tilda)
##plot line from after the transformation
y_pred <- X%*%B_initial
plot(y ~ y_pred ,xlim = range(y_pred), ylim = range(-.05, 1.5))
x_line <- seq(min(y_pred), max(y_pred), length = 100)
y_line <- pf(t(X), B_initial)
lines(x_line, sort(y_line), col = "red", lty = 1, lwd = 2)

##confusion matrix function
library(caret)
confmat <- function(y_actual, y_pred, cutoff = 0.5){
  #convert predicted values to 0 or 1 based on threshold
  y_bin <- ifelse(y_pred >= cutoff, 1,0)
  cmat <- confusionMatrix(factor(y), factor(y_bin))
  metrics <- cmat$byClass
  prevalence <- metrics["Prevalence"]
  accuracy <- cmat$overall["Accuracy"]
  sensitivity <- metrics["Sensitivity"]
  specificity <- metrics["Specificity"]
  true_p <- cmat$byClass["Pos Pred Value"] * cmat$byClass["Precision"]
  false_p <- cmat$byClass["Neg Pred Value"] * cmat$byClass["Precision"]
  false_discovery_rate <- as.numeric(false_p / (false_p + true_p))
  formatted_fdr <- sprintf("%.7f", false_discovery_rate)
  diagnostic_odds_ratio <- metrics["Diagnostic Odds Ratio"]
  
  out <- invisible(list(prevalence, accuracy, sensitivity, specificity, formatted_fdr, diagnostic_odds_ratio))
  return(out)

}

confmat(y, y_pred)
metrics <- cmat$byClass
metrics["Prevalence"]


