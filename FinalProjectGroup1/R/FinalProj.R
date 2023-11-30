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


# @description Function models the relationship between the independent variables and the probability of a particular outcome. The function then list the predicted, actual, and initial values.
#' @param resp y \code{vector} of length n.
#' @param pred X \code{matrix} n number rows beta plus 1 number of columns with first column being entirely ones.
#' @return A \list{numeric} giving beta coefficients estimates, a vector of predicted y-values using predicted estimates, and original y-values and x-matrix.
#' @author Elena Gagliano
#' @author Helen Wu
#' @author Max Van Horn
#' @examples
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' x3 <- rnorm(100)
#' X <- cbind(1,x1, x2, x3) #1 for the intercept term
#' beta0 <- -1
#' beta1 <- 2
#' beta2 <- -0.5
#' beta3 <- 1.5
#' true_beta <- matrix(c(beta0, beta1, beta2, beta3), ncol = 1)
#' p <- 1/(1 + exp(-(X%*%true_beta)))
#' y <- rbinom(100,1,p)
#' B_initial <- solve(t(X)%*%X)%*%t(X)%*%y

log_reg <- function(X, y){
  colnames(X)[1] <- "Intercept"
  B_initial <- solve(t(X)%*%X)%*%t(X)%*%y
  result <- optim(par = B_initial, fn = beta_ls, x = X, y = y)
  parameters <- result$par
  y_pred <- X%*%parameters
  out <- list("betas" = parameters, "y_pred" = y_pred, "y_actual" = y, "X" = X)
  class(out) <- "my_b"

  return(out)

}


#' @title Logistic Regression Plots
#'
#' @description This function uses logistic regression computations to plot a fitted logistic curve in the data
#' @param output log_reg function
#' @return A graph{numeric} plotting predicted y-values against actual y-values and a logistic regression curve fitted using transformed x values and beta coefficient estimates
#' @author Elena Gagliano
#' @author Helen Wu
#' @author Max Van Horn
#' @examples
#' set.seed(123)
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' x3 <- rnorm(100)
#' X <- cbind(1,x1, x2, x3) #1 for the intercept term
#' beta0 <- -1
#' beta1 <- 2
#' beta2 <- -0.5
#' beta3 <- 1.5
#' true_beta <- matrix(c(beta0, beta1, beta2, beta3), ncol = 1)
#' p <- 1/(1 + exp(-(X%*%true_beta)))
#' y <- rbinom(100,1,p)
#' B_initial <- solve(t(X)%*%X)%*%t(X)%*%y
plot.my_b <- function(obj){
  y_pred <- obj$y_pred
  y_actual <- obj$y_actual
  X <- obj$X
  parameters <- obj$parameters

  plot(y_actual ~ y_pred ,xlim = range(y_pred), ylim = range(-.05, 1.5))

  x_line <- seq(min(y_pred), max(y_pred), length = length(y_pred))
  y_line <- pf(t(X), obj$betas)
  lines(x_line, sort(y_line), col = "red", lty = 1, lwd = 2)
}

test <- log_reg(X, y)
plot(test)


## OUTLINE for bootstrap confidence intervals -- will need to be updated once logistic regression is complete

#' @title Bootstrap
#'
#' @description This function resamples the data and lets the user choose the significance level α to obtain for the 1−α confidence intervals for β, and the number of bootstraps which by default is 20.
#' @param resp A \code{vector} of dimension n.
#' @param pred A \code{matrix} containing predictors.
#' @param beta A \code{vector} containing coefficients.
#' @param norm A \code{character} defining loss to use (default is `L2`).
#' @return A \code{numeric} giving value of loss at \code{beta}
#' @author Elena Gagliano
#' @author Helen Wu
#' @author Max Van Horn
#' @importFrom stats #put package and then fxn from package you rely on
#' @export
#' @examples
bootstrap_conf_intervals <- function(X, y, alpha = 0.05, n_bootstraps = 20) {
  n <- nrow(X)
  beta_bootstraps <- matrix(NA, ncol = n_bootstraps, nrow = length(log_reg(X, y)$betas))

  for (i in 1:n_bootstraps) {
    # Sample with replacement
    indices <- sample(1:n, replace = TRUE)
    X_bootstrap <- X[indices, ]
    y_bootstrap <- y[indices]

    # Run logistic regression on the bootstrap sample
    beta_bootstraps[, i] <- log_reg(X_bootstrap, y_bootstrap)$betas
  }

  # Calculate confidence intervals
  lower <- apply(beta_bootstraps, 1, function(row) quantile(row, alpha / 2))
  upper <- apply(beta_bootstraps, 1, function(row) quantile(row, 1 - alpha / 2))

  intervals <- data.frame(lower = lower, upper = upper)
  return(intervals)
}


# Bootstrap confidence intervals
intervals <- bootstrap_conf_intervals(X, y)

# Print the resulting intervals
print(intervals)


##confusion matrix function
#' @title Confusion Matrix
#'
#' @description This function creates the confusion matrix that measures the accuracy of the function
#' @param resp A \code{vector} of dimension n.
#' @param pred A \code{matrix} containing predictors.
#' @param beta A \code{vector} containing coefficients.
#' @param norm A \code{character} defining loss to use (default is `L2`).
#' @return A \code{numeric} giving value of loss at \code{beta}
#' @author Elena Gagliano
#' @author Helen Wu
#' @author Max Van Horn
#' @importFrom stats #put package and then fxn from package you rely on
#' @export
#' @examples
library(caret)
y_pred <- log_reg(X,y)$y_pred
confmat <- function(y_actual, y_pred, cutoff = 0.5){
  #convert predicted values to 0 or 1 based on threshold
  y_bin <- ifelse(y_pred >= cutoff, 1,0)
  cmat <- confusionMatrix(factor(y_actual), factor(y_bin))
  metrics <- cmat$byClass
  prevalence <- as.numeric(metrics["Prevalence"])
  accuracy <- as.numeric(cmat$overall["Accuracy"])
  sensitivity <- as.numeric(metrics["Sensitivity"])
  specificity <- as.numeric(metrics["Specificity"])
  true_p <- cmat$byClass["Pos Pred Value"] * cmat$byClass["Precision"]
  false_p <- cmat$byClass["Neg Pred Value"] * cmat$byClass["Precision"]
  false_discovery_rate <- as.numeric(false_p / (false_p + true_p))
  diagnostic_odds_ratio <- as.numeric(metrics["Sensitivity"]/(1 - metrics["Specificity"]))

  out <- invisible(list("Prevalence" = prevalence, "Accuracy" = accuracy, "Sensitivity" = sensitivity, "Specificity" = specificity, "False Discovery Rate" = false_discovery_rate, "Diagnostic Odds Ratio" = diagnostic_odds_ratio))
  return(out)

}

confmat(y, y_pred)


#plot metrics cut off grid


# Function to plot metrics over a grid of cutoff values
#' @title Confusion Matrix Plots
#'
#' @description This function plots a grid that includes all of the metrics measured by the confusion matrix
#' @param resp A \code{vector} of dimension n.
#' @param pred A \code{matrix} containing predictors.
#' @param beta A \code{vector} containing coefficients.
#' @param norm A \code{character} defining loss to use (default is `L2`).
#' @return A \code{numeric} giving value of loss at \code{beta}
#' @author Elena Gagliano
#' @author Helen Wu
#' @author Max Van Horn
#' @importFrom stats #put package and then fxn from package you rely on
#' @export
#' @examples
plot_metrics <- function(X, y, cutoff_values = seq(0.1, 0.9, by = 0.1)) {
  result <- log_reg(X, y)$betas
  predicted_probs <- pf(t(X), result)
  metrics_matrix <- matrix(NA, nrow = length(cutoff_values), ncol = 7,
                           dimnames = list(NULL, c("Cutoff", "Prevalence", "Accuracy", "Sensitivity", "Specificity", "False Discovery Rate", "Diagnostic Odds Ratio")))

  # Calculate metrics for each cutoff
  for (i in seq_along(cutoff_values)) {
    metrics <- confmat(y, predicted_probs, cutoff_values[i])
    metrics_matrix[i, ] <- c(cutoff_values[i], metrics[[1]], metrics[[2]], metrics[[3]], metrics[[4]], metrics[[5]], metrics[[6]])
  }

  # Plot metrics
  par(mfrow = c(3, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  metrics_names <- c("Prevalence", "Accuracy", "Sensitivity", "Specificity", "False Discovery Rate", "Diagnostic Odds Ratio")
  for (i in seq_along(metrics_names)) {
    plot(metrics_matrix[, "Cutoff"], metrics_matrix[, i + 1], type = "l", col = i + 1,
         xlab = "Cutoff", ylab = metrics_names[i], main = paste("Metric vs. Cutoff"))

  }
}



