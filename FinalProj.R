
pf <- function(x, beta){
  out <- 1/(1 + exp(-t(x)%*%beta))
  return(out)
}





#need a matrix with first column as predictors, and second column as response as input
beta_ls <- function(matrix, beta){
  sum <- 0
  for(i in 1:nrow(matrix)){
    xi <- as.matrix(unlist(matrix[i,1:4]), ncol = 1) #format it so dimensions work out later
    yi <- matrix[i,5]
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

test <- beta_ls(data,B_initial)
result <- optim(par = B_initial, fn = beta_ls, matrix = data)
