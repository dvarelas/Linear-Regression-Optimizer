## lr optim 

rest_mle <- function(y,
                     x,
                     initialValues = NULL,
                     restrictionsL = NULL,
                     restrictionsU = NULL,
                     maxit = 5000,
                     trace = TRUE) {
  
  # loglikelihood of a linear regression model
  loglik <- function(y,x,beta) {
    b0 <- beta[1]
    beta <- beta[-1]
    n <- length(x)
    sigma <- sd(y-(b0 + x%*%beta))
    ll <- -(n/2)*log(2*pi)-(n/2)*log(sigma^2) - (1/(2*sigma^2))*sum((y-(b0 + x%*%beta))^2)
    return(-ll)
  }
  
  # define initial values
  if (is.null(initialValues)) {
    ls <- lm(y~x)
    betaInit <- as.numeric(ls$coefficients)
  } else {
    betaInit <- initialValues
  }
  
  ## define restrictions for the parameters
  # Lower bound
  if (is.null(restrictionsL)) {
    lower <- rep(-1e8, length(betaInit))
  } else {
    lower <- restrictionsL
  }
  
  # Upper bound
  if (is.null(restrictionsU)) {
    upper <- rep(1e8, length(betaInit))
  } else {
    upper <- restrictionsU
  }
  
  trace_level <- ifelse(trace == TRUE, 3, 0)
  
  x <- as.matrix(x)
  
  rmle <- optim(betaInit,
                loglik, x = x, y = y, 
                lower = lower, upper = upper, method = "L-BFGS-B",
                control = list(trace = trace_level, maxit = maxit))
  
  beta <- rmle$par
  
  return(beta)
}


