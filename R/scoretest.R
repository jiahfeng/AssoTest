#' score test with missing data for continuous and binary outcomes
#'
#' @param x vector of covariates
#' @param y outcome of interest. Quantitative for family="gaussian". For family="binomial" should be either a factor with two levels
#' @param s a covariate of interets, allowed to be missing
#' @param w potentially high-dimensional vector of covariates, if null, mean model
#' @param r indicator of whether s is observed
#' @param family response type. See above
#' @importFrom  glmnet glmnet
#' @importFrom  stats BIC lm coef glm pnorm var logLik pchisq cor
#' @export
scoretest <- function(x , y , s , w = NULL , r, family = c("gaussian", "binomial")){
  family <- match.arg(family)
  y <- drop(y)
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  n <- nrowy
  dimx <- dim(x)
  nrowx <- as.integer(dimx[1])
  if (nrowy != nrowx)
    stop(paste("number of observations in y (", nrowy,
               ") not equal to the number of rows of x (", nrowx, ")", sep = ""))

  alphahat <- glm(y~x, family = family)$coefficient
  x_IC <- cbind(1,x)
  if (family=="binomial") {
    lg_1 <- exp(x_IC%*%alphahat)/(1+exp(x_IC%*%alphahat))
    lg_2 <- exp(x_IC%*%alphahat)/(1+exp(x_IC%*%alphahat))^2
    A <- y - lg_1
    A_prime <- -c(lg_2) * x_IC
  } else {
    A <- y-x_IC%*%alphahat
    A_prime <- -x_IC
  }

  s <-drop(s)
  dims <- dim(s)
  nrows <- ifelse(is.null(dims), length(s), dims[1])
  r <- drop(r)
  if (length(r) != nrows)
    stop(paste("length of r (", length(r),
               ") not equal to the length of s (", nrows, ")", sep = ""))
  sobs <- s[r==1]
  if (any(is.na(sobs)))
    stop("observed s should not contain NA value")

  if(is.null(w)) {
    wtilde <- matrix(1, n, 1)
  } else {
    wtilde <- cbind(1,w)
  }
  sfit <- lm(s[r==1]~0+wtilde[r==1,])
  gammahat <- sfit$coefficient
  residual <- s[r==1] - as.matrix(wtilde[r==1,]) %*% gammahat
  total <- sum((s[r==1] - mean(s[r==1]))^2)
  res <- sum(residual^2)
  rsq <- 1 - res/total
  logL <- logLik(sfit)[1]
  s[r==0] <- 9999 ###change NA value to 9999
  I_alpha <- solve(-t(A_prime) %*% x_IC)
  I_gamma <- solve(t(r*wtilde) %*% wtilde)
  I_alphabeta <- t(A_prime) %*% (r*s + (1-r)*(wtilde%*%gammahat))
  I_gammabeta <- t((1-r)*wtilde) %*% A

  er <- rep(mean(A[r==1]), n)
  er[r==0] <- mean(A[r==0])
  expand_SS <- A * (r*s+(1-r)*(wtilde%*%gammahat))
  SS <- sum(expand_SS) / sqrt(n)
  expand_score1 <- r * (A + wtilde%*%I_gamma%*%I_gammabeta) * (s - wtilde%*%gammahat)
  expand_score2 <- A * (wtilde%*%gammahat + x_IC%*%I_alpha%*%I_alphabeta)
  estimated_variance <- var(expand_score1 + expand_score2)
  pv <- pnorm(abs(SS), mean=0, sd=sqrt(estimated_variance), lower.tail=FALSE)*2

  object <- list(scorestatistic = SS, pv = pv, estimated_variance = estimated_variance, rsquared = rsq)
  object
}
