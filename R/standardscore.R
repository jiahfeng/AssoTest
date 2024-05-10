#' score test with missing data
#'
#' @param x vector of covariates
#' @param y outcome of interest. Quantitative for family="gaussian". For family="binomial" should be either a factor with two levels
#' @param s a covariate of interets, allowed to be missing
#' @param family response type. See above
#' @export
completecase <- function(x , y , s, family = c("gaussian", "binomial")){
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
    A_prime <- -c(lg_2)
  } else {
    A <- y - x_IC%*%alphahat
    A_prime <- -1
  }

  s_impute <-drop(s)
  diff_alpha <- matrix(rep(0,ncol(x_IC)), ncol(x_IC), 1 )
  diff_beta <- t(s_impute) %*% A
  XandS <- cbind(x_IC, s_impute)
  Hessian <-t(XandS) %*%  (A_prime*XandS)
  U <- rbind(diff_alpha, diff_beta)
  SS <- t(U) %*% solve(-Hessian, U)

  pv <- pchisq(abs(SS), df=1, lower.tail = FALSE)

  estimated_variance <- mean(Hessian/n)



  object <- list(scorestatistic = SS, pv = pv, estimated_variance = estimated_variance)
  object
}
