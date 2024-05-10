#' BIC model selection
#'
#' @param s a covariate of interets, allowed to be missing
#' @param x vector of covariates
#' @param w potentially high-dimensional vector of covariates
#' @param r indicator of whether s is observed
#' @param lambda A user supplied lambda sequence
#'
#' @importFrom  stats cor
#' @export
modelselect <- function(s, x, w, r, lambda = NULL){
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
  np <- dim(w)
  if (np[1] != nrows)
    stop(paste("number of rows of s(", nrows,
               ") not equal to the number of rows of w (", np[1], ")", sep = ""))
  dimw <- np[2]
  if (is.null(np) | (np[2] <= 1))
  stop("W should be a matrix with 2 or more columns")
  wobs <- w[r==1,]
  x <- cbind(1, x)
  xobs <- x[r==1,]
  gamma_x <- coef(lm(sobs ~ 0 + xobs))
  s_residual <- s - x%*%gamma_x
  rho <- cor(s_residual[r==1], wobs)

  all.BIC <- sapply(lambda, function(x){
    M <- which(abs(rho)>x)
    wtilde <- cbind(xobs, wobs[,M])
    slinear <- lm(sobs ~ 0 + wtilde)
    BIC_slinear <- BIC(slinear)
    return(list(BIC_slinear, M))
  })
  BIC_value <- unlist(all.BIC[1,])
  BIC_value[which(is.infinite(BIC_value))] <- NA
  M_value <- all.BIC[2,]
  index_lambda <- which.min(BIC_value)
  M_BIC <- unlist(M_value[index_lambda])
  while (length(M_BIC) > 40) {
    BIC_value[index_lambda] <- NA
    index_lambda <- which.min(BIC_value)
    M_BIC <- M_value[[index_lambda]]
  }

  wselect <- w[,M_BIC]
  wtilde <- cbind(x, wselect)
  sfit <- lm(s[r==1]~0+wtilde[r==1,])
  gammahat <- sfit$coefficient
  residual <- s[r==1] - as.matrix(wtilde[r==1,]) %*% gammahat
  total <- sum((s[r==1] - mean(s[r==1]))^2)
  res <- sum(residual^2)
  rsq <- 1 - res/total
  logL <- logLik(sfit)[1]
  s_impute <- wtilde%*%gammahat
  s_impute[r==1] <- s[r==1]


  object <- list(bic = BIC_value[index_lambda], lambda = lambda[index_lambda],  model = M_BIC, s_impute = s_impute)
  object
}


