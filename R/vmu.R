## Vector unfolding model
## X ... fixed coordinates columns

vmu <- function (delta, ndim = 2, center = TRUE, scale = FALSE, col.coord = NULL) { 
  m <- dim(delta)[2]
  N <- dim(delta)[1]
  P1 <- t(scale(t(delta), center = center, scale = scale))     ## centering
  X <- col.coord
  if (is.null(X)) {
    S <- svd(P1)
    X <- (m-1)^(1/2) * S$v[,1:ndim]
    Y <- (m-1)^(-1/2) * S$u[,1:ndim] %*% diag(S$d)[1:ndim,1:ndim]
    var <- sum(S$d[1:ndim]^2)/sum(S$d^2)
  } else {
    if (nrow(X) != ncol(delta)) {stop("Number of rows in X must match number of columns in delta.")}
    if (ncol(X) != ndim) {stop("Number of columns in X must match ndim.")}
    Y <-  t(ginv((t(X) %*% X)) %*% t(X) %*% t(P1))
    Y <- (m-1)^(-1/2) * Y 
    var <- cor(c(X%*%t(Y)), c(t(P1)))^2
  }
    
  rownames(X) <- colnames(delta)
  rownames(Y) <- rownames(delta)
  colnames(X) <- colnames(Y) <- paste0("D", 1:ndim)
    
  result <- list(conf.row = Y, conf.col = X, VAF = var, ndim = ndim, call = match.call()) 
  class(result) <- "vmu"
  return(result)
}


