gravity <- function(X, lambda = 1)
{
  Y <- X
  Y[Y > 0] <- 1
  C <- t(Y) %*% Y                               ## co-occurrence matrix
  diag(C) <- 0
  rs <- rowSums(C)
  cs <- colSums(C)
  grav.disty <- sqrt((rs %*% t(cs))/C)          ## gravity equation
  rownames(grav.disty) <- colnames(grav.disty)
  W <- matrix(1, ncol = ncol(grav.disty), nrow = nrow(grav.disty))  ## weight matrix
  W[grav.disty == Inf] <- 0               
  grav.disty[grav.disty == Inf] <- NA  
  grav.disty <- as.dist(grav.disty^lambda)
  diag(C) <- NA                                 ## blank out diagonals (since they don't correspond to word frequencies due to binarization)
  
  result <- list(gravdiss = grav.disty, weightmat = as.dist(W), co.occ = C)
  result         
}