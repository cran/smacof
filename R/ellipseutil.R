## ------------- Utility functions for pseudo confidence ellipses ---------------------


## normalize delta such that SS off-diagonal elements = 1
normDeltaOff <- function(delta, w) {
  if (!is.list(delta)) {
    delta <- list(delta)
    w <- list(w)
  }
  delta <- lapply(delta, as.matrix)
  w <- lapply(w, as.matrix)
  
  m <- length(delta)  
  s <- 0
  for (k in 1:m) {
    s <- s + sum (w[[k]] * delta[[k]]^2)
  }
  s <- sqrt (4 / s)
  for (k in 1:m) {
    delta[[k]] <- s * delta[[k]]
  }
  return(delta)
}


## additional iterations X <-> T_k iterations for INDSCAL
smacofINDSCAL <- function (w, delta, b, x, itmax = 100, eps = 1e-6) {
  n <- nrow (x)
  p <- ncol (x)
  m <- length (w)
  itel <- 1
  sold <- Inf
  repeat {
    sa <- 0.0
    v <- matrix (0, n * p, n * p)
    u <- matrix (0, n * p, n * p)
    for (k in 1:m) {
      cmat <- tcrossprod (b[[k]])
      xmat <- x %*% b[[k]]
      dmat <- as.matrix (dist (xmat))
      vmat <- -w[[k]]
      diag(vmat) <- -rowSums(vmat)
      bmat <- -delta[[k]] * ifelse (dmat == 0, 0, 1 / dmat)
      diag(bmat) <- -rowSums(bmat)
      sa <- sa + sum (w[[k]] * (delta[[k]] - dmat) ^ 2) / 2.0
      v <- v + kronecker (cmat, vmat)
      u <- u + kronecker (cmat, bmat)
    }
    x <- matrix (ginv (v) %*% u %*% as.vector (x), n, p)
    sb <- 0.0
    for (k in 1:m) {
      xmat <- x %*% b[[k]]
      dmat <- as.matrix (dist (xmat))
      vmat <- -w[[k]]
      diag(vmat) <- -rowSums(vmat)
      hmat <- crossprod (x, vmat %*% x)
      bmat <- -delta[[k]] * ifelse (dmat == 0, 0, 1 / dmat)
      diag(bmat) <- -rowSums(bmat)
      gmat <- crossprod (x, bmat %*% x)
      sb <- sb + sum (w[[k]] * (delta[[k]] - dmat) ^ 2) / 2.0
      kmat <- ginv (hmat) %*% gmat
      b[[k]] <- diag (drop (kmat %*% diag (b[[k]])))
    }
   
    if ((itel == itmax) || (sold - sb < eps))
      break
    sold <- sb
    itel <- itel + 1
  }
  return (list (x = x, b = b, s = sb))
}



## --------------- additional idioscal iterations ------------------
smacofIDIOSCAL <- function (w, delta, b, x, itmax = 100, eps = 1e-6) {
    n <- nrow (x)
    p <- ncol (x)
    m <- length (w)
    itel <- 1
    sold <- Inf
    repeat {
      sa <- 0.0
      v <- matrix (0, n * p, n * p)
      u <- matrix (0, n * p, n * p)
      for (k in 1:m) {
        cmat <- tcrossprod (b[[k]])
        xmat <- x %*% b[[k]]
        dmat <- as.matrix (dist (xmat))
        vmat <- -w[[k]]
        diag(vmat) <- -rowSums(vmat)
        bmat <- -delta[[k]] * ifelse (dmat == 0, 0, 1 / dmat)
        diag(bmat) <- -rowSums(bmat)
        sa <-
          sa + sum (w[[k]] * (delta[[k]] - dmat) ^ 2) / 2.0
        v <- v + kronecker (cmat, vmat)
        u <- u + kronecker (cmat, bmat)
      }
      x <- matrix (ginv (v) %*% u %*% as.vector (x), n, p)
      sb <- 0.0
      for (k in 1:m) {
        xmat <- x %*% b[[k]]
        dmat <- as.matrix (dist (xmat))
        vmat <- -w[[k]]
        diag(vmat) <- -rowSums(vmat)
        hmat <- crossprod (x, vmat %*% x)
        bmat <- -delta[[k]] * ifelse (dmat == 0, 0, 1 / dmat)
        diag(bmat) <- -rowSums(bmat)
        gmat <- crossprod (x, bmat %*% x)
        sb <-
          sb + sum (w[[k]] * (delta[[k]] - dmat) ^ 2) / 2.0
        kmat <- ginv (hmat) %*% gmat
        for (j in 1:p) {
          b[[k]][, j] <- kmat %*% b[[k]][, j]
        }
      }
      if ((itel == itmax) || (sold - sb < eps))
        break
      sold <- sb
      itel <- itel + 1
    }
    return (list (x = x, b = b, s = sb))
  }



## --------------- ordinary MDS or identity restriction ------------
smacofNODIFF <- function (w, delta, x, itmax = 100, eps = 1e-6) {
  n <- nrow (x)
  p <- ncol (x)
  m <- length (w)
  itel <- 1
  sold <- Inf
  repeat {
    snew <- 0.0
    v <- matrix (0, n * p, n * p)
    u <- matrix (0, n * p, n * p)
    for (k in 1:m) {
      dmat <- as.matrix (dist (x))
      vmat <- -w[[k]]
      diag(vmat) <- -rowSums(vmat)
      bmat <- -delta[[k]] * ifelse (dmat == 0, 0, 1 / dmat)
      diag(bmat) <- -rowSums(bmat)
      snew <- snew + sum (w[[k]] * (delta[[k]] - dmat) ^ 2) / 2.0
      v <- v + kronecker (diag (p), vmat)
      u <- u + kronecker (diag (p), bmat)
    }
    x <- matrix (ginv (v) %*% u %*% as.vector (x), n, p)
    if ((itel == itmax) || (sold - snew < eps))
      break
    sold <- snew
    itel <- itel + 1
  }
  return (list (x = x, s = snew))
}




## ------------------------ derivatives --------------------------
smacofDerivativesX <- function (w, delta, b, x) {
  n <- nrow (x)
  p <- ncol (x)
  m <- length (w)
  s <- 0.0
  v <- matrix (0, n * p, n * p)
  u <- matrix (0, n * p, n * p)
  r <- matrix (0, n * p, n * p)
  for (k in 1:m) {
    cmat <- tcrossprod (b[[k]])
    xmat <- x %*% b[[k]]
    dmat <- as.matrix (dist (xmat))
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        ei <- ifelse (i == 1:n, 1, 0)
        ej <- ifelse (j == 1:n, 1, 0)
        amat <- outer (ei - ej, ei - ej)
        kmat <- kronecker (cmat, amat)
        kx <- drop (kmat %*% as.vector(x))
        kxkx <- outer (kx, kx)
        s <- s + w[[k]][i, j] * (delta[[k]][i, j] - dmat [i, j]) ^ 2
        v <- v + w[[k]][i, j] * kmat
        u <- u + w[[k]][i, j] * (delta[[k]][i, j] / dmat[i, j]) * kmat
        r <- r + w[[k]][i, j] * (delta[[k]][i, j] / (dmat [i, j] ^ 3)) * kxkx
      }
    }
  }
  g <- drop ((v - u) %*% as.vector(x))
  h <- v - u + r
  return (list(s = s, g = 0, h = h))
}







