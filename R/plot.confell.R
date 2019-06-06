# plot method for confidence ellipses

plot.confell <- function(x, eps = 0.05, plot.dim = c(1,2), col = 1, label.conf = list(label = TRUE, pos = 3, cex = 0.8), 
                         ell = list(lty = 1, lwd = 1, col = 1), main, xlab, ylab, xlim, ylim, asp = 1, type = "p", pch = 20, ...)
                            
                            
{
# x ... object of class confeff
  s <- x1 <- plot.dim[1]
  t <- y1 <- plot.dim[2]
  
  if (missing(main)) main <- paste("Pseudo-Confidence Ellipses") else main <- main
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  
  if (missing(xlim)) xlim <- range(x$X[,x1])*1.1
  if (missing(ylim)) ylim <- range(x$X[,y1])*1.1
  
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  
  if (is.null(ell$lty)) ell$lty <- 1
  if (is.null(ell$lwd)) ell$lwd <- 1
  if (is.null(ell$col)) ell$col <- 1
  
  n <- nrow(x$X)
  
  ## -------
  plot(x$X[,x1], x$X[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
  if (label.conf$label) {
    if (label.conf$pos == 5) {
      thigmophobe.labels(x$X[,x1], x$X[,y1], labels = rownames(x$X), cex = label.conf$cex, text.pos = NULL, col = col)  
    } else {
      text(x$X[,x1], x$X[,y1], labels = rownames(x$X), cex = label.conf$cex, pos = label.conf$pos, col = col)
    }
  }
  
  ## add ellipsoids
  z <- cbind (sin (seq(-pi, pi, length = 100)), cos(seq(-pi, pi, length = 100))) * sqrt(2 * eps * x$s)
  for (i in 1:n) {
    ii <- c((s - 1) * n + i, (t - 1) * n + i)
    y <- x$X[i, c(s, t)]
    amat <- x$h[ii, ii]
    heig <- eigen(amat)
    zi <- z %*% diag (1 /sqrt(heig$values))
    zi <- tcrossprod(zi, heig$vectors)
    zi <- zi + matrix (y, 100, 2, byrow = TRUE)
    lines(x = zi[, 1], y = zi[, 2], lty = ell$lty, lwd = ell$lwd, col = ell$col)
  }
}



  


  

 
