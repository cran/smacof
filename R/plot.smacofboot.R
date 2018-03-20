# plot method for bootstrap smacof

plot.smacofboot <- function(x, plot.dim = c(1,2), col = 1, label.conf = list(label = TRUE, pos = 3, cex = 0.8), 
                            ell = list(lty = 1, lwd = 1), main, xlab, ylab, xlim, ylim, asp = 1, type = "p", pch = 20, ...)
                            
                            
{
# x ... object of class smacofboot
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  if (missing(main)) main <- paste("MDS Bootstrap Plot") else main <- main
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  
  if (is.null(ell$lty)) ell$lty <- 1
  if (is.null(ell$lwd)) ell$lwd <- 1
  
  n <- x$nobj
  aus <- list()
  for (i in 1:n) aus[[i]] <- ellipse(x$cov[[i]][plot.dim, plot.dim], level = 1-x$alpha, type = 'l', centre = x$conf[i,plot.dim])
  limcoor <- rbind(do.call(rbind, aus), x$conf[,plot.dim])
  
  if (missing(xlim)) xlim <- range(limcoor[,1])*1.1
  if (missing(ylim)) ylim <- range(limcoor[,2])*1.1
  
  
  plot(x$conf[,x1], x$conf[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
  if (label.conf$label) {
    if (label.conf$pos == 5) {
      thigmophobe.labels(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), 
                         cex = label.conf$cex, text.pos = NULL, 
                         col = col)  
    } else {
      text(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), 
           cex = label.conf$cex, pos = label.conf$pos, 
           col = col)
      
    }
  }
  
  if (length(col) != n) col <- rep(col[1], n)
  for (i in 1:n) lines(aus[[i]][,1], aus[[i]][,2], lty = ell$lty, lwd = ell$lwd, col = col[i]) 
  names(aus) <- rownames(x$conf)
  invisible(aus)
}


  


  

 
