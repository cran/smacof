# plot method for Jackknife smacof

plot.smacofJK <- function(x, plot.dim = c(1,2), hclpar = list(c = 50, l = 70), col.p, col.l, 
                          plot.lines = TRUE, legend = FALSE, inset = c(-0.2, 0), cex.legend = 0.7,
                          main, xlab, ylab, xlim, ylim, asp = 1, ...)
{
# x ... object of class smacofJK
  n <- x$nobj
  if (missing(col.p)) {
    hclcolors <- rainbow_hcl(n, c = hclpar[[1]], l = hclpar[[2]])
    hclcolors1 <- rainbow_hcl(n, c = hclpar[[1]], l = hclpar[[2]]+20)
  } else {
    hclcolors <- col.p
    hclcolors1 <- rep(col.l, n)
  }
  if (legend) hclcolors1 <- hclcolors
  
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  if (missing(main)) main <- paste("SMACOF Jackknife Plot") else main <- main
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  
  #if (missing(type)) type <- "n" else type <- type
      
  x0 <- x$smacof.conf
  y0 <- x$comparison.conf
  yy <- x$jackknife.conf
  xcoor <- cbind(apply(yy, 3, function(xc) xc[,x1]), y0[,x1], x0[,x1]) 
  ycoor <- cbind(apply(yy, 3, function(yc) yc[,y1]), y0[,y1], x0[,y1])
  
  if (missing(xlim)) xlim <- range(xcoor)
  if (missing(ylim)) ylim <- range(ycoor)
  
  if (legend) op <- par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
  plot(x0[, 1:2], type = "n", xlab = xlab, ylab = ylab, main = main, col = hclcolors, xlim = xlim, ylim = ylim, asp = asp, ...)
  points(y0[, 1:2], col = hclcolors1, cex = 0.6, pch = 16)
	
  if (plot.lines) {
	 for (i in 1:n) {
	   for (j in 1:n) {
	    points(yy[j, 1, i], yy[j, 2, i], cex = 0.2, col = hclcolors1[j], pch = 20)
	    lines(matrix (c (y0[j, 1:2], yy[j, 1:2, i]), 2, 2, byrow=TRUE), col = hclcolors1[j], lty = 1, lwd = 0.5)
	   }
	 }
  }
  
  if (!legend) {
    text(x0[, 1:2], rownames(x$smacof.conf), cex = 0.8, col = hclcolors) 
  } else {
    legend("topright", legend = rownames(x$smacof.conf), pch = 16, col = hclcolors1, 
           title = "Objects", cex = cex.legend, inset = inset)
    par(op)
  }
}


  


  

 
