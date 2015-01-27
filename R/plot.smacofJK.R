# plot method for Jackknife smacof

plot.smacofJK <- function(x, plot.dim = c(1,2), hclpar = list(start = 30, end = 200), plot.lines = TRUE, main, xlab, ylab, xlim, ylim, ...)
{
# x ... object of class smacofJK
  n <- x$nobj
  hclcolors <- rainbow_hcl(n, start = hclpar[[1]], end = hclpar[[2]])

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
  
  plot(x0[, 1:2], type = "n", xlab = xlab, ylab = ylab, main = main, col = hclcolors, xlim = xlim, ylim = ylim, ...)
  text(x0[, 1:2], rownames(x$smacof.conf), cex = 0.8, col = hclcolors)
  points(y0[, 1:2], col = hclcolors, cex = 0.6, pch = 16)
	
  if (plot.lines) {
	 for (i in 1:n) {
	   for (j in 1:n) {
	    #text (yy[j, 1, i], yy[j, 2, i], as.character(i), cex = 0.5, pos = ppos, col = hclcolors[j])
      points(yy[j, 1, i], yy[j, 2, i], cex = 0.2, col = hclcolors[j], pch = 20)
	    lines(matrix (c (y0[j, 1:2], yy[j, 1:2, i]), 2, 2, byrow=TRUE), col = hclcolors[j], lty = 1, lwd = 0.5)
	   }
	 }
  }
}


  


  

 
