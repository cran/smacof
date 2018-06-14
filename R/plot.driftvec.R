plot.driftvec <- function(x, adjust = 1, main, xlim, ylim, xlab = "Dimension 1", ylab = "Dimension 2", pch = 20, asp = 1, 
                          col.conf = "black", col.drift = "lightgray", 
                          label.conf = list(label = TRUE, pos = 3, col = "black", cex = 0.8), ...) {
  
  if (missing(main)) main <- "Drift Vectors" else main <- main
  
  conf <- x$fitsym$conf
  X <- x$driftcoor
  coor <- X - conf    ## center
  angle <- atan2(coor[,2], coor[,1])
  length <- sqrt(coor[,1]^2 + coor[,2]^2)   ## length 
  length1 <- length/adjust                  ## stretch/squeeze
  
  x1 <- length1 * cos(angle) + conf[,1]     ## move back
  y1 <- length1 * sin(angle) + conf[,2]
  
  
  if (missing(xlim)) xlim <- range(c(x$fitsym$conf[,1], x1))*1.1
  if (missing(ylim)) ylim <- range(c(x$fitsym$conf[,2], y1))*1.1
  
  plot(x$fitsym$conf, main = main, asp = asp, col = col.conf, pch = pch, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
  if (label.conf$label) text(x$fitsym$conf[,1], x$fitsym$conf[,2], labels = rownames(x$fitsym$conf), 
                             cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
  arrows(conf[,1], conf[,2], x1, y1, length = 0.10, lwd = 0.5, col = col.drift)  
  
}

