plot.mdsbi <- function(x, vecscale = NULL, plot.dim = c(1,2), sphere = TRUE, col = 1, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), 
                       vec.conf = list(col = 1, cex = 0.8, length = 0.1), 
                       identify = FALSE, type = "p", pch = 20, 
                       asp = 1, main, xlab, ylab, xlim, ylim, ...) {
  
  ## argument checks
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$col)) label.conf$col <- 1
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  if (identify) label.conf$label <- FALSE
  
  if (is.null(vec.conf$length)) vec.conf$length <- 0.1
  if (is.null(vec.conf$col)) vec.conf$col <- 1
  if (is.null(vec.conf$cex)) vec.conf$cex <- 0.8
  
  if (missing(main)) main <- paste("Configuration Plot") else main <- main
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  
  circle <- function(x, y, r, ...) {
    ang <- seq(0, 2*pi, length = 100)
    xx <- x + r * cos(ang)
    yy <- y + r * sin(ang)
    polygon(xx, yy, ...)
  }
  
  x$conf <- x$model$X           ## configurations
  
  if (is.null(vecscale) && ncol(x$coef) == 1) vecscale <- 1
  if (is.null(vecscale)) vecscale <- vecscale(x$coef)  ## scale
  x$coef <- vecscale * x$coef
  
  temp <- rbind(x$conf, t(x$coef))
  if (missing(xlim)) xlim <- range(temp[,x1])*1.1
  if (missing(ylim)) ylim <- range(temp[,y1])*1.1
  
  ## configuration plot
  plot(x$conf[,x1], x$conf[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
  if (label.conf$label) {
    if (label.conf$pos[1] == 5) {
      thigmophobe.labels(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), 
                         cex = label.conf$cex, text.pos = NULL, 
                         col = label.conf$col)  
    } else {
      text(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), 
           cex = label.conf$cex, pos = label.conf$pos, 
           col = label.conf$col)
      
    }
  }
 
  if ((any(class(x) == "smacofSP")) && (sphere)) {
    radius <- sqrt(x$conf[2,1]^2 + x$conf[2,2]^2)
    circle(0, 0, radius, lty = 2, border = "lightgray")
  }
  if (identify) {
    identify(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), cex = label.conf$cex, pos = label.conf$pos, 
             col = label.conf$col)
  }
  
  ## add vectors
  abline(h = 0, v = 0, lty = 2, col = "darkgray")
  
  n <- ncol(x$coef)
  xycoor <- t(x$coef[c(x1, y1), ])
  posvec <- apply(xycoor, 1, sign)[2,] + 2     
  for (i in 1:n) {
    arrows(0, 0, x$coef[x1, i], x$coef[y1, i], length = vec.conf$length, col = vec.conf$col)
    text(x$coef[x1, i], x$coef[y1, i], labels = colnames(x$coef)[i], col = vec.conf$col, pos = posvec[i], cex = vec.conf$cex)
  }
}

## the code below is copied from https://raw.githubusercontent.com/friendly/candisc/refs/heads/master/R/vecscale.R
## (C) by Michael Friendly; LICENSE: GPL-2/3
## curl -sL https://raw.githubusercontent.com/friendly/candisc/refs/heads/master/R/vecscale.R | md5sum | cut -d ' ' -f 1
## -> db6ac7dcd414e1bee9dcd088e6281da4

#' Scale vectors to fill the current plot
#' 
#' Calculates a scale factor so that a collection of vectors nearly fills the
#' current plot, that is, the longest vector does not extend beyond the plot
#' region.
#' 
#' 
#' @param vectors a two-column matrix giving the end points of a collection of
#'        vectors
#' @param bbox the bounding box of the containing plot region within which the
#'        vectors are to be plotted## 
#' @param origin origin of the vectors
#' @param factor maximum length of the rescaled vectors relative to the maximum
#'        possible
#' @return scale factor, the multiplier of the vectors
#' @author Michael Friendly
#' @seealso \code{\link{vectors}}
#' @keywords manip
#' @examples
#' 
#' bbox <- matrix(c(-3, 3, -2, 2), 2, 2)
#' colnames(bbox) <- c("x","y")
#' rownames(bbox) <- c("min", "max")
#' bbox
#' 
#' vecs <- matrix( runif(10, -1, 1), 5, 2)
#' 
#' plot(bbox)
#' arrows(0, 0, vecs[,1], vecs[,2], angle=10, col="red")
#' (s <- vecscale(vecs))
#' arrows(0, 0, s*vecs[,1], s*vecs[,2], angle=10)
#' 
vecscale <- function(vectors, 
		bbox=matrix(par("usr"), 2, 2),
		origin=c(0, 0), factor=0.95) {	
	scale <- c(sapply(bbox[,1] - origin[1], function(dist) dist/vectors[,1]), 
			sapply(bbox[,2] - origin[2], function(dist) dist/vectors[,2])) 
	scale <- factor * min(scale[scale > 0])
	scale
}
