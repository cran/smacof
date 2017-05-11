# plot method for all smacof objects

plot.smacof <- function(x, plot.type = "confplot", plot.dim = c(1,2), sphere = TRUE, bubscale = 1, 
                        col = 1, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), 
                        hull.conf = list(hull = FALSE, col = 1, lwd = 1, ind = NULL),
                        shepard.x = NULL, identify = FALSE, type = "p", pch = 20, cex = 0.5, asp = 1, 
                        main, xlab, ylab, xlim, ylim, col.hist = NULL, ...)

# x ... object of class smacof
# plot.type ... types available: "confplot", "Shepard", "resplot", "histogram"
# sphere ... if TRUE, sphere is drawn for spherical smacof
  
{
  ## --- check type args:
  plot.type <- match.arg(plot.type, c("confplot", "Shepard", "resplot","bubbleplot", "stressplot", "histogram"), several.ok = FALSE)
  
  ## --- check label lists
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$col)) label.conf$col <- 1
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  if (identify) label.conf$label <- FALSE
  
  ## --- check hull list
  if (is.null(hull.conf$hull)) hull.conf$hull <- FALSE
  if (is.null(hull.conf$col)) hull.conf$col <- 1
  if (is.null(hull.conf$lwd)) hull.conf$lwd <- 1
  if (is.null(hull.conf$ind)) hull.conf$ind <- NULL
  if (hull.conf$hull && is.null(hull.conf$ind)) stop("Index vector for hulls needs to be specified!")
  
  #--------------- utility function for circle drawing -----------------
  circle <- function(x, y, r, ...) {
    ang <- seq(0, 2*pi, length = 100)
    xx <- x + r * cos(ang)
    yy <- y + r * sin(ang)
    polygon(xx, yy, ...)
  }
  #------------ end utility functions ----------------

  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  #----------------- configuration plot ---------------------
  if (plot.type == "confplot") {
    
    
    #if (type == "n") label.conf$pos <- NULL
      
    if (missing(main)) main <- paste("Configuration Plot") else main <- main
    if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$conf[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$conf[,y1])*1.1
    
    #if (missing(type)) type <- "n" else type <- type
    #if (identify) type <- "p"
    
    plot(x$conf[,x1], x$conf[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
         xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
    if (label.conf$label) {
      if (label.conf$pos == 5) {
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
    
    if (hull.conf$hull) {
      ind <- hull.conf$ind
      n <- dim(x$conf)[1]
      M <- as.data.frame(x$conf)
      XX <- cbind(M, ind) 
      X.sort <- XX[order(ind), ]
      xx <- yy <- NULL
      k <- 0
      for (i in 1:n) { 
        v <- X.sort$ind[i+1]
        if (i==n) v="$"
        if (X.sort$ind[i] == v ) { 
          k<-k+1 
        } else { 
          von <- i-k 
          xx <- X.sort[von:i, 1] 
          yy <- X.sort[von:i, 2]
          hpts <- chull(x = xx, y = yy)
          hpts <- c(hpts, hpts[1])
          lines(xx[hpts], yy[hpts], col = hull.conf$col, lwd = hull.conf$lwd) 
          k<-0 
        } 
      } 
    }
    
  }

  #---------------- Shepard diagram ------------------
  if (plot.type == "Shepard") {
    if (missing(main)) main <- paste("Shepard Diagram") else main <- main
    if (missing(xlab)) {
      if (is.null(shepard.x)) xlab <- "Dissimilarities" else xlab <- "Proximities"
    } else xlab <- xlab
    
    if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab
    
    notmiss <- as.vector(x$weightmat > 0)
    if (is.null(shepard.x)) {
      xcoor <- as.vector(x$delta)
    } else {
      xcoor <- as.vector(as.dist(shepard.x))
    }  
    
    if (missing(xlim)) xlim <- range(xcoor[notmiss], na.rm = TRUE)
    if (missing(ylim)) ylim <- range(as.vector(x$confdist)[notmiss])
    
    plot(xcoor[notmiss], as.vector(x$confdist)[notmiss], main = main, type = "p", pch = pch, cex = cex,
         xlab = xlab, ylab = ylab, col = "darkgray", xlim = xlim, ylim = ylim, ...)
    notmiss.iord <- notmiss[x$iord]
    points((xcoor[x$iord])[notmiss.iord], (as.vector(x$dhat[x$iord]))[notmiss.iord], type = "b", pch = pch, cex = cex)
    
    #points(as.vector(x$delta[x$iord]),as.vector(x$dhat[x$iord]), type = "b", pch = 20, cex = .5)
 }

  #--------------- Residual plot --------------------
  if (plot.type == "resplot") {
    if (missing(main)) main <- paste("Residual plot") else main <- main
    if (missing(xlab)) xlab <- "Normalized Dissimilarities (d-hats)" else xlab <- xlab
    if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab
    #resmat <- residuals(x)
 
    if (missing(xlim)) xlim <- range(as.vector(x$dhat))
    if (missing(ylim)) ylim <- range(as.vector(x$confdist))
    
    plot(as.vector(x$dhat), as.vector(x$confdist), main = main, type = "p", col = "darkgray",
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,...)
    abline(0, 1)
  }

  #----------------------- Stress decomposition -----------------
  if (plot.type == "stressplot") {
    if (missing(main)) main <- paste("Stress Decomposition Chart") else main <- main
    if (missing(xlab)) xlab <- "Objects" else xlab <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab
 
    spp.perc <- sort(x$spp, decreasing = TRUE)
    xaxlab <- names(spp.perc)

    if (missing(xlim)) xlim1 <- c(1,length(spp.perc)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- (range(spp.perc)*c(0.9, 1.1)) else ylim1 <- ylim
    
    op <- par(mar = c(3, 4, 4, 2))
    plot(1:length(spp.perc), spp.perc, xaxt = "n", type = "p", xlab = " ", ylab = ylab, main = main, xlim = xlim1, ylim = ylim1, ...)
    mtext(xlab, side = 1, padj = 2)
    text(1:length(spp.perc), spp.perc, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(spp.perc)) lines(c(i,i), c(spp.perc[i],0), col = "lightgray", lty = 2)                            
    par(op)
  }

  #------------------------------ bubble plot -------------------------
  if (plot.type == "bubbleplot")
  {

    if (missing(main)) main <- paste("Bubble Plot") else main <- main
    if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$conf[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$conf[,y1])*1.1
    
    spp.perc <- x$spp
    bubsize <- spp.perc/length(spp.perc)*(bubscale + 3)
    
    
    plot(x$conf, cex = bubsize, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, asp = asp, ...)
    xylabels <- x$conf
    ysigns <- sign(x$conf[,y1])
    xylabels[,2] <- (abs(x$conf[,y1])-(x$conf[,y1]*(bubsize/50)))*ysigns 
    text(xylabels, rownames(x$conf), pos = 3,cex = 0.7)
  }  
  
  #------------------------------ histogram plot -------------------------
  if (plot.type == "histogram")
  {
    if (missing(main)) main <- paste("Weighted histogram") else main <- main
    if (missing(xlab)) xlab <- paste("Dissimilarity") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Frequency") else ylab <- ylab

    #if (missing(xlim)) xlim <- range(breaks)
    
    wtd.hist(x$delta, weight = x$weightmat, main = main, xlab = xlab, ylab = ylab, col = col.hist, ...) 
  }

}
