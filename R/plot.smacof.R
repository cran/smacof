# plot method for all smacof objects

plot.smacof <- function(x, plot.type = "confplot", plot.dim = c(1,2), sphere = TRUE, bubscale = 3, 
                        col = 1, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), identify = FALSE, 
                        type = "p", pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, ...)

# x ... object of class smacof
# plot.type ... types available: "confplot", "Shepard", "resplot"
# sphere ... if TRUE, sphere is drawn for spherical smacof
  
{
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
  
  if ((any(class(x) == "smacofSP")) && (x$alg == "dual")) {             #remove first column
     x$obsdiss <- as.dist(as.matrix(x$obsdiss1)[,-1][-1,])
     x$confdiss <- as.dist(as.matrix(x$confdiss)[,-1][-1,])
  }   
  
  if (type == "n") label.conf$pos <- NULL

  #----------------- configuration plot ---------------------
  if (plot.type == "confplot") {
    
    if (missing(main)) main <- paste("Configuration Plot") else main <- main
    if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$conf[,x1])
    if (missing(ylim)) ylim <- range(x$conf[,y1])
    
    #if (missing(type)) type <- "n" else type <- type
    #if (identify) type <- "p"
    
    plot(x$conf[,x1], x$conf[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
         xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
    if (label.conf[[1]]) text(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), 
                              cex = label.conf$cex, pos = label.conf$pos, 
                              col = label.conf$col)
      
    if ((any(class(x) == "smacofSP")) && (sphere)) {
      if (x$alg == "primal") {                     
        radius <- sqrt(x$conf[2,1]^2 + x$conf[2,2]^2)
        circle(0, 0, radius, lty = 2, border = "lightgray")
      }
    }
    
    if (identify) {
      identify(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), cex = label.conf$cex, pos = label.conf$cex, 
               col = label.conf$col)
    }
    
  }

  #---------------- Shepard diagram ------------------
  if (plot.type == "Shepard") {
    if (missing(main)) main <- paste("Shepard Diagram") else main <- main
    if (missing(xlab)) xlab <- "Dissimilarities" else xlab <- xlab
    if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab

    if (missing(xlim)) xlim <- range(as.vector(x$delta))
    if (missing(ylim)) ylim <- range(as.vector(x$confdiss))

    plot(as.vector(x$delta), as.vector(x$confdiss), main = main, type = "p", pch = 20, cex = .5,
         xlab = xlab, ylab = ylab, col = "darkgray", xlim = xlim, ylim = ylim, ...)

    points(as.vector(x$delta[x$iord]),as.vector(x$obsdiss[x$iord]), type = "b", pch = 20, cex = .5)
    
#     if (x$type == "ordinal") {
#       isofit <- isoreg(as.vector(x$delta), as.vector(x$confdiss))  #isotonic regression
#       points(sort(isofit$x), isofit$yf, type = "b", pch = 20)
#     } else {
#       regfit <- lsfit(as.vector(x$delta), as.vector(x$confdiss))   #linear regression
#       abline(regfit, lwd = 0.5)
#     }
  }

  #--------------- Residual plot --------------------
  if (plot.type == "resplot") {
    if (missing(main)) main <- paste("Residual plot") else main <- main
    if (missing(xlab)) xlab <- "Normalized Dissimilarities (d-hats)" else xlab <- xlab
    if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab
    #resmat <- residuals(x)

    if (missing(xlim)) xlim <- range(as.vector(x$obsdiss))
    if (missing(ylim)) ylim <- range(as.vector(x$confdiss))
    
    plot(as.vector(x$obsdiss), as.vector(x$confdiss), main = main, type = "p", col = "darkgray",
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,...)
    abline(0, 1)
  }

  #----------------------- Stress decomposition -----------------
  if (plot.type == "stressplot") {
    if (missing(main)) main <- paste("Stress Decomposition Chart") else main <- main
    if (missing(xlab)) xlab <- "Objects" else xlab <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab

    spp.perc <- sort((x$spp/sum(x$spp)*100), decreasing = TRUE)
    xaxlab <- names(spp.perc)

    if (missing(xlim)) xlim1 <- c(1,length(spp.perc)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(spp.perc) else ylim1 <- ylim
    
    plot(1:length(spp.perc), spp.perc, xaxt = "n", type = "p",
         xlab = xlab, ylab = ylab, main = main, xlim = xlim1, ylim = ylim1, ...)
    text(1:length(spp.perc), spp.perc, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(spp.perc)) lines(c(i,i), c(spp.perc[i],0), col = "lightgray", lty = 2)
                                  
  }

  #------------------------------ bubble plot -------------------------
  if (plot.type == "bubbleplot")
  {

    if (missing(main)) main <- paste("Bubble Plot") else main <- main
    if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$conf[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$conf[,y1])*1.1
    
    spp.perc <- x$spp/sum(x$spp)*100
    bubsize <- (max(spp.perc) - spp.perc + 1)/length(spp.perc)*bubscale
    
    plot(x$conf, cex = bubsize, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,...)
    xylabels <- x$conf
    ysigns <- sign(x$conf[,y1])
    xylabels[,2] <- (abs(x$conf[,y1])-(x$conf[,y1]*(bubsize/50)))*ysigns 
    text(xylabels, rownames(x$conf), pos = 1,cex = 0.7)
  }  

}
