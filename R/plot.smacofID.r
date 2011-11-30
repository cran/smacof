# plot method for all smacof objects

plot.smacofID <- function(x, plot.type = "confplot", plot.dim = c(1,2), bubscale = 5, 
                          label.conf = list(label = TRUE, pos = 1, col = 1), identify = FALSE, 
                          type, main, xlab, ylab, xlim, ylim, ...)

# x ... object of class smacofID
# plot.type ... types available: "confplot", "Shepard", "resplot"
# Shepard plot and resplot are performed over sum of distances
  
{
  
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  
  if (plot.type == "confplot") {
    if (missing(main)) main <- paste("Group Configurations") else main <- main
    if (missing(xlab)) xlab <- paste("Configurations D", x1,sep = "") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Configurations D", y1,sep = "") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$gspace[,x1])
    if (missing(ylim)) ylim <- range(x$gspace[,y1])
    
    if (missing(type)) type <- "n" else type <- type
    if (identify) type <- "p"
    
    ppos <- label.conf[[2]]
    if (type == "n") ppos <- NULL

    plot(x$gspace[,x1], x$gspace[,y1], main = main, type = type, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
    
    if (!identify) {
      if (label.conf[[1]]) text(x$gspace[,x1], x$gspace[,y1], labels = rownames(x$gspace), cex = 0.8, pos = ppos, col = label.conf[[3]])    
    } else {
      identify(x$gspace[,x1], x$gspace[,y1], labels = rownames(x$gspace), cex = 0.8)
    }
  }

  #---------------- Shepard diagram ------------------
   if (plot.type == "Shepard") {
    if (missing(main)) main <- paste("Shepard Diagram") else main <- main
    if (missing(xlab)) xlab <- "Aggregated Observed Dissimilarities" else xlab <- xlab
    if (missing(ylab)) ylab <- "Aggregated Configuration Distances" else ylab <- ylab

    delta <- sumList(x$delta)
    confdiss <- sumList(x$confdiss)
    
    if (missing(xlim)) xlim <- range(as.vector(delta))
    if (missing(ylim)) ylim <- range(as.vector(confdiss))

    plot(as.vector(delta), as.vector(confdiss), main = main, type = "p", pch = 1,
         xlab = xlab, ylab = ylab, col = "darkgray", xlim = xlim, ylim = ylim, ...)

    if (!x$metric) {
      isofit <- isoreg(as.vector(delta), as.vector(confdiss))  #isotonic regression
      points(sort(isofit$x), isofit$yf, type = "b", pch = 16)
    } else {
      regfit <- lsfit(as.vector(delta), as.vector(confdiss))   #linear regression
      abline(regfit, lwd = 0.5)
    }
  }

  #--------------- Residual plot -------------------- 
  if (plot.type == "resplot") {
    if (missing(main)) main <- paste("Residual plot") else main <- main
    if (missing(xlab)) xlab <- "Aggregated Normalized Dissimilarities (d-hats)" else xlab <- xlab
    if (missing(ylab)) ylab <- "Aggregated Configuration Distances" else ylab <- ylab
    obsdiss <- sumList(x$obsdiss)
    confdiss <- sumList(x$confdiss)

    if (missing(xlim)) xlim <- range(as.vector(obsdiss))
    if (missing(ylim)) ylim <- range(as.vector(confdiss))

    plot(as.vector(obsdiss), as.vector(confdiss), main = main, type = "p", col = "darkgray", xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,...)
    abline(lsfit(as.vector(obsdiss), as.vector(confdiss)))

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
    if (missing(xlab)) xlab <- paste("Configurations D", x1,sep = "") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Configurations D", y1,sep = "") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$gspace[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$gspace[,y1])*1.1
    
    spp.perc <- x$spp/sum(x$spp)*100
    bubsize <- (max(spp.perc) - spp.perc + 1)/length(spp.perc)*bubscale
    
    plot(x$gspace, cex = bubsize, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
    xylabels <- x$gspace
    ysigns <- sign(x$gspace[,y1])
    xylabels[,2] <- (abs(x$gspace[,y1])-(x$gspace[,y1]*(bubsize/50)))*ysigns 
    text(xylabels, rownames(x$gspace), pos = 1,cex = 0.7)
  }  

  

  
}
