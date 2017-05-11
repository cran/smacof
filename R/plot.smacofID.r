# plot method for all smacof objects

plot.smacofID <- function(x, plot.type = "confplot", plot.dim = c(1,2), bubscale = 1, col = 1, 
                          label.conf = list(label = TRUE, pos = 3, col = 1), identify = FALSE, 
                          type = "p", pch = 20,  cex = 0.5, asp = 1, plot.array,
                          main, xlab, ylab, xlim, ylim, ...)

# x ... object of class smacofID
# plot.type ... types available: "confplot", "bubbleplot", "stressplot", "Shepard"
# Shepard plot and resplot are performed over sum of distances
  
{
  ## --- check type args:
  plot.type <- match.arg(plot.type, c("confplot", "Shepard", "resplot","bubbleplot", "stressplot"), several.ok = FALSE)
  
  
  ## --- check label lists
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$col)) label.conf$col <- 1
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  if (identify) label.conf$label <- FALSE
  
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  
  if (plot.type == "confplot") {
    if (missing(main)) main <- paste("Group Configuration") else main <- main
    if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab

    if (missing(xlim)) xlim <- range(x$gspace[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$gspace[,y1])*1.1
  
        
    plot(x$gspace[,x1], x$gspace[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
         xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
    if (label.conf$label) text(x$gspace[,x1], x$gspace[,y1], labels = rownames(x$gspace), 
                               cex = label.conf$cex, pos = label.conf$pos, 
                               col = label.conf$col)
    
    if (identify) {
       identify(x$gspace[,x1], x$gspace[,y1], labels = rownames(x$gspace), cex = 0.8)
    }
  }

  #---------------- Shepard diagram ------------------
  if (plot.type == "Shepard") {
    if (missing(main)) main <- paste("Shepard Diagram (Summed Distances)") else main <- main
    if (missing(xlab)) xlab <- "Observed Dissimilarities" else xlab <- xlab
    if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab
    
    if(missing(plot.array)) use_individual_distances <- FALSE else use_individual_distances <- TRUE
    
    if (use_individual_distances) {
      
      # do Shepard plots of each subject
      nvars <- length(x$delta)
      
      if (length(plot.array) < 2){
        npanv <- plot.array[1]
        npanh <- plot.array[1]
      } else {
        npanv <- plot.array[1]
        npanh <- plot.array[2]
      }
      
      if (npanv == 0 | npanh == 0) {
        npanv <- ceiling(sqrt(nvars)) 
        npanh <- floor(sqrt(nvars))
        if (npanv * npanh < nvars) npanv <- npanv + 1
      }
      
      # prevent invalid dimensions
      npanv <- max(npanv, 1)
      npanh <- max(npanh, 1)
      
      if (npanv == 1 && npanh == 1) parop <- FALSE else parop <- TRUE
      if (parop) op <- par(mfrow = c(npanv, npanh))
      
      if (is.null(names(x$conf))) namevec <- 1:nvars else namevec <- names(x$conf)
      for (i in 1:nvars) {
        main <- paste("Shepard Diagram", namevec[i])
        notmiss <- as.vector(x$weightmat[[i]] > 0)
        xcoor <- (as.vector(x$delta)[[i]])[notmiss]
        ycoor <- (as.vector(x$confdiss)[[i]])[notmiss]
        xlim <- range(xcoor)
        ylim <- range(ycoor)
        plot(xcoor, ycoor, main = main, type = "p", pch = pch, cex = cex,
             xlab = xlab, ylab = ylab, col = "darkgray", xlim = xlim, ylim = ylim)
        iord <- order(xcoor)
        points(xcoor[iord], x$dhat[[i]][iord], type = "b", pch = pch, cex = cex)
      }
      # restore old parameter setting
      if (parop) on.exit(par(op))
      
    } else {
      
      # Make a Shepard plot of the distance data summed across subjects
      delta <- sumList(x$delta)
      confdiss <- sumList(x$confdist)
      weightdiss <- sumList(x$weightmat)
      notmiss <- as.vector(weightdiss > 0)
      xcoor <- as.vector(delta)[notmiss]
      ycoor <- as.vector(confdiss)[notmiss]
      if (missing(xlim)) xlim <- range(xcoor)
      if (missing(ylim)) ylim <- range(ycoor)
      
      plot(xcoor, ycoor, main = main, type = "p", pch = 1,
           xlab = xlab, ylab = ylab, col = "darkgray", xlim = xlim, ylim = ylim, ...)
      
      if (x$type == "ordinal") {
        isofit <- isoreg(xcoor, ycoor)  #isotonic regression
        points(sort(isofit$x), isofit$yf, type = "b", pch = 16)
      } else {
        regfit <- lsfit(xcoor, ycoor)   #linear regression
        abline(regfit, lwd = 0.5)
      }
    }
  }
  
  #--------------- Residual plot -------------------- 
  if (plot.type == "resplot") {
      if (missing(main)) main <- paste("Residual plot") else main <- main
      if (missing(xlab)) xlab <- "Aggregated Normalized Dissimilarities (d-hats)" else xlab <- xlab
      if (missing(ylab)) ylab <- "Aggregated Configuration Distances" else ylab <- ylab
      obsdiss <- sumList(x$dhat)
      confdist <- sumList(x$confdist)
      
      if (missing(xlim)) xlim <- range(as.vector(obsdiss))
      if (missing(ylim)) ylim <- range(as.vector(confdist))
      
      plot(as.vector(obsdiss), as.vector(confdist), main = main, type = "p", col = "darkgray", xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,...)
      abline(lsfit(as.vector(obsdiss), as.vector(confdist)))
      
  }
  
  
  
  #----------------------- Stress decomposition -----------------
  if (plot.type == "stressplot") {
    if (missing(main)) main <- paste("Stress Decomposition Chart") else main <- main
    if (missing(xlab)) xlab <- "Objects" else xlab <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab

    spp.perc <- sort(x$spp, decreasing = TRUE)
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

    if (missing(xlim)) xlim <- range(x$gspace[,x1])*1.1
    if (missing(ylim)) ylim <- range(x$gspace[,y1])*1.1
    
    spp.perc <- x$spp/sum(x$spp)*100
    bubsize <- spp.perc/length(spp.perc)*(bubscale + 3)
    
    plot(x$gspace, cex = bubsize, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
    xylabels <- x$gspace
    ysigns <- sign(x$gspace[,y1])
    xylabels[,2] <- (abs(x$gspace[,y1])-(x$gspace[,y1]*(bubsize/50)))*ysigns 
    text(xylabels, rownames(x$gspace), pos = 3,cex = 0.7)    
  }  

}
