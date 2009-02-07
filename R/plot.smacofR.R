# plot method for rectangular smacof

plot.smacofR <- function(x, plot.type = "confplot", joint = FALSE, plot.dim = c(1,2), 
                         main, xlab, ylab, xlim, ylim, ...)

# x ... object of class smacofR
# plot.type ... types available: "confplot", "Shepard", "stressplot", "resplot"
# joint ... if TRUE, row and column configurations in 1 plot
  
{
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  #--------------------------- configuration plot -----------------------
  if (plot.type == "confplot") {
    if (!joint) {                                                                                 #separate devices
      if (missing(main)) main1 <- paste("Configuration Plot - Columns") else main1 <- main        #plot column configurations
      if (missing(xlab)) xlab1 <- paste("Column Configurations D", x1,sep = "") else xlab1 <- xlab
      if (missing(xlab)) ylab1 <- paste("Column Configurations D", y1,sep = "") else ylab1 <- ylab

      fullconf <- rbind(x$conf.col[,c(x1,y1)],x$conf.row[,c(x1,y1)])
      if (missing(xlim)) xlim <- range(fullconf)
      if (missing(ylim)) ylim <- range(fullconf)

      plot(x$conf.col[,x1], x$conf.col[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = "n", xlim = xlim, ylim = ylim,...)
      text(x$conf.col[,x1], x$conf.col[,y1], labels = rownames(x$conf.col), col = "BLUE")

      if (missing(main)) main1 <- paste("Configuration Plot - Rows") else main1 <- main       #plot row configurations
      if (missing(xlab)) xlab1 <- paste("Column Configurations D", x1,sep = "") else xlab1 <- xlab
      if (missing(ylab)) ylab1 <- paste("Column Configurations D", y1,sep = "") else ylab1 <- ylab
      dev.new()
      plot(x$conf.row[,x1], x$conf.row[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = "n", xlim = xlim, ylim = ylim, ...)
      text(x$conf.row[,x1], x$conf.row[,y1], labels = rownames(x$conf.row), col = "RED")
      
    } else {                                                                                     #joint plot
      if (missing(main)) main1 <- paste("Joint Configuration Plot") else main1 <- main        
      if (missing(xlab)) xlab1 <- paste("Configurations D", x1,sep = "") else xlab1 <- xlab  
      if (missing(ylab)) ylab1 <- paste("Configurations D", y1,sep = "") else ylab1 <- ylab

      fullconf <- rbind(x$conf.col[,c(x1,y1)],x$conf.row[,c(x1,y1)])
      if (missing(xlim)) xlim <- range(fullconf)
      if (missing(ylim)) ylim <- range(fullconf)
      
      plot(x$conf.col[,x1], x$conf.col[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = "n", xlim = xlim, ylim = ylim, ...)
      text(x$conf.col[,x1], x$conf.col[,y1], labels = rownames(x$conf.col), col = "BLUE")
      points(x$conf.row[,x1], x$conf.row[,y1], main = main1, xlab = xlab1, ylab = ylab1, type = "n")
      text(x$conf.row[,x1], x$conf.row[,y1], labels = rownames(x$conf.row), col = "RED")
    }
  }

  #---------------- Shepard diagram ------------------
  if (plot.type == "Shepard") {
    if (missing(main)) main <- paste("Shepard Diagram") else main <- main
    if (missing(xlab)) xlab <- "Dissimilarities" else xlab <- xlab
    if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab

    ocdiss <- c(as.vector(x$obsdiss), as.vector(x$confdiss))
    if (missing(xlim)) xlim <- range(ocdiss)
    if (missing(ylim)) ylim <- range(ocdiss)

    plot(as.vector(x$obsdiss), as.vector(x$confdiss), main = main, type = "p", pch = 1,
         xlab = xlab, ylab = ylab, col = "lightgray", xlim = xlim, ylim = ylim, ...)

    if (!x$metric) {
      isofit <- isoreg(as.vector(x$obsdiss), as.vector(x$confdiss))  #isotonic regression
      points(sort(isofit$x), isofit$yf, type = "b", pch = 16)
    } else {
      abline(0,1, lty = 2)
    }
   
  }

  #--------------- Residual plot --------------------
  if (plot.type == "resplot") {
    if (missing(main)) main <- paste("Residual plot") else main <- main
    if (missing(xlab)) xlab <- "Configuration Distances" else xlab <- xlab
    if (missing(ylab)) ylab <- "Residuals" else ylab <- ylab
    resmat <- residuals(x)

    fullres <- c(as.vector(x$confdiss), as.vector(resmat)) 
    if (missing(xlim)) xlim <- range(fullres)
    if (missing(ylim)) ylim <- range(fullres)
    
    plot(as.vector(x$confdiss), as.vector(resmat), main = main, type = "p",
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
    abline(h = 0, col = "lightgray", lty = 2)  
  }

  #----------------------- Stress decomposition -----------------
  if (plot.type == "stressplot") {
    if (missing(main)) main1 <- paste("Stress Decomposition Chart - Rows") else main1 <- main
    if (missing(main)) main2 <- paste("Stress Decomposition Chart - Columns") else main2 <- main
    
    if (missing(xlab)) xlab1 <- "Row Objects" else xlab1 <- xlab
    if (missing(xlab)) xlab2 <- "Column Objects" else xlab2 <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab
    stress.ri <- ((as.matrix(x$obsdiss) - as.matrix(x$confdiss))^2) 
    
    # row-wise
    stress.r <- rowSums(stress.ri)
    decomp.stress <- stress.r/(sum(stress.r))*100
    sdecomp.stress <- sort(decomp.stress, decreasing = TRUE)
    xaxlab <- names(sdecomp.stress)

    if (missing(xlim)) xlim1 <- c(1,length(decomp.stress)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(sdecomp.stress) else ylim1 <- ylim
    
    plot(1:length(decomp.stress), sdecomp.stress, xaxt = "n", type = "p",
         xlab = xlab1, ylab = ylab, main = main1, xlim = xlim1, ylim = ylim1, ...)
    text(1:length(decomp.stress), sdecomp.stress, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(sdecomp.stress)) lines(c(i,i), c(sdecomp.stress[i],0), col = "lightgray", lty = 2)                               
 
    #column-wise
    dev.new()
    stress.c <- colSums(stress.ri)
    decomp.stress <- stress.c/(sum(stress.c))*100
    sdecomp.stress <- sort(decomp.stress, decreasing = TRUE)
    xaxlab <- names(sdecomp.stress)

    if (missing(xlim)) xlim1 <- c(1,length(decomp.stress)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(sdecomp.stress) else ylim1 <- ylim

    plot(1:length(decomp.stress), sdecomp.stress, xaxt = "n", type = "p",
         xlab = xlab2, ylab = ylab, main = main2, xlim = xlim1, ylim = ylim1, ...)
    text(1:length(decomp.stress), sdecomp.stress, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(sdecomp.stress)) lines(c(i,i), c(sdecomp.stress[i],0), col = "lightgray", lty = 2)                               
 
  }

}
