# plot method for all smacof objects

plot.smacof <- function(x, plot.type = "confplot", plot.dim = c(1,2), sphere = TRUE, 
                        main, xlab, ylab, xlim, ylim, ...)

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
  
  if (x$model == "Spherical SMACOF (dual)") {             #remove first column
     x$obsdiss <- as.dist(as.matrix(x$obsdiss1)[,-1][-1,])
     x$confdiss <- as.dist(as.matrix(x$confdiss)[,-1][-1,])
  }   

  #----------------- configuration plot ---------------------
  if (plot.type == "confplot") {
    if (missing(main)) main <- paste("Configuration Plot") else main <- main
    if (missing(xlab)) xlab <- paste("Configurations D", x1,sep = "") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Configurations D", y1,sep = "") else ylab <- ylab

    #fullconf <- rbind(x$conf[,x1],x$conf[,y1])
    #if (missing(xlim)) xlim <- range(fullconf)
    #if (missing(ylim)) ylim <- range(fullconf)
    if (missing(xlim)) xlim <- range(x$conf[,x1])
    if (missing(ylim)) ylim <- range(x$conf[,y1])
    
    plot(x$conf[,x1], x$conf[,y1], main = main, type = "n", xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
    if ((any(class(x) == "smacofSP")) && (sphere)) {
      if (x$model == "Spherical SMACOF (dual)") {                     #dual smacof centered around first configuration row
        radius <- sqrt((abs(x$conf[2,1])+abs(x$conf[1,1]))^2 + (abs(x$conf[2,2])+abs(x$conf[1,2]))^2)    #sphere radius dual
        circle(x$conf[1,1], x$conf[1,2], radius, lty = 2, border = "lightgray")
      } else {
        radius <- sqrt(x$conf[2,1]^2 + x$conf[2,2]^2)
        circle(0, 0, radius, lty = 2, border = "lightgray")
      }
    }
      text(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), cex = 0.8)    
  }

  #---------------- Shepard diagram ------------------

  if (plot.type == "Shepard") {
    if (missing(main)) main <- paste("Shepard Diagram") else main <- main
    if (missing(xlab)) xlab <- "Dissimilarities" else xlab <- xlab
    if (missing(ylab)) ylab <- "Configuration Distances" else ylab <- ylab

    #ocdiss <- c(as.vector(x$obsdiss), as.vector(x$confdiss))
    #if (missing(xlim)) xlim <- range(ocdiss)
    #if (missing(ylim)) ylim <- range(ocdiss)
    if (missing(xlim)) xlim <- range(as.vector(x$obsdiss))
    if (missing(ylim)) ylim <- range(as.vector(x$confdiss))

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

    if (missing(xlim)) xlim <- range(as.vector(x$confdiss))
    if (missing(ylim)) ylim <- range(as.vector(resmat[lower.tri(resmat)]))
    
    plot(as.vector(x$confdiss), as.vector(resmat[lower.tri(resmat)]), main = main, type = "p",
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
    abline(h = 0, col = "lightgray", lty = 2)  
  }

  #----------------------- Stress decomposition -----------------
 
  if (plot.type == "stressplot") {
    if (missing(main)) main <- paste("Stress Decomposition Chart") else main <- main
    if (missing(xlab)) xlab <- "Objects" else xlab <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab
    stress.ri <- ((as.matrix(x$obsdiss) - as.matrix(x$confdiss))^2)   #sorted decomposed stress values
    stress.r <- rowSums(stress.ri)
    decomp.stress <- stress.r/(sum(stress.r))*100
    sdecomp.stress <- sort(decomp.stress, decreasing = TRUE)
    xaxlab <- names(sdecomp.stress)

    if (missing(xlim)) xlim1 <- c(1,length(decomp.stress)) else xlim1 <- xlim
    if (missing(ylim)) ylim1 <- range(sdecomp.stress) else ylim1 <- ylim
    
    plot(1:length(decomp.stress), sdecomp.stress, xaxt = "n", type = "p",
         xlab = xlab, ylab = ylab, main = main, xlim = xlim1, ylim = ylim1, ...)
    text(1:length(decomp.stress), sdecomp.stress, labels = xaxlab, pos = 3, cex = 0.8)
    for (i in 1:length(sdecomp.stress)) lines(c(i,i), c(sdecomp.stress[i],0), col = "lightgray", lty = 2)
                                  
  }

  #if (plot.type == "smearing") {
  #  delta.r <- as.matrix(x$confdiss)[1,]
  #  bw <- npregbw(formula=delta.r~x$conf[,1]+x$conf[,2], tol=.1, ftol=.1)
  #  model <- npreg(bws = bw)

    #predict(model, newdata)
    #x... sequence(min(x$conf[,1],max(x$conf[,1]))
    #y... sequence(min(x$conf[,2],max(x$conf[,2]))

}
