# plot method for all smacof objects

plot.smacofID <- function(x, plot.type = "confplot", plot.dim = c(1,2), main, xlab, ylab, ...)

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
    plot(x$gspace[,x1], x$gspace[,y1], main = main, type = "n", xlab = xlab, ylab = ylab, ...)
    text(x$gspace[,x1], x$gspace[,y1], labels = rownames(x$gspace), cex = 0.8)
  }

  #---------------- Shepard diagram ------------------
  if (plot.type == "Shepard") {
    if (missing(main)) main <- paste("Shepard Diagram") else main <- main
    if (missing(xlab)) xlab <- "Sum Observed Distances" else xlab <- xlab
    if (missing(ylab)) ylab <- "Sum Configuration Distances" else ylab <- ylab
    obsdiss <- sumList(x$obsdiss)
    confdiss <- sumList(x$confdiss)
    isofit <- isoreg(as.vector(obsdiss), as.vector(confdiss))  #isotonic regression
    plot(as.vector(obsdiss), as.vector(confdiss), main = main, type = "p", pch = 1,
         xlab = xlab, ylab = ylab, col = "lightgray", ...)
    points(sort(isofit$x), isofit$yf, type = "b", pch = 16)
  }

  #--------------- Residual plot --------------------
  if (plot.type == "resplot") {
    if (missing(main)) main <- paste("Residual plot") else main <- main
    if (missing(xlab)) xlab <- "Sum Configuration Distances" else xlab <- xlab
    if (missing(ylab)) ylab <- "Residuals" else ylab <- ylab
    obsdiss <- sumList(x$obsdiss)
    confdiss <- sumList(x$confdiss)
    resmat <- as.matrix(obsdiss - confdiss)
    plot(as.vector(confdiss), as.vector(resmat[lower.tri(resmat)]), main = main,
         xlab = xlab, ylab = ylab, ...)
    abline(h = 0, col = "lightgray", lty = 2)  
  }

  if (plot.type == "stressplot") {
    if (missing(main)) main <- paste("Stress Decomposition Chart") else main <- main
    if (missing(xlab)) xlab <- "Objects" else xlab <- xlab
    if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab
    obsdiss <- sumList(x$obsdiss)
    confdiss <- sumList(x$confdiss)
    stress.ri <- ((as.matrix(obsdiss) - as.matrix(confdiss))^2)   #sorted decomposed stress values
    stress.r <- rowSums(stress.ri)
    decomp.stress <- stress.r/(sum(stress.r))*100
    sdecomp.stress <- sort(decomp.stress, decreasing = TRUE)
    xaxlab <- names(sdecomp.stress)
    plot(1:length(decomp.stress), sdecomp.stress, xaxt = "n", type = "p",
         xlab = xlab, ylab = ylab, main = main, ...)
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
