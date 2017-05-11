# plot method for IC exploration

plot.icexplore <- function(x, main = "IC Plot", xlab = "Dimension 1", ylab = "Dimension 2", cex.scale = 10, xlim, ylim, 
                           cex, asp = 1, ...)
                            
                            #hclpar = list(c = 50, l = 70), col.p, col.l, plot.lines = TRUE, , ...)
{
# x ... object of class smacofboot
  
  if (missing(xlim)) xlim <- range(x$conf[,1])*1.1
  if (missing(ylim)) ylim <- range(x$conf[,2])*1.1
  if (missing(cex)) cex <- (cex.scale*x$stressvec/2)^2
  
  textplot(x$conf[,1], x$conf[,2], words = rownames(x$conf), cex = cex, 
           xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab, ...)
}


  


  

 
