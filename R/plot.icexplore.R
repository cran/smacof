# plot method for IC exploration

plot.icexplore <- function(x, main = "IC Plot", xlab = "Dimension 1", ylab = "Dimension 2", cex.scale = 10, 
                           cex, asp = 1, ...)
                            
                           
{
# x ... object of class smacofboot
  
  #if (missing(xlim)) xlim <- range(x$conf[,1])*1.1
  #if (missing(ylim)) ylim <- range(x$conf[,2])*1.1
  if (missing(cex)) cex <- (cex.scale*x$stressvec/2)^2
  
  ## color shading
  cols <- rev(c("#F7F7F7", "#D9D9D9", "#BDBDBD", "#969696", "#636363", "#252525"))
  pal <- colorRampPalette(cols)
  stress_order <- findInterval(cex, sort(cex))

  textplot(x$conf[,1], x$conf[,2], words = rownames(x$conf), cex = cex, col = pal(nrow(x$conf))[stress_order], 
           asp = asp, main = main, xlab = xlab, ylab = ylab)
}


  


  

 
