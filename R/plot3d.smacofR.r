`plot3d.smacofR` <-
function(x, plot.dim = c(1,2,3), joint = FALSE, xlab, ylab, zlab, 
         col, main, bgpng = NULL, ax.grid = TRUE, sphere.rgl = FALSE,...)
{
#S3 plot method for objects of class "smacof"
#plot.dim ... vector of length 3 with dimensions to be plotted against

 if (x$ndim < 3) stop("No 3D plots can be drawn for ndim < 3 !") 
 if (length(plot.dim) !=  3) stop("plot.dim must be of length 3!")
 pd1 <- plot.dim[1]
 pd2 <- plot.dim[2]
 pd3 <- plot.dim[3]
 if (pd3 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")

 if (is.null(bgpng)) {
   texture1 <- NULL
 } else {
   texture1 <- system.file(paste("textures/",bgpng,sep=""), package = "rgl")
 }

#------------------ configuration plot ---------------
  if (!joint) {                                                                                 #separate devices
      if (missing(main)) main1 <- paste("Configuration Plot - Columns") else main1 <- main        #plot column configurations
      if (missing(xlab)) xlab1 <- paste("Column Configurations D", pd1,sep = "") else xlab1 <- xlab
      if (missing(ylab)) ylab1 <- paste("Column Configurations D", pd2,sep = "") else ylab1 <- ylab
      if (missing(zlab)) zlab1 <- paste("Column Configurations D", pd3,sep = "") else zlab1 <- zlab
      if (missing(col)) col1 <- "BLUE"     
      
      open3d()
      bg3d(sphere = sphere.rgl, texture = texture1, back = "filled", color = "white")
      text3d(x$conf.col[,pd1], , x$conf.col[,pd3], texts = rownames(x$conf.col), col = col1, 
             alpha = 1, textenvmap = TRUE, lit = TRUE, ...)
      axes3d(c('x','y','z'), labels = TRUE, color = "black", alpha = 1)
      title3d(xlab = xlab1, ylab = ylab1, zlab = zlab1, main = main1, color = "black", alpha = 1)
  
      if (missing(main)) main1 <- paste("Configuration Plot - Rows") else main1 <- main        #plot column configurations
      if (missing(xlab)) xlab1 <- paste("Row Configurations D", pd1,sep = "") else xlab1 <- xlab
      if (missing(ylab)) ylab1 <- paste("Row Configurations D", pd2,sep = "") else ylab1 <- ylab
      if (missing(zlab)) zlab1 <- paste("Row Configurations D", pd3,sep = "") else zlab1 <- zlab
      if (missing(col)) col1 <- "RED"     
      
      open3d()
      bg3d(sphere = sphere.rgl, texture = texture1, back = "filled", color = "white")
      text3d(x$conf.row[,pd1], x$conf.row[,pd2], x$conf.row[,pd3], texts = rownames(x$conf.row), col = col1, 
      alpha = 1, textenvmap = TRUE, lit = TRUE, ...)
      axes3d(c('x','y','z'), labels = TRUE, color = "black", alpha = 1)
      title3d(xlab = xlab1, ylab = ylab1, zlab = zlab1, main = main1, color = "black", alpha = 1)
         
    } else {

      if (missing(main)) main1 <- paste("Joint Configuration Plot") else main1 <- main        
      if (missing(xlab)) xlab1 <- paste("Configurations D", pd1,sep = "") else xlab1 <- xlab  
      if (missing(ylab)) ylab1 <- paste("Configurations D", pd2,sep = "") else ylab1 <- ylab
      if (missing(zlab)) zlab1 <- paste("Configurations D", pd3,sep = "") else zlab1 <- ylab
      
      fullconf <- rbind(x$conf.col[,c(pd1,pd2,pd3)],x$conf.row[,c(pd1,pd2,pd3)])
      xlim1 <- range(fullconf[,1])
      ylim1 <- range(fullconf[,2])
      zlim1 <- range(fullconf[,3])
      
      open3d()
      bg3d(sphere = sphere.rgl, texture = texture1, back = "filled", color = "white")
      text3d(x$conf.col[,pd1], x$conf.col[,pd2], x$conf.col[,pd3], texts = rownames(x$conf.col), col = "BLUE", 
      alpha = 1, textenvmap = TRUE, lit = TRUE, ...)
      text3d(x$conf.row[,pd1], x$conf.row[,pd2], x$conf.row[,pd3], texts = rownames(x$conf.row), col = "RED", 
      alpha = 1, textenvmap = TRUE, lit = TRUE, ...)
      axes3d(c('x','y','z'), labels = TRUE, color = "black", alpha = 1)
      title3d(xlab = xlab1, ylab = ylab1, zlab = zlab1, main = main1, color = "black", alpha = 1)
   }
    
}

