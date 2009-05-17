`plot3d.smacofID` <-
function(x, plot.dim = c(1,2,3), xlab, ylab, zlab, 
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

 if (missing(xlab)) xlab <- paste("Dimension",pd1)
 if (missing(ylab)) ylab <- paste("Dimension",pd2)
 if (missing(zlab)) zlab <- paste("Dimension",pd3)

 if (is.null(bgpng)) {
   texture1 <- NULL
 } else {
   texture1 <- system.file(paste("textures/",bgpng,sep=""), package = "rgl")
 }

#------------------ configuration plot ---------------
 
 if (missing(col)) col1 <- "blue" else col1 <- col
 if (missing(main)) main1 <- paste("Group Configurations")  else main1 <- main
  x1 <- x$gspace[,pd1]
  y1 <- x$gspace[,pd2]
  z1 <- x$gspace[,pd3]

  open3d()
  bg3d(sphere = sphere.rgl, texture = texture1, back = "filled", color = "white")
    
  text3d(x1, y1, z1, texts = names(x1), col = col1, ...)
  axes3d(c('x','y','z'), labels = TRUE, color = "black")
  title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "black")   


}

