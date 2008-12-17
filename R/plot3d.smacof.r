`plot3d.smacof` <-
function(x, plot.dim = c(1,2,3), sphere = TRUE, xlab, ylab, zlab, 
         col, main, bgpng = "particle.png", ax.grid = TRUE, sphere.rgl = TRUE,...)
{
#S3 plot method for objects of class "smacof"
#plot.dim ... vector of length 3 with dimensions to be plotted against

 if (x$ndim < 3) stop("No 3D plots can be drawn for ndim < 3 !") 
 if (length(plot.dim) !=  3) stop("plot.dim must be of length 3!")
 pd1 <- plot.dim[1]
 pd2 <- plot.dim[2]
 pd3 <- plot.dim[3]
 if (pd3 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")

 x1 <- x$conf[,pd1]
 y1 <- x$conf[,pd2]
 z1 <- x$conf[,pd3]

 if (missing(xlab)) xlab <- paste("Dimension",pd1)
 if (missing(ylab)) ylab <- paste("Dimension",pd2)
 if (missing(zlab)) zlab <- paste("Dimension",pd3)

 if (is.null(bgpng)) {
   texture1 <- NULL
 } else {
   texture1 <- system.file(paste("textures/",bgpng,sep=""), package = "rgl")
 }

#------------------ configuration plot ---------------
 if (missing(main)) main1 <- "Configuration Plot"  else main1 <- main
 if (missing(col)) col1 <- "blue" else col1 <- col

 rgl.open()
 rgl.bg(sphere = sphere.rgl, texture = texture1, back = "filled", color = "white")

 
 if ((any(class(x) == "smacofSP")) && (sphere)) {
   if (x$model == "Spherical SMACOF (dual)") {                     #dual smacof centered around first configuration row
     a.x1 <- abs(x1)
     a.y1 <- abs(y1)
     a.z1 <- abs(z1)
     radius.sphere <- sqrt(((a.x1[2]+a.x1[1])^2) + ((a.y1[2]+a.y1[1])^2) + ((a.z1[2]+a.z1[1])^2)) 
     rgl.spheres(x1[1], y1[1], z1[1], radius = radius.sphere, col = "white", alpha = 0.8, back = "cull", front = "line")
   } else {
     radius.sphere <- sqrt(x1^2 + y1^2 + z1^2)                           #Pythagoras 3D
     rgl.spheres(0,0,0, radius = radius.sphere, col = "white", alpha = 0.8, back = "cull", front = "line")
   }
 }
 text3d(x1, y1, z1, texts = rownames(x$conf), col = col1, alpha = 1, ...)
 axes3d(c('x','y','z'), labels = TRUE, color = "black", alpha = 1)
 title3d(xlab = xlab, ylab = ylab, zlab = zlab, main = main1, color = "black", alpha = 1)
}

