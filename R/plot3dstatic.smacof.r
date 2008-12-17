plot3dstatic.smacof <- function(x, plot.dim = c(1,2,3), main, xlab, ylab, zlab, col, ...)
{
#produces static 3D-scatterplot

options(locatorBell = FALSE)
if (x$ndim < 3) stop("No 3D plots can be drawn for ndim < 3 !")
if (length(plot.dim) !=  3) stop("plot.dim must be of length 3!")
pd1 <- plot.dim[1]
pd2 <- plot.dim[2]
pd3 <- plot.dim[3]
if (pd3 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")

if (missing(xlab)) xlab <- paste("Dimension",pd1)
if (missing(ylab)) ylab <- paste("Dimension",pd2)
if (missing(zlab)) zlab <- paste("Dimension",pd3)

x1 <- x$conf[,pd1]
y1 <- x$conf[,pd2]
z1 <- x$conf[,pd3]

#----------------------------- configuration plot -------------------------------------
if (missing(main)) main1 <- "Configuration Plot" else main1 <- main
if (missing(col)) col1 <- "blue" else col1 <- col

pr <- scatterplot3d(x1, y1, z1, type = "n", main = main1, xlab = xlab, ylab = ylab,
                    zlab = zlab, ...)
text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$conf), col = col1)

}
