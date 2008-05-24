plot3dstatic.smacofR <- function(x, plot.dim = c(1,2,3), main, xlab, ylab, zlab, col, joint = FALSE, ...)
{
#produces static 3D-scatterplot

options(locatorBell = FALSE)
if (x$ndim < 3) stop("No 3D plots can be drawn for ndim < 3 !")
if (length(plot.dim) !=  3) stop("plot.dim must be of length 3!")
pd1 <- plot.dim[1]
pd2 <- plot.dim[2]
pd3 <- plot.dim[3]
if (pd3 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")

#----------------------------- configuration plot -------------------------------------
if (!joint) {                                                                                 #separate devices
      if (missing(main)) main1 <- paste("Configuration Plot - Columns") else main1 <- main        #plot column configurations
      if (missing(xlab)) xlab1 <- paste("Column Configurations D", pd1,sep = "") else xlab1 <- xlab
      if (missing(ylab)) ylab1 <- paste("Column Configurations D", pd2,sep = "") else ylab1 <- ylab
      if (missing(zlab)) zlab1 <- paste("Column Configurations D", pd3,sep = "") else zlab1 <- zlab
      if (missing(col)) col1 <- "BLUE"
      x1 <- x$conf.col[,pd1]
      y1 <- x$conf.col[,pd2]
      z1 <- x$conf.col[,pd3]

      pr <- scatterplot3d(x1, y1, z1, type = "n", main = main1, xlab = xlab1, ylab = ylab1,
                    zlab = zlab1, ...)
      text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$conf.col), col = col1)
      
      get(getOption("device"))()
      
      if (missing(main)) main1 <- paste("Configuration Plot - Rows") else main1 <- main        #plot column configurations
      if (missing(xlab)) xlab1 <- paste("Row Configurations D", pd1,sep = "") else xlab1 <- xlab
      if (missing(ylab)) ylab1 <- paste("Row Configurations D", pd2,sep = "") else ylab1 <- ylab
      if (missing(zlab)) zlab1 <- paste("Row Configurations D", pd3,sep = "") else zlab1 <- zlab
      if (missing(col)) col1 <- "RED"
      
      x1 <- x$conf.row[,pd1]
      y1 <- x$conf.row[,pd2]
      z1 <- x$conf.row[,pd3]

      pr <- scatterplot3d(x1, y1, z1, type = "n", main = main1, xlab = xlab1, ylab = ylab1,
                    zlab = zlab1, ...)
      text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$conf.row), col = col1)

  } else {

      if (missing(main)) main1 <- paste("Joint Configuration Plot") else main1 <- main
      if (missing(xlab)) xlab1 <- paste("Configurations D", pd1,sep = "") else xlab1 <- xlab
      if (missing(ylab)) ylab1 <- paste("Configurations D", pd2,sep = "") else ylab1 <- ylab
      if (missing(zlab)) zlab1 <- paste("Configurations D", pd3,sep = "") else zlab1 <- ylab

      nc <- nrow(x$conf.col)
      nr <- nrow(x$conf.row)
      fullconf <- rbind(x$conf.col[,c(pd1,pd2,pd3)],x$conf.row[,c(pd1,pd2,pd3)])
      xlim1 <- range(fullconf[,1])
      ylim1 <- range(fullconf[,2])
      zlim1 <- range(fullconf[,3])

      x1 <- fullconf[,pd1]
      y1 <- fullconf[,pd2]
      z1 <- fullconf[,pd3]
      
      pr <- scatterplot3d(x1, y1, z1, type = "n", main = main1, xlab = xlab1, ylab = ylab1,
                    zlab = zlab1, xlim = xlim1, ylim = ylim1, zlim = zlim1, ...)

      text(pr$xyz.convert(x1[1:nc], y1[1:nc], z1[1:nc]), labels = rownames(x$conf.col), col = "BLUE")
      text(pr$xyz.convert(x1[(nc+1):(nc+nr)], y1[(nc+1):(nc+nr)], z1[(nc+1):(nc+nr)]), labels = rownames(x$conf.row), col = "RED")
  }
}



