svm_mdsplot <- function(mds_object, svm_object, class, legend1 = TRUE, legend2 = TRUE, 
                        inset = c(-0.2, 0.5), plot.dim = c(1,2), by = 0.01, 
                        main, xlab, ylab, xlim, ylim, ...) 
{
  X <- mds_object$conf[,plot.dim]
  X1 <- mds_object$conf
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  if (missing(main)) main <- paste("MDS-SVM Configuration Plot") else main <- main
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  if (missing(xlim)) xlim <- range(X[, 1])*1.1 else xlim <- xlim
  if (missing(ylim)) ylim <- range(X[, 2])*1.1 else ylim <- ylim
  
  nc <- svm_object$nclasses
  minmax <- range(c(X))*2.2      
  seqs <- seq(minmax[1], minmax[2], by = by)
  dnames <- attr(svm_object$terms, "term.labels")[plot.dim]
  
  xgrid <- expand.grid(list(seqs, seqs))
  colnames(xgrid) <- dnames
  
  ## ---- modify xgrid for >2D solutions
  if (ncol(X1) > 2) {
    col_ind <- which(colnames(X1) %in% colnames(xgrid))
    xgrid1 <- as.data.frame(matrix(0, ncol = ncol(X1), nrow = nrow(xgrid)))
    colnames(xgrid1) <- colnames(X1)
    xgrid1[, col_ind] <- xgrid
    xgrid <- xgrid1
  }
  lev <- predict(svm_object, newdata = xgrid)
  levmat <- matrix(as.numeric(lev), sqrt(length(lev)))
  
  op <- par(pty = "s", mar = c(5.1, 4.1, 4.1, 9.1), xpd = TRUE)
  image(seqs, seqs, levmat, xlab = xlab, ylab = ylab, main = main, col = rainbow_hcl(nc), asp = 1,
        xlim = xlim, ylim = ylim)
  points(X, cex = 0.7, pch = (1:nc)[as.numeric(class)])      
  text(X, labels = rownames(X), pos = 3, cex = 0.7)
  if (legend1) legend("topright", inset = c(inset[1], 0), legend = svm_object$levels, pch = 15, title = "Regions", col = rainbow_hcl(nc), cex = 0.7)
  if (legend2) legend("topright", inset = c(inset[1], inset[2]), legend = svm_object$levels, title = "Classes", pch = 1:nc, cex = 0.7)
  par(op)
}
