\name{plot.smacof}
\alias{plot.smacof}
\alias{plot.smacofR}
\alias{plot.smacofID}


\title{2D SMACOF plots}
\description{These methods provide various 2D plots for SMACOF models.
}
\usage{
\method{plot}{smacof}(x, plot.type = "confplot", plot.dim = c(1,2), sphere = TRUE, 
                      bubscale = 1, col = 1, label.conf = list(label = TRUE, pos = 3, 
                      col = 1, cex = 0.8), hull.conf = list(hull = FALSE, col = 1, 
                      lwd = 1, ind = NULL), shepard.x = NULL, identify = FALSE, 
                      type = "p", pch = 20, cex = 0.5, asp = 1, main, xlab, ylab, 
                      xlim, ylim, col.hist = NULL, ...)

\method{plot}{smacofR}(x, plot.type = "confplot", what = c("both", "columns", "rows"), 
                       plot.dim = c(1,2), col.rows = hcl(0), col.columns = hcl(240), 
                       label.conf.rows = list(label = TRUE, pos = 3, 
                       col = hcl(0, l = 50), cex = 0.8), 
                       label.conf.columns = list(label = TRUE, pos = 3, 
                       col = hcl(240, l = 50), cex = 0.8),  
                       shepard.x = NULL, col.dhat = NULL, type = "p", pch = 20,
                       cex = 0.5, asp = 1, main, xlab, ylab, xlim, ylim, ...)

\method{plot}{smacofID}(x, plot.type = "confplot", plot.dim = c(1,2), bubscale = 1, 
                        col = 1, label.conf = list(label = TRUE, pos = 3, col = 1, 
                        cex = 0.8), identify = FALSE, type = "p", pch = 20,  cex = 0.5, 
                        asp = 1, plot.array, main, xlab, ylab, xlim, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"smacof"}, \code{"smacofR"}, and \code{"smacofID"} (see details)}
  \item{plot.type}{String indicating which type of plot to be produced: \code{"confplot"}, \code{"resplot"} 
  \code{"Shepard"}, \code{"stressplot"}, \code{"bubbleplot"} \code{"histogram"} (see details)}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{type}{What type of plot should be drawn (see also \code{\link[graphics:plot.default]{plot}}).}
  \item{pch}{Plot symbol.}
  \item{cex}{Symbol size.}
  \item{asp}{Aspect ratio.}
  \item{col}{Point color.}
  \item{sphere}{In case of spherical smacof, whether sphere should be plotted or not.}
  \item{bubscale}{Scaling factor (size) for the bubble plot.}
  \item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color). If \code{pos = 5} labels are placed away from the nearest point.}
  \item{hull.conf}{Option to add convex hulls to a configuration plot. Hull index needs to be provided.}
  \item{shepard.x}{Shepard plot only: original data (e.g. correlation matrix) can be provided for plotting on x-axis.}
  \item{identify}{If \code{TRUE}, the \code{identify()} function is called internally that allows to add configuration labels by mouse click.}
  \item{what}{For unfolding only: Whether row coordinates, column coordinates, or both should be plotted.}
  \item{col.rows}{Row colors in unfolding configuration plot.} 
  \item{col.columns}{Column colors in unfolding configuration plot.} 
  \item{col.dhat}{Shepard plot only: color specification of the dhats. For row conditional transformations in unfolding a vector of the length of the number of rows should be specified.}
  \item{label.conf.rows}{List with arguments for plotting the labels of the row configurations in an unfolding configuration plot (logical value whether to plot labels or not, label position, label color).} 
  \item{label.conf.columns}{List with arguments for plotting the labels of the columns configurations in an unfolding configuration plot (logical value whether to plot labels or not, label position, label color).}
  \item{col.hist}{Color of the borders of the histogram.}
  \item{plot.array}{Array arrangements of plots for individual difference models (see details).}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics:plot.default]{plot}} for detailed information.}
}

\details{\code{mds()} and \code{smacofSym()} create an object of class \code{"smacof"}, \code{unfolding()}, \code{prefscal()}, and \code{smacofRect()} produce \code{"smacofR"}, and \code{smacofIndDiff()} generates \code{"smacofID"}.

Plot description:
 
- Configuration plot (\code{plot.type = "confplot"}): Plots the MDS configuration.

- Residual plot (\code{plot.type = "resplot"}): Plots the disparities (d-hats) distances against  
the fitted distances. 

- Shepard diagram (\code{plot.type = "Shepard"}): Diagram with the observed dissimilarities against the fitted distances including (isotonic) regression line.

- Stress decomposition plot (\code{plot.type = "stressplot"}): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding smacof function to percentages (sum is 100). The higher the contribution, the worse the fit. 

- Bubble plot (\code{plot.type = "bubbleplot"}, not available for rectangular SMACOF): Combines the configuration plot with the point stress contribution. The larger the bubbles, the worse the fit. 

- Histogram (\code{plot.type = "histogram"}: gives a weighted histogram of the dissimilarities. For optional arguments, see \code{\link[weights]{wtd.hist}}.

For \code{smacofIndDiff()} the residual plot, Shepard diagram, and stress plot are based on the sum of the residuals across individuals/ways. The configuration plot represents the group stimulus space (i.e., joint configuration). If \code{plot.array} is not specified, it produces a Shepard plot of the distances summed across subjects, if \code{plot.array = 0} it produces a sqrt(nsubjects) times sqrt(nsubjects) array of graph panels, if \code{plot.array = 3} it produces 3x3 arrays of graph panels, if \code{plot.array = c(2, 3)} it produces 2x3 arrays of graph panels, and if \code{plot.array = c(3, 2, 5)} produces 3x2 arrays of panels (only the first two values are used).
}

\seealso{\code{\link{plot.procr}}}
\examples{

## 2D plots for simple MDS
data(trading)
res <- mds(trading)
plot(res, plot.type = "confplot")
plot(res, plot.type = "confplot", label.conf = list(pos = 5)) ## avoid overlapping labels
plot(res, plot.type = "Shepard")
plot(res, plot.type = "stressplot")
plot(res, plot.type = "resplot")
plot(res, plot.type = "bubbleplot")
plot(res, plot.type = "histogram")

## Add convex hulls to configuration plot
r <- cor(PVQ40, use = "pairwise.complete.obs")
diss <- sim2diss(r, method = "corr") 
res <- mds(delta = diss, type = "ordinal") 
codes <- substring(colnames(PVQ40), 1, 2)  ## supplementary variable
plot(res, hull.conf = list(hull = TRUE, ind = codes, col = "coral1", lwd = 2))

## Shepard plots
ekmanD <- sim2diss(ekman)
fit1 <- mds(ekmanD, type = "ordinal")
plot(fit1, plot.type = "Shepard")
plot(fit1, plot.type = "Shepard", shepard.x = ekman)  ## original data on x-axis

## Joint configuration plot and row/column stressplots for unfolding
data(breakfast)
res <- unfolding(breakfast)
plot(res, plot.type = "confplot")
plot(res, plot.type = "stressplot")
}

\keyword{ hplot }
