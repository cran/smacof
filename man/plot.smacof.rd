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
                      col = 1, cex = 0.8), shepard.x = NULL, identify = FALSE, 
                      type = "p", pch = 20, asp = 1, main, xlab, ylab, 
                      xlim, ylim, col.hist = NULL, ...)

\method{plot}{smacofR}(x, plot.type = "confplot", joint = TRUE, plot.dim = c(1,2), 
                       col.rows = hcl(0), col.columns = hcl(240), 
                       label.conf.rows = list(label = TRUE, pos = 3, 
                       col = hcl(0, l = 50), cex = 0.8), 
                       label.conf.columns = list(label = TRUE, pos = 3, 
                       col = hcl(240, l = 50), cex = 0.8), type = "p", pch = 20, 
                       asp = 1, main, xlab, ylab, xlim, ylim, ...)

\method{plot}{smacofID}(x, plot.type = "confplot", plot.dim = c(1,2), bubscale = 1, 
                        col = 1, label.conf = list(label = TRUE, pos = 3, col = 1), 
                        identify = FALSE, type = "p", pch = 20,  asp = 1, 
                        main, xlab, ylab, xlim, ylim, ...)
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
  \item{type}{What type of plot should be drawn (see also \code{\link[graphics]{plot}}).}
  \item{pch}{Plot symbol.}
  \item{asp}{Aspect ratio.}
  \item{col}{Point color.}
  \item{sphere}{In case of spherical smacof, whether sphere should be plotted or not.}
  \item{bubscale}{Scaling factor (size) for the bubble plot.}
  \item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color).}
  \item{shepard.x}{Shepard plot only: original data (e.g. correlation matrix) can be provided for plotting on x-axis.}
  \item{identify}{If \code{TRUE}, the \code{identify()} function is called internally that allows to add configuration labels by mouse click.}
  \item{joint}{If \code{TRUE}, the configurations are plotted jointly in rectangular smacof.}
  \item{col.rows}{Row colors in rectangular configuration plot.} 
  \item{col.columns}{Column colors in rectangular configuration plot.} 
  \item{label.conf.rows}{List with arguments for plotting the labels of the row configurations in a rectangular configuration plot (logical value whether to plot labels or not, label position, label color).} 
  \item{label.conf.columns}{List with arguments for plotting the labels of the columns configurations in a rectangular configuration plot (logical value whether to plot labels or not, label position, label color).}
  \item{col.hist}{Color of the borders of the histogram.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{plot}} for detailed information.}
}

\details{\code{smacofSym()} creates object of class \code{"smacof"}, whereas \code{smacofRect()} produces 
\code{"smacofR"} and \code{smacofIndDiff()} generates \code{"smacofID"}.

Plot description:
 
- Configuration plot (\code{plot.type = "confplot"}): Plots the MDS configurations.

- Residual plot (\code{plot.type = "resplot"}): Plots the normalized dissimilarities (d-hats) distances against  
the fitted distances. 

- Shepard diagram (\code{plot.type = "Shepard"}): Diagram with the observed dissimilarities against the fitted distances including
(isotonic) regression line.

- Stress decomposition plot (\code{plot.type = "stressplot"}): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding 
smacof function to percentages (sum is 100). The higher the contribution, the worse the fit. 

- Bubble plot (\code{plot.type = "bubbleplot"}, not available for rectangular SMACOF): Combines the configuration plot with the point stress contribution. The larger the bubbles, the worse the fit. 

- Histogram (\code{plot.type = "histogram"}: gives a weighted histogram of the dissimilarities. For optional arguments, see \code{\link[weights]{wtd.hist}}.

For \code{smacofIndDiff()} the residual plot, Shepard diagram, and stress plot are based on the sum of the residuals across individuals/ways. The configuration plot represents the group stimulus space (i.e., joint configurations).

}

\seealso{\code{\link{plot.procr}}}
\examples{

## 2D plots for basic SMACOF
data(trading)
res <- smacofSym(trading)
plot(res, plot.type = "confplot")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "stressplot")
plot(res, plot.type = "resplot")
plot(res, plot.type = "bubbleplot")

data(wish)
wish.dist <- sim2diss(wish, method = 7)
res <- smacofSym(wish.dist, type = "ordinal")
res
plot(res, ylim = c(-0.6, 1))    ## aspect ratio = 1 (default)
plot(res, label.conf = list(FALSE), asp = NA) ## asp not 1, no labels
plot(res, type = "n")           ## labels only (at configuration coordinates)

## Shepard plots
ekmanD <- sim2diss(ekman)
fit1 <- mds(ekmanD, type = "ordinal")
plot(fit1, plot.type = "Shepard")
plot(fit1, plot.type = "Shepard", shepard.x = ekman)  ## original data

## Joint configuration plot and row/column stressplots for rectangular SMACOF
data(breakfast)
res <- smacofRect(breakfast)
plot(res, plot.type = "confplot")

plot(res, plot.type = "stressplot")

plot(res, type = "p", pch = 25)

plot(res, type = "p", pch = 25, col.columns = 3, 
label.conf.columns = list(label = TRUE, pos = 3, col = 3), col.rows = 8, 
label.conf.rows = list(label = TRUE, pos = 3, col = 8))

plot(res, joint = TRUE, type = "p", pch = 25, col.columns = 4, 
label.conf.columns = list(label = TRUE, pos = 3, col = 4), 
col.rows = 8, label.conf.rows = list(label = TRUE, pos = 3, col = 8))
}

\keyword{ hplot }
