\name{plot.smacof}
\alias{plot.smacof}
\alias{plot.smacofR}
\alias{plot.smacofID}


\title{2D SMACOF plots}
\description{These methods provide various 2D plots for SMACOF models.
}
\usage{
\method{plot}{smacof}(x, plot.type = "confplot", plot.dim = c(1,2), sphere = TRUE, 
bubscale = 3, label.conf = list(label = TRUE, pos = 1, col = 1), 
identify = FALSE, type, main, xlab, ylab, xlim, ylim, ...)

\method{plot}{smacofR}(x, plot.type = "confplot", joint = FALSE, plot.dim = c(1,2), 
col.rows = "red", col.columns = "blue", 
label.conf.rows = list(label = TRUE, pos = 1, col = "red"), 
label.conf.columns = list(label = TRUE, pos = 1, col = "blue"), 
type, main, xlab, ylab, xlim, ylim, ...)

\method{plot}{smacofID}(x, plot.type = "confplot", plot.dim = c(1,2), bubscale = 5, 
label.conf = list(label = TRUE, pos = 1, col = 1), identify = FALSE, type, main, xlab, 
ylab, xlim, ylim, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"smacof"}, \code{"smacofR"}, and \code{"smacofID"} (see details)}
  \item{plot.type}{String indicating which type of plot to be produced: \code{"confplot"}, \code{"resplot"} 
  \code{"Shepard"}, \code{"stressplot"} (see details)}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{type}{What type of plot should be drawn (see also \code{\link[graphics]{plot}}).}
  \item{sphere}{In case of spherical smacof, whether sphere should be plotted or not.}
  \item{bubscale}{Scaling factor (size) for the bubble plot.}
  \item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color).}
  \item{identify}{If \code{TRUE}, the \code{identify()} function is called internally that allows to add configuration labels by mouse click.}
  \item{joint}{If \code{TRUE}, the configurations are plotted jointly in rectangular smacof.}
  \item{col.rows}{Row colors in rectangular configuration plot.} 
  \item{col.columns}{Column colors in rectangular configuration plot.} 
  \item{label.conf.rows}{List with arguments for plotting the labels of the row configurations in a rectangular configuration plot (logical value whether to plot labels or not, label position, label color).} 
  \item{label.conf.columns}{List with arguments for plotting the labels of the columns configurations in a rectangular configuration plot (logical value whether to plot labels or not, label position, label color).}
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

- Stress decomposition plot (\code{plot.type = "stressplot"}): Plots the stress contribution in of each observation. The higher the contribution, the worse the fit. 

- Bubble plot (\code{plot.type = "bubbleplot"}, not available for rectangular SMACOF): Combines the configuration plot with the point stress contribution. The larger the bubbles, the better the fit. 

For \code{smacofIndDiff()} the residual plot, Shepard diagram, and stress plot are based on the sum of the residuals across individuals/ways. The configuration plot represents the group stimulus space (i.e., joint configurations).

}

\seealso{\code{\link{plot3d.smacof}}}
\examples{

## 2D plots for spherical SMACOF
data(trading)
res <- smacofSym(trading)
plot(res, plot.type = "confplot")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "stressplot")
plot(res, plot.type = "resplot")
plot(res, plot.type = "bubbleplot")


## Joint configuration plot and row/column stressplots for rectangular SMACOF
data(breakfast)
res <- smacofRect(breakfast)
plot(res, plot.type = "confplot", joint = TRUE)

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
