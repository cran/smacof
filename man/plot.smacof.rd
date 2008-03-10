\name{plot.smacof}
\alias{plot.smacof}
\alias{plot.smacofR}
\alias{plot.smacofID}


\title{2D SMACOF plots}
\description{These methods provide various 2D plots for SMACOF models.
}
\usage{
\method{plot}{smacof}(x, plot.type = "confplot", plot.dim = c(1,2), sphere = TRUE, 
main, xlab, ylab, xlim, ylim, ...)

\method{plot}{smacofR}(x, plot.type = "confplot", joint = FALSE, plot.dim = c(1,2), 
main, xlab, ylab, xlim, ylim, ...)

\method{plot}{smacofID}(x, plot.type = "confplot", plot.dim = c(1,2), main, xlab, ylab, xlim, ylim, ...)

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
  \item{sphere}{In case of spherical smacof, whether sphere should be plotted or not.}
  \item{joint}{If \code{TRUE}, the configurations are plotted jointly in rectangular smacof.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{plot}}} in package \code{scatterplot3d} for detailed information.
}

\details{\code{smacofSym()} creates object of class \code{"smacof"}, whereas \code{smacofRect()} produces 
\code{"smacofR"} and \code{smacofIndDiff()} generates \code{"smacofID"}.

Plot description:
 
- Configuration plot (\code{plot.type = "confplot"}): Plots the MDS configurations.

- Residual plot (\code{plot.type = "resplot"}): Plots the configuration distances against  
the corresponding residuals. 

- Shepard diagram (\code{plot.type = "Shepard"}): Diagram with the observed against the fitted distances including
isotonic regression line.

- Stress decomposition plot (\code{plot.type = "stressplot"}): Plots the stress contribution in of each observation.

For \code{smacofIndDiff()} the residual plot, Shepard diagram, and stress plot are based on the sum of the residuals across individuals/ways. 
The configuration plot represents the group stimulus space (i.e., joint configurations).

}

\seealso{\code{\link{plot3d.smacof}}}
\examples{

## 2D plots for spherical SMACOF
data(trading)
res <- smacofSym(trading)
plot(res, plot.type = "confplot")
plot(res, plot.type = "Shepard")
plot(res, plot.type = "stressplot")

## Joint configuration plot and row/column stressplots for rectangular SMACOF
data(breakfast)
res <- smacofRect(breakfast)
plot(res, plot.type = "confplot", joint = TRUE)
plot(res, plot.type = "stressplot")

}

\keyword{ hplot }
