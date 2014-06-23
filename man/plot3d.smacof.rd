\name{plot3d.smacof}
\alias{plot3d.smacof}
\alias{plot3d.smacofR}
\alias{plot3d.smacofID}
\alias{plot3dstatic}
\alias{plot3dstatic.smacof}
\alias{plot3dstatic.smacofR}
\alias{plot3dstatic.smacofID}



\title{3D SMACOF plots}
\description{These methods produce static and dynamic 3D configuration plots for SMACOF models.
}
\usage{
\method{plot3d}{smacof}(x, plot.dim = c(1,2,3), sphere = FALSE, xlab, ylab, zlab, 
col, main, bgpng = NULL, ax.grid = TRUE, sphere.rgl = FALSE,...)

\method{plot3d}{smacofR}(x, plot.dim = c(1,2,3), joint = FALSE, xlab, ylab, zlab, 
col, main, bgpng = NULL, ax.grid = TRUE, sphere.rgl = FALSE,...)

\method{plot3d}{smacofID}(x, plot.dim = c(1,2,3), xlab, ylab, zlab, 
col, main, bgpng = NULL, ax.grid = TRUE, sphere.rgl = FALSE,...)

\method{plot3dstatic}{smacof}(x, plot.dim = c(1,2,3), main, xlab, ylab, zlab, col, ...)

\method{plot3dstatic}{smacofR}(x, plot.dim = c(1,2,3), main, xlab, ylab, zlab, col, joint = FALSE, ...)

\method{plot3dstatic}{smacofID}(x, plot.dim = c(1,2,3), main, xlab, ylab, zlab, col, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{"smacof"}, \code{"smacofR"}, and \code{"smacofID"} (see details)}
  \item{plot.dim}{Vector of length 3 with dimensions to be plotted.}
  \item{sphere}{Spherical SMACOF: Whether sphere should be plotted or not.}
  \item{joint}{Rectangular SMACOF: If \code{TRUE}, the configurations are plotted jointly.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{zlab}{Label of z-axis.}
  \item{col}{Color of the text labels.}
  \item{bgpng}{Background image from rgl library; \code{NULL} for white background}
  \item{ax.grid}{If \code{TRUE}, axes grid is plotted.}
  \item{sphere.rgl}{If \code{TRUE}, rgl sphere (background) is plotted.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{plot}}} in package \code{scatterplot3d} for detailed information.
}

\details{\code{smacofSym()} creates object of class \code{"smacof"}, whereas \code{smacofRect()} produces 
\code{"smacofR"} and \code{smacofIndDiff()} generates \code{"smacofID"}.

For \code{smacofIndDiff()} the configuration plot represents the group stimulus space (i.e., joint configurations).
}
\seealso{\code{\link{plot.smacof}}}
\examples{

## 3D plot for spherical SMACOF
data(trading)
res <- smacofSphere(trading, ndim = 3, verbose = FALSE, itmax = 5000)
plot3d(res, plot.type = "confplot")
plot3dstatic(res)

## Group stimulus space for rectangular SMACOF
data(breakfast)
res <- smacofRect(breakfast, ndim = 3)
plot3d(res, joint = TRUE)

}

\keyword{ hplot }
