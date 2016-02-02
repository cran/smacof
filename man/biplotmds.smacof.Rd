\name{biplotmds}
\alias{biplotmds}
\alias{biplotmds.smacof}
\alias{biplotmds.smacofID}
\alias{plot.mdsbi}


\title{MDS Biplots}
\description{Regresses external variables on a MDS configuration which results in a MDS biplot.
}
\usage{
\method{biplotmds}{smacof}(object, extvar, scale = TRUE)
\method{plot}{mdsbi}(x, plot.dim = c(1,2), sphere = TRUE, col = 1, 
                     label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), 
                     vec.conf = list(col = 1, cex = 0.8, length = 0.1), 
                     identify = FALSE, type = "p", pch = 20, 
                     asp = 1, main, xlab, ylab, xlim, ylim, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"smacof"} or \code{"smacofID"}}
  \item{extvar}{Data frame with external variables}
  \item{scale}{If \code{TRUE}, external variables are standardized internally}
     
  \item{x}{Object of class \code{"mdsbi"}}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{pch}{Plot symbol.}
  \item{asp}{Aspect ratio.}
  \item{col}{Point color.}
  \item{type}{What type of plot should be drawn.}
  \item{sphere}{In case of spherical smacof, whether sphere should be plotted or not.}
  \item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color).}
  \item{vec.conf}{List with arguments for arrows and arrow labels of the external variables.}
  \item{identify}{If \code{TRUE}, the \code{identify()} function is called internally that allows to add configuration labels by mouse click.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{plot}} for detailed information.}
}

\details{If a model for individual differences is provided, the external variables are regressed on the group stimulus space configurations. 
}

\value{   
Returns an object belonging to classes \code{"mlm"} and \code{"mdsbi"}. See \code{\link{lm}} for details.

  \item{R2vec}{Vector containing the single R2 values.}
}

\references{
Greenacre, M. (2010). Biplots in Practice. Fundacion BBVA, Bilbao, Spain
}

\author{Patrick Mair}

\examples{
## morsecode data with external scales
res <- mds(morse)
fitbi <- biplotmds(res, morsescales[,2:3])
plot(fitbi, main = "MDS Biplot")
}

\keyword{models}
