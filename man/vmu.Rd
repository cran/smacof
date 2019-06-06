\name{vmu}
\alias{vmu}
\alias{print.vmu}
\alias{plot.vmu}

\title{Vector Model of Unfolding}
\description{
Computes the metric vector model of unfolding (VMU) on rectangular input data (preferences, ratings) with the individuals (rows) represented as vectors in the biplot. There is also the option to fix the column coordinates. 
}
\usage{
vmu(delta, ndim = 2, center = TRUE, scale = FALSE, col.coord = NULL)

\method{plot}{vmu}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Data frame or matrix of preferences, ratings, dissimilarities}
  \item{ndim}{Number of dimensions}
  \item{center}{If \code{TRUE} input data are centered row-wise.}
  \item{scale}{If \code{TRUE} input data are scaled row-wise.}
  \item{col.coord}{Optional fixed coordinates for the column objects in \code{delta}.}
  \item{x}{Object of class \code{"vmu"}.}
  \item{...}{Additional arguments passed to \code{biplot} in \code{stats}.}
}

\value{
  \item{conf.row}{Row coordinates}
  \item{conf.col}{Column coordinates}
  \item{VAF}{variance accounted for}
}

\references{
Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.

Borg, I., Groenen, P. J. F., & Mair, P. (2018). Applied Multidimensional Scaling and Unfolding (2nd ed.). Springer.

Tucker, L. R. (1960). Intra-individual and inter-individual multidimensionality. In H. Gulliksen & S. Messick (Eds.), Psychological scaling: Theory and applications (pp. 155-167). Wiley.
}

\author{Ingwer Borg and Patrick Mair}

\seealso{\code{\link{biplot}}, \code{\link{unfolding}}}

\examples{
## VMU on portrait value questionnaire ratings
fit_vmu <- vmu(PVQ40agg)         ## fit 2D VMU
fit_vmu
plot(fit_vmu, cex = c(1, 0.7))   ## call biplot from stats

## VMU with fixed column coordinates (circular)
tuv <- matrix(0, nrow = 10, ncol = 2)
alpha <- -360/10
for (i in 1:10){
  alpha <- alpha+360/10
  tuv[i, 1]<- cos(alpha*pi/180)
  tuv[i, 2] <- sin(alpha*pi/180) 
}
fit_vmu2 <- vmu(PVQ40agg, col.coord = tuv)  ## fit 2D circular VMU
fit_vmu2
plot(fit_vmu2, cex = c(1, 0.7))
}

\keyword{multivariate}
