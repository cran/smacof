\name{fitCircle}
\alias{fitCircle}

\title{Fitting circle into point configuration}
\description{
Utility function for fitting a circle into 2D point configurations.
}
\usage{
fitCircle(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector with x-coordinates}
  \item{y}{Vector with y-coordinates}
}

\value{
  \item{cx}{x-coordinate center}
  \item{cy}{y-coordinate center}
  \item{radius}{circle radius}
}

\references{
Pratt, V. (1987). Direct least-squares fitting of algebraic surfaces. Computer Graphics, Vol. 21, pages 145-152.
}

\author{
Retrieved on R help: \url{http://r.789695.n4.nabble.com/Fit-circle-with-R-td4635157.html}.
}

\examples{
## Dataset on Schwartz values:
require(plotrix)
valsD <- 1 - cor(indvalues)
fit <- mds(valsD)
plot(fit, main = "MDS Value Circle")
circle <- fitCircle(fit$conf[,1], fit$conf[,2])
draw.circle(circle$cx, circle$cy, radius = circle$radius, border = "gray") 
}

\keyword{models}
