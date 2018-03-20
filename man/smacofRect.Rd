\name{unfolding}
\alias{unfolding}
\alias{smacofRect}
\alias{prefscal}

\title{Nonmetric unfolding}
\description{Variant of smacof for rectangular matrices (typically ratings, preferences) that allows for nonmetric transformations. Also known as nonmetric unfolding.
}
\usage{
unfolding(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
         conditionality = c("matrix", "row"), lambda = 0.5, omega = 0.1, 
         circle = c("none", "row", "column"), weightmat = NULL, init = NULL, 
         ties = c("primary", "secondary"), verbose = FALSE, relax = TRUE, 
         itmax = 10000, eps = 1e-6, spline.degree = 2, spline.intKnots = 2)

smacofRect(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
         conditionality = c("matrix", "row"), lambda = 0.5, omega = 0.1, 
         circle = c("none", "row", "column"), weightmat = NULL, init = NULL, 
         ties = c("primary", "secondary"), verbose = FALSE, relax = TRUE, 
         itmax = 10000, eps = 1e-6, spline.degree = 2, spline.intKnots = 2)

prefscal(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
         conditionality = c("matrix", "row"), lambda = 0.5, omega = 0.1, 
         circle = c("none", "row", "column"), weightmat = NULL, init = NULL, 
         ties = c("primary", "secondary"), verbose = FALSE, relax = TRUE, 
         itmax = 10000, eps = 1e-6, spline.degree = 2, spline.intKnots = 2)
}

\arguments{
  \item{delta}{Data frame or matrix of preferences, ratings, dissimilarities.}
  \item{ndim}{Number of dimensions}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"}, or \code{"mspline"}}
  \item{conditionality}{A single transformations are applied for the entire matrix \code{"matrix"} or for each row separately \code{"row"}}
  \item{lambda}{Penalty strength balancing the loss contribution of stress and the penalty}
  \item{omega}{Penalty width determines for what values of the variation coefficient the penalty should become active.}
  \item{circle}{If \code{"column"}, the column configurations are restricted to be on a circle, if \code{"row"}, row configurations are on a circle, if \code{"none"}, there are no restrictions on row and column configurations}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Matrix with starting values for configurations (optional)}
  \item{ties}{Tie specification for \code{ordinal} transformations: \code{primary} unties the ties and \code{secondary} keeps the ties tied}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization after 100 iterations. It tends to reduce the number of iterations by a factor 2.}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for an \code{"mspline"} transformation}
  \item{spline.intKnots}{Number of interior knots of the spline for a \code{"mspline"} transformation}
}
\details{Unfolding tries to match a rectangular matrix \code{delta} of dissimilarities between row and column objects by Euclidean distances between row and column points. Badness of fit is measured by raw Stress as the sum of squared differences between \code{delta} and the Euclidean distances. Instead of dissimilarities optimal transformations (dhats) can be found. The dhats should be a function of the original \code{delta} restricted to be \code{"ratio"}, \code{"interval"}, \code{"ordinal"}, or \code{"mspline"}. These transformations can be the same for the entire matrix (\code{conditionality = "matrix"}) of data or different per row (\code{conditionality = "row"}). To avoid a degenerate solution with all dhats and distances equal to 1, the prefscal penalty is is used. A penalty is added based on the variation coefficient of the dhats (mean dhat divided by the standard deviation of the dhats). The penalty width (\code{omega}) weights the penalty and determines from what value of the variation coefficient of the dhats the penalty should become active. The penalty strength (\code{lambda}) is needed to ensure that the penalty can be strong enough. 

Creates an object of class \code{smacofR}. 
}

\value{
  \item{obsdiss}{Observed dissimilarities, corresponds to \code{delta}}
  \item{confdist}{Configuration dissimilarities}
  \item{dhat}{Matrix with optimal transformation of size \code{delta}}
  \item{iord}{List of size 1 for matrix conditional and size \code{nrow(delta)} for row conditional with the index that orders the dhats. Needed for the Shepard plot}
  \item{conf.row}{Matrix of final row configurations}
  \item{conf.col}{Matrix of final column configurations}
  \item{stress}{Final, normalized stress value}
  \item{pstress}{Penalized stress value (the criterion that is minimized)}
  \item{spp.row}{Stress per point, rows}
  \item{spp.col}{Stress per point, columns}
  \item{congvec}{Vector of congruency coefficients}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nind}{Number of individuals (rows)}
  \item{trans}{Transformation}
  \item{conditionality}{Conditionality of the transformation}
  \item{nobj}{Number of objects (columns)}
}
\references{de Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 

Busing, F. M. T. A., Groenen, P. J. F., & Heiser, W. J. (2005). Avoiding degeneracy in multidimensional unfolding by penalizing on the coefficient of variation. Psychometrika, 70, 71-98.
}

\author{Patrick Groenen, Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{plot.smacof}}, \code{\link{smacofConstraint}}, \code{\link{smacofSym}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere}}}

\examples{
## Ratio unfolding
res <- unfolding(breakfast)
res

## various configuration plots
plot(res)
plot(res, type = "p", pch = 25)
plot(res, type = "p", pch = 25, col.columns = 3, 
     label.conf.columns = list(label = TRUE, pos = 3, col = 3), 
     col.rows = 8, label.conf.rows = list(label = TRUE, pos = 3, col = 8))

## Shepard plot
plot(res, "Shepard")

## Stress decomposition chart
plot(res, "stressplot")

\dontrun{
## Ordinal unfolding row conditional
## Note that ordinal unfolding may need many iterations (several thousands)
res <- unfolding(breakfast, type = "ordinal", conditionality = "row", omega = 0.1, itmax = 3000)
res
## Shepard plot
plot(res, "Shepard")
plot(res)
}

}

\keyword{models}
