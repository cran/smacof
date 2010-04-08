\name{smacofConstraint}
\alias{smacofConstraint}

\title{SMACOF Constraint}
\description{
SMACOF with constraints on the configuration
}
\usage{
smacofConstraint(delta, constraint = "linear", external, ndim = 2, weightmat = NULL, startconf = NULL,
metric = TRUE, ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{constraint}{Type of constraint: \code{"linear"}, \code{"unique"}, \code{"diagonal"}, or a user-specified function (see details)}
  \item{external}{Data frame or matrix with external covariates, or list for simplex and circumplex (see details)}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{startconf}{Optional matrix with starting values for configurations (see details)} 
  \item{metric}{If \code{FALSE} non-metric MDS is performed}
  \item{ties}{Tie specification for non-metric MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}
\details{The user can specify a function with the following arguments: configuration matrix with starting values, 
matrix V (based on the weight matrix, see package vignette), external scale matrix. 
The function must return a matrix of resulting configurations. 

A matrix with starting configurations can be specifiied. For \code{constraint = "linear"} it has to be of dimension (n x p). For \code{constraint = "unique"} it is typically of the form X = (Y|D) with D as (n x n) diagonal matrix and Y (n x p) configuration matrix. Hence X is of dimension (n x (n + p)). For \code{constraint = "diagonal"} it is typically of dimension  (n x q) where q is the number of columns of the external scale matrix (and thus number of dimensions). If \code{constraint} is user-specified the specification of \code{startconf} is mandatory.

The argument \code{external} allows for the specification of a covariate data frame (or matrix) of dimension (n x q). Alternatively, for simplex fitting the user can specify a list of the following structure: \code{external = list("simplex", dim2)} with \code{dim2} denoting the dimension of the simplex with dim2 < n. For a circumplex the list has to be of the following form: \code{external = list("circumplex", dim2, k1, k2)} with 1 <= k1 <= k2 <= n (see also examples section). k1 and k2 denote the circumplex width.
}
\value{
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{confdiss}{Configuration dissimilarities}
  \item{conf}{Matrix of final configurations}
  \item{stress.m}{stress value for metric MDS}
  \item{stress.nm}{stress value for non-metric MDS (if computed)}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{de Leeuw, J. \& Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 

de Leeuw, J., \& Heiser, W. (1980). Multidimensional scaling with restrictions on the configurations. In P. R. Krishnaiah (eds.), Multivariate Analysis V, pp. 501-522. North-Holland. 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}}
\examples{

## SMACOF with linear configuration constraints
data(kinshipdelta)
data(kinshipscales)
res.lin1 <- smacofConstraint(kinshipdelta, constraint = "linear", external = kinshipscales)

## SMACOF with unique constraints
res.unique <- smacofConstraint(kinshipdelta, constraint = "unique", external = kinshipscales)

## SMACOF with diagonal constraints
res.diag <- smacofConstraint(kinshipdelta, constraint = "diagonal", external = kinshipscales)

## Fitting a simplex with q = 14 (i.e., n-1), diagonal constraints
res.simp <- smacofConstraint(kinshipdelta, constraint = "diagonal", external = list("simplex", 14))

## Fitting a circumplex with q = 10, k1 = 2, k2 = 8, diagonal constraints
res.circ <- smacofConstraint(kinshipdelta, constraint = "diagonal", external = list("circumplex", 10, 2, 8))

}

\keyword{models}
