\name{smacofConstraint}
\alias{smacofConstraint}

\title{SMACOF Constraint}
\description{
SMACOF with constraints on the configurations.
}
\usage{
smacofConstraint(delta, constraint = "linear", external, ndim = 2, 
weightmat = NULL, init = NULL, metric = TRUE, ties = "primary", verbose = FALSE, 
modulus = 1, itmax = 1000, eps = 1e-6)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{constraint}{Type of constraint: \code{"linear"}, \code{"unique"}, \code{"diagonal"}, or a user-specified function (see details)}
  \item{external}{Data frame or matrix with external covariates, or list for simplex and circumplex (see details)}
  \item{ndim}{Number of dimensions}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Optional matrix with starting values for configurations} 
  \item{metric}{If \code{FALSE} non-metric MDS is performed}
  \item{ties}{Tie specification for non-metric MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
}
\details{The argument \code{external} is mandatory and allows for the specification of a covariate data frame (or matrix) of dimension (n x q). Alternatively, for simplex fitting the user can specify a list of the following structure: \code{external = list("simplex", dim2)} with \code{dim2} denoting the dimension of the simplex with dim2 < n. For a circumplex fitting, the list has to be of the following form: \code{external = list("circumplex", dim2, k1, k2)} with 1 <= k1 <= k2 <= n (see also examples section). k1 and k2 denote the circumplex width.

In constraint smacof, the configuration matrix is subject of a constraint based on the external scales (predictors). This constraint can be specified using the \code{constraint} argument. We provide the following standard setting: 

For \code{constraint = "linear"} the configurations X are decomposed linearly, i.e. X = ZC where Z is the known predictor matrix specified using \code{external}. 

The same for \code{constraint = "diagonal"} where X needs to be of dimension  (n x q) where q is the number of columns of the external scale matrix (and thus number of dimensions). Here, C is restricted to be diagonal. 

For \code{constraint = "unique"} we get the Bentler-Weeks uniqueness model. Hence X is of dimension (n x (n + p)). This implies that we fit a certain number of dimensions p and, in addition we extract n additional dimensions where each object is scored on a separate dimension. More technical details can be found in the corresponding JSS article (reference see below).

In addition, the user can specify his own constraint function with the following arguments: configuration matrix with starting values (\code{init}) (mandatory in this case), matrix V (\code{weightmat}; based on the weight matrix, see package vignette), external scale matrix (\code{external}). The function must return a matrix of resulting configurations. 


}
\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{confdiss}{Configuration dissimilarities}
  \item{conf}{Matrix of final configurations}
  \item{stress.m}{stress value for metric MDS}
  \item{stress.nm}{stress value for non-metric MDS (if computed)}
  \item{spp}{Stress per point}
  \item{ndim}{Number of dimensions}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{de Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 

de Leeuw, J., & Heiser, W. (1980). Multidimensional scaling with restrictions on the configurations. In P. R. Krishnaiah (eds.), Multivariate Analysis V, pp. 501-522. North-Holland. 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere.primal}}, \code{\link{smacofSphere.dual}}}
\examples{

## SMACOF with linear configuration constraints
data(kinshipdelta)
data(kinshipscales)
res.lin1 <- smacofConstraint(kinshipdelta, constraint = "linear", external = kinshipscales)

## X = ZC decomposition
Z <- as.matrix(kinshipscales)     ## matrix Z with external constraints
X <- res.lin1$conf                ## resulting configurations X
C <- solve(t(Z)\%*\%Z)\%*\%t(Z)\%*\%X   ## compute C out of X = ZC
C
Z%*%C                             ## check: should be equal to X   

## SMACOF with diagonal constraints
res.diag <- smacofConstraint(kinshipdelta, constraint = "diagonal", 
external = kinshipscales, ndim = 3)
X <- res.diag$conf                ## resulting configurations X
C <- solve(t(Z)\%*\%Z)\%*\%t(Z)\%*\%X   ## compute C out of X = ZC
round(C, 3)
Z\%*\%C                             ## check: should be equal to X   

## SMACOF with unique constraints (Bentler-Weeks model)
res.unique <- smacofConstraint(kinshipdelta, constraint = "unique", 
external = kinshipscales)

## Fitting a simplex with q = 4 (i.e., n-1), diagonal constraints
res.simp <- smacofConstraint(kinshipdelta, constraint = "diagonal", 
external = list("simplex", 4), ndim = 3)
plot3d(res.simp)

## Fitting a circumplex with q = 3, k1 = 1, k2 = 2, diagonal constraints
res.circ <- smacofConstraint(kinshipdelta, constraint = "diagonal", 
external = list("circumplex", 3, 1, 2), ndim = 3)
plot3d(res.circ)


}

\keyword{models}
