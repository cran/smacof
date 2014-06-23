\name{smacofConstraint}
\alias{smacofConstraint}

\title{SMACOF Constraint}
\description{
SMACOF with internal constraints on the configurations.
}
\usage{
smacofConstraint(delta, constraint = "linear", external, ndim = 2, 
                 type = c("ratio", "interval", "ordinal", "mspline"), weightmat = NULL,
                 init = NULL, ties = "primary", verbose = FALSE, modulus = 1, 
                 itmax = 1000, eps = 1e-6, spline.intKnots = 4, spline.degree = 2, 
                 constraint.type = c("ratio", "interval", "ordinal", "spline", 
                 "mspline"), constraint.ties = "primary", 
                 constraint.spline.intKnots = 2, constraint.spline.degree = 2)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{Either a symmetric dissimilarity matrix or an object of class \code{"dist"}}
  \item{constraint}{Type of constraint: \code{"linear"}, \code{"unique"}, \code{"diagonal"}, or a user-specified function (see details)}
  \item{external}{Data frame or matrix with external covariates, or list for simplex and circumplex (see details)}
  \item{ndim}{Number of dimensions}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"}, or \code{"mspline"} (nonmetric MDS)}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Optional matrix with starting values for configurations} 
  \item{ties}{Tie specification for non-metric MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
  \item{constraint.type}{Transformation for \code{external} covariates: \code{"ratio"},
                 \code{"interval"}, \code{"ordinal"}, \code{"spline"}, \code{"spline"}, or         
                 \code{"mspline"})}
  \item{constraint.ties}{Tie specification for \code{external} covariates with 
                \code{constraint.type = "ordinal"}: \code{"primary"}, 
                \code{"secondary"}, or \code{"tertiary"}}
  \item{constraint.spline.intKnots}{Number of interior knots for \code{external} covariates
               with \code{constraint.type = "spline"} or \code{"mspline"}}
  \item{constraint.spline.degree}{Degree of the spline for \code{external} covariates
               with \code{constraint.type = "spline"} or \code{"mspline"}}
  
}
\details{The argument \code{external} is mandatory and allows for the specification of a covariate data frame (or matrix) of dimension (n x q). Alternatively, for simplex fitting the user can specify a list of the following structure: \code{external = list("simplex", dim2)} with \code{dim2} denoting the dimension of the simplex with dim2 < n. For a circumplex fitting, the list has to be of the following form: \code{external = list("circumplex", dim2, k1, k2)} with 1 <= k1 <= k2 <= n (see also examples section). k1 and k2 denote the circumplex width.

In constraint smacof, the configuration matrix is subject of a constraint based on the external scales (predictors). This constraint can be specified using the \code{constraint} argument. We provide the following standard setting: 

For \code{constraint = "linear"} the configurations X are decomposed linearly, i.e. X = ZC where Z is the known predictor matrix specified using \code{external}. 

The same for \code{constraint = "diagonal"} where X needs to be of dimension  (n x q) where q is the number of columns of the external scale matrix (and thus number of dimensions). Here, C is restricted to be diagonal. 

For \code{constraint = "linear"} or \code{"diagonal"}, the external covariates Z can be optimally transformed as specified by \code{constraint.type}. Choosing the number of covariates equal to the number of dimensions together with \code{constraint.type = "ordinal"}, \code{constraint.ties = "primary"} will effectively restrict the configuration to parallel regions defined by the categories of the covariates. Note that missing values of the covariates are estimated by the model.

For \code{constraint = "unique"} we get the Bentler-Weeks uniqueness model. Hence X is of dimension (n x (n + p)). This implies that we fit a certain number of dimensions p and, in addition we extract n additional dimensions where each object is scored on a separate dimension. More technical details can be found in the corresponding JSS article (reference see below).

In addition, the user can specify his own constraint function with the following arguments: configuration matrix with starting values (\code{init}) (mandatory in this case), matrix V (\code{weightmat}; based on the weight matrix, see package vignette), external scale matrix (\code{external}). The function must return a matrix of resulting configurations. 


}
\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{confdiss}{Configuration dissimilarities}
  \item{conf}{Matrix of final configurations}
  \item{C}{Matrix with restrictions}
  \item{stress}{Stress-1 value}
  \item{spp}{Stress per point}
  \item{ndim}{Number of dimensions}
  \item{extvars}{List for each external covariate with a list of class \code{"optscal"}}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{de Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/} 

de Leeuw, J., & Heiser, W. (1980). Multidimensional scaling with restrictions on the configurations. In P. R. Krishnaiah (eds.), Multivariate Analysis V, pp. 501-522. North-Holland. 
}
\author{Jan de Leeuw and Patrick Mair}

\seealso{\code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere}}}
\examples{

## SMACOF with linear configuration constraints
data(kinshipdelta)
data(kinshipscales)
res.lin1 <- smacofConstraint(kinshipdelta, constraint = "linear", external = kinshipscales)

## X = ZC decomposition
X <- res.lin1$conf                ## resulting configurations X
Z <- as.matrix(kinshipscales) 
C <- res.lin1$C
C
Z\%*\%C                             ## check: should be equal to X   

## SMACOF with diagonal constraints
res.diag <- smacofConstraint(kinshipdelta, constraint = "diagonal", 
external = kinshipscales, ndim = 3)
X <- res.diag$conf                ## resulting configurations X
C <- res.diag$C
Z\%*\%C                             ## check: should be equal to X   

## SMACOF with parallel regional constraints 
data(morse)
data(morsescales)
res.unc <- smacofSym(morse,  type = "ordinal", ties = "primary", ndim = 2)
res.parreg <- smacofConstraint(morse, type = "ordinal", ties = "primary", 
constraint = "linear", external = as.data.frame(morsescales[,2:3]), 
constraint.type = "ordinal", constraint.ties = "primary", ndim = 2, 
init = res.unc$conf)
X <- res.parreg$conf                ## resulting configurations X
C <- res.parreg$C                   ## Matrix of weights
Z <- res.parreg$external            ## optimally transformed external variables
Z\%*\%C                               ## check: should be equal to X   


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
