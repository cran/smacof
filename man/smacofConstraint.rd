\name{smacofConstraint}
\alias{smacofConstraint}

\title{SMACOF Constraint}
\description{
SMACOF with internal constraints on the configurations.
}
\usage{
smacofConstraint(delta, constraint = "unrestricted", external, ndim = 2, 
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
  \item{constraint}{Type of constraint: \code{"unrestricted"}, \code{"unique"}, \code{"diagonal"}, or a user-specified function (see details)}
  \item{external}{Data frame or matrix with external covariates, or list for simplex and circumplex (see details)}
  \item{ndim}{Number of dimensions}
  \item{type}{MDS type: \code{"interval"}, \code{"ratio"}, \code{"ordinal"} (nonmetric MDS), or \code{"mspline"}}
  \item{weightmat}{Optional matrix with dissimilarity weights}
  \item{init}{Optional matrix with starting values for configurations. If \code{NULL} random starts are used (see details).} 
  \item{ties}{Tie specification for non-metric MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{spline.degree}{Degree of the spline for \code{"mspline"} MDS type}
  \item{spline.intKnots}{Number of interior knots of the spline for \code{"mspline"} MDS type}
  \item{constraint.type}{Transformation for \code{external} covariates: \code{"ratio"},
                 \code{"interval"}, \code{"ordinal"}, \code{"spline"}, or         
                 \code{"mspline"})}
  \item{constraint.ties}{Tie specification for \code{external} covariates with 
                \code{constraint.type = "ordinal"}: \code{"primary"}, 
                \code{"secondary"}, or \code{"tertiary"}}
  \item{constraint.spline.intKnots}{Number of interior knots for \code{external} covariates
               with \code{constraint.type = "spline"} or \code{"mspline"}}
  \item{constraint.spline.degree}{Degree of the spline for \code{external} covariates
               with \code{constraint.type = "spline"} or \code{"mspline"}}
  
}
\details{The argument \code{external} is mandatory to specify and requires a data frame (or matrix) of dimension (n x q). Alternatively, for simplex fitting the user can specify a list of the following structure: \code{external = list("simplex", dim2)} with \code{dim2} denoting the dimension of the simplex with dim2 < n. For a circumplex fitting, the list has to be of the following form: \code{external = list("circumplex", dim2, k1, k2)} with \eqn{1 \leq k1 \leq k2 \leq n} (see also examples section). k1 and k2 denote the circumplex width.

In constraint smacof, the configuration matrix \eqn{X} is subject to a constraint based on the external scales (predictors \eqn{Z} specified using \code{external}) of the following linear form: \eqn{X = ZC}. The type of constraint in \eqn{C} can be specified using the \code{constraint} argument. We provide the following standard setting: 

For \code{constraint = "unrestricted"}, \eqn{C} is unrestricted. Note that \code{"linear"} still works as well for backward compatibility.  

The same for \code{constraint = "diagonal"} where \eqn{X} needs to be of dimension \eqn{(n x q)} where \eqn{q} is the number of columns of the external scale matrix (and thus number of dimensions). Here, \eqn{C} is restricted to be diagonal. 

For \code{constraint = "unrestricted"} or \code{"diagonal"}, the external covariates \eqn{Z} can be optimally transformed as specified by \code{constraint.type}. Choosing the number of covariates equal to the number of dimensions together with \code{constraint.type = "ordinal"}, \code{constraint.ties = "primary"} will effectively restrict the configuration to parallel regions defined by the categories of the covariates. Note that missing values of the covariates are estimated by the model.

For \code{constraint = "unique"} we get the Bentler-Weeks uniqueness model. Hence \eqn{X} is of dimension \eqn{(n x (n + p))}. This implies that we fit a certain number of dimensions p and, in addition we extract n additional dimensions where each object is scored on a separate dimension. More technical details can be found in the corresponding JSS article (reference see below).

In addition, the user can specify his own constraint function with the following arguments: configuration matrix with starting values (\code{init}) (mandatory in this case), matrix \eqn{V} (\code{weightmat}; based on the weight matrix, see package vignette), external scale matrix (\code{external}). The function must return a matrix of resulting configurations. 

If no starting configuration is provided, a random starting solution is used. In most applications, this is not a good idea in order to find a well fitting model. The user can fit an exploratory MDS using \code{mds()} first, and use the resulting configurations as starting configuration for \code{smacofConstraint()}. Alternatively, if the user has starting configurations determined by some underlying theory, they can be used as well. 

}
\value{
  \item{delta}{Observed dissimilarities}
  \item{obsdiss}{Observed dissimilarities, normalized}
  \item{confdist}{Configuration dissimilarities}
  \item{conf}{Matrix of final configurations}
  \item{C}{Matrix with restrictions}
  \item{stress}{Stress-1 value}
  \item{spp}{Stress per point}
  \item{resmat}{Matrix with squared residuals}
  \item{rss}{Residual sum-of-squares}
  \item{weightmat}{Weight matrix}
  \item{ndim}{Number of dimensions}
  \item{extvars}{List for each external covariate with a list of class \code{"optscal"}}
  \item{init}{Starting configuration}
  \item{model}{Type of smacof model}
  \item{niter}{Number of iterations}
  \item{nobj}{Number of objects}
}
\references{
De Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{https://www.jstatsoft.org/v31/i03/} 

Mair, P., Groenen, P. J. F., & De Leeuw, J. (2020). More on multidimensional scaling and unfolding in
R: smacof version 2. Journal of Statistical Software, Forthcoming.

De Leeuw, J., & Heiser, W. (1980). Multidimensional scaling with restrictions on the configurations. In P. R. Krishnaiah (eds.), Multivariate Analysis V, pp. 501-522. North-Holland. 

Borg, I., & Lingoes, J. C. (1980). A model and algorithm for multidimensional scaling with external constraints on the distances. Psychometrika, 45, 25-38.
}


\seealso{\code{\link{smacofSym}}, \code{\link{smacofRect}}, \code{\link{smacofIndDiff}}, \code{\link{smacofSphere}}}
\examples{


## theoretical grid restrictions (rectangles; keep covariate ties tied)
fit.rect1 <- mds(rectangles, type = "ordinal", init = rect_constr) 
fit.rect2 <- smacofConstraint(rectangles, type = "ordinal", ties = "secondary",
                        constraint = "diagonal", init = fit.rect1$conf, 
                        external = rect_constr, constraint.type = "ordinal")
plot(fit.rect2)

## regional restrictions morse code data (signal length, strength)
fitMorse1 <- mds(morse, type = "ordinal")
fitMorse1
fitMorse2 <- smacofConstraint(morse, type = "ordinal", constraint = "unrestricted",
                              external = morsescales[,2:3], 
                              constraint.type = "ordinal", 
                              init = fitMorse1$conf)
fitMorse2
plot(fitMorse2)

## facial expression data I (axial restriction, C diagonal)
Delta <- FaceExp
attr(Delta, "Labels") <- NULL            
fitFace <- mds(Delta, type = "ordinal")   ## starting solution
Z <- FaceScale[, c(1,3)]                  ## external variables
fitFaceC1 <- smacofConstraint(Delta, type = "ordinal", 
  constraint = "diagonal", external = Z, constraint.type = "ordinal", 
  init = fitFace$conf)
fitFaceC1$C 
plot(fitFaceC1, xlab = "Pleasant-Unpleasant", ylab = "Tension-Sleep", 
  main = "Face Expression (Diagonal Restriction)")

## facial expression data II (C unrestricted)
fitFaceC3 <- smacofConstraint(Delta, type = "ordinal", 
  constraint = "unrestricted", external = Z, constraint.type = "ordinal", 
  init = fitFace$conf)
fitFaceC3$C   
plot(fitFaceC3, main = "Face Expression (C Unrestricted, Ordinal Transformation)")
}

\keyword{multivariate}
