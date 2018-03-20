\name{transform}
\alias{transform}
\alias{transPrep}

\title{Internal Dissimilarity Transformation}
\description{Utility functions for optimal scaling calls (used internally)
}
\usage{
transPrep(x, trans = "ordinals", spline.intKnots = 4, spline.degree = 2, missing = "none")
transform(Target, x, w = rep(1,length(x$x)), normq = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{Target}{unconstrained vector of target values}
   \item{x}{object of type optScal}
   \item{w}{vector non-negative weights}
   \item{normq}{sum of squares normalization}
   \item{trans}{type of transformation}
   \item{spline.intKnots}{interior spline knots}
   \item{spline.degree}{spline degree}
   \item{missing}{missing treatment}
}  
