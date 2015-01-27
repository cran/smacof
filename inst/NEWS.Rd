\name{NEWS}
\title{News for Package \pkg{smacof}}
\encoding{UTF-8}

\section{Changes in Version 1.6-0}{
 \itemize{
    \item Bug in \code{smacofIndDiff} with initial starting configurations fixed. Stress is now normalized to Stress-1.
    \item Several wrapper functions added: \code{mds}, \code{unfolding}, \code{indscal}, \code{idioscal}. 
}
}
\section{Changes in Version 1.5-0}{
 \itemize{
    \item Recoding of transformations:
    \itemize{
       \item Ordinal transformations are coded in C++ and therefore much faster.
       \item Added monotone spline transformations (\code{type = mspline} with number of interior knots \code{spline.intKnots} and degree of the spline \code{spline.degree = 2}.       
    }
    \item \code{smacofConstraint} allows for transformations of the external variables through \code{constraint.type} (\code{"ratio", "interval", "ordinal", "spline", "spline",} or \code{"mspline"}). Missing values (\code{NA}) are estimated by the model. Choosing the number of external variables equal to the number of dimensions together with \code{constraint.type = "ordinal", constraint.ties = "primary"} will effectively restrict the configuration to parallel regions defined by the categories of the external variables.
    \item Stress values reported by \code{verbose = TRUE} during iterations of \code{smacofSym}, \code{smacofConstraint}, and \code{smacofRect} are normalized Stress values.
 }
}

\section{Changes in Version 1.4-0}{
  \itemize{
    \item Gravity model for dissimilarity computation based on co-occurrences added: see \code{gravity()}. 
    \item The \code{circle} argument in \code{smacofRect()} allows for circular restrictions on the configurations.
    \item All smacof functions now return Stress-1 (before it was squared Stress-1); the value return by the function is now called \code{stress}. 
    \item The \code{metric} argument is deprecated. Now \code{type} should be used (either \code{"ordinal"}, \code{"interval"}, or \code{"ratio"})
    for pretty much all smacof functions. 
    \item \code{smacofSphere.primal()} and \code{smacofSphere.dual()} are now merged into a single function \code{smacofSphere()} with a corresponding \code{algorithm} argument. 
    \item \code{smacofConstraint()} now returns weight matrix \code{C}.
    \item Plot updates:
    \itemize{
       \item For unfolding (\code{smacofRect()}) configurations: \code{joint = TRUE} as default. 
       \item Configuration plot looks a bit fancier now (labels and points as default), also aspect ratio = 1 as default.
       \item A note in plot help file added that in the stress plot the SPP from a various smacof functions is rescaled to percentages (which sum up to 100). 
    }
    \item Datasets:
    \itemize{
       \item Labels for rectangle dataset now provided in English.
       \item Dataset on Republican statements added (document-term matrix).  
    }
  }
}
