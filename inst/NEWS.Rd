\name{NEWS}
\title{News for Package \pkg{smacof}}
\encoding{UTF-8}

\section{Changes in Version 1.8-13}{
  \itemize{
    \item Aspect ratio 1 added as default for the bubble plot option in \code{plot.smacof()}. 
}  
}


\section{Changes in Version 1.8-12}{
  \itemize{
    \item Function dissWeight added to create weights as a function of the dissimilarities. 
    Removed dissimilarities with weights 0 from Shepard plot. multistart.smacof added. 
}  
}


\section{Changes in Version 1.8-11}{
  \itemize{
    \item In mds, smacofSym, smacofConstraint, smacofIndDiff, and smacofSphere the order of computations in update changed. This update improves the speed of computations for large data when using weightmat.
}  
}


\section{Changes in Version 1.8-10}{
  \itemize{
    \item Hue change in jackknife plot.
}  
}


\section{Changes in Version 1.8-9}{
  \itemize{
    \item Some fixes in the permutation tests.
}  
}

\section{Changes in Version 1.8-6}{
  \itemize{
    \item Permutation tests based on permutations of the original data matrix added.
}  
}

\section{Changes in Version 1.8-5}{
  \itemize{
    \item \code{biplot} function for creating MDS biplots added.
}  
}

\section{Changes in Version 1.8-3}{
  \itemize{
    \item Optimization in \code{smacofConstraint} improved (monotone regression in two directions).  
    \item Some plot cosmetics for \code{smacofIndDiff} and \code{Procrustes} calls. 
    \item Transformations in \code{smacofIndDiff} use C implementation (mspline can be done as well). 
}  
}


\section{Changes in Version 1.8-0}{
 \itemize{
    \item Drift vector model added.
    \item Datasets added. 
    \item \code{smacofConstraint} allows now for data frame and matrices for external variables, sanity checks added.
    \item Verbose output now says "stress (raw)".
    \item init default argument is now \code{"torgerson"} (before it was \code{NULL} which used Torgerson scaling anyway).
    \item Bug for optimal scaling on external variables in \code{smacofConstraint} removed. 
    \item SPP-plot for unfolding models fixed. 
}
}

\section{Changes in Version 1.7-0}{
 \itemize{
    \item Some datasets added. 
    \item Stress-per-point (spp) is now returned as percentage stress contribution including the squared error matrix (\code{resmat}). 
    \item If there are missing input dissimilarities, the corresponding weight matrix is created automatically.
    \item Sanity check for negative input dissimilarities. 
    \item Some plot cosmetrics.
    \item Output value structure for dual algorithm in \code{smacofSphere} changed, circle added in configuration plot. 
    \item Function \code{randomstress} for computing (average) stress values based on random input dissimilarities. 
    \item Unidimensional scaling added (all permutations): \code{uniscale}. 
    \item Procrustes added: \code{Procrustes} including print and plot methods.    
    \item Function \code{torgerson} (for classical scaling) is now exported. 
    \item Bug in \code{smacofConstraint} for external missing values removed.
}
}
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
       \item Ordinal transformations are coded in C.
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
