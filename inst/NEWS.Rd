\name{NEWS}
\title{News for Package \pkg{smacof}}
\encoding{UTF-8}

\section{Changes in Version 2.1-5}{
  \itemize{
    \item JSS reference updated.
    \item Principal axis transformation argument added to mds().
    \item Correlation between X and Y-hat added to Procrustes().
    \item sim2diss() fixed for missing similarities.
    }
} 


\section{Changes in Version 2.1-4}{
  \itemize{
    \item bugfix in unfolding() for user-specified starting solution.
    }
} 


\section{Changes in Version 2.1-3}{
  \itemize{
    \item stress0() is now an object of class "smacof".
    }
} 


\section{Changes in Version 2.1-2}{
  \itemize{
    \item Bug in smacofConstraint() for constraint.ties fixed. 
    }
} 

\section{Changes in Version 2.1-1}{
  \itemize{
    \item svm_mdsplot() updated.
    \item Nutshell vignettes added.
    }
}  

\section{Changes in Version 2.1-0}{
  \itemize{
    \item Alienation coefficient added to Procrustes().
    \item Procrustes() returns pairwise distances.
    \item jackknife() renamed to jackmds().
    \item Legend option for jakcknife added to plot.smacofJK().
    \item In smacofConstraint() "linear" constraint changed to "unrestricted".
    \item In unfolding() "matrix" changed to "unconditional". 
    \item omega default in unfolding() set to 1. 
    \item Guttman datasets added (Guttman1965, Guttman1991). 
    \item Second rectangle dataset added (area-shape). 
    \item svm_mdsplot() added for facet plotting based on SVM fit.
    \item Argument "scale" in biplotmds() added. 
    }
}


\section{Changes in Version 2.0-0}{
  \itemize{
    \item Color shading added to plot.icexplore().
    \item Fixed column coordinates version of vmu() implemented. 
    \item Vignette added.
    \item Help files for indvalues and PVQ40 data updated.
    \item method.dat = "rows" as default in permtest.smacofR(). 
    }
}

\section{Changes in Version 1.13-1}{
  \itemize{
    \item Slight acceleration in unfolding when using transformations.
    }
}

\section{Changes in Version 1.13-0}{
  \itemize{
    \item Pseudo-confidence ellipses added via confEllipse().
    \item Inverse MDS removed because it was buggy. Will be added later on again. 
    }
}

\section{Changes in Version 1.12-1}{
  \itemize{
    \item Extra sanity check added for fixed coordinates in unfolding added.
    \item Parallelization is added for rowconditional unfolding.
    }
}


\section{Changes in Version 1.12-0}{
  \itemize{
    \item External unfolding added.
    \item In case of NA's, proper weight matrix is created in unfolding.
    \item Vector unfolding model added. 
    }
}

\section{Changes in Version 1.11-0}{
  \itemize{
    \item Circular unfolding for all unfolding versions implemented.
    }
}

\section{Changes in Version 1.10-8}{
  \itemize{
    \item variable names in morsescales corrected.
    \item adjust argument in plot.driftvec added.
    }
}

\section{Changes in Version 1.10-7}{
  \itemize{
    \item sanity check in smacofConstraint for diagonal C. 
    \item stability coefficient for bootstrap added. 
    }
}

\section{Changes in Version 1.10-6}{
  \itemize{
    \item print.imds() for inverseMDS() added.  
    \item Gravity computation slightly modified.
    \item Vaziri-Xu dataset added. 
    \item Facial expression data added.  
    \item Legend arguments for Procrustes plots added.
    }
}


\section{Changes in Version 1.10-5}{
  \itemize{
    \item some fixes in the output structure of inverseMDS().  
    \item biplotmds() and plot.mdsbi() fixed for single covariate vectors.
    \item Additional conversion options added to sim2diss().
    }
}

\section{Changes in Version 1.10-4}{
  \itemize{
    \item verbose output in uniscale(). 
    \item transform() exported. 
    }
}


\section{Changes in Version 1.10-3}{
  \itemize{
    \item plot.smacofboot() updated.
    }
}

\section{Changes in Version 1.10-2}{
  \itemize{
    \item unfolding() updated with relaxed update that often reduces the iterations by a factor 2).
    \item unfolding() updated with an overal rescaling to the coordinates, distances, and d-hats.
    }
}


\section{Changes in Version 1.10-1}{
  \itemize{
    \item smacofRect() updated to do nonlinear unfolding added matrix and row conditional).
    \item prefscal() is an alias for unfolding() and smacofRect()
    \item plot.smacofR() updated to show row transformations by smacofR and show column labels.
    \item mspline transformation in mds() and smacofSym() now includes an intercept.
    }
}

\section{Changes in Version 1.10-0}{
  \itemize{
    \item smacofRect updated for nonlinear unfolding added (only matrix conditional).
    \item plot.smacofR updated to show transformations by smacofR.
}
}


\section{Changes in Version 1.9-7}{
  \itemize{
    \item Procrustes help file updated.
    \item Bug in SPP fixed. 
}
}


\section{Changes in Version 1.9-7}{
  \itemize{
    \item Some help files fixed.
}
}

\section{Changes in Version 1.9-6}{
  \itemize{
    \item Registration for C function calls fixed.
}
}

\section{Changes in Version 1.9-5}{
  \itemize{
    \item Plot option for summed distances for \code{smacofIndDiff} models added. 
}
}

\section{Changes in Version 1.9-4}{
  \itemize{
    \item Fast optimal scaling in \code{smacofSphere} including msplines. 
    \item Primal algorithm in \code{smacofSphere} accelerated.
    \item Bug for smacofSP Shepard plots fixed. 
    \item \code{cex} argument for unfolding configuration plot activated.
    \item Thigmophobe label option in \code{plot.mdsbi} added. 
    \item Multiple panels for \code{plot.smacofID} improved (thanks to Lew Harvey for the fix). 
}
}


\section{Changes in Version 1.9-3}{
  \itemize{
    \item \code{type = "ratio"} for all models restored.
}
}

\section{Changes in Version 1.9-2}{
  \itemize{
    \item Shepard plots for IndDiff-models and unfolding fixed. For IndDiff, a Shepard plot for each individual is produced.
    \item New \code{scale} argument added for scaling the arrow endpoints (coefficients) in \code{plot.mdsbi}.
}
}

\section{Changes in Version 1.9-1}{
  \itemize{
    \item Canadian newspaper dataset added (word co-occurrences). 
    \item \code{stress0()} function added to compute stress for a particular (starting) configuration.
    \item Bug in \code{smacofIndDiff} with missing input dissimilarities fixed. 
    \item \code{icExplore()} function incl. print/plot S3 added for initial configuration exploration.
    \item \code{KIPT} dataset added.
}
}

\section{Changes in Version 1.9-0}{
  \itemize{
    \item Crime dataset added. 
    \item \code{smacofConstraint()} updated to allow for regional restrictions (diagonal case). 
    \item \code{type = "interval"} as new default in \code{mds()}. 
    \item \code{confdist} instead of \code{confdiss} as output value in \code{mds()}.
    \item Labels in \code{plot.smacof()} are now positioned using \code{thigmophobe.labels()}. 
    \item Rectangle data corrected, help file updated. 
    \item Convex hull option added to \code{plot.smacof}.
    \item Age and gender added as attributes to \code{PVQ40} dataset.
    \item New \code{what} argument in \code{plot.smacofR} for row/column/joint plots added. 
    \item \code{bootmds()} incl. plot and print S3 added for bootstrapping MDS solutions.
    \item Help file for \code{perception} data updated.
}  
}



\section{Changes in Version 1.8-16}{
  \itemize{
    \item Bug in \code{smacofConstraint()} removed: now returns proper  matrix \code{C} when \code{constraint = "diagonal"}. 
}  
}


\section{Changes in Version 1.8-15}{
  \itemize{
    \item \code{plot.smacofJK()} improved. 
}  
}

\section{Changes in Version 1.8-14}{
  \itemize{
    \item Bug fix in \code{plot.smacofID()}. 
}  
}

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
