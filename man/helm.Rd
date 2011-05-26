\name{helm}
\alias{helm}
\docType{data}
\title{Helm's color data}
\description{This dataset contains dissimilarity data for individual difference scaling from an experiment carried out by Helm (1959). 
}
\usage{data(helm)}
\format{List containing objects of class \code{dist}
}
\details{A throughout description of the experiment can be found in Borg and Groenen (2005, p. 451) with the corresponding Table 21.1. containing distance estimates for color pairs. There were 14 subjects that rated the similarity of colors, 2 of whom replicated the experiment. 10 subjects have a normal color vision (labelled by N1 to N10 in our list object), 4 of them are red-green deficient in varying degrees. In this dataset we give the dissimilarity matrices for each of the subjects, including the replications. They are organized as a list of length 16 suited for \code{smacofIndDiff} computations.

The authors thank Michael Friendly and Phil Spector for data preparation. 
}

\references{Helm, C. E. (1959). A multidimensional ratio scaling analysis of color relations. Technical Report, Princeton University and Educational Testing Service. Princeton, NJ. 

Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling: Theory and Applications (2nd edition). New York: Springer. 
}
\examples{
data(helm)
length(helm)
helm
}
\keyword{datasets}
