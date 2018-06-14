\name{VaziriXu}
\alias{VaziriXu}
\docType{data}
\title{Visual Object Representations}
\description{Contains two similarity matrices related to an experiment on visual object represenations. 

similarities (correlations) of 54 OCP (see O'Reilly, Chatman, and Caldwell, 1991) items. The last three columns contain the facet assigned by Bilsky and Jehn (2002) as well as the external variables for regional restrictions. 
}
\usage{data(VaziriXu)}
\format{List of two similarity matrices from two experimental conditions: artificial and real object categories.

V1-V4: early visual areas

LO: lateral occipital regions

VOT: ventral occipitotemporal regions

V3A, V3B, IPS0-IPS4: topographic regions along the intraparietal sulcus

Inferior IPS: inferior intraparietal sulcus

Superior IPS: superior intraparietal sulcus
}

\references{Vaziri-Pashkam M, Xu Y (2018). An information-driven two-pathway characterization of occipito-temporal and posterior parietal visual object representations. Cerebral Cortex, forthcoming.
}
\examples{
vx1 <- sim2diss(VaziriXu[[1]])
vx2 <- sim2diss(VaziriXu[[2]])
}
\keyword{datasets}
