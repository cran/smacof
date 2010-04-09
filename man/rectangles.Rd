\name{rectangles}
\alias{rectangles}
\alias{rect_constr}
\docType{data}
\title{Rectangles}
\description{These data are based on an experiment by Borg and Leutner (1983). They constructed rectangles based on a specific grid structure. The size of the rectangles is given in the constraints matrix. In total, we have 16 rectangles. 21 subjects had to rate twice the similarity of each pair of rectangles on a scale from 0 (identical) to 9 (very different). The dist object contains the average ratings which are dissimilarities. 
}

\usage{data(rectangles)
data(rect_constr)
}

\format{The rectangles are object of class \code{dist}, the constraints are given as matrix
}
\details{Note that the labels are given in German. 
}

\references{
Borg, I., & Leutner, D. (1983). Dimensional models for the perception of rectangles. Perception and Psychophysics, 34, 257-269.

Borg, I., Groenen, P. J. F., & Mair, P. (2010). Multidimensionale Skalierung. Muenchen: Hampp Verlag.

}
\examples{
data(rectangles)
data(rect_constr)
}
\keyword{datasets}
