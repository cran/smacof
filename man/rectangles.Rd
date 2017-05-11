\name{rectangles}
\alias{rectangles}
\alias{rect_constr}
\docType{data}
\title{Rectangles}
\description{These data are based on an experiment by Borg and Leutner (1983). They constructed rectangles on the basis of the grid design (see \code{rect_constr}). Each point in this grid defines a rectangle. Rectangle 16, for example, had a width of 4.25 cm and a height of
1.25 cm; rectangle 4 was 3.00 cm wide and 2.75 cm tall. A total of 21 persons rated (twice) the similarity of each pair of these 16 rectangles. 
on a 10-point scale ranging from from 0 = equal/identical to 9 = very different. The means of these ratings over persons and replications are given in \code{rectangles}. 
}

\usage{data(rectangles)
data(rect_constr)
}

\format{The rectangles are object of class \code{dist}, the constraints are given as matrix
}

\references{
Borg, I., & Leutner, D. (1983). Dimensional models for the perception of rectangles. Perception and Psychophysics, 34, 257-269.

Borg, I., Groenen, P. J. F., & Mair, P. (2017). Applied Multidimensional Scaling and Unfolding. New York, Springer.

}
\examples{
rectangles
rect_constr
}
\keyword{datasets}
