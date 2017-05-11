\name{PVQ40}
\alias{PVQ40}
\alias{PVQ40agg}
\docType{data}
\title{Portrait Value Questionnaire}
\description{The PVQ40 (Schwartz et al., 1999) consists of 40 items, each a short portrait of one person. For example, to
measure power, the PVQ includes two portraits (male/female versions): It is important to him to be rich. He wants to have a lot of money and expensive things. It is important to him to get respect from others. He wants people to do
what he says. Respondent indicates on 6-point bipolar rating scale (1 ... not at all like me, 6 ... very much like me) the degree to which the description also fits him-/herself. Gender and age of the participants are added as attributes.
}
\usage{data(PVQ40)}
\format{PVQ40 data of 151 adults from various states in the USA:

sd1-sd4: self-direction

po1-po3: power

un1-un6: universalism

ac1-ac4: achievement

se1-se5: security

st1-st3: stimulation

co1-co4: conformity

tr1-tr4: tradition

he1-he3: hedonism

be1-be4: benevolence

Age and Gender are added as attributes.
}

\details{
\code{PVQ40agg} is an aggregated version of \code{PVQ40} where the item scores belonging to the same value are averaged. 
}

\references{Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}
\examples{
str(PVQ40)
head(PVQ40)
attr(PVQ40, "Gender")
attr(PVQ40, "Age")
str(PVQ40agg)
}
\keyword{datasets}
