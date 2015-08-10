\name{PVQ40}
\alias{PVQ40}
\docType{data}
\title{Portrait Value Questionnaire}
\description{The PVQ40 (Schwartz et al., 1999) consists of 40 items, each a short portrait of one person. For example, to
measure power, the PVQ includes two portraits (male/female versions): It is important to him to be rich. He wants to have a lot of money and expensive things. It is important to him to get respect from others. He wants people to do
what he says.
Respondent indicates on 6-point bipolar rating scale (1 ... not at all like me, 6 ... very much like me) the degree to which the
description also fits him-/herself:
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
}

\references{Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}
\examples{
pvqdiss <- sim2diss(cor(PVQ40, use = "complete.obs"))
fit <- mds(pvqdiss, type = "ordinal")
plot(fit)
plot(fit, plot.type = "Shepard")

}
\keyword{datasets}
