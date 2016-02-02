\name{Procrustes}
\alias{Procrustes}
\alias{plot.procr}

\title{Procrustean Similarity Transformations}
\description{Solves the Procrustean problem of fitting one (MDS) configuration (testee) to another (target) MDS configuration.
}

\usage{
Procrustes(X, Y)

\method{plot}{procr}(x, plot.type = "jointplot", plot.dim = c(1,2), main, xlab, ylab, 
                     xlim, ylim, asp = 1, pch = 20, col.X = "cadetblue", col.Y = "gray", 
                     col.Yhat = "coral1", 
                     label.conf = list(label = TRUE, pos = 3, cex = 0.8), 
                     arrows = TRUE, length = 0.10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Target configuration}
  \item{Y}{Testee configuration}
  \item{x}{Object of class \code{procr}}
  \item{plot.type}{Either \code{"jointplot"} or \code{"transplot"}}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{pch}{Plot symbol.}
  \item{asp}{Aspect ratio.}
  \item{col.X}{Color target configuration.}
  \item{col.Y}{Color testess configuration.}
  \item{col.Yhat}{Color transformed configuration.}
  \item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color).}
  \item{length}{length of the edges of the arrow head (in inches).}
  \item{arrows}{For \code{"transplot"} only, whether arrows should be plotted or not}
  \item{...}{Additional plot arguments.}
}

\details{\code{Y} is going to be modified by finding an optimal dilation factor, an optimal translation and rotation for \code{Y}
such that it is a similar as possible to \code{X}. \code{X} remains untouched.  
}

\value{
Returns an object of class \code{procr} with:

  \item{X}{Input target configuration}
  \item{Y}{Input testee configuration}
  \item{Yhat}{Procrustes transformed (fitted) configuration}
  \item{translation}{Translation vector}
  \item{dilation}{Dilation factor}
  \item{rotation}{Rotation-reflection matrix}
  \item{confdistX}{Configuration distances X}
  \item{confdistY}{Configuration distances Y}
  \item{confdistYhat}{Configuration distances of fitted configuration}
  \item{congcoef}{Congruence coefficient}
}
\references{Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.
}

\examples{
## artificial example:
X <- matrix(c(1, -1, -1, 1, 2, 2, -2, -2), ncol = 2)
Y <- matrix(c(0.07, 0.93, 1.93, 1.07, 2.62, 3.12, 1.38, 0.88), ncol = 2)
op <- par(mfrow = c(1,2))
plot(X[,1], X[,2], xlim = c(-3, 3), ylim = c(-2, 3.5), asp = 1, xlab = "", ylab = "")
rect(-1, -2, 1, 2)
points(Y[,1], Y[,2], xlim = c(-3, 3), col = "gray")
polygon(Y[,1], Y[,2], border = "gray")

fitp <- Procrustes(X, Y)
plot(fitp$Yhat[,1], fitp$Yhat[,2], col = "red", xlim = c(-3, 3), ylim = c(-2, 3.5), 
     asp = 1, xlab = "", ylab = "")
polygon(fitp$Yhat[,1], fitp$Yhat[,2], border = "red")
par(op)

## MDS example:
eastD <- sim2diss(EW_eng$east)
attr(eastD, "Labels") <- abbreviate(attr(eastD, "Labels"))
fit.east <- mds(eastD, type = "ordinal")
westD <- sim2diss(EW_eng$west)
attr(westD, "Labels") <- abbreviate(attr(westD, "Labels"))
fit.west <- mds(westD, type = "ordinal", init = torgerson(eastD))

fit.proc <- Procrustes(fit.east$conf, fit.west$conf)
fit.proc

## Configuration plots; Procrustes plots.
plot(fit.east, main = "MDS East Germany")   ## MDS plot East Germany
plot(fit.west, main = "MDS West Germany")   ## MDS plot West Germany

## Procrustes configurations (X and Yhat)
plot(fit.proc, ylim = c(-1, 1))     
legend("topleft", legend = c("East Germany", "West Germany"), 
       col = c("cadetblue", "coral1"), pch = 19)

## Procrustes transformations (Y and Yhat)
plot(fit.proc, plot.type = "transplot", length = 0.05, ylim = c(-1,1))  
legend("topleft", legend = c("West Germany (untransformed)", "West Germany (transformed)"), 
       col = c("gray", "coral1"), pch = 19)

}

\keyword{models}
 