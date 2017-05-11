\name{biplotmds}
\alias{biplotmds}
\alias{biplotmds.smacof}
\alias{biplotmds.smacofID}
\alias{plot.mdsbi}


\title{MDS Biplots}
\description{Regresses external variables on a MDS configuration which results in a MDS biplot.
}
\usage{
\method{biplotmds}{smacof}(object, extvar, scale = TRUE)
\method{plot}{mdsbi}(x, vecscale = NULL, plot.dim = c(1,2), sphere = TRUE, col = 1, 
                     label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), 
                     vec.conf = list(col = 1, cex = 0.8, length = 0.1), 
                     identify = FALSE, type = "p", pch = 20, 
                     asp = 1, main, xlab, ylab, xlim, ylim, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{"smacof"} or \code{"smacofID"}.}
  \item{extvar}{Data frame with external variables.}
  \item{scale}{If \code{TRUE}, external variables are standardized internally.}
  
  \item{x}{Object of class \code{"mdsbi"}.}
  \item{vecscale}{Scaling factor for regression coefficients, either a single number or \code{NULL} (see details).}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{pch}{Plot symbol.}
  \item{asp}{Aspect ratio.}
  \item{col}{Point color.}
  \item{type}{What type of plot should be drawn.}
  \item{sphere}{In case of spherical smacof, whether sphere should be plotted or not.}
  \item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color).}
  \item{vec.conf}{List with arguments for arrows and arrow labels of the external variables.}
  \item{identify}{If \code{TRUE}, the \code{identify()} function is called internally that allows to add configuration labels by mouse click.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{plot}} for detailed information.}
}

\details{If a model for individual differences is provided, the external variables are regressed on the group stimulus space configurations. In the biplot only the relative length of the vectors and their direction matters. Using the scale argument the user can control for the relative length of the vectors. If \code{vecscale = NULL}, the \code{vecscale()} function from the \pkg{candisc} package is used which tries to automatically calculate the scale factor so that the vectors approximately fill the same space as the configuration. 
}

\value{   
Returns an object belonging to classes \code{"mlm"} and \code{"mdsbi"}. See \code{\link{lm}} for details.

  \item{R2vec}{Vector containing the R2 values.}
}

\references{
Greenacre, M. (2010). Biplots in Practice. Fundacion BBVA, Bilbao, Spain
}

\seealso{\code{\link{plot.smacof}}}

\examples{
## morse code data with external scales
res <- mds(morse)
fitbi <- biplotmds(res, morsescales[,2:3])
plot(fitbi, main = "MDS Biplot", vecscale = 0.5)

## wish data with external economic development factor
diss <- sim2diss(wish, method = 7)
res <- mds(diss, type = "ordinal")
ecdev <- data.frame(ecdev = c(3,1,3,3,8,3,7,9,4,7,10,6))
fitbi <- biplotmds(res, ecdev)
plot(fitbi, main = "MDS Biplot", vecscale = 1)
plot(fitbi, main = "MDS Biplot", vecscale = 0.5, xlim = c(-1, 1), 
vec.conf = list(col = "red", length = 0.05))

## Ekman's color data (by Michael Friendly)
require(colorspace)
wavelengths <- attr(ekman, "Labels")
colors <- c("#2600F0", "#0028FF", "#0092FF", "#00B2FF", "#00FFFF", "#00FF61", "#77FF00", "#B3FF00",
            "#FFF200", "#FFBE00", "#FF9B00", "#FF5700", "#F60000", "#D60000")
ekmanD <- sim2diss(ekman)
res <- mds(ekmanD, type = "ordinal")
RGB <- t(col2rgb(colors)) / 255
HCL <- as(hex2RGB(colors), "polarLUV")
HCL <- slot(HCL, "coords")
fit <- biplotmds(res, cbind(RGB, HCL))
fit
plot(fit, vecscale = 0.5, cex = 6, col = colors,
	   label.conf=list(cex = 1, pos = ifelse(wavelengths < 490, 4, 3)),
	   vec.conf = list(cex = 1.2), main = "Ekman configuration and color properties" )
}

\keyword{hplot}
