\name{svm_mdsplot}
\alias{svm_mdsplot}

\title{Support Vector Machine MDS}
\description{Plots 2D MDS configuration including facets as determined by an SVM. 
}
\usage{
svm_mdsplot(mds_object, svm_object, class, legend1 = TRUE, legend2 = TRUE, 
            inset = c(-0.2, 0.5), plot.dim = c(1,2), by = 0.01, 
            main, xlab, ylab, xlim, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mds_object}{Object of class \code{"smacofB"}, i.e., an MDS solution from \code{smacofSym()} or \code{smacofConstraint}.}
  \item{svm_object}{Object of class \code{"svm"}, i.e., an SVM solution from \code{svm} or \code{tune.svm}}
  \item{class}{Vector of class assignments (facets) for each object.}
  \item{legend1}{If \code{TRUE}, facet legend is added.}
  \item{legend2}{If \code{TRUE}, class legend is added.}
  \item{inset}{Inset distance from the margins for both legends as a fraction of the plot region when legend is placed by keyword.}
  \item{plot.dim}{Vector with dimensions to be plotted.}
  \item{by}{Scaling factor for resolution (the smaller, the higher the resolution).}
  \item{main}{Plot title.}
  \item{xlab}{Label of x-axis.}
  \item{ylab}{Label of y-axis.}
  \item{xlim}{Scale x-axis.}
  \item{ylim}{Scale y-axis.}
  \item{\dots}{Further plot arguments passed: see \code{\link[graphics]{image}} for detailed information.}
}

\details{Using the SVM implementation of \code{e1071} one can determine facets in an MDS configuration based on an SVM fit. This function plots the resulting facets on top of the 2D MDS configuration. Note that this function is work in progress. 
}

\seealso{\code{\link[e1071]{svm}}, \code{\link[e1071:tune.wrapper]{tune.svm}}}

\examples{
## Guttman intelligence data
Delta <- sim2diss(Guttman1965[[1]])
class <- Guttman1965[[2]]

## ordinal MDS fit
mds_gut <- mds(Delta, ndim = 2, type = "ordinal")  
mds_gut
cols <- rainbow_hcl(4)[as.numeric(class)]
plot(mds_gut, col = cols, label.conf = list(col = cols))
legend("bottomright", legend = levels(class), cex = 0.7, col = rainbow_hcl(4), pch = 19)

## radial SVM fit
X <- mds_gut$conf                    ## extract configuration
dat <- data.frame(class = class, X)  ## merge with class vector
costvec <- 2^seq(-4, 4)              ## tuning parameter grid
gamma <- seq(0.01, 0.5, 10)

set.seed(111)
svm_gut <- tune.svm(class ~ D1 + D2, data = dat, kernel = "radial", 
                    cross = 10, cost = costvec)$best.model
svm_gut
preds <- predict(svm_gut, data = dat)  ## predicted classes
table(obs = class, pred = preds)       ## confusion matrix

svm_mdsplot(mds_gut, svm_gut, dat$class, inset = c(-0.3, 0.5))
}

\keyword{ models }
\keyword{ hplot }
