# converts similarity matrix into dissimilarities

sim2diss <- function(s, method = "corr", to.dist = FALSE)
{
  # s... similarity matrix (not necessarily symmetric, nor quadratic)
  # method ... various methods provided
  # to.dist ... if TRUE, it creates an object of class "dist", if FALSE a matrix.
  s <- as.matrix(s)
  EPS <- .Machine$double.eps
  
  if (!is.numeric(method)) method <- match.arg(method, c("corr", "reverse", "reciprocal", "ranks", 
                                                          "exp", "Gaussian", "cooccurrence", "gravity", 
                                                          "confusion", "transition", "membership", "probability"), 
                                                several.ok = FALSE)
  
  if (method == "corr") {
    if (any(s < -1) || any(s > 1)) stop( "Correlations expected for correlation transform." )
    dissmat <- sqrt(1-s)
  } 
  if (method == "reverse") dissmat <- max(s, na.rm = TRUE) + min(s, na.rm = TRUE) - s
  if (method == "reciprocal") {
    s[s == 0] <- NA
    dissmat <- 1/s
  }
  if (method == "ranks") {
    dissmat <- matrix(rank(-s), dim(s))
    colnames(dissmat) <- colnames(s)
  }  
  if (method == "exp") dissmat <- -log((EPS+s)/(EPS+max(s)))
  if (method == "Gaussian") dissmat <- sqrt(-log((EPS+s)/(EPS+max(s))))
  if (method == "cooccurrence") { 
    rsum <- rowSums(s, na.rm = TRUE)
    csum <- colSums(s, na.rm = TRUE)
    tsum <- sum(s, na.rm = TRUE)
    s <- (tsum*s)/(rsum%*%t(csum))
    dissmat <- (1/(1+s))
  }  
  if (method == "gravity") {
    s[s == 0] <- NA
    rsum <- rowSums(s, na.rm = TRUE)
    csum <- colSums(s, na.rm = TRUE)
    tsum <- sum(s, na.rm = TRUE)
    s <- (rsum%*%t(csum))/(tsum*s)
    dissmat <- sqrt(s)
  }
  if (method == "confusion") {
    if (any(s < 0) || any(s > 1)) stop( "Proportions expected for confusion transform!" )
    dissmat <- 1-s
  }
  if (method == "transition") {
    if (any(s < 0)) stop( "Frequencies expected for transition transform." )
    s[s == 0] <- NA
    dissmat <- 1/sqrt(s)
  }
  if (method == "membership") dissmat <- 1-s
  if (method == "probability") {
    if (any(s < 0) || any(s > 1)) stop( "Probabilities expected for probability transform." )
    s[s == 0] <- NA
    dissmat <- 1/sqrt(asin(s))
  }
  
  if (is.numeric(method)) dissmat <- method - s

  if (to.dist) dissmat <- as.dist(dissmat)

  return(dissmat)
}
  
                       
                    
 
