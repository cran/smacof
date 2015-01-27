indscal <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal"), 
                   weightmat = NULL, init = NULL, ties = "primary", 
                   verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6) {
  
  smacofIndDiff(delta = delta, ndim = ndim, type = type, constraint = "indscal",
                weightmat = weightmat, init = init, ties = ties, verbose = verbose, modulus = modulus, 
                itmax = itmax, eps = eps)
}