#' Computes the Lower and Upper Correlation Bounds
#'
#' \code{ComputeCorrGpois} computes the lower and upper correlation bounds of pairwise
#' correlations between any pair of generalized Poisson variables using the Generate, Sort,
#' and Correlate (GSC) algorithm described in Demirtas and Hedeker (2011).
#'
#' @param theta.vec rate parameters in the generalized Poisson distribution. It is assumed that the
#'  length of the vector is at least two, and each value has to be a positive number.
#' @param lambda.vec dispersion parameters in the generalized Poisson distribution. It is assumed that the length
#'  of the vector is at least two. All lambda values have to be less than 1. 
#'  For lambda < 0, lambda must be greater than or equal to -theta/4.
#' @param verbose logical variable that determines whether to display the traces. Default is set to TRUE.
#' @return Lower and upper correlation bounds.
#' @examples
#'  \donttest{
#'  ComputeCorrGpois(c(3, 2, 5, 4), c(0.3, 0.2, 0.5, 0.6), verbose = TRUE)
#'  ComputeCorrGpois(c(4, 5), c(-0.45, -0.11), verbose = TRUE)}
#' @references
#'  Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds.
#'  \emph{The American Statistician}, \bold{65(2)}, 104-109.
#' @export
ComputeCorrGpois = function(theta.vec, lambda.vec, verbose = TRUE) {
  no.gpois = length(theta.vec)
  samples = 1e+05
  u = matrix(NA, nrow = no.gpois, ncol = samples)
  for (i in 1:no.gpois){
    u[i,] = GenUniGpois(theta.vec[i], lambda.vec[i], samples, method = "Inversion", details = FALSE)$data
  }
  maxmat = minmat = diag(NA, no.gpois)
  errorCount = 0
  for (i in 1:no.gpois){
    for (j in 1:no.gpois){
      if (i != j) {
        maxcor = cor(sort(u[i,]), sort(u[j,]))
        mincor = cor(sort(u[i,]), sort(u[j,], decreasing = TRUE))
        minmat[i, j] = mincor
        maxmat[i, j] = maxcor
        if (verbose == TRUE) {
          cat(".")
        }
      }
    }
  }
  if (verbose == TRUE) {
    cat("\n")
  }
  return(list(min = minmat, max = maxmat))
}

