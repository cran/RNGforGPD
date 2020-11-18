#' @include ComputeCorrGpois.R
#' @importFrom corpcor is.positive.definite
NULL


#' Validates Pairwise Correlations
#'
#' \code{ValidCorrGpois} checks the validity of the values of pairwise correlations including
#' positive definiteness, symmetry, and correctness of the dimensions.
#'
#' @param corMat a positive definite target correlation matrix whose entries are within the valid correlation bounds.
#' @param theta.vec rate parameters in the generalized Poisson distribution. It is assumed that the
#'  length of the vector is at least two, and each value has to be a positive number.
#' @param lambda.vec dispersion parameters in the generalized Poisson distribution. It is assumed that the length
#'  of the vector is at least two. All lambda values have to be less than 1. For lambda < 0, lambda must be greater than or equal to -theta/4.
#' @param verbose logical variable that determines whether to display the traces. Default is set to TRUE.
#' @return TRUE or FALSE.
#' @examples
#' \donttest{
#'  ValidCorrGpois(matrix(c(1, 0.9, 0.9, 1), byrow = TRUE, nrow = 2), 
#'                 c(0.5, 0.5), c(0.1, 0.105), verbose = TRUE)
#'  ValidCorrGpois(matrix(c(1, 0.9, 0.9, 1), byrow = TRUE, nrow = 2), 
#'                 c(3, 2), c(-0.3, -0.2), verbose = TRUE)}
#' @references
#'  Amatya, A. and Demirtas, H. (2017). PoisNor: An R package for generation of multivariate data with
#'  Poisson and normal marginals. \emph{Communications in Statistics - Simulation and Computation},
#'  \bold{46(3)}, 2241-2253.
#'  
#'  Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds.
#'  \emph{The American Statistician}, \bold{65(2)}, 104-109. 
#' @export
ValidCorrGpois = function (corMat, theta.vec, lambda.vec, verbose = TRUE) {
  no.gpois = length(theta.vec)
  errorCount = 0
  if (ncol(corMat) != (no.gpois)) {
    stop("Dimension of correlation matrix does not match the number of variables!\n")
  }
  if (is.positive.definite(corMat) == FALSE) {
    stop("Specified correlation matrix is not positive definite! \n")
  }
  if (isSymmetric(corMat) == FALSE) {
    stop("Specified correlation matrix is not symmetric! \n")
  }
  if (sum(corMat > 1) > 0) {
    stop("Correlation values cannot be greater than 1! \n")
  }
  if (sum(corMat < (-1)) > 0) {
    stop("Correlation values cannot be less than -1! \n")
  }
  if (sum(diag(corMat) != 1) > 0) {
    stop("All diagonal elements of the correlation matrix must be 1! \n")
  }
  maxcor = ComputeCorrGpois(theta.vec, lambda.vec, verbose)$max
  mincor = ComputeCorrGpois(theta.vec, lambda.vec, verbose)$min
  for (i in 1:nrow(corMat)) {
    for (j in 1:i) {
      if (errorCount == 0) {
        if (verbose == TRUE) {
        cat(".")
        }
      }
      if (i != j) {
        if (corMat[i, j] > maxcor[i, j] | corMat[i, j] <  mincor[i, j]) {
          cat("\n corMat[", i, ",", j, "] must be between ", round(mincor[i,j], 3), " and ", round(maxcor[i,j], 3), "\n")
          errorCount = errorCount + 1
        }
      }
    }
  }
  if (verbose == TRUE) {
  cat("\n")
  }
  if (errorCount > 0) {
    stop("Range violation occurred in the target correlation matrix!\n")
  }
  return(TRUE)
}


