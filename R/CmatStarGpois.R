#' @include ValidCorrGpois.R
#' @include CorrNNGpois.R
#' @importFrom corpcor is.positive.definite
#' @importFrom Matrix nearPD
NULL


#' Computes Intermediate Correlation Matrix
#'
#' \code{CmatStarGpois} computes an intermediate correlation matrix that will be used to obtain
#' the target correlation matrix using the inverse CDF transformation method in \code{GenMVGpois}.
#' If the intermediate correlation martrix is not positive definite, the nearest positive definite
#' matrix is used.
#'
#' @param corMat target correlation matrix.
#' @param theta.vec theta.vec rate parameters in the generalized Poisson distribution. It is assumed that the
#'  length of the vector is at least two, and each value has to be a positive number.
#' @param lambda.vec dispersion parameters in the generalized Poisson distribution. It is assumed that the length
#'  of the vector is at least two. All lambda values have to be < 1. For lambda < 0, lambda must be >= -theta/4.
#' @return intermediate correlation matrix.
#' @examples
#' \donttest{
#'  lambda.vec = c(-0.2, 0.2, -0.3)
#'  theta.vec = c(1, 3, 4)
#'  M = c(0.352, 0.265, 0.342)
#'  N = diag(3)
#'  N[lower.tri(N)] = M
#'  TV = N + t(N)
#'  diag(TV) = 1
#'  cstar = CmatStarGpois(TV, theta.vec, lambda.vec)
#'  cstar}
#' @references
#'  Yahav, I. and Shmueli, G. (2012). On generating multivariate Poisson data in management science applications.
#'  \emph{Applied Stochastic Models in Business and Industry}, \bold{28(1)}, 91-102.
#' @export
CmatStarGpois = function(corMat, theta.vec, lambda.vec){
  no.gpois = length(theta.vec)
  if (ValidCorrGpois(corMat, theta.vec, lambda.vec)) {
    corMat.star = diag(nrow(corMat))
    for (i in 1:nrow(corMat)) {
      for (j in 1:no.gpois) {
        if (i != j) {
          corMat.star[i,j] = CorrNNGpois(c(theta.vec[i], theta.vec[j]),
                                          c(lambda.vec[i], lambda.vec[j]),
                                          corMat[i, j])
        }
        cat(".")
      }
    }
  }
  cat("\n")
  if (!is.positive.definite(corMat.star)) {
    warning("Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
    corMat.star = as.matrix(nearPD(corMat.star, corr = TRUE, keepDiag = TRUE)$mat)
  }
  corMat.star = (corMat.star + t(corMat.star))/2
  return(corMat.star)
}

