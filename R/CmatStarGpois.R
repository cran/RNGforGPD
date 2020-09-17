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
#' @param theta.vec rate parameters in the generalized Poisson distribution. It is assumed that the
#'  length of the vector is at least two, and each value has to be a positive number.
#' @param lambda.vec dispersion parameters in the generalized Poisson distribution. It is assumed that the length
#'  of the vector is at least two. All lambda values have to be less than 1.
#'  For lambda < 0, lambda must be greater than or equal to -theta/4.
#' @param verbose logical variable that determines whether to display the traces. Default is set to TRUE.
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
#'  cstar = CmatStarGpois(TV, theta.vec, lambda.vec, verbose = TRUE)
#'  cstar}
#' @references
#'  Yahav, I. and Shmueli, G. (2012). On generating multivariate Poisson data in management science applications.
#'  \emph{Applied Stochastic Models in Business and Industry}, \bold{28(1)}, 91-102.
#' @export
CmatStarGpois = function(corMat, theta.vec, lambda.vec, verbose = TRUE) {
  no.gpois = length(theta.vec)
  if (ValidCorrGpois(corMat, theta.vec, lambda.vec, verbose)) {
    corMat.star = diag(nrow(corMat))
    # lower matrix index
    g <- expand.grid(row = 1:nrow(corMat), col = 1:ncol(corMat))
    g.lower.tri <- g[lower.tri(corMat, diag = TRUE),]
    corMat.lower.index <- g.lower.tri[-which(g.lower.tri[,1]==g.lower.tri[,2]),]
    for (i in 1:nrow(corMat.lower.index)) {
      i.temp <- corMat.lower.index[i,1]
      j.temp <- corMat.lower.index[i,2]
      corMat.star[i.temp,j.temp] <- CorrNNGpois(c(theta.vec[i.temp], theta.vec[j.temp]),
                                                c(lambda.vec[i.temp], lambda.vec[j.temp]),
                                                corMat[i.temp,j.temp])
      if (verbose == TRUE) {
        cat(".")
      }
    }
    # upper matrix index
    g.upper.tri <- g[upper.tri(corMat, diag = TRUE),]
    corMat.upper.index <- g.upper.tri[-which(g.upper.tri[,1]==g.upper.tri[,2]),]
    for (i in 1:nrow(corMat.upper.index)) {
      i.temp <- corMat.upper.index[i,1]
      j.temp <- corMat.upper.index[i,2]
      sym.index <- intersect(which(corMat.lower.index[,2]==i.temp), which(corMat.lower.index[,1]==j.temp))
      corMat.star[i.temp, j.temp] <- corMat.star[corMat.lower.index[sym.index,1], corMat.lower.index[sym.index,2]]
    }
  }
  if (verbose == TRUE) {
    cat("\n")
  }
  if (!is.positive.definite(corMat.star)) {
    warning("Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
    corMat.star = as.matrix(nearPD(corMat.star, corr = TRUE, keepDiag = TRUE)$mat)
  }
  corMat.star = (corMat.star + t(corMat.star))/2
  return(corMat.star)
}

