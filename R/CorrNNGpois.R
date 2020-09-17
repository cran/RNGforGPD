#' @include GenUniGpois.R
#' @importFrom VGAM vglm

NULL


#' Adjusts the Target Correlation
#'
#' \code{CorrNNGpois} adjusts the actual/realized correlation to the target correlation bounds for
#' a pair of generalized Poisson variables.
#'
#' @param theta.vec rate parameters in the generalized Poisson distribution. It is assumed that the
#'  length of the vector is at least two, and each value has to be a positive number.
#' @param lambda.vec dispersion parameters in the generalized Poisson distribution. It is assumed that the length
#'  of the vector is at least two. All lambda values have to be less than 1. 
#'  For lambda < 0, lambda must be greater than or equal to -theta/4.
#' @param r desired target correlation.
#' @return the adjusted target correlation.
#' @examples
#' \donttest{
#'  CorrNNGpois(c(0.1, 10), c(0.1, 0.2), 0.5)
#'  CorrNNGpois(c(0.1, 10), c(-0.01, -0.02), 0.5)
#'  CorrNNGpois(c(4, 2.3), c(-0.32,-0.3), 0.7)
#'  CorrNNGpois(c(14, 10), c(-0.8, -0.3), 0.99)}
#' @references
#'  Yahav, I. and Shmueli, G.(2012), On generating multivariate Poisson data in management science applications.
#'  \emph{Applied Stochastic Models in Business and Industry}, \bold{28(1)}, 91-102.
#' @export
CorrNNGpois = function(theta.vec, lambda.vec, r) {
  if(abs(r) > 1) stop("The desired correlation r has to be within (-1,1)")
  samples = 1e+05
  u = GenUniGpois(theta.vec[1], lambda.vec[1], samples, method = "Inversion", details = FALSE)$data
  v = GenUniGpois(theta.vec[2], lambda.vec[2], samples, method = "Inversion", details = FALSE)$data
  maxcor = cor(sort(u), sort(v))
  mincor = cor(sort(u), sort(v, decreasing = TRUE))
  a = -maxcor * mincor/(maxcor + mincor)
  b = log((maxcor + a)/a)
  c = -a
  corrected = log((r + a)/a)/b
  if( abs(corrected) > 1 ) cat("The actual correlation, ",corrected,", is not feasible!", sep = "")
  else{ return(corrected) }
}


