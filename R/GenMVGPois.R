#' @include GenUniGpois.R
#' @include ValidCorrGpois.R
#' @include ComputeCorrGpois.R
#' @include CorrNNGpois.R
#' @include CmatStarGpois.R
#' @include QuantileGpois.R
#' @importFrom mvtnorm rmvnorm
#' @importFrom corpcor is.positive.definite
#' @importFrom stats pnorm
#' @import VGAM
NULL


#' Generates Data from Multivariate Generalized Poisson Distribution
#'
#' \code{GenMVGpois} simulates a sample of size \emph{sample.size} from a set of multivariate generalized
#' Poisson variables with correlation matrix \emph{cmat.star} and pre-specified marginals.
#'
#' @param sample.size desired sample size (number of rows) for the multivariate generalized Poisson data
#' @param no.gpois dimension of the multivariate generalized Poisson distribution.
#' @param cmat.star intermediate correlation matrix.
#' @param theta.vec rate parameters in the generalized Poisson distribution. It is assumed that the
#'  length of the vector is at least two, and each value has to be a positive number.
#' @param lambda.vec dispersion parameters in the generalized Poisson distribution. It is assumed that the length
#'  of the vector is at least two. All lambda values have to be < 1. For lambda < 0, lambda must be >= -theta/4.
#' @param details index of whether to display the specified and empirical values of parameters. Default is set to be TRUE.
#' @return data that follow multivariate generalized Poisson distribution.
#' @examples
#' \donttest{
#'  sample.size = 10000; no.gpois = 3
#'  lambda.vec = c(0.2, 0.2, 0.3); theta.vec = c(1, 3, 4)
#'  M = c(0.352, 0.265, 0.342); N = diag(3); N[lower.tri(N)] = M
#'  TV = N + t(N); diag(TV) = 1
#'  cstar = CmatStarGpois(TV, theta.vec, lambda.vec)
#'  data = GenMVGpois(sample.size, no.gpois, cstar, theta.vec, lambda.vec, details = FALSE)
#'  apply(data, 2, mean) # empirical means
#'  theta.vec / (1 - lambda.vec) # theoretical means
#'  apply(data, 2, var) # empirical variances
#'  theta.vec / (1 - lambda.vec)^3 # theoretical variances
#'  cor(data) # empirical correlation matrix
#'  TV # specified correlation matrix}
#' @references
#'  Amatya, A. and Demirtas, H. (2015). Simultaneous generation of multivariate mixed data with Poisson
#'  and normal marginals. \emph{Journal of Statistical Computation and Simulation}, \bold{85(15)}, 3129-3139.
#'
#'  Amatya, A. and Demirtas, H. (2017). PoisNor: An R package for generation of multivariate data with
#'  Poisson and normal marginals. \emph{Communications in Statistics - Simulation and Computation},
#'  \bold{46(3)}, 2241-2253.
#'  
#'  Demirtas, H. (2017). On accurate and precise generation of generalized
#'  Poisson variates. \emph{Communications in Statistics - Simulation and Computation},
#'  \bold{46(1)}, 489-499.
#'
#'  Yahav, I. and Shmueli, G.(2012). On generating multivariate Poisson data in management science applications.
#'  \emph{Applied Stochastic Models in Business and Industry}, \bold{28(1)}, 91-102.
#' @export
GenMVGpois = function(sample.size, no.gpois, cmat.star, theta.vec, lambda.vec, details = TRUE) {
  if (sample.size == 1 & details == TRUE) {
    stop("Parameter estimates are exactly the same as input of parameters in the case when n = 1!")
  }
  if (length(theta.vec) != no.gpois) {
    stop("Dimension of the theta vector does not match the number of generalized Poisson variables!\n")
  }
  if (length(lambda.vec) != no.gpois) {
    stop("Dimension of the lambda vector does not match the number of generalized Poisson variables!\n")
  }
  if (nrow(cmat.star) != no.gpois) {
    stop("Dimension of cmat.star and number of variables do not match!\n")
  }
  XX = rmvnorm(sample.size, rep(0, no.gpois), cmat.star)
  YY = NULL
  for (i in 1:no.gpois) {
    UU = pnorm(XX[, i])
    XXgpois = QuantileGpois(UU, theta.vec[i],lambda.vec[i], FALSE)
    YY = cbind(YY, XXgpois)
  }
  colnames(YY) = NULL
  if (sample.size != 1) {
    model = vglm(YY ~ 1, genpoisson(zero = 1))
    coef = matrix(Coef(model), 2, no.gpois)
    emp.theta = round(coef[2, ], 6)
    emp.lambda = round(coef[1, ], 6)
    emp.corr = cor(YY)
    if (details == TRUE) {
      my.res = list(data = YY,
                  specified.theta = theta.vec,
                  empirical.theta = emp.theta,
                  specified.lambda = lambda.vec,
                  empirical.lambda = emp.lambda,
                  specified.corr = cmat.star,
                  empirical.corr = emp.corr)
      print(my.res[2:6])
      return(my.res)
    }
  }
  if (details == FALSE) {
    return(YY)
  }
}

