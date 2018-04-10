#' Generates Unvariate and Multivariate Generalized Poisson Variables
#'
#' @description
#'  This package is about generating unvariate and multivariate data that follow the generalized
#'  Poisson distributions.There are seven functions in the package: \code{\link{GenUniGpois}} generates univariate the generalized Poisson variables;
#' \code{\link{GenMVGpois}} generates the multivariate generalized Poisson variables;
#' \code{\link{ValidCorrGpois}} checks the validity of the values of pairwise correlations;
#' \code{\link{ComputeCorrGpois}} computes the lower and upper correlation bounds of a pairwise correlation between a pair of generalized Poisson variables;
#' \code{\link{CorrNNGpois}} adjusts the target correlation for a pair of generalized Poisson variables;
#' \code{\link{QuantileGpois}} computes the quantile of a given generalized Poisson distribution;
#' \code{\link{CmatStarGpois}} computes an intermediate correlation matrix. To learn more about this package please refer to both the reference manual and the vignette file.
#'
#' @docType package
#'
#' @name RNGforGPD-package
#'
#' @details
#' \tabular{ll}{Package: \tab RNGforGPD\cr
#'              Type: \tab Package\cr
#'              Version: \tab 1.0\cr
#'              Date: \tab 2018-04-09\cr
#'              License: \tab GPL-2 | GPL-3}
#' @author
#'  Hesen Li, Ruizhe Chen, Hai Nguyen, Yu-Che Chung, Ran Gao, Hakan Demirtas
#'
#'  Maintainer: Ruizhe Chen <rchen18@uic.edu>
#'
#' @references
#'  Demirtas, H. (2017). On accurate and precise generation of generalized
#'  Poisson variates. \emph{Communications in Statistics - Simulation and Computation},
#'  \bold{46(1)}, 489-499.
#'
#'  Yahav, I. and Shmueli, G. (2012). On generating multivariate Poisson data in management science applications.
#'  \emph{Applied Stochastic Models in Business and Industry}, \bold{28(1)}, 91-102.
#'
#'  Amatya, A. and Demirtas, H. (2015). Simultaneous generation of multivariate mixed data with Poisson
#'  and normal marginals. \emph{Journal of Statistical Computation and Simulation}, \bold{85(15)}, 3129-3139.
#'
#'  Amatya, A. and Demirtas, H. (2017). PoisNor: An R package for generation of multivariate data with
#'  Poisson and normal marginals. \emph{Communications in Statistics - Simulation and Computation},
#'  \bold{46(3)}, 2241-2253.
#'
#'  Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper correlation bounds.
#'  \emph{The American Statistician}, \bold{65(2)}, 104-109.
#'
#'
NULL



