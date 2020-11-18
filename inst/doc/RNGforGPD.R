## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
) 

## ---- echo = FALSE, include = FALSE-------------------------------------------
library(RNGforGPD)
library(corpcor)
library(mvtnorm)
library(Matrix)

## -----------------------------------------------------------------------------
GenUniGpois(2, 0.9, 5000, method = "Branching")
GenUniGpois(5, -0.4, 5000, method = "Inversion")
GenUniGpois(12, 0.5, 5000, method = "Normal-Approximation")
data = GenUniGpois(3, 0.9, 10, method = "Build-Up", details = FALSE)
data
data = GenUniGpois(10, 0.4, 10, method = "Chop-Down", details = FALSE)
data

## -----------------------------------------------------------------------------
set.seed(3406)
ComputeCorrGpois(c(3, 2, 5, 4), c(0.3, 0.2, 0.5, 0.6), verbose = FALSE)
ComputeCorrGpois(c(4, 5), c(-0.45, -0.11), verbose = FALSE)

## -----------------------------------------------------------------------------
ValidCorrGpois(matrix(c(1, 0.9, 0.9, 1), byrow = TRUE, nrow = 2), c(0.5, 0.5), c(0.1, 0.105), verbose = TRUE)
ValidCorrGpois(matrix(c(1, 0.9, 0.9, 1), byrow = TRUE, nrow = 2), c(3, 2), c(-0.3, -0.2), verbose = TRUE)

## -----------------------------------------------------------------------------
QuantileGpois(0.98, 1, -0.2, details = TRUE)
QuantileGpois(0.80, 2, 0.025, details = FALSE)

## -----------------------------------------------------------------------------
set.seed(3406)
CorrNNGpois(c(0.1,10), c(0.1, 0.2),0.5)
lambda.vec = c(-0.2, 0.2, -0.3)
theta.vec = c(1, 3, 4)
M = c(0.352, 0.265, 0.342)
N = diag(3)
N[lower.tri(N)] = M
TV = N + t(N)
diag(TV) = 1
cstar = CmatStarGpois(TV, theta.vec, lambda.vec, verbose = TRUE)
cstar

## -----------------------------------------------------------------------------
set.seed(3406)
lambda.vec = c(-0.2, 0.2, -0.3)
theta.vec = c(1, 3, 4)
M = c(0.352, 0.265, 0.342)
N = diag(3)
N[lower.tri(N)] = M
TV = N + t(N)
diag(TV) = 1
cstar = CmatStarGpois(TV, theta.vec, lambda.vec, verbose = TRUE)
sample.size = 10000; no.gpois = 3
data = GenMVGpois(sample.size, no.gpois, cstar, theta.vec, lambda.vec, details = FALSE)
apply(data, 2, mean) # empirical means
theta.vec / (1 - lambda.vec) # theoretical means
apply(data, 2, var) # empirical variances
theta.vec / (1 - lambda.vec)^3 # theoretical variances
cor(data) # empirical correlation matrix
TV # specified correlation matrix

