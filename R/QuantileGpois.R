#' Computes Quantiles
#'
#' \code{QuantileGpois} computes the quantile for the generalized Poisson distribution
#'  for specified values of percentile, lambda and theta parameters.
#'
#' @param p percentile of the generalized Poisson distribution.
#' @param theta the rate parameter in the generalized Poisson distribution. It has to be a positive number.
#' @param lambda the dispersion parameter in the generalized Poisson distribution.
#'  It has to be < 1. For lambda < 0, lambda must be >= -theta/4.
#' @param details show the detailed information of probability and cumulative probability.
#'  Default is set as FALSE.
#' @return quantile of the specified distribution if the parameter \emph{details} is set as FALSE.
#'  quantile and the detailed information of probability and cumulative probability if
#'  the parameter \emph{details} is set as TRUE.
#' @examples
#'  QuantileGpois(0.98,1,-0.2,details = TRUE)
#'  QuantileGpois(0.80,2,0.025,details = FALSE)
#' @references
#'  Demirtas, H. (2017). On accurate and precise generation of generalized
#'  Poisson variates. \emph{Communications in Statistics - Simulation and Computation},
#'   \bold{46(1)}, 489-499.
#' @export
QuantileGpois = function(p, theta, lambda, details = FALSE){
  # Check if parameters are valid
  if(theta <= 0) stop("Theta has to be greater than 0!")
  if(lambda > 1) stop("Lambda has to be less than 1!")

  # since m >= 4, check lower bound for lambda when labda is negative
  if (lambda < 0 & lambda < (-theta)/4) stop (paste("For lambda<0, lambda must be >= -theta/4 = ", (-theta)/4, "!", sep=""))
  
  # Determine m
  m = numeric(1)
  if (lambda < 0) {
    mod = floor(-theta/lambda)
    if (-theta-mod * lambda == 0) m = mod
    else m = mod + 1
  }
  p.in = p
  upper = max(p.in)
  s = numeric(10000)
  q = numeric(length(p.in))
  w = exp(-lambda)
  
  #P(X = 0)
  p = exp(-theta)
  # s is for cumulative probability
  s[1] = p
  if (details) cat("x = 0, P(X = x) =",p, ",P(X <= x) =", s[1], "\n")
  i=1
  s[1] = as.numeric(p)
  if (details) cat("x = 0, P(X = x) =",p , ", P(X <= x) =", s[1], "\n")
  i = 1
  while (s[i] < upper) {
    #P(X = x)
    if (lambda < 0 & i > m) break
    else{
      p = theta * (theta + lambda * i)^(i - 1) * exp(-theta-lambda*i) / factorial(i)
      if(is.infinite(p)) {cat("error: p goes to infinite", "\n"); break  }
      if (i == 10000) {
        temp = numeric(10000)
        s = c(s, temp)
      }
        s[i + 1] = s[i] + p
        if (details) cat("x =", i, ", P(X = x) =", p, ", P(X <= x) =", s[i+1], "\n")
        i=i+1
     }
  }
  
  # For lambda < 0 to eliminate the truncation error ###
  if (lambda < 0) {
    Fm = s[i]
    s[1:i] = s[1:i] * Fm^(-1)
    if (details){
      cat("When lambda is negative, we need to account for truncation error\n")
      cat("The adjusted CDF are:", s[1:i])
    }
  }
  
  # quantile for x
  if(!is.infinite(p)) {for (j in 1:length(p.in)) {
    i = 1
    while (p.in[j] > s[i]) {
      i = i + 1
    }
    q[j] = i - 1
  }
    return(q)
  }
  else {
    for (j in 1:length(p.in)) {
      i = 1
      while (p.in[j] > s[i]) {
        i = i + 1
      }
      q[j] = i - 1
    }
    return(q)
  }
}

