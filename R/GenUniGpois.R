#' Generates Univariate Generalized Poisson Variates
#'
#' \code{GenUniGpois} generates univariate random variables from the generalized Poisson
#' distribution using one of the five methods including Inversion, Branching, Normal-Approximation, Build-Up, and Chop-Down.
#'
#' @param theta the rate parameter in the generalized Poisson distribution. It has to be a positive number.
#' @param lambda the dispersion parameter in the generalized Poisson distribution.
#'  It has to be less than 1. For lambda < 0, lambda must be greater than or equal to -theta/4.
#' @param n number of data points that is to be generated.
#' @param details index to indicate whether to print out the estimates of parameters. Default is set to TRUE.
#' @param method index to specify one of the five methods for generating univariate GPD variable: 
#' "Inversion", "Branching", "Normal-Approximation", "Build-Up" or "Chop-Down".
#' @details
#'  All five methods come from Demirtas (2017).
#'  When lambda equals to 0, it is the ordinary Poisson distribution, so there is no need to specify the method.
#'  "Branching" only works when lambda is positive.
#'  When theta is less than 10, the "Normal-Approximation" may not be reliable.
#' @return A list that includes generated data, specified and empirical values of theta and lambda, and the specified method.
#' @examples
#' GenUniGpois(5, -0.4, 100, method = "Inversion")
#' GenUniGpois(2, 0.9, 100, method = "Branching")
#' GenUniGpois(12, 0.5, 100, method = "Normal-Approximation")
#' data <- GenUniGpois(3, 0.9, 10000, method = "Build-Up", details = FALSE)
#' data <- GenUniGpois(10, 0.4, 10, method = "Chop-Down", details = FALSE)
#' @references
#'  Demirtas, H. (2017). On accurate and precise generation of generalized
#'  Poisson variates. \emph{Communications in Statistics - Simulation and Computation},
#'   \bold{46(1)}, 489-499.
#' @export
GenUniGpois <- function(theta, lambda, n, details = TRUE, method) {
  # Check if parameters are valid
  tol <- .Machine$double.eps^0.5
  if (n <= 0 | abs(n - round(n)) > tol) {
    stop("n must be a positive integer.")
  }
  if(n == 1) stop("Parameter estimates are exactly the same as input of parameters in the case when n = 1!")
  if(theta <= 0) stop("theta has to be greater than 0!")
  if(lambda >= 1) stop("lambda has to be less than 1!")
  if (lambda < 0 & lambda < (-theta) / 4) stop(paste("For lambda < 0, lambda must be greater than or equal to -theta/4 which is ", (-theta)/4, "!", sep = ""))
  myset = numeric(n)
  if (lambda == 0) { # 1. rpois function
    myset = rpois(n, theta)
    if (n != 1) {
      
        params = MOM.genpois(myset)
        emp.theta = round(params$emp.theta, 6)
        emp.lambda = round(params$emp.lambda, 6)
      
    }
  } else if (method == "Inversion"){ # 2. Inversion Method (when lambda is negative)
    w = exp(-lambda)
    for (i in 1:n){
      mys = exp(-theta)
      myp = mys
      x = 0
      u = runif(1)
      while (u > mys){
        x = x + 1
        myc = theta - lambda + lambda * x
        myp = w * myc * (1 + lambda/myc)^(x-1) * myp * x^(-1)
        mys = mys + myp
      }
      myset[i] = x
    }
    if (n != 1) {
      
        params = MOM.genpois(myset)
        emp.theta = round(params$emp.theta, 6)
        emp.lambda = round(params$emp.lambda, 6)

    }
  } else if (method == "Branching") { # 3. Branching Method (when lambda is positive)
    if (lambda < 0) stop("Lambda should be greater than 0!")
    for (i in 1:n){
      index = 0
      y = rpois(1,theta)
      x = y
      if (y <= 0) myset[i] = x
      else while(index == 0){
        z = rpois(1, lambda * y)
        x = x + z
        y = z
        index = 1 * (y <= 0)
      }
      myset[i] = x
    }
    if (n != 1) {

        params = MOM.genpois(myset)
        emp.theta = round(params$emp.theta, 6)
        emp.lambda = round(params$emp.lambda, 6)
        
    }
  } else if (method == "Normal-Approximation") { # 4. Normal-Approximation Method (when theta is large enough which is bigger than 10)
    if (theta < 10) warning("Normal approximation may not be reliable for theta less than 10!")
    mym = theta / (1 - lambda)
    myv = sqrt(theta * (1 - lambda)^(-3))
    y = rnorm(n)
    x = floor(mym + myv * y + 0.5)
    x[x < 0] = 0
    myset = x
    if (n != 1) {
      
        params = MOM.genpois(myset)
        emp.theta = round(params$emp.theta, 6)
        emp.lambda = round(params$emp.lambda, 6)

    }
  } else if (method == "Build-Up") { # 5. Build-Up Method
    mynumx = numeric(n)
    t = exp(-theta)
    for (i in 1:n) {
      u = runif(1)
      x = 0
      px = t
      s = px
      while (u > s) {
        x = x + 1
        px = theta * (theta + lambda * x)^(x - 1) * exp(-theta - lambda * x) / factorial(x)
        s = s + px
      }
      myset[i] = x
    }
    if (n != 1) {

        params = MOM.genpois(myset)
        emp.theta = round(params$emp.theta, 6)
        emp.lambda = round(params$emp.lambda, 6)

    }
  } else if (method == "Chop-Down") { # 6. Chop-Down Method
    mynump = numeric(n)
    t = exp(-theta)
    for (i in 1:n) {
      u = runif(1)
      x = 0
      px = t
      while (u > px) {
        u = u - px
        x = x + 1
        px = theta * (theta + lambda * x)^(x - 1) * exp(-theta - lambda * x) / factorial(x)
      }
      myset[i] = x
    }
    if (n != 1) {

        params = MOM.genpois(myset)
        emp.theta = round(params$emp.theta, 6)
        emp.lambda = round(params$emp.lambda, 6)

    }
  }
  if (details == TRUE) {
    print(paste("Specified theta is ", theta, ", empirical theta is ",
                emp.theta, ", specified lambda is ", lambda, ", empirical lambda is ", emp.lambda, ".", sep = ""))
  }
  return(invisible(list(data = myset,
                        specified.theta = theta,
                        empirical.theta = as.numeric(emp.theta),
                        specified.lambda = lambda,
                        empirical.lambda = as.numeric(emp.lambda),
                        method = method)))
}


