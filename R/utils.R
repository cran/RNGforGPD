MOM.genpois <- function(data){
  
  m <- mean(data); v <- var(data)
  
  emp.theta <- sqrt(m^3 / v)
  emp.lambda <- 1 - sqrt(m / v)
  
  return(list(
    emp.theta = emp.theta,
    emp.lambda = emp.lambda
  ))
  
}