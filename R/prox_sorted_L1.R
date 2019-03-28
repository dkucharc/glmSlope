#' Compute the prox for the sorted L1 norm
#' @param  x an input vector
#' @param lambda vectors of \eqn{\lambda}'s
#' 
#' @return The prox estimations for the sorted L1 norm
#' @references M. Bogdan et al. (2015) \emph{SLOPE--Adaptive variable selection via convex optimization}, \url{http://dx.doi.org/10.1214/15-AOAS842} 
#' @export
prox_sorted_L1 <- function(x, lambda) {
  if (is.complex(x))
  {
    sign = Arg(x)
    x = Mod(x)
  }
  else
  {
    sign = sign(x)
    x = abs(x)
  }
  
  s = sort(x, decreasing = TRUE, index.return=TRUE)
  s$ix <- s$ix - 1
  
  result <- prox_sorted_L1_C(s$x, lambda, s$ix)
  
  result <- result * sign
  
  return(result)
}