#' Compute the prox for the sorted L1 norm
#' @param  x an input vector
#' @param lambda vectors of $\lamda$'s
#' 
#' @return The prox estimations for the sorted L1 norm
#' @examples 
#' proxSortedL1(x, lambda)
proxSortedL1 = function (x, lambda) {
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
  
  s = sort(x, decreasing=TRUE, index.return=TRUE)
  s$ix <- s$ix - 1
  
  result <- prox_sorted_L1_C(s$x, lambda, s$ix)
  
  result <- result * sign
  
  return(result)
}