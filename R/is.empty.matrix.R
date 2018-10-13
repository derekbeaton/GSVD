#' @export
#'
#' @title \code{is.empty.matrix}: test if a matrix contains all 0s.
#'
#' @description \code{is.empty.matrix} takes a matrix and tests if it contains all 0s.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix contains all 0s, FALSE if the matrix does not.


is.empty.matrix <- function(x,tol=.Machine$double.eps){

  x <- as.matrix(x)
  x[abs(x) < tol] <- 0

  if(sum(abs(x))==0){
    return(TRUE)
  }else{
    return(FALSE)
  }

}
