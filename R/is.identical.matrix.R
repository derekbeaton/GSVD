#' @export
#'
#' @title \code{is.identical.matrix}: test if a matrix contains all identical values.
#'
#' @description \code{is.identical.matrix} takes a matrix and tests if it contains all identical values.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#' @param round.digits number of decimal places to round to (via \code{\link{round}} function).
#'
#' @return A boolean. TRUE if the matrix contains all identical values, FALSE if the matrix does not.

is.identical.matrix <- function(x,tol=.Machine$double.eps, round.digits = 12){

  x <- as.matrix(x)
  x[abs(x) < tol] <- 0
  x <- round(x,digits=round.digits)

  if(length(unique(c(x)))==1){
    return(TRUE)
  }else{
    return(FALSE)
  }

}
