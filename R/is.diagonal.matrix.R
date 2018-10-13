#' @export
#'
#' @title \code{is.diagonal.matrix}: test if a matrix is a diagonal matrix.
#'
#' @description \code{is.diagonal.matrix} takes a matrix and tests if it is a diagonal matrix.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix is a diagonal matrix, FALSE if the matrix is not.


  ## I stole this from somewhere... but I don't remember where
    ## I should find where I stole this from.
    ## I believe this: https://stackoverflow.com/questions/11057639/identifying-if-only-the-diagonal-is-non-zero
is.diagonal.matrix <- function(x,tol=.Machine$double.eps){
  if(is.null(dim(x))){
    stop("is.diagonal.matrix: X is not a matrix.")
  }
  x[ x^2 < tol ] <- 0
  return(all(x[lower.tri(x)] == 0, x[upper.tri(x)] == 0))
}
