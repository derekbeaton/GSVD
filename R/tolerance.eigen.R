#' @export
#'
#' @title \code{tolerance.eigen}: An eigenvalue decomposition to truncate potentially spurious (near machine precision) components.
#'
#' @description \code{tolerance.eigen} eliminates likely spurious components: any eigenvalue (squared singular value) below a tolerance level is elminated.
#'    The (likely) spurious eigen values and vectors are then eliminated from \code{$vectors} and \code{$values}.
#'    The use of a real positive value for \code{tol} will eliminate any small valued components.
#'    With \code{tol}, \code{tolerance.eigen} will stop if any singular values are complex or negative.
#'
#' @param x A data matrix of size for input to the eigen value decomposition (\code{\link{eigen}})
#' @param tol Default is \code{sqrt(.Machine$double.eps)}. A tolerance level for eliminating near machine precision components.
#' Use of this parameter causes \code{tolerance.eigen} to stop if negative or complex eigen values are detected.
#' The use of \code{tol < 0}, \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf}, or \code{NULL} passes through to \code{\link{eigen}}.
#' @param ... Further arguments to \code{\link{eigen}}. See \code{\link{eigen}}.
#'
#' @return A list with two elements (like \code{eigen}):
#'  \item{values}{ A vector containing the eigen values of x > \code{tol}.}
#'  \item{vectors}{ A matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension \code{min(c(ncol(x), nv, length(d))}.}
#'
#' @seealso \code{\link{eigen}}
#'
#' @author Derek Beaton
#' @keywords multivariate, diagonalization, eigen

tolerance.eigen <- function(x, tol = sqrt(.Machine$double.eps), ...) {

  eigen_res <- eigen(x, ...)

  # if tolerance is any of these values, just do nothing; send back the EVD results as is.
  if( (is.null(tol) | is.infinite(tol) | is.na(tol) | is.nan(tol) | tol < 0) ){

    return(eigen_res)

  }
  ## once you go past this point you *want* the tolerance features.

  if(any(unlist(lapply(eigen_res$values,is.complex)))){
    stop("tolerance.eigen: eigen values ($values) are complex.")
  }
  # if( (any(abs(eigen_res$values) > tol) ) & (any(sign(eigen_res$values) != 1)) ){

  # if( (any(abs(eigen_res$values) < tol) ) ){
  if( any( (abs(eigen_res$values) > tol) & (sign(eigen_res$values)==-1) ) ){
    stop("tolerance.eigen: eigen values ($values) are negative with a magnitude above 'tol'.")
  }

  evs.to.keep <- which(!(eigen_res$values < tol))
  if(length(evs.to.keep)==0){
    stop("tolerance.eigen: All eigen values were below 'tol'")
  }

  eigen_res$values <- eigen_res$values[evs.to.keep]
    ## this would happen if only.values=TRUE
  if(!is.null(eigen_res$vectors)){
    eigen_res$vectors <- eigen_res$vectors[,evs.to.keep]
    rownames(eigen_res$vectors) <- colnames(x)

    ## force consistent directions as best as possible:
    if( sign(eigen_res$vectors[1]) == -1 ){
      eigen_res$vectors <- eigen_res$vectors * -1
    }
  }

  class(eigen_res) <- c("list", "GSVD", "eigen")
  return(eigen_res)
}
