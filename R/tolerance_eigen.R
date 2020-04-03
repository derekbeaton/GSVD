#' @export
#'
#' @title \code{tolerance_eigen}: An eigenvalue decomposition to truncate potentially spurious (near machine precision) components.
#'
#' @description \code{tolerance_eigen} eliminates likely spurious components: any eigenvalue (squared singular value) below a tolerance level is elminated.
#'    The (likely) spurious eigen values and vectors are then eliminated from \code{$vectors} and \code{$values}.
#'    The use of a real positive value for \code{tol} will eliminate any small valued components.
#'    With \code{tol}, \code{tolerance_eigen} will stop if any singular values are complex or negative.
#'
#' @param x A data matrix of size for input to the eigen value decomposition (\code{\link{eigen}})
#' @param tol Default is \code{sqrt(.Machine$double.eps)}. A tolerance level for eliminating near machine precision components.
#' Use of this parameter causes \code{tolerance_eigen} to stop if negative or complex eigen values are detected.
#' The use of \code{tol < 0}, \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf}, or \code{NULL} passes through to \code{\link{eigen}}.
#' @param ... Further arguments to \code{\link{eigen}}. See \code{\link{eigen}}.
#'
#' @return A list with two elements (like \code{eigen}):
#'  \item{values}{ A vector containing the eigen values of x > \code{tol}.}
#'  \item{vectors}{ A matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension \code{min(c(ncol(x), nv, length(d))}.}
#'
#' @seealso \code{\link{eigen}}
#'
#' @examples
#'  data(wine)
#'  cor_X <- cor(as.matrix(wine$objective))
#'  s_asis <- tolerance_eigen(cor_X)
#'  s_.Machine <- tolerance_eigen(cor_X, tol=.Machine$double.eps)
#'  s_000001 <- tolerance_eigen(cor_X, tol=.000001)
#'
#' @author Derek Beaton
#' @keywords multivariate

tolerance_eigen <- function(x, tol = sqrt(.Machine$double.eps), ...) {

  ### this probably needs a try-catch
  eigen_res <- eigen(x, ...)

  # if tolerance is any of these values, just do nothing; send back the EVD results as is.
  if( (is.null(tol) | is.infinite(tol) | is.na(tol) | is.nan(tol) | tol < 0) ){

    return(eigen_res)

  }
  ## once you go past this point you *want* the tolerance features.

  if(any(unlist(lapply(eigen_res$values,is.complex)))){
    stop("tolerance_eigen: eigen values ($values) are complex.")
  }
  # if( (any(abs(eigen_res$values) > tol) ) & (any(sign(eigen_res$values) != 1)) ){

  # if( (any(abs(eigen_res$values) < tol) ) ){
  if( any( (abs(eigen_res$values) > tol) & (sign(eigen_res$values)==-1) ) ){
    stop("tolerance_eigen: eigen values ($values) are negative with a magnitude above 'tol'.")
  }

  evs.to.keep <- which(!(eigen_res$values < tol))
  if(length(evs.to.keep)==0){
    stop("tolerance_eigen: All eigen values were below 'tol'")
  }

  eigen_res$values <- eigen_res$values[evs.to.keep]

    ## this would happen if only.values=TRUE
  if(!is.null(eigen_res$vectors)){
    eigen_res$vectors <- eigen_res$vectors[,evs.to.keep]
    rownames(eigen_res$vectors) <- colnames(x)

    ## new way inspired by FactoMineR but with some changes
    vector_signs <- ifelse(colSums(eigen_res$vectors) < 0, -1, 1)
    eigen_res$vectors <- t(t(eigen_res$vectors) * vector_signs)

  }

  # class(eigen_res) <- c("eigen", "GSVD", "list")
  # return(eigen_res)
  return(eigen_res)
}
