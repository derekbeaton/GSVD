#' @export
#'
#' @title \code{tolerance.svd}: An SVD to truncate potentially spurious (near machine precision) components.
#'
#' @description \code{tolerance.svd} eliminates likely spurious components: any eigenvalue (squared singular value) below a tolerance level is elminated.
#'    The (likely) spurious singular values and vectors are then eliminated from \code{$u}, \code{$d}, and \code{$v}.
#'    Additionally, all values in \code{abs($u)} or \code{abs($v)} that fall below the \code{tol} are set to 0.
#'    The use of a real positive value for \code{tol} will eliminate any small valued components.
#'    With \code{tol}, \code{tolerance.svd} will stop if any singular values are complex or negative.
#'
#' @param x A data matrix of size for input to the singular value decomposition (\code{\link{svd}})
#' @param nu The number of left singular vectors to be computed. Default is \code{min(dim(x))}
#' @param nv The number of right singular vectors to be computed. Default is \code{min(dim(x))}
#' @param tol Default is \code{.Machine$double.eps}. A tolerance level for eliminating near machine precision components.
#' Use of this parameter causes \code{tolerance.svd} to stop if negative or complex singular values are detected.
#' The use of \code{tol < 0}, \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf}, or \code{NULL} passes through to \code{\link{svd}}.
#'
#' @return A list with three elements (like \code{svd}):
#'  \item{d}{ A vector containing the singular values of x > \code{tol}.}
#'  \item{u}{ A matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension \code{min(c(nrow(x), nu, length(d))}.}
#'  \item{v}{ A matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension \code{min(c(ncol(x), nv, length(d))}.}
#'
#' @seealso \code{\link{svd}}
#'
#' @examples
#'  data(wine)
#'  X <- scale(as.matrix(wine$objective))
#'  s_asis <- tolerance.svd(X)
#'  s_sqrt.Machine <- tolerance.svd(X,tol=sqrt(.Machine$double.eps))
#'  s_000001 <- tolerance.svd(X,tol=.000001)
#'
#' @author Derek Beaton
#' @keywords multivariate, diagonalization, eigen

tolerance.svd <- function(x, nu=min(dim(x)), nv=min(dim(x)), tol = .Machine$double.eps) {

  ## the R SVD is much faster/happier when there are more rows than columns in a matrix
    ## however, even though a transpose can speed up the SVD, there is a slow down to then set the U and V back to where it was
    ## so I will remove this for now. I just need to keep it in mind.

  # x.dims <- dim(x)
  # x.is.transposed <- F
  # if( (x.dims[1]*10) < x.dims[2]){ # * 10 to make it worth the transpose.
  #   x.is.transposed <- T
  #   x <- t(x)
  # }

  ## nu and nv are pass through values.
  svd.res <- svd(x, nu = nu, nv = nv)

  # if tolerance is any of these values, just do nothing; send back the SVD results as is.
  if( (is.null(tol) | is.infinite(tol) | is.na(tol) | is.nan(tol) | tol < 0) ){

    return(svd.res)

  }
  ## once you go past this point you *want* the tolerance features.


  if(any(unlist(lapply(svd.res$d,is.complex)))){
    stop("tolerance.svd: Singular values ($d) are complex.")
  }
  # if( (any(abs(svd.res$d) > tol) ) & (any(sign(svd.res$d) != 1)) ){
  if( any( (svd.res$d^2 < tol) & (sign(svd.res$d)==-1) ) ){
    stop("tolerance.svd: Singular values ($d) are negative with a magnitude above 'tol'.")
  }

  svs.to.keep <- which(!(svd.res$d^2 < tol))
  if(length(svs.to.keep)==0){
    stop("tolerance.svd: All (squared) singular values were below 'tol'")
  }

  svd.res$d <- svd.res$d[svs.to.keep]

  ## are these checks necessary? problably...
  if(nu >= length(svs.to.keep)){
    svd.res$u <- as.matrix(svd.res$u[,svs.to.keep])
  }else{
    svd.res$u <- as.matrix(svd.res$u[,1:nu])
  }

  if(nv >= length(svs.to.keep)){
    svd.res$v <- as.matrix(svd.res$v[,svs.to.keep])
  }else{
    svd.res$v <- as.matrix(svd.res$v[,1:nv])
  }

  rownames(svd.res$u) <- rownames(x)
  rownames(svd.res$v) <- colnames(x)

  ## force consistent directions as best as possible:
  if( sign(svd.res$u[1]) == -1 ){
    svd.res$u <- svd.res$u * -1
    svd.res$v <- svd.res$v * -1
  }

  return(svd.res)
}
