#' @export
#'
#' @title \code{matrix.generalized.inverse}: psuedo inverse and rebuild of a matrix
#'
#' @description \code{matrix.generalized.inverse} takes in a matrix and will compute the psuedo-inverse via the singular value decomposition.
#'  Additionally, the psuedo-inverse can be computed for a lower rank estimate of the matrix.
#'
#' @param x data matrix to compute the pseudo-inverse of
#' @param k the number of components to retain in order to build a lower rank estimate of \code{x}
#' @param ... parameters to pass through to \code{\link{tolerance.svd}}
#'
#' @return The (possibly lower rank) generalized inverse of \code{x}
#'
#' @seealso \code{\link{tolerance.svd}} and \code{\link{matrix.exponent}}
#'
#' @examples
#'  data(wine)
#'  X <- as.matrix(wine$objective)
#'  X.inv <- matrix.generalized.inverse(X)
#'  X.inv %*% X ## is approximately an identity.
#'
#' @author Derek Beaton
#'
#' @keywords multivariate, diagonalization, eigen, pseudo-inverse, Moore-Penrose, generalized inverse

#matrix.generalized.inverse <- mgi <- m.g.i <- function(x, k=0, ...){
matrix.generalized.inverse <- function(x, k=0, ...){

  ##stolen from MASS::ginv()
  if (length(dim(x)) > 2L || !(is.numeric(x) || is.complex(x)))
    stop("matrix.generalized.inverse: 'x' must be a numeric or complex matrix")
  if (!is.matrix(x))
    x <- as.matrix(x)

  k <- round(k)
  if(k <= 0){
    k <- min(nrow(x),ncol(x))
  }

  ## I need to test which is faster...

  ## the special cases:

  ## is diagonal
  if(is.diagonal.matrix(x)){
    return( diag( diag(x)^(-1) ) )

  }
  ## is vector
  if( any(dim(x)==1) ){
    return( x^(-1) )
  }

  #res <- tolerance.svd(x,...)  ## just go with the defaults of this or allow pass through?
  #comp.ret <- 1:min(length(res$d),k)
  #return( sweep(res$v[,comp.ret],2,res$d[comp.ret],"/") %*% t(res$u[,comp.ret]) )

  res <- tolerance.svd(x, nu = k, nv = k,...)
  if(k > length(res$d)){
    k <- length(res$d)
  }
  return( sweep(res$v,2,res$d[1:k],"/") %*% t(res$u) )
}
