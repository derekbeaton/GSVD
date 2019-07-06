
#' @export
#'
#' @title \code{matrix.exponent}: raise matrix to a power and rebuild lower rank version
#'
#' @description \code{matrix.exponent} takes in a matrix and will compute raise that matrix to some arbitrary power via the singular value decomposition.
#'  Additionally, the matrix can be computed for a lower rank estimate of the matrix.
#'
#' @param x data matrix
#' @param power the power to raise \code{x} by (e.g., 2 is squared)
#' @param k the number of components to retain in order to build a lower rank estimate of \code{x}
#' @param ... parameters to pass through to \code{\link{tolerance.svd}}
#'
#' @return The (possibly lower rank) raised to an arbitrary \code{power} version of \code{x}
#'
#' @seealso \code{\link{tolerance.svd}} and \code{\link{matrix.generalized.inverse}}
#'
#' @examples
#'  data(wine)
#'  X <- as.matrix(wine$objective)
#'  X.power_1 <- matrix.exponent(X)
#'  X / X.power_1
#'
#'  ## other examples.
#'  X.power_2 <- matrix.exponent(X,power=2)
#'  X.power_negative.1.div.2 <- matrix.exponent(X,power=-1/2)
#'
#'  X.power_negative.1 <- matrix.exponent(X,power=-1)
#'  X.power_negative.1 / (X %^% -1)
#'
#' @author Derek Beaton
#'
#' @keywords multivariate, diagonalization, eigen

#matrix.exponent <- me <- m.e <- function(x, power = 1, k = 0, ...){
matrix.exponent <- function(x, power = 1, k = 0, ...){

  ##stolen from MASS::ginv()
  if (length(dim(x)) > 2L || !(is.numeric(x) || is.complex(x)))
    stop("matrix.exponent: 'x' must be a numeric or complex matrix")
  if (!is.matrix(x))
    x <- as.matrix(x)

  k <- round(k)
  if(k<=0){
    k <- min(nrow(x),ncol(x))
  }

  ## should be tested for speed.

  #res <- tolerance.svd(x,...)
  #comp.ret <- 1:min(length(res$d),k)
  #return( (res$u[,comp.ret] * matrix(res$d[comp.ret]^power,nrow(res$u[,comp.ret]),ncol(res$u[,comp.ret]),byrow=T)) %*% t(res$v[,comp.ret]) )


  ## the special cases:
  ## power = 0
  if(power==0){
    x <- diag(1,nrow(x),ncol(x))
    attributes(x)$message.to.user = "https://www.youtube.com/watch?v=9w1y-kMPNcM"
    return( x )
  }
  ## is diagonal
  if(is.diagonal.matrix(x)){
    return( diag( diag(x)^power ) )

  }
  ## is vector
  if( any(dim(x)==1) ){
    return( x^power )
  }

  res <- tolerance.svd(x, nu = k, nv = k, ...)
  if(k > length(res$d)){
    k <- length(res$d)
  }
  return( sweep(res$u,2,res$d[1:k]^power,"*") %*% t(res$v) )

}



#' @export
#'
#' @title Matrix exponentiation
#'
#' @description takes in a matrix and will compute raise that matrix to some arbitrary power via the singular value decomposition.
#'  Additionally, the matrix can be computed for a lower rank estimate of the matrix.
#'
#' @param x data matrix
#' @param power the power to raise \code{x} by (e.g., 2 is squared)
#'
#' @return \code{x} raised to an arbitrary \code{power}
#'
#' @seealso \code{\link{matrix.exponent}}
#'
#' @examples
#'  data(wine)
#'  X <- as.matrix(wine$objective)
#'  X %^% 2 # power of 2
#'  X %^% -1 # (generalized) inverse
#'
#' @author Derek Beaton
#'
#' @keywords multivariate, diagonalization, eigen
#'
`%^%` <- function(x,power){
  matrix.exponent(x,power=power)
}


