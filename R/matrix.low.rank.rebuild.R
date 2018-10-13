#' @export
#'
#' @title \code{matrix.low.rank.rebuild}: raise matrix to a power and rebuild lower rank version
#'
#' @description \code{matrix.low.rank.rebuild} takes in a matrix and will rebuild a lower rank estimate.
#'
#' @param x data matrix
#' @param k a flexible parameter in order to build a lower rank estimate of \code{x}. \\ If the value is greater than 1, k is the number of consecutive components (from first) to retain. \\ If 0 < k < 1, this function rebuilds a matrix with the components that explain k*100 percent of the variance. \\ k can also be a vector of numbers that correspond to arbitrary components (e.g., k=c(1,3,5) rebuilds \code{x} with the first, third, and fifth components only).
#' @param ... parameters to pass through to \code{\link{tolerance.svd}}.
#'
#' @return A low rank version of \code{x}.
#'
#' @seealso \code{\link{tolerance.svd}}, \code{\link{matrix.exponent}}, and \code{\link{matrix.generalized.inverse}}
#'
#' @author Derek Beaton
#'
#' @keywords multivariate, diagonalization, eigen, low rank, rank


## TODO: rebuild with gsvd() not tolerance.svd() [though technically it's just a pass through at that point...]
#### the use of gsvd() will require putting the weights and whatnot back in.

matrix.low.rank.rebuild <- function(x, k = 0, ...){

  ## quick tests for escape
  if( !is.numeric(k) ){
    stop("k is not numeric")
  }else{
    if( is.infinite(k) | is.nan(k) ){
      stop("k is Inf or NaN")
    }
  }

  if( length(k)==0 ){
    stop("k is of length 0")
  }

  ##stolen from MASS::ginv()
  if (length(dim(x)) > 2L || !(is.numeric(x) || is.complex(x)))
    stop("matrix.low.rank.rebuild: 'x' must be a numeric or complex matrix")
  if (!is.matrix(x))
    x <- as.matrix(x)


  if( length(k) > 1 ){

    ## this allows for a rebuild of arbitrary components...
    comps.to.rebuild <- unique( round(ifelse(k < 1, 1, k)) )
    res <- tolerance.svd(x, nu = max(comps.to.rebuild), nv = max(comps.to.rebuild), ...)
    return( sweep(res$u[,comps.to.rebuild],2,res$d[comps.to.rebuild],"*") %*% t(res$v[,comps.to.rebuild]) )


  }else { #if(length(k) == 1){
    if( k <= 0 ){
      stop("k is a negative or 0 value")
    }

    ## this allows for a rebuild of some cumulative percentage
    if(k > 0 & k < 1 ){
      res <- tolerance.svd(x)
      comps.to.rebuild <- 1:(max(which( cumsum(res$tau) <= k))+1)
      return( sweep(res$u[,comps.to.rebuild],2,res$d[comps.to.rebuild],"*") %*% t(res$v[,comps.to.rebuild]) )

    ## this rebuilds based on first K components where K = round(k)
    }else {#if( (k > 0) ){
      res <- tolerance.svd(x, nu = round(k), nv = round(k) )
      return( sweep(res$u,2,res$d[1:round(k)],"*") %*% t(res$v) )

    }
  }


}
