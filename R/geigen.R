#' @export
#'
#' @title Generalized Eigenvalue decomposition
#'
#' @description
#' \code{geigen} takes in right (\code{RW}) constraints (usually diagonal matrices, but any positive semi-definite matrix) that are applied to the data (\code{DAT})
#'   Right constraints are used for the orthogonality conditions.
#'
#' @param DATA a data matrix to decompose
#' @param RW \bold{R}ight \bold{W}eights -- the constraints applied to the right side of the matrix.
#' @param k total number of components to return though the full variance will still be returned (see \code{d.orig}). If 0, the full set of components are returned.
#' @param tol default is .Machine$double.eps. A parameter with two roles: A tolerance level for (1) eliminating (tiny variance or negative or imaginary) components and (2) converting all values < tol to 0 in \code{u} and \code{v}.
#'
#' @return A list with nine elements:
#' \item{values_orig}{A vector containing the eigen values > \code{tol}.}
#' \item{tau}{A vector that contains the (original) explained variance per component (eigenvalues derived from \code{$d.orig}.}
#' \item{values}{A vector containing the eigen values of x > \code{tol}. Length is \code{min(length(d.orig), k)}}
#' \item{vectors}{Right (columns) singular vectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{generalized_vectors}{Right (columns) generalized singular vectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{scores}{Right (columns) component scores. Dimensions are \code{ncol(DAT)} by k.}
#'
#' @seealso \code{\link{tolerance.svd}}, \code{\link{gsvd}} and \code{\link{svd}}
#'
#' @examples
#'
#'
#' @author Derek Beaton
#' @keywords multivariate, diagonalization, eigen


## I need to figure out the names... I should homogenize the names across gsvd & geigen

geigen <- function(DATA, RW, k = 0, tol=.Machine$double.eps){

  DAT <- as.matrix(DAT)
  # DAT[abs(DAT) < tol] <- 0
  RW.is.vector <- RW.is.missing <- F ##asuming everything is a matrix.

  ### These are here out of convenience for the tests below. They started to get too long.
  if( !missing(RW) ){
    if(is.empty.matrix(RW)){
      stop("geigen: RW is empty (i.e., all 0s")
    }
  }

  if( missing(RW) ){
    RW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(RW) ) {
      RW.is.vector <- T
    }else if(!RW.is.vector){

      if( is.identity.matrix(RW) ){
        RW.is.missing <- T
        warning("gsvd: RW was an identity matrix. RW will not be used in the GSVD.")
      }else if( is.diagonal.matrix(RW) ){

        RW <- diag(RW)

        if( length(RW) != DAT.dims[2] ){
          stop("gsvd:length(RW) does not equal ncol(DAT)")
        }else{
          RW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(RW) != ncol(RW) | nrow(RW) != DAT.dims[2] ){
        stop("gsvd:nrow(RW) does not equal ncol(RW) or ncol(DAT)")
      }
    }
  }

    ### the sqrt here may not be correct...
  if(!RW.is.missing){
    if( RW.is.vector ){  ## replace with sweep
      DAT <- sweep(DAT,2,sqrt(RW),"*")
    }else{
      DAT <- DAT %*% (RW %^% (1/2))
    }
  }


}
