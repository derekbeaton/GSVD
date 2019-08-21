#' @export
#'
#' @title Generalized eigenvalue decomposition
#'
#' @description
#' \code{geigen} takes in constraints (\code{W}), (usually diagonal matrices, but any positive semi-definite matrix) that are applied to the data (\code{DAT}).
#'   Constraints are used for the orthogonality conditions.
#'
#' @param DAT a square, symmetric data matrix to decompose
#' @param W \bold{W}eights -- the constraints applied to the matrix and thus the eigen vectors.
#' @param k total number of components to return though the full variance will still be returned (see \code{d.orig}). If 0, the full set of components are returned.
#' @param tol default is \code{sqrt(.Machine$double.eps)}. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components.
#' @param symmetric if \code{DAT} is symmetric, set as TRUE. See \code{\link{eigen}}.
#'
#' @return A list with eight elements:
#' \item{d.orig}{A vector containing the singular values of DAT above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of DAT above the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}.}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values of DAT}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values of DAT}
#' \item{v}{Eigenvectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{q}{Generalized eigenvectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{fj}{Component scores. Dimensions are \code{ncol(DAT)} by k.}
#'
#' @seealso \code{\link{tolerance.eigen}}, \code{\link{tolerance.svd}}, \code{\link{gsvd}} and \code{\link{svd}}
#'
#' @examples
#' Observed <- authors/sum(authors)
#' row.w <- rowSums(Observed)
#' row.W <- diag(1/row.w)
#' col.w <- colSums(Observed)
#' col.W <- diag(1/col.w)
#' Expected <- row.w %o% col.w
#' Deviations <- Observed - Expected
#' ca.res_gsvd <- gsvd(Deviations,row.W,col.W)
#' ca.res_geigen <- geigen(t(Deviations) %*% diag(1/row.w) %*% Deviations, col.W)
#'
#' @author Derek Beaton
#' @keywords multivariate, diagonalization, eigen


### Should I do some PSD checks?


geigen <- function(DAT, W, k = 0, tol= sqrt(.Machine$double.eps), symmetric){

  # preliminaries
  DAT.dims <- dim(DAT)
  if(length(DAT.dims)!=2){
    stop("gsvd: DAT must have dim length of 2 (i.e., rows and columns)")
  }
  DAT <- as.matrix(DAT)

  W.is.vector <- W.is.missing <- F ##asuming everything is a matrix.

  ### These are here out of convenience for the tests below. They started to get too long.
  if( !missing(W) ){
    if(is.empty.matrix(W)){
      stop("geigen: W is empty (i.e., all 0s")
    }
  }

  if( missing(W) ){
    W.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(W) ) {
      W.is.vector <- T
    }else if(!W.is.vector){

      if( is.identity.matrix(W) ){
        W.is.missing <- T
        warning("gsvd: W was an identity matrix. W will not be used in the GSVD.")
      }else if( is.diagonal.matrix(W) ){

        W <- diag(W)

        if( length(W) != DAT.dims[2] ){
          stop("gsvd:length(W) does not equal ncol(DAT)")
        }else{
          W.is.vector <- T  #now it's a vector
        }

      }else if( nrow(W) != ncol(W) | nrow(W) != DAT.dims[2] ){
        stop("gsvd:nrow(W) does not equal ncol(W) or ncol(DAT)")
      }
    }
  }

  if(!W.is.missing){
    if( W.is.vector ){  ## replace with sweep
      DAT <- sweep(sweep(DAT, 2, sqrt(W), "*"), 1, sqrt(W), "*")
    }else{
      DAT <- (W %^% (1/2)) %*% DAT %*% (W %^% (1/2))
    }
  }

  if(k<=0){
    k <- min(nrow(DAT),ncol(DAT))
  }

  if(missing(symmetric)){
    symmetric <- isSymmetric(DAT)
  }

  res <- tolerance.eigen(DAT, tol=tol, symmetric=symmetric)


  res$l.orig <- res$values
    res$values <- NULL #ew
  res$d.orig <- sqrt(res$l.orig)
  res$tau <- (res$l.orig/sum(res$l.orig)) * 100

  components.to.return <- min(length(res$d.orig),k) #a safety check

  res$d <- res$d.orig[1:components.to.return]
  res$l <- res$l.orig[1:components.to.return]
  res$v <- as.matrix(res$vectors[,1:components.to.return])
    res$vectors <- NULL #ew

  if(!W.is.missing){
    if(W.is.vector){
      res$q <- sweep(res$v,1,1/sqrt(W),"*")
      res$fj <- sweep(sweep(res$q,1,W,"*"),2,res$d,"*")
    }else{
      res$q <- (W %^% (-1/2)) %*% res$v
      res$fj <- sweep((W %*% res$q),2,res$d,"*")
    }
  }else{
    res$q <- res$v
    res$fj <- sweep(res$q,2,res$d,"*")
  }

  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(DAT)

  class(res) <- c("list", "GSVD", "geigen")
  return(res)

}
