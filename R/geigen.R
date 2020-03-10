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


### I have made changes, and these should really have formal tests


geigen <- function(DAT, W, k = 0, tol= sqrt(.Machine$double.eps), symmetric){

  # preliminaries
  DAT.dims <- dim(DAT)
  if(length(DAT.dims)!=2){
    stop("gsvd: DAT must have dim length of 2 (i.e., rows and columns)")
  }
  ## should check square-ness here.
  if(DAT.dims[1] != DAT.dims[2]){
    stop("gsvd: DAT must be square (i.e., have the same number of rows and columns)")
  }

  DAT <- as.matrix(DAT)

  W.is.vector <- is.vector(W)
  W.is.missing <- missing(W)


  ### These are here out of convenience for the tests below. They started to get too long.
  if( !W.is.missing ){
    if(is.empty.matrix(W)){
      stop("geigen: W is empty (i.e., all 0s")
    }
  }


  if(!W.is.vector){

    if( nrow(W) != ncol(W) | nrow(W) != DAT.dims[2] ){
      stop("gsvd:nrow(W) does not equal ncol(W) or ncol(DAT)")
    }

    if( is.identity.matrix(W) ){

      W.is.missing <- T

    }else if( is.diagonal.matrix(W) ){

      W <- diag(W)
      W.is.vector <- T  #now it's a vector

      # if( length(W) != DAT.dims[2] ){ ## I should be able to get rid of this because of above
      #   stop("gsvd:length(W) does not equal ncol(DAT)")
      # }else{
      #   W.is.vector <- T  #now it's a vector
      # }
    }
  }

  if(!W.is.missing){
    if( W.is.vector ){
      # DAT <- sweep(sweep(DAT, 2, sqrt(W), "*"), 1, sqrt(W), "*")

      sqrt_W <- sqrt(W)
      ## this assumes square & symmetric; but I need to account for when it isn't by transposing back or moving the division
      DAT <- t(DAT / sqrt_W) / sqrt_W
    }else{
      sqrt_W <- W %^% (1/2)
      DAT <- sqrt_W %*% DAT %*% sqrt_W
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
    res$values <- NULL #this is actually OK
  res$d.orig <- sqrt(res$l.orig)
  res$tau <- (res$l.orig/sum(res$l.orig)) * 100

  components.to.return <- min(length(res$d.orig),k) #a safety check

  res$d <- res$d.orig[1:components.to.return]
  res$l <- res$l.orig[1:components.to.return]
  res$v <- as.matrix(res$vectors[,1:components.to.return])
    res$vectors <- NULL

  if(!W.is.missing){
    if(W.is.vector){
      res$q <- sweep(res$v,1,1/sqrt(W),"*") ## can replace this and also use sqrt_W
      res$fj <- sweep(sweep(res$q,1,W,"*"),2,res$d,"*") ## can replace this
    }else{
      res$q <- (W %^% (-1/2)) %*% res$v
      res$fj <- sweep((W %*% res$q),2,res$d,"*") ## can replace part of this
    }
  }else{
    res$q <- res$v
    res$fj <- sweep(res$q,2,res$d,"*") ## can replace part of this
  }

  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(DAT)

  class(res) <- c("list", "GSVD", "geigen")
  return(res)

}
