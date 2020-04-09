#' @export
#'
#' @title Generalized eigenvalue decomposition
#'
#' @description
#' \code{geigen} takes in constraints (\code{W}), (usually diagonal matrices, but any positive semi-definite matrix) that are applied to the data (\code{X}).
#'   Constraints are used for the orthogonality conditions.
#'
#' @param X a square, symmetric data matrix to decompose
#' @param W \bold{W}eights -- the constraints applied to the matrix and thus the eigen vectors.
#' @param k total number of components to return though the full variance will still be returned (see \code{d.orig}). If 0, the full set of components are returned.
#' @param tol default is \code{sqrt(.Machine$double.eps)}. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components.
#' @param symmetric if \code{X} is symmetric, set as TRUE. See \code{\link{eigen}}.
#'
#' @return A list with eight elements:
#' \item{d.orig}{A vector containing the singular values of X above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of X above the tolerance threshold (\code{tol}).}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values of X}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values of X}
#' \item{v}{Eigenvectors. Dimensions are \code{ncol(X)} by k.}
#' \item{q}{Generalized eigenvectors. Dimensions are \code{ncol(X)} by k.}
#' \item{fj}{Component scores. Dimensions are \code{ncol(X)} by k.}
#'
#' @seealso \code{\link{tolerance_eigen}}, \code{\link{gsvd}} and \code{\link{gplssvd}}
#'
#' @examples
#'
#' ## (Metric) Multidimensional Scaling
#' data(wine, package="GSVD")
#' D <- as.matrix(dist(wine$objective))^2
#' masses <- rep(1/nrow(D), nrow(D))
#' Xi <- matrix(-masses, length(masses), length(masses))
#' diag(Xi) <- (1-masses)
#' mds.res_geigen <- geigen((-D / (nrow(D) * 2)), Xi)
#'
#' ## Principal components analysis: "covariance"
#' cov_X <- as.matrix(cov(wine$objective))
#' cov_pca.res_geigen <- geigen(cov_X)
#'
#' ## Principal components analysis: "correlation"
#' cor_X <- as.matrix(cor(wine$objective))
#' cor_pca.res_geigen <- geigen(cor_X)
#'
#' @author Derek Beaton
#' @keywords multivariate

geigen <- function(X, W, k = 0, tol= sqrt(.Machine$double.eps), symmetric){

  # preliminaries
  X_dimensions <- dim(X)
  ## stolen from MASS::ginv()
  if (length(X_dimensions) > 2 || !(is.numeric(X) || is.complex(X))){
    stop("geigen: 'X' must be a numeric or complex matrix")
  }
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }

  # check square-ness here.
  if(X_dimensions[1] != X_dimensions[2]){
    stop("geigen: X must be square (i.e., have the same number of rows and columns)")
  }

  # a few things about W for stopping conditions
  W_is_missing <- missing(W)
  if(!W_is_missing){

    W_is_vector <- is.vector(W)

    if(!W_is_vector){

      if( nrow(W) != ncol(W) | nrow(W) != X_dimensions[2] ){
        stop("geigen: nrow(W) does not equal ncol(W) or ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(W)){
        stop("geigen: W is empty (i.e., all 0s")
      }
    }

    if(W_is_vector){
      if(length(W)!=X_dimensions[1]){
        stop("geigen: length(W) does not equal nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(all(abs(W)<=tol)){
        stop("geigen: W is empty (i.e., all 0s")
      }
    }
  }

  ## convenience checks *could* be removed* if problematic
  # convenience checks & conversions; these are meant to minimize W's memory footprint
  if(!W_is_missing){
    if( !W_is_vector) {
      # if( is_identity_matrix(W) ){
      #   W_is_missing <- T
      #   W <- substitute() # neat! this makes it go missing
      # }

      if( !W_is_vector & is_diagonal_matrix(W) ){
        W <- diag(W)
        W_is_vector <- T  # now it's a vector
      }
    }

    if( W_is_vector & all(W==1) ){
      W_is_missing <- T
      W <- substitute() # neat! this makes it go missing
    }
  }

  # this manipulates X as needed
  if(!W_is_missing){
    if( W_is_vector ){

      sqrt_W <- sqrt(W)
      X <- t(t(X * sqrt_W) * sqrt_W)

    }else{

      W <- as.matrix(W)
      # sqrt_W <- W %^% (1/2)
      sqrt_W <- sqrt_psd_matrix(W)
      X <- sqrt_W %*% X %*% sqrt_W

    }
  }

  # all the decomposition things
  if(k<=0){
    k <- min(X_dimensions)
  }

  if(missing(symmetric)){
    symmetric <- isSymmetric(X)
  }

  res <- tolerance_eigen(X, tol=tol, symmetric=symmetric)

  res$l.orig <- res$values
    res$values <- NULL
  res$d.orig <- sqrt(res$l.orig)
  # res$tau <- (res$l.orig/sum(res$l.orig)) * 100

  components_to_return <- min(length(res$d.orig),k)

  res$d <- res$d.orig[1:components_to_return]
  res$l <- res$l.orig[1:components_to_return]
  res$v <- res$vectors[,1:components_to_return, drop = FALSE]
    res$vectors <- NULL




  # make scores according to weights
  if(!W_is_missing){
    if(W_is_vector){

      res$q <- res$v / sqrt_W
      res$fj <- t(t(res$q * W) * res$d)

    }else{

      # res$q <- (W %^% (-1/2)) %*% res$v
      res$q <- invsqrt_psd_matrix(W) %*% res$v
      res$fj <- t(t(W %*% res$q) * res$d)

    }
  }else{

    res$q <- res$v
    res$fj <- t(t(res$q) * res$d)

  }

  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(X)

  class(res) <- c("geigen", "GSVD", "list")
  return(res)

}
