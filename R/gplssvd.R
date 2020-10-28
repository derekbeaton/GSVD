#' @export
#'
#' @title Generalized partial least squares singular value decomposition
#'
#' @description
#' \code{gplssvd} takes in left (\code{XLW}, \code{YLW}) and right (\code{XRW}, \code{YRW}) constraints (usually diagonal matrices, but any positive semi-definite matrix is fine) that are applied to each data matrix (\code{X} and \code{Y}), respectively.
#'   Right constraints for each data matrix are used for the orthogonality conditions.
#'
#' @param X a data matrix
#' @param Y a data matrix
#' @param XLW \bold{X}'s \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the \code{X} matrix and thus left singular vectors.
#' @param YLW \bold{Y}'s \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the \code{Y} matrix and thus left singular vectors.
#' @param XRW \bold{X}'s \bold{R}ight \bold{W}eights -- the constraints applied to the right side (columns) of the \code{X} matrix and thus right singular vectors.
#' @param YRW \bold{Y}'s \bold{R}ight \bold{W}eights -- the constraints applied to the right side (columns) of the \code{Y} matrix and thus right singular vectors.
#' @param k total number of components to return though the full variance will still be returned (see \code{d_full}). If 0, the full set of components are returned.
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components.
#'
#' @return A list with thirteen elements:
#' \item{d_full}{A vector containing the singular values of DAT above the tolerance threshold (based on eigenvalues).}
#' \item{l_full}{A vector containing the eigen values of DAT above the tolerance threshold (\code{tol}).}
#' \item{d}{A vector of length \code{min(length(d_full), k)} containing the retained singular values}
#' \item{l}{A vector of length \code{min(length(l_full), k)} containing the retained eigen values}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(DAT)} by k.}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{v}{Right (columns) singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fj}{Right (columns) component scores.}
#' \item{ly}{Latent variable scores for rows of \code{Y}}
#'
#' @seealso \code{\link{tolerance_svd}}, \code{\link{geigen}} and \code{\link{gsvd}}
#'
#' @examples
#'
#'  # Three "two-table" technique examples
#'  data(wine)
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'
#'  ## Partial least squares (correlation)
#'  pls.res <- gplssvd(X, Y)
#'
#'
#'  ## Canonical correlation analysis (CCA)
#'  ### NOTE:
#'  #### This is not "traditional" CCA because of the generalized inverse.
#'  #### However the results are the same as standard CCA when data are not rank deficient.
#'  #### and this particular version uses tricks to minimize memory & computation
#'  cca.res <- gplssvd(
#'      X = MASS::ginv(t(X)),
#'      Y = MASS::ginv(t(Y)),
#'      XRW=crossprod(X),
#'      YRW=crossprod(Y)
#'  )
#'
#'
#'  ## Reduced rank regression (RRR) a.k.a. redundancy analysis (RDA)
#'  ### NOTE:
#'  #### This is not "traditional" RRR because of the generalized inverse.
#'  #### However the results are the same as standard RRR when data are not rank deficient.
#'  rrr.res <- gplssvd(
#'      X = MASS::ginv(t(X)),
#'      Y = Y,
#'      XRW=crossprod(X)
#'  )
#'
#' @author Derek Beaton
#' @keywords multivariate

gplssvd <- function(X, Y, XLW, YLW, XRW, YRW, k = 0, tol = .Machine$double.eps){

  # preliminaries
  X_dimensions <- dim(X)
  ## stolen from MASS::ginv()
  if (length(X_dimensions) > 2 || !(is.numeric(X) || is.complex(X))){
    stop("gplssvd: 'X' must be a numeric or complex matrix")
  }
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }

  # preliminaries
  Y_dimensions <- dim(Y)
  ## stolen from MASS::ginv()
  if (length(Y_dimensions) > 2 || !(is.numeric(Y) || is.complex(Y))){
    stop("gplssvd: 'Y' must be a numeric or complex matrix")
  }
  if (!is.matrix(Y)){
    Y <- as.matrix(Y)
  }

  # check that row dimensions match
  if(X_dimensions[1] != Y_dimensions[1]){
    stop("gplssvd: X and Y must have the same number of rows")
  }

  # a few things about XLW for stopping conditions
  XLW_is_missing <- missing(XLW)
  if(!XLW_is_missing){

    XLW_is_vector <- is.vector(XLW)

    if(!XLW_is_vector){

      if( nrow(XLW) != ncol(XLW) | nrow(XLW) != X_dimensions[1] ){
        stop("gplssvd: nrow(XLW) does not equal ncol(XLW) or nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(XLW)){
        stop("gplssvd: XLW is empty (i.e., all 0s")
      }
    }

    if(XLW_is_vector){
      if(length(XLW)!=X_dimensions[1]){
        stop("gplssvd: length(XLW) does not equal nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(XLW)<=tol)){
      if(!are_all_values_positive(XLW)){
        stop("gplssvd: XLW is not strictly positive values")
      }
    }
  }

  # a few things about XRW for stopping conditions
  XRW_is_missing <- missing(XRW)
  if(!XRW_is_missing){

    XRW_is_vector <- is.vector(XRW)

    if(!XRW_is_vector){

      if( nrow(XRW) != ncol(XRW) | nrow(XRW) != X_dimensions[2] ){
        stop("gplssvd: nrow(XRW) does not equal ncol(XRW) or ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(XRW)){
        stop("gplssvd: XRW is empty (i.e., all 0s")
      }
    }

    if(XRW_is_vector){
      if(length(XRW)!=X_dimensions[2]){
        stop("gplssvd: length(XRW) does not equal ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(XRW)<=tol)){
      if(!are_all_values_positive(XRW)){
        stop("gplssvd: XRW is not strictly positive values")
      }
    }
  }

  # a few things about YLW for stopping conditions
  YLW_is_missing <- missing(YLW)
  if(!YLW_is_missing){

    YLW_is_vector <- is.vector(YLW)

    if(!YLW_is_vector){

      if( nrow(YLW) != ncol(YLW) | nrow(YLW) != Y_dimensions[1] ){
        stop("gplssvd: nrow(YLW) does not equal ncol(YLW) or nrow(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(YLW)){
        stop("gplssvd: YLW is empty (i.e., all 0s")
      }
    }

    if(YLW_is_vector){
      if(length(YLW)!=Y_dimensions[1]){
        stop("gplssvd: length(YLW) does not equal nrow(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(YLW)<=tol)){
      if(!are_all_values_positive(YLW)){
        stop("gplssvd: YLW is not strictly positive values")
      }
    }
  }

  # a few things about YRW for stopping conditions
  YRW_is_missing <- missing(YRW)
  if(!YRW_is_missing){

    YRW_is_vector <- is.vector(YRW)

    if(!YRW_is_vector){

      if( nrow(YRW) != ncol(YRW) | nrow(YRW) != Y_dimensions[2] ){
        stop("gplssvd: nrow(YRW) does not equal ncol(YRW) or ncol(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(YRW)){
        stop("gplssvd: YRW is empty (i.e., all 0s")
      }
    }

    if(YRW_is_vector){
      if(length(YRW)!=Y_dimensions[2]){
        stop("gplssvd: length(YRW) does not equal ncol(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(YRW)<=tol)){
      if(!are_all_values_positive(YRW)){
        stop("gplssvd: YRW is not strictly positive values")
      }
    }
  }

  ## convenience checks *could* be removed* if problematic
  # convenience checks & conversions; these are meant to minimize XLW's memory footprint
  if(!XLW_is_missing){

    ## this is a matrix call; THIS MAKES IT GO MISSING AND NEXT CHECKS FAIL.
      ### technically, the following two checks actually handle this.
      ### it is a diagonal, and then it's a vector that's all 1s
    if( !XLW_is_vector ){

      # if(is_identity_matrix(XLW) ){
      #   XLW_is_missing <- T
      #   XLW <- substitute() # neat! this makes it go missing
      # }

      ## this is a matrix call
      if( is_diagonal_matrix(XLW) ){
        XLW <- diag(XLW)
        XLW_is_vector <- T  # now it's a vector
      }
    }

    ## this is a vector call
    if( XLW_is_vector & all(XLW==1) ){
      XLW_is_missing <- T
      XLW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize XRW's memory footprint
  if(!XRW_is_missing){
    if( !XRW_is_vector ){

      # if( is_identity_matrix(XRW) ){
      #   XRW_is_missing <- T
      #   XRW <- substitute() # neat! this makes it go missing
      # }

      if( is_diagonal_matrix(XRW) ){
        XRW <- diag(XRW)
        XRW_is_vector <- T  # now it's a vector
      }
    }

    if( XRW_is_vector & all(XRW==1) ){
      XRW_is_missing <- T
      XRW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize YLW's memory footprint
  if(!YLW_is_missing){
    if( !YLW_is_vector ){

      # if( is_identity_matrix(YLW) ){
      #   YLW_is_missing <- T
      #   YLW <- substitute() # neat! this makes it go missing
      # }

      if( is_diagonal_matrix(YLW) ){
        YLW <- diag(YLW)
        YLW_is_vector <- T  # now it's a vector
      }

    }

    if( YLW_is_vector & all(YLW==1) ){
      YLW_is_missing <- T
      YLW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize YRW's memory footprint
  if(!YRW_is_missing){
    if( !YRW_is_vector ){

      # if( is_identity_matrix(YRW) ){
      #   YRW_is_missing <- T
      #   YRW <- substitute() # neat! this makes it go missing
      # }

      if( is_diagonal_matrix(YRW) ){
        YRW <- diag(YRW)
        YRW_is_vector <- T  # now it's a vector
      }

    }

    if( YRW_is_vector & all(YRW==1) ){
      YRW_is_missing <- T
      YRW <- substitute() # neat! this makes it go missing
    }
  }


  # this manipulates X as needed based on XLW
  if(!XLW_is_missing){

    if( XLW_is_vector ){
      # X <- sweep(X,1,sqrt(XLW),"*") ## replace the sweep with * & t()
      X <- X * sqrt(XLW)
    }else{
      XLW <- as.matrix(XLW)
      # X <- (XLW %^% (1/2)) %*% X
      X <- sqrt_psd_matrix(XLW) %*% X
    }

  }
  # this manipulates X as needed based on XRW
  if(!XRW_is_missing){

    if( XRW_is_vector ){
      sqrt_XRW <- sqrt(XRW)
      # X <- sweep(X,2, sqrt_XRW,"*") ## replace the sweep with * & t()
      X <- t(t(X) * sqrt_XRW)
    }else{
      XRW <- as.matrix(XRW)
      # X <- X %*% (XRW %^% (1/2))
      X <- X %*% sqrt_psd_matrix(XRW)
    }

  }
  # this manipulates Y as needed based on YLW
  if(!YLW_is_missing){

    if( YLW_is_vector ){
      # Y <- sweep(Y,1,sqrt(YLW),"*")  ## replace the sweep with * & t()
      Y <- Y * sqrt(YLW)
    }else{
      YLW <- as.matrix(YLW)
      Y <- sqrt_psd_matrix(YLW) %*% Y
    }

  }
  # this manipulates Y as needed based on YRW
  if(!YRW_is_missing){

    if( YRW_is_vector ){
      sqrt_YRW <- sqrt(YRW)
      # Y <- sweep(Y,2,sqrt_YRW,"*")  ## replace the sweep with * & t()
      Y <- t(t(Y) * sqrt_YRW)
    }else{
      YRW <- as.matrix(YRW)
      # Y <- Y %*% (YRW %^% (1/2))
      Y <- Y %*% sqrt_psd_matrix(YRW)
    }

  }


  # all the decomposition things
  if(k<=0){
    k <- min(X_dimensions, Y_dimensions)
  }

  res <- tolerance_svd( t(X) %*% Y, nu=k, nv=k, tol=tol)
  res$d_full <- res$d
  res$l_full <- res$d_full^2
  # res$tau <- (res$l_full/sum(res$l_full)) * 100
  components.to.return <- min(length(res$d_full),k) #a safety check
  res$d <- res$d_full[1:components.to.return]
  res$l <- res$d^2
  res$u <- res$u[,1:components.to.return, drop = FALSE]
  res$v <- res$v[,1:components.to.return, drop = FALSE]

  res$lx <- X %*% res$u
  res$ly <- Y %*% res$v


  # make scores according to weights
  if(!XRW_is_missing){
    if(XRW_is_vector){

      # res$p <- sweep(res$u,1,1/sqrt_XRW,"*")
      res$p <- res$u / sqrt_XRW
      # res$fi <- sweep(sweep(res$p,1,XRW,"*"),2,res$d,"*")
      res$fi <- t(t(res$p * XRW) * res$d)

    }else{

      # res$p <- (XRW %^% (-1/2)) %*% res$u
      res$p <- invsqrt_psd_matrix(XRW) %*% res$u
      # res$fi <- sweep((XRW %*% res$p),2,res$d,"*")
      res$fi <- t(t(XRW %*% res$p) * res$d)

    }
  }else{

    res$p <- res$u
    # res$fi <- sweep(res$p,2,res$d,"*")
    res$fi <- t(t(res$p) * res$d)

  }

  if(!YRW_is_missing){
    if(YRW_is_vector){

      # res$q <- sweep(res$v,1,1/sqrt_YRW,"*")
      res$q <- res$v /sqrt_YRW
      # res$fj <- sweep(sweep(res$q,1,YRW,"*"),2,res$d,"*")
      res$fj <- t(t(res$q * YRW) * res$d)

    }else{

      # res$q <- (YRW %^% (-1/2)) %*% res$v
      res$q <- invsqrt_psd_matrix(YRW) %*% res$v
      # res$fj <- sweep((YRW %*% res$q),2,res$d,"*")
      res$fj <- t(t(YRW %*% res$q) * res$d)

    }
  }else{

    res$q <- res$v
    # res$fj <- sweep(res$q,2,res$d,"*")
    res$fj <- t(t(res$q) * res$d)

  }

  rownames(res$fi) <- rownames(res$u) <- rownames(res$p) <- colnames(X)
  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(Y)

  rownames(res$lx) <- rownames(X)
  rownames(res$ly) <- rownames(Y)

  class(res) <- c("gplssvd", "GSVD", "list")
  return(res)

}
