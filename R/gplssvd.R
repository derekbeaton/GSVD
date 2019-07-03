## this provides the basis of all GPLS algorithms. So this would move over into gplssvd as the companion to it.
### but it possibly needs a new name...

### with just a small amount of preprocs, this needs to produce:
### PLS, CCA, RRR, and PLS-CA

gplssvd <- function(X, Y, XLW=rep(1, nrow(ZX)), YLW=rep(1, nrow(ZY)), XRW=rep(1, nrow(ZX)), YRW=rep(1, nrow(ZY)), k = 0, tol = .Machine$double.eps){


  ## this and GSVD and GEIGEN need a bit of cleaning up.

  # preliminaries
  X.dims <- dim(X)
  if(length(X.dims)!=2){
    stop("gplssvd: X must have dim length of 2 (i.e., rows and columns)")
  }
  X <- as.matrix(X)

  Y.dims <- dim(Y)
  if(length(Y.dims)!=2){
    stop("gplssvd: Y must have dim length of 2 (i.e., rows and columns)")
  }
  Y <- as.matrix(Y)

  if(X.dims[1] != Y.dims[1]){
    stop("gplssvd: X and Y must have the same number of rows")
  }

  XRW.is.vector <- XLW.is.vector <- XRW.is.missing <- XLW.is.missing <- F -> YRW.is.vector -> YLW.is.vector -> YRW.is.missing -> YLW.is.missing

  ### These are here out of convenience for the tests below. They started to get too long.
  if( !missing(XLW) ){
    if(is.empty.matrix(XLW)){
      stop("gplssvd: XLW is empty (i.e., all 0s")
    }
  }
  if( !missing(XRW) ){
    if(is.empty.matrix(XRW)){
      stop("gplssvd: XRW is empty (i.e., all 0s")
    }
  }
  if( !missing(YLW) ){
    if(is.empty.matrix(YLW)){
      stop("gplssvd: YLW is empty (i.e., all 0s")
    }
  }
  if( !missing(YRW) ){
    if(is.empty.matrix(YRW)){
      stop("gplssvd: YRW is empty (i.e., all 0s")
    }
  }



  ### oof.
  # check if XLW and XRW are missing, if they are vectors, or if they are diagonal matrices.
  if( missing(XLW) ){
    XLW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(XLW) ) {
      XLW.is.vector <- T
    }else if(!XLW.is.vector){

      if( is.identity.matrix(XLW) ){
        XLW.is.missing <- T
        warning("gplssvd: XLW was an identity matrix. XLW will not be used in the gplssvd.")
      }else if( is.diagonal.matrix(XLW) ){

        XLW <- diag(XLW)

        if( length(XLW) != X.dims[1] ){
          stop("gplssvd:length(XLW) does not equal nrow(X)")
        }else{
          XLW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(XLW) != ncol(XLW) | nrow(XLW) != X.dims[1] ){
        stop("gplssvd:nrow(XLW) does not equal ncol(XLW) or nrow(X)")
      }
    }
  }


  if( missing(XRW) ){
    XRW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(XRW) ) {
      XRW.is.vector <- T
    }else if(!XRW.is.vector){

      if( is.identity.matrix(XRW) ){
        XRW.is.missing <- T
        warning("gplssvd: XRW was an identity matrix. XRW will not be used in the gplssvd.")
      }else if( is.diagonal.matrix(XRW) ){

        XRW <- diag(XRW)

        if( length(XRW) != X.dims[2] ){
          stop("gplssvd:length(XRW) does not equal ncol(X)")
        }else{
          XRW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(XRW) != ncol(XRW) | nrow(XRW) != X.dims[2] ){
        stop("gplssvd:nrow(XRW) does not equal ncol(XRW) or ncol(X)")
      }
    }
  }

  # check if YLW and YRW are missing, if they are vectors, or if they are diagonal matrices.
  if( missing(YLW) ){
    YLW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(YLW) ) {
      YLW.is.vector <- T
    }else if(!YLW.is.vector){

      if( is.identity.matrix(YLW) ){
        YLW.is.missing <- T
        warning("gplssvd: YLW was an identity matrix. YLW will not be used in the gplssvd.")
      }else if( is.diagonal.matrix(YLW) ){

        YLW <- diag(YLW)

        if( length(YLW) != X.dims[1] ){
          stop("gplssvd:length(YLW) does not equal nrow(Y)")
        }else{
          YLW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(YLW) != ncol(YLW) | nrow(YLW) != X.dims[1] ){
        stop("gplssvd:nrow(YLW) does not equal ncol(YLW) or nrow(Y)")
      }
    }
  }


  if( missing(YRW) ){
    YRW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(YRW) ) {
      YRW.is.vector <- T
    }else if(!YRW.is.vector){

      if( is.identity.matrix(YRW) ){
        YRW.is.missing <- T
        warning("gplssvd: YRW was an identity matrix. YRW will not be used in the gplssvd.")
      }else if( is.diagonal.matrix(YRW) ){

        YRW <- diag(YRW)

        if( length(YRW) != X.dims[2] ){
          stop("gplssvd:length(YRW) does not equal ncol(Y)")
        }else{
          YRW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(YRW) != ncol(YRW) | nrow(YRW) != X.dims[2] ){
        stop("gplssvd:nrow(YRW) does not equal ncol(YRW) or ncol(Y)")
      }
    }
  }



  if(!XLW.is.missing){
    if( XLW.is.vector ){  ## replace with sweep
      X <- sweep(X,1,sqrt(XLW),"*")
    }else{
      X <- (XLW %^% (1/2)) %*% X
    }
  }

  if(!XRW.is.missing){
    if( XRW.is.vector ){  ## replace with sweep
      X <- sweep(X,2,sqrt(XRW),"*")
    }else{
      X <- X %*% (XRW %^% (1/2))
    }
  }


  if(!YLW.is.missing){
    if( YLW.is.vector ){  ## replace with sweep
      Y <- sweep(Y,1,sqrt(YLW),"*")
    }else{
      Y <- (YLW %^% (1/2)) %*% Y
    }
  }

  if(!YRW.is.missing){
    if( YRW.is.vector ){  ## replace with sweep
      Y <- sweep(Y,2,sqrt(YRW),"*")
    }else{
      Y <- Y %*% (YRW %^% (1/2))
    }
  }


  if(k<=0){
    k <- min(X.dims, Y.dims)
  }

  res <- tolerance.svd(t(X) %*% Y, nu=k, nv=k, tol=tol)
  res$d.orig <- res$d
  res$l.orig <- res$d.orig^2
  res$tau <- (res$l.orig/sum(res$l.orig)) * 100
  components.to.return <- min(length(res$d.orig),k) #a safety check
  ## u and v should already be k vectors but: be safe.
  res$d <- res$d.orig[1:components.to.return]
  res$l <- res$d^2
  res$u <- as.matrix(res$u[,1:components.to.return])
  res$v <- as.matrix(res$v[,1:components.to.return])

  res$lx <- X %*% res$u
  res$ly <- Y %*% res$v



  ## the logic here should match the one from above
  if(!XRW.is.missing){
    if(XRW.is.vector){
      res$p <- sweep(res$u,1,1/sqrt(XRW),"*")
      res$fi <- sweep(sweep(res$p,1,XRW,"*"),2,res$d,"*")
    }else{
      res$p <- (XRW %^% (-1/2)) %*% res$u
      res$fi <- sweep((XRW %*% res$p),2,res$d,"*")
    }
  }else{
    res$p <- res$u
    res$fi <- sweep(res$p,2,res$d,"*")
  }

  if(!YRW.is.missing){
    if(YRW.is.vector){
      res$q <- sweep(res$v,1,1/sqrt(YRW),"*")
      res$fj <- sweep(sweep(res$q,1,YRW,"*"),2,res$d,"*")
    }else{
      res$q <- (YRW %^% (-1/2)) %*% res$v
      res$fj <- sweep((YRW %*% res$q),2,res$d,"*")
    }
  }else{
    res$q <- res$v
    res$fj <- sweep(res$q,2,res$d,"*")
  }

  rownames(res$fi) <- rownames(res$u) <- rownames(res$p) <- colnames(X)
  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(Y)

  rownames(res$lx) <- rownames(X)
  rownames(res$ly) <- rownames(Y)

  class(res) <- c("list", "GSVD", "gplssvd")
  return(res)

}
