#' GSVD
#'
#' @description The generalized singular value decomposition (GSVD) generalizes the standard SVD (see \code{\link{svd}}) procedure through addition of (optional) constraints applied to the rows and/or columns of a matrix.
#' @details While dedicated to the GSVD, this package also includes some nice features that are helpful for matrix analyses. For example there are tests to check if matrices are empty (\code{\link{is.empty.matrix}}) or identity (\code{\link{is.identity.matrix}}); also included are operations such as generalized inverse (\code{\link{matrix.generalized.inverse}}) and matrix exponents (\code{\link{matrix.exponent}} or \code{\link{\%^\%}}).
#' @seealso \code{\link{gsvd}}, \code{\link{geigen}}, \code{\link{tolerance.svd}}, \code{\link{\%^\%}}
#' @examples
#'  ## an example of correspondence analysis.
#'  data(authors)
#'  Observed <- authors/sum(authors)
#'  row.w <- rowSums(Observed)
#'    row.W <- diag(1/row.w)
#'  col.w <- colSums(Observed)
#'    col.W <- diag(1/col.w)
#'  Expected <- row.w %o% col.w
#'  Deviations <- Observed - Expected
#'  ca.res <- gsvd(Deviations,row.W,col.W)
#'
#'
#'  # several examples of principal component analysis
#'  data(wine)
#'  wine.objective <- wine$objective
#'  ## "covariance" PCA
#'  cov.pca.data <- scale(wine.objective,scale=FALSE)
#'  cov.pca.res <- gsvd(cov.pca.data)
#'  ## "correlation" PCA
#'  cor.pca.data <- scale(wine.objective,scale=TRUE)
#'  cor.pca.res <- gsvd(cor.pca.data)
#'  ## "correlation" PCA with GSVD constraints
#'  cor.pca.res2 <- gsvd(cov.pca.data,RW=1/apply(wine.objective,2,var))
#'
#'
#'  # Three "two-table" technique examples
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'  R <- t(X) %*% Y
#'
#'  ## an example of partial least squares (correlation)
#'  pls.res <- gsvd(R)
#'
#'
#'  ## an example of canonical correlation analysis (CCA)
#'  ### NOTE:
#'  #### This is not "traditional" CCA because of the generalized inverse.
#'  #### However results are the same as standard CCA when data are not rank deficient.
#'  cca.res <- gsvd(
#'      DAT=(crossprod(X) %^% -1) %*% R %*% (crossprod(Y) %^% -1),
#'      LW=crossprod(X),
#'      RW=crossprod(Y)
#'  )
#'  cca.res$lx <- (X %*% cca.res$p)
#'  cca.res$ly <- (Y %*% cca.res$q)
#'
#'  ## an example of reduced rank regression (RRR) a.k.a. redundancy analysis (RDA)
#'  ### NOTE:
#'  #### This is not "traditional" RRR because of the generalized inverse.
#'  #### However the results are the same as standard RRR when data are not rank deficient.
#'  rrr.res <- gsvd(
#'      DAT=(crossprod(X) %^% -1) %*% R,
#'      LW=crossprod(X)
#'  )
#'  rrr.res$lx <- (X %*% rrr.res$p)
#'  rrr.res$ly <- (Y %*% rrr.res$q)
#'
#' @references
#'   Abdi, H. (2007). Singular Value Decomposition (SVD) and Generalized Singular Value Decomposition (GSVD). In N.J. Salkind (Ed.): \emph{Encyclopedia of Measurement and Statistics}.Thousand Oaks (CA): Sage. pp. 907-912.\cr
#'   Beaton, D., Fatt, C. R. C., & Abdi, H. (2014). An ExPosition of multivariate analysis with the singular value decomposition in R. \emph{Computational Statistics & Data Analysis}, 72, 176-189.
#'   Greenacre, M. (1984). Theory and applications of correspondence analysis. \emph{Academic Press.}\cr
#'   Holmes, S. (2008). Multivariate data analysis: the French way. \emph{In Probability and statistics: Essays in honor of David A. Freedman (pp. 219-233)}. Institute of Mathematical Statistics. \cr
#'   Holmes, S., & Josse, J. (2017). Discussion of “50 Years of Data Science”. \emph{Journal of Computational and Graphical Statistics}, 26(4), 768-769. \cr
#'   Husson, F., Josse, J., & Saporta, G. (2016). Jan de Leeuw and the French school of data analysis. \emph{Journal of Statistical Software}, 73(6), 1-17. \cr
#'   Lebart  L.,  Morineau  A., & Warwick  K.  (1984). Multivariate  Descriptive  Statistical  Analysis. \emph{J. Wiley}, New York. \cr
#'   Yanai, H., Takeuchi, K., & Takane, Y. (2011). Projection Matrices, Generalized Inverse Matrices, and Singular Value Decomposition. \emph{Springer-Verlag, New-York.}\cr
#'
#' @keywords multivariate svd genearlized matrix decomposition variance component orthogonal
#'
"_PACKAGE"

#' @export
#'
#' @title Generalized SVD
#'
#' @description
#' \code{gsvd} takes in left (\code{LW}) and right (\code{RW}) constraints (usually diagonal matrices, but any positive semi-definite matrix is fine) that are applied to the data (\code{DAT})
#'   Left and right constraints are used for the orthogonality conditions.
#'
#' @param DAT a data matrix to decompose
#' @param LW \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the matrix and thus left singular vectors.
#' @param RW \bold{R}ight \bold{W}eights -- the constraints applied to the right side (columns) of the matrix and thus right singular vectors.
#' @param k total number of components to return though the full variance will still be returned (see \code{d.orig}). If 0, the full set of components are returned.
#' @param tol default is .Machine$double.eps. A parameter with two roles: A tolerance level for (1) eliminating (tiny variance or negative or imaginary) components and (2) converting all values < tol to 0 in \code{u} and \code{v}.
#'
#' @return A list with nine elements:
#' \item{d.orig}{A vector containing the singular values of DAT above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of DATabove the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}.}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values of DAT}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values of DAT}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(DAT)} by k.}
#' \item{p}{Left (rows) generalized singular vectors. Dimensions are \code{nrow(DAT)} by k.}
#' \item{fi}{Left (rows) component scores. Dimensions are \code{nrow(DAT)} by k.}
#' \item{v}{Right (columns) singular vectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{q}{Right (columns) generalized singular vectors. Dimensions are \code{ncol(DAT)} by k.}
#' \item{fj}{Right (columns) component scores. Dimensions are \code{ncol(DAT)} by k.}
#'
#' @seealso \code{\link{tolerance.svd}} and \code{\link{svd}}
#'
#' @examples
#'
#'  ## an example of correspondence analysis.
#'  data(authors)
#'  Observed <- authors/sum(authors)
#'  row.w <- rowSums(Observed)
#'    row.W <- diag(1/row.w)
#'  col.w <- colSums(Observed)
#'    col.W <- diag(1/col.w)
#'  Expected <- row.w %o% col.w
#'  Deviations <- Observed - Expected
#'  ca.res <- gsvd(Deviations,row.W,col.W)
#'
#'
#'  # several examples of principal component analysis
#'  data(wine)
#'  wine.objective <- wine$objective
#'  ## "covariance" PCA
#'  cov.pca.data <- scale(wine.objective,scale=FALSE)
#'  cov.pca.res <- gsvd(cov.pca.data)
#'  ## "correlation" PCA
#'  cor.pca.data <- scale(wine.objective,scale=TRUE)
#'  cor.pca.res <- gsvd(cor.pca.data)
#'  ## "correlation" PCA with GSVD constraints
#'  cor.pca.res2 <- gsvd(cov.pca.data,RW=1/apply(wine.objective,2,var))
#'
#'
#'  # Three "two-table" technique examples
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'  R <- t(X) %*% Y
#'
#'  ## an example of partial least squares (correlation)
#'  pls.res <- gsvd(R)
#'
#'
#'  ## an example of canonical correlation analysis (CCA)
#'  ### NOTE:
#'  #### This is not "traditional" CCA because of the generalized inverse.
#'  #### However the results are the same as standard CCA when data are not rank deficient.
#'  cca.res <- gsvd(
#'      DAT=(crossprod(X) %^% -1) %*% R %*% (crossprod(Y) %^% -1),
#'      LW=crossprod(X),
#'      RW=crossprod(Y)
#'  )
#'  cca.res$lx <- (X %*% cca.res$p)
#'  cca.res$ly <- (Y %*% cca.res$q)
#'
#'
#'  ## an example of reduced rank regression (RRR) a.k.a. redundancy analysis (RDA)
#'  ### NOTE:
#'  #### This is not "traditional" RRR because of the generalized inverse.
#'  #### However the results are the same as standard RRR when data are not rank deficient.
#'  rrr.res <- gsvd(
#'      DAT=(crossprod(X) %^% -1) %*% R,
#'      LW=crossprod(X)
#'  )
#'  rrr.res$lx <- (X %*% rrr.res$p)
#'  rrr.res$ly <- (Y %*% rrr.res$q)
#'
#' @author Derek Beaton
#' @keywords multivariate, diagonalization, eigen


### I need to do some PSD checks.

gsvd <- function(DAT, LW, RW, k = 0, tol = .Machine$double.eps){

  # preliminaries
  DAT.dims <- dim(DAT)
  if(length(DAT.dims)!=2){
    stop("gsvd: DAT must have dim length of 2 (i.e., rows and columns)")
  }
  DAT <- as.matrix(DAT)
  RW.is.vector <- LW.is.vector <- RW.is.missing <- LW.is.missing <- F ##asuming everything is a matrix.

  ### These are here out of convenience for the tests below. They started to get too long.
  if( !missing(LW) ){
    if(is.empty.matrix(LW)){
      stop("gsvd: LW is empty (i.e., all 0s")
    }
  }
  if( !missing(RW) ){
    if(is.empty.matrix(RW)){
      stop("gsvd: RW is empty (i.e., all 0s")
    }
  }

  # check if LW and RW are missing, if they are vectors, or if they are diagonal matrices.
  if( missing(LW) ){
    LW.is.missing <- T
  }else{ # it's here and we have to check!

    if ( is.vector(LW) ) {
      LW.is.vector <- T
    }else if(!LW.is.vector){

      if( is.identity.matrix(LW) ){
        LW.is.missing <- T
        warning("gsvd: LW was an identity matrix. LW will not be used in the GSVD.")
      }else if( is.diagonal.matrix(LW) ){

        LW <- diag(LW)

        if( length(LW) != DAT.dims[1] ){
          stop("gsvd:length(LW) does not equal nrow(DAT)")
        }else{
          LW.is.vector <- T  #now it's a vector
        }

      }else if( nrow(LW) != ncol(LW) | nrow(LW) != DAT.dims[1] ){
        stop("gsvd:nrow(LW) does not equal ncol(LW) or nrow(DAT)")
      }
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



  if(!LW.is.missing){
    if( LW.is.vector ){  ## replace with sweep
      DAT <- sweep(DAT,1,sqrt(LW),"*")
    }else{
      DAT <- (LW %^% (1/2)) %*% DAT
    }
  }

  if(!RW.is.missing){
    if( RW.is.vector ){  ## replace with sweep
      DAT <- sweep(DAT,2,sqrt(RW),"*")
    }else{
      DAT <- DAT %*% (RW %^% (1/2))
    }
  }


  if(k<=0){
    k <- min(nrow(DAT),ncol(DAT))
  }
  res <- tolerance.svd(DAT,nu=k,nv=k,tol=tol)
  res$d.orig <- res$d
  res$l.orig <- res$d.orig^2
  res$tau <- (res$l.orig/sum(res$l.orig)) * 100
  components.to.return <- min(length(res$d.orig),k) #a safety check
  ## u and v should already be k vectors but: be safe.
  res$d <- res$d.orig[1:components.to.return]
  res$l <- res$d^2
  res$u <- as.matrix(res$u[,1:components.to.return])
  res$v <- as.matrix(res$v[,1:components.to.return])



  ## the logic here should match the one from above
  if(!LW.is.missing){
    if(LW.is.vector){
      res$p <- sweep(res$u,1,1/sqrt(LW),"*")
      res$fi <- sweep(sweep(res$p,1,LW,"*"),2,res$d,"*")
    }else{
      res$p <- (LW %^% (-1/2)) %*% res$u
      res$fi <- sweep((LW %*% res$p),2,res$d,"*")
    }
  }else{
    res$p <- res$u
    res$fi <- sweep(res$p,2,res$d,"*")
  }

  if(!RW.is.missing){
    if(RW.is.vector){
      res$q <- sweep(res$v,1,1/sqrt(RW),"*")
      res$fj <- sweep(sweep(res$q,1,RW,"*"),2,res$d,"*")
    }else{
      res$q <- (RW %^% (-1/2)) %*% res$v
      res$fj <- sweep((RW %*% res$q),2,res$d,"*")
    }
  }else{
    res$q <- res$v
    res$fj <- sweep(res$q,2,res$d,"*")
  }


  rownames(res$fi) <- rownames(res$u) <- rownames(res$p) <- rownames(DAT)
  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(DAT)

  return(res)
  #return(list(fi = fi, fj = fj, p = p, q = q, u = res$u, v = res$v, d = d, d.orig = d.orig, tau = tau))
}
