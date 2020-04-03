#' GSVD
#'
#' @description The generalized singular value decomposition (GSVD) generalizes the standard SVD (see \code{\link{svd}}) procedure through addition of (optional) constraints applied to the rows and/or columns of a matrix.
#' @details A package specifically designed for the generalized eigen decomposition, generalized singular value decomposition, and generalized partial least squares-singular value decomposition.
#' Each decomposition allows for the use of weights (or constraints) applied to the columns and/or rows of each data matrix.
#' This package provides these decompositions as the core for a very large list of standard statistical---and typically multivariate---approaches, including but not limited to principal components analysis, metric multidimensional scaling, standard and multiple correspondence analysis, partial least squares, canonical correlation analysis, and multivariate linear regression/reduced rank regression/redundancy analysis.
#' @seealso \code{\link{gsvd}}, \code{\link{geigen}}, \code{\link{gplssvd}}
#'
#' @references
#'   Abdi, H. (2007). Singular Value Decomposition (SVD) and Generalized Singular Value Decomposition (GSVD). In N.J. Salkind (Ed.): \emph{Encyclopedia of Measurement and Statistics}.Thousand Oaks (CA): Sage. pp. 907-912.\cr
#'   Beaton, D., Fatt, C. R. C., & Abdi, H. (2014). An ExPosition of multivariate analysis with the singular value decomposition in R. \emph{Computational Statistics & Data Analysis}, 72, 176-189. \cr
#'   Cherry, S. (1996). Singular Value Decomposition Analysis and Canonical Correlation Analysis. Journal of Climate, 9(9), 2003–2009. Retrieved from JSTOR. \cr
#'   Gower, J. C., Gardner-Lubbe, S., & Le Roux, N. J. DATA ANALYSIS: GOOD BUT…. Statistica Applicata - Italian Journal of Applied Statistics, 29. \cr
#'   Greenacre, M. (1984). Theory and applications of correspondence analysis. \emph{Academic Press.}\cr
#'   Holmes, S. (2008). Multivariate data analysis: the French way. \emph{In Probability and statistics: Essays in honor of David A. Freedman (pp. 219-233)}. Institute of Mathematical Statistics. \cr
#'   Holmes, S., & Josse, J. (2017). Discussion of “50 Years of Data Science”. \emph{Journal of Computational and Graphical Statistics}, 26(4), 768-769. \cr
#'   Husson, F., Josse, J., & Saporta, G. (2016). Jan de Leeuw and the French school of data analysis. \emph{Journal of Statistical Software}, 73(6), 1-17. \cr
#'   Jolliffe, I. T., & Cadima, J. (2016). Principal component analysis: A review and recent developments. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 374(2065).\cr
#'   Lebart  L.,  Morineau  A., & Warwick  K.  (1984). Multivariate  Descriptive  Statistical  Analysis. \emph{J. Wiley}, New York. \cr
#'   Nguyen, L. H., & Holmes, S. (2019). Ten quick tips for effective dimensionality reduction. PLOS Computational Biology, 15(6), e1006907. https://doi.org/10.1371/journal.pcbi.1006907 \cr
#'   Yanai, H., Takeuchi, K., & Takane, Y. (2011). Projection Matrices, Generalized Inverse Matrices, and Singular Value Decomposition. \emph{Springer-Verlag, New-York.}\cr
#'
#' @keywords internal
#'
"_PACKAGE"

#' @export
#'
#' @title Generalized singular value decomposition
#'
#' @description
#' \code{gsvd} takes in left (\code{LW}) and right (\code{RW}) constraints (usually diagonal matrices, but any positive semi-definite matrix is fine) that are applied to the data (\code{X}).
#'   Left and right constraints are used for the orthogonality conditions.
#'
#' @param X a data matrix to decompose
#' @param LW \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the matrix and thus left singular vectors.
#' @param RW \bold{R}ight \bold{W}eights -- the constraints applied to the right side (columns) of the matrix and thus right singular vectors.
#' @param k total number of components to return though the full variance will still be returned (see \code{d.orig}). If 0, the full set of components are returned.
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components.
#'
#' @return A list with eleven elements:
#' \item{d.orig}{A vector containing the singular values of X above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of X above the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}).}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values of X}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values of X}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(X)} by k.}
#' \item{p}{Left (rows) generalized singular vectors. Dimensions are \code{nrow(X)} by k.}
#' \item{fi}{Left (rows) component scores. Dimensions are \code{nrow(X)} by k.}
#' \item{v}{Right (columns) singular vectors. Dimensions are \code{ncol(X)} by k.}
#' \item{q}{Right (columns) generalized singular vectors. Dimensions are \code{ncol(X)} by k.}
#' \item{fj}{Right (columns) component scores. Dimensions are \code{ncol(X)} by k.}
#'
#' @seealso \code{\link{tolerance_svd}}, \code{\link{geigen}} and \code{\link{gplssvd}}
#'
#' @examples
#'
#'  data(wine, package="GSVD")
#'  wine.objective <- wine$objective
#'  ## Principal components analysis: "covariance"
#'  ## "covariance" PCA
#'  cov.pca.data <- scale(wine.objective,scale=FALSE)
#'  cov.pca.res <- gsvd(cov.pca.data)
#'
#'  ## Principal components analysis: "correlation"
#'  cor.pca.data <- scale(wine.objective,scale=TRUE)
#'  cor.pca.res <- gsvd(cor.pca.data)
#'
#'  ## Principal components analysis: "correlation" with the covariance matrix and constraints
#'  cor.pca.res2 <- gsvd(cov.pca.data,RW=1/apply(wine.objective,2,var))
#'
#'
#'  ## Correspondence analysis
#'  data(authors, package="GSVD")
#'  Observed <- authors/sum(authors)
#'  row.w <- rowSums(Observed)
#'  col.w <- colSums(Observed)
#'  Expected <- row.w %o% col.w
#'  Deviations <- Observed - Expected
#'  ca.res <- gsvd(Deviations,1/row.w,1/col.w)
#'
#'
#'  ## Multiple correspondence analysis
#'  data("snps.druguse", package="GSVD")
#'  X <- model.matrix(~ .,
#'      data=snps.druguse$DATA1,
#'      contrasts.arg = lapply(snps.druguse$DATA1, contrasts, contrasts=FALSE))[,-1]
#'  Observed <- X/sum(X)
#'  row.w <- rowSums(Observed)
#'  col.w <- colSums(Observed)
#'  Expected <- row.w %o% col.w
#'  Deviations <- Observed - Expected
#'  ca.res <- gsvd(Deviations,1/row.w,1/col.w)
#'
#' @author Derek Beaton
#' @keywords multivariate


gsvd <- function(X, LW, RW, k = 0, tol = .Machine$double.eps){

  # preliminaries
  X_dimensions <- dim(X)
  ## stolen from MASS::ginv()
  if (length(X_dimensions) > 2 || !(is.numeric(X) || is.complex(X))){
    stop("gsvd: 'X' must be a numeric or complex matrix")
  }
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }


  # a few things about LW for stopping conditions
  LW_is_missing <- missing(LW)
  if(!LW_is_missing){

    LW_is_vector <- is.vector(LW)

    if(!LW_is_vector){

      if( nrow(LW) != ncol(LW) | nrow(LW) != X_dimensions[1] ){
        stop("gsvd: nrow(LW) does not equal ncol(LW) or nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(LW)){
        stop("gsvd: LW is empty (i.e., all 0s")
      }
    }

    if(LW_is_vector){
      if(length(LW)!=X_dimensions[1]){
        stop("gsvd: length(LW) does not equal nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(all(abs(LW)<=tol)){
        stop("gsvd: LW is empty (i.e., all 0s")
      }
    }
  }

  # a few things about RW for stopping conditions
  RW_is_missing <- missing(RW)
  if(!RW_is_missing){

    RW_is_vector <- is.vector(RW)

    if(!RW_is_vector){

      if( nrow(RW) != ncol(RW) | nrow(RW) != X_dimensions[2] ){
        stop("gsvd: nrow(RW) does not equal ncol(RW) or ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(RW)){
        stop("gsvd: RW is empty (i.e., all 0s")
      }
    }

    if(RW_is_vector){
      if(length(RW)!=X_dimensions[2]){
        stop("gsvd: length(RW) does not equal ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(all(abs(RW)<=tol)){
        stop("gsvd: RW is empty (i.e., all 0s")
      }
    }
  }


  ## convenience checks *could* be removed* if problematic
  # convenience checks & conversions; these are meant to minimize LW's memory footprint
  if(!LW_is_missing){
    if( !LW_is_vector){

      # if( is_identity_matrix(LW) ){
      #   LW_is_missing <- T
      #   LW <- substitute() # neat! this makes it go missing
      # }

      if( is_diagonal_matrix(LW) ){
        LW <- diag(LW)
        LW_is_vector <- T  # now it's a vector
      }
    }

    if( LW_is_vector & all(LW==1) ){
      LW_is_missing <- T
      LW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize RW's memory footprint
  if(!RW_is_missing){
    if( !RW_is_vector ){

      # if( is_identity_matrix(RW) ){
      #   RW_is_missing <- T
      #   RW <- substitute() # neat! this makes it go missing
      # }

      if( !RW_is_vector & is_diagonal_matrix(RW) ){
        RW <- diag(RW)
        RW_is_vector <- T  # now it's a vector
      }
    }

    if( RW_is_vector & all(RW==1) ){
      RW_is_missing <- T
      RW <- substitute() # neat! this makes it go missing
    }
  }


  # this manipulates X as needed based on XLW
  if(!LW_is_missing){

    if( LW_is_vector ){
      sqrt_LW <- sqrt(LW)
      # X <- sweep(X,1,sqrt_LW,"*") ## replace the sweep with * & t()
      X <- X * sqrt_LW
    }else{
      LW <- as.matrix(LW)
      # X <- (LW %^% (1/2)) %*% X
      X <- sqrt_psd_matrix(LW) %*% X
    }

  }
  # this manipulates X as needed based on XRW
  if(!RW_is_missing){

    if( RW_is_vector ){
      sqrt_RW <- sqrt(RW)
      # X <- sweep(X,2, sqrt_RW,"*") ## replace the sweep with * & t()
      X <- t(t(X) * sqrt_RW)
    }else{
      RW <- as.matrix(RW)
      X <- X %*% sqrt_psd_matrix(RW)
    }

  }

  # all the decomposition things
  if(k<=0){
    k <- min(X_dimensions)
  }

  res <- tolerance_svd(X,nu=k,nv=k,tol=tol)
  res$d.orig <- res$d
  res$l.orig <- res$d.orig^2
  res$tau <- (res$l.orig/sum(res$l.orig)) * 100
  components.to.return <- min(length(res$d.orig),k) #a safety check
  res$d <- res$d.orig[1:components.to.return]
  res$l <- res$d^2
  res$u <- res$u[,1:components.to.return, drop = FALSE]
  res$v <- res$v[,1:components.to.return, drop = FALSE]



  # make scores according to weights
  if(!LW_is_missing){
    if(LW_is_vector){

      # res$p <- sweep(res$u,1,1/sqrt_LW,"*")
      res$p <- res$u / sqrt_LW
      # res$fi <- sweep(sweep(res$p,1,LW,"*"),2,res$d,"*")
      res$fi <- t(t(res$p * LW) * res$d)

    }else{

      # res$p <- (LW %^% (-1/2)) %*% res$u
      res$p <- invsqrt_psd_matrix(LW) %*% res$u
      # res$fi <- sweep((LW %*% res$p),2,res$d,"*")
      res$fi <- t(t(LW %*% res$p) * res$d)

    }
  }else{

    res$p <- res$u
    # res$fi <- sweep(res$p,2,res$d,"*")
    res$fi <- t(t(res$p) * res$d)

  }

  if(!RW_is_missing){
    if(RW_is_vector){

      # res$q <- sweep(res$v,1,1/sqrt_RW,"*")
      res$q <- res$v / sqrt_RW
      # res$fj <- sweep(sweep(res$q,1,RW,"*"),2,res$d,"*")
      res$fj <- t(t(res$q * RW) * res$d)

    }else{

      res$q <- invsqrt_psd_matrix(RW) %*% res$v
      # res$fj <- sweep((RW %*% res$q),2,res$d,"*")
      res$fj <- t(t(RW %*% res$q) * res$d)

    }
  }else{

    res$q <- res$v
    # res$fj <- sweep(res$q,2,res$d,"*")
    res$fj <- t(t(res$q)  * res$d)

  }


  rownames(res$fi) <- rownames(res$u) <- rownames(res$p) <- rownames(X)
  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(X)

  class(res) <- c("gsvd", "GSVD", "list")
  return(res)

}
