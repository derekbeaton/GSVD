% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsvd.R
\name{gsvd}
\alias{gsvd}
\title{Generalized singular value decomposition}
\usage{
gsvd(X, LW, RW, k = 0, tol = .Machine$double.eps)
}
\arguments{
\item{X}{a data matrix to decompose}

\item{LW}{\bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the matrix and thus left singular vectors.}

\item{RW}{\bold{R}ight \bold{W}eights -- the constraints applied to the right side (columns) of the matrix and thus right singular vectors.}

\item{k}{total number of components to return though the full variance will still be returned (see \code{d_full}). If 0, the full set of components are returned.}

\item{tol}{default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components.}
}
\value{
A list with eleven elements:
\item{d_full}{A vector containing the singular values of X above the tolerance threshold (based on eigenvalues).}
\item{l_full}{A vector containing the eigen values of X above the tolerance threshold (\code{tol}).}
\item{d}{A vector of length \code{min(length(d_full), k)} containing the retained singular values of X}
\item{l}{A vector of length \code{min(length(l_full), k)} containing the retained eigen values of X}
\item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(X)} by k.}
\item{p}{Left (rows) generalized singular vectors. Dimensions are \code{nrow(X)} by k.}
\item{fi}{Left (rows) component scores. Dimensions are \code{nrow(X)} by k.}
\item{v}{Right (columns) singular vectors. Dimensions are \code{ncol(X)} by k.}
\item{q}{Right (columns) generalized singular vectors. Dimensions are \code{ncol(X)} by k.}
\item{fj}{Right (columns) component scores. Dimensions are \code{ncol(X)} by k.}
}
\description{
\code{gsvd} takes in left (\code{LW}) and right (\code{RW}) constraints (usually diagonal matrices, but any positive semi-definite matrix is fine) that are applied to the data (\code{X}).
  Left and right constraints are used for the orthogonality conditions.
}
\examples{

 data(wine, package="GSVD")
 wine.objective <- wine$objective
 ## Principal components analysis: "covariance"
 ## "covariance" PCA
 cov.pca.data <- scale(wine.objective,scale=FALSE)
 cov.pca.res <- gsvd(cov.pca.data)

 ## Principal components analysis: "correlation"
 cor.pca.data <- scale(wine.objective,scale=TRUE)
 cor.pca.res <- gsvd(cor.pca.data)

 ## Principal components analysis: "correlation" with the covariance matrix and constraints
 cor.pca.res2 <- gsvd(cov.pca.data,RW=1/apply(wine.objective,2,var))


 ## Correspondence analysis
 data(authors, package="GSVD")
 Observed <- authors/sum(authors)
 row.w <- rowSums(Observed)
 col.w <- colSums(Observed)
 Expected <- row.w \%o\% col.w
 Deviations <- Observed - Expected
 ca.res <- gsvd(Deviations,1/row.w,1/col.w)


 ## Multiple correspondence analysis
 data("snps.druguse", package="GSVD")
 X <- model.matrix(~ .,
     data=snps.druguse$DATA1,
     contrasts.arg = lapply(snps.druguse$DATA1, contrasts, contrasts=FALSE))[,-1]
 Observed <- X/sum(X)
 row.w <- rowSums(Observed)
 col.w <- colSums(Observed)
 Expected <- row.w \%o\% col.w
 Deviations <- Observed - Expected
 ca.res <- gsvd(Deviations,1/row.w,1/col.w)

}
\seealso{
\code{\link{tolerance_svd}}, \code{\link{geigen}} and \code{\link{gplssvd}}
}
\author{
Derek Beaton
}
\keyword{multivariate}
