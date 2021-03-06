% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tolerance_svd.R
\name{tolerance_svd}
\alias{tolerance_svd}
\title{\code{tolerance_svd}: An SVD to truncate potentially spurious (near machine precision) components.}
\usage{
tolerance_svd(x, nu = min(dim(x)), nv = min(dim(x)), tol = .Machine$double.eps)
}
\arguments{
\item{x}{A data matrix of size for input to the singular value decomposition (\code{\link{svd}})}

\item{nu}{The number of left singular vectors to be computed. Default is \code{min(dim(x))}}

\item{nv}{The number of right singular vectors to be computed. Default is \code{min(dim(x))}}

\item{tol}{Default is \code{.Machine$double.eps}. A tolerance level for eliminating near machine precision components.
Use of this parameter causes \code{tolerance_svd} to stop if negative or complex singular values are detected.
The use of \code{tol < 0}, \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf}, or \code{NULL} passes through to \code{\link{svd}}.}
}
\value{
A list with three elements (like \code{svd}):
 \item{d}{ A vector containing the singular values of x > \code{tol}.}
 \item{u}{ A matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension \code{min(c(nrow(x), nu, length(d))}.}
 \item{v}{ A matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension \code{min(c(ncol(x), nv, length(d))}.}
}
\description{
\code{tolerance_svd} eliminates likely spurious components: any eigenvalue (squared singular value) below a tolerance level is elminated.
   The (likely) spurious singular values and vectors are then eliminated from \code{$u}, \code{$d}, and \code{$v}.
   Additionally, all values in \code{abs($u)} or \code{abs($v)} that fall below the \code{tol} are set to 0.
   The use of a real positive value for \code{tol} will eliminate any small valued components.
   With \code{tol}, \code{tolerance_svd} will stop if any singular values are complex or negative.
}
\examples{
 data(wine)
 X <- scale(as.matrix(wine$objective))
 s_asis <- tolerance_svd(X)
 s_.Machine <- tolerance_svd(X, tol= .Machine$double.eps)
 s_000001 <- tolerance_svd(X, tol=.000001)

}
\seealso{
\code{\link{svd}}
}
\author{
Derek Beaton
}
\keyword{multivariate}
