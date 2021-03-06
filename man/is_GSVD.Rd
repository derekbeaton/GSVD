% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_methods_helpers.R
\name{is_GSVD}
\alias{is_GSVD}
\title{is_GSVD}
\usage{
is_GSVD(x)
}
\arguments{
\item{x}{object to test}
}
\value{
boolean. \code{TRUE} if the object is of class GSVD, FALSE otherwise.
}
\description{
Tests if the \code{x} object is of class type "GSVD"
}
\details{
The three primary functions in the \code{GSVD} package produce an inherited (hierarchical) class structure where all of them are of type "GSVD". Those functions are \code{\link{geigen}}, \code{\link{gsvd}}, and \code{\link{gplssvd}}.
}
\seealso{
\code{\link{inherits}}
}
