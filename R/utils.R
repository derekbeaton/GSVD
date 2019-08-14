#' @export
#'
#' @title \code{is.diagonal.matrix}: test if a matrix is a diagonal matrix.
#'
#' @description \code{is.diagonal.matrix} takes a matrix and tests if it is a diagonal matrix.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix is a diagonal matrix, FALSE if the matrix is not.


## I stole this from somewhere... but I don't remember where
## I should find where I stole this from.
## I believe this: https://stackoverflow.com/questions/11057639/identifying-if-only-the-diagonal-is-non-zero
is.diagonal.matrix <- function(x,tol=.Machine$double.eps){
  if(is.null(dim(x))){
    stop("is.diagonal.matrix: X is not a matrix.")
  }
  x[ x^2 < tol ] <- 0
  return(all(x[lower.tri(x)] == 0, x[upper.tri(x)] == 0))
}


#' @export
#'
#' @title \code{is.empty.matrix}: test if a matrix contains all 0s.
#'
#' @description \code{is.empty.matrix} takes a matrix and tests if it contains all 0s.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix contains all 0s, FALSE if the matrix does not.


is.empty.matrix <- function(x,tol=.Machine$double.eps){

  x <- as.matrix(x)
  x[abs(x) < tol] <- 0

  if(sum(abs(x))==0){
    return(TRUE)
  }else{
    return(FALSE)
  }

}


#' @export
#'
#' @title \code{is.identical.matrix}: test if a matrix contains all identical values.
#'
#' @description \code{is.identical.matrix} takes a matrix and tests if it contains all identical values.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#' @param round.digits number of decimal places to round to (via \code{\link{round}} function).
#'
#' @return A boolean. TRUE if the matrix contains all identical values, FALSE if the matrix does not.

is.identical.matrix <- function(x,tol=.Machine$double.eps, round.digits = 12){

  x <- as.matrix(x)
  x[abs(x) < tol] <- 0
  x <- round(x,digits=round.digits)

  if(length(unique(c(x)))==1){
    return(TRUE)
  }else{
    return(FALSE)
  }

}


#' @export
#'
#' @title \code{is.identity.matrix}: test if a matrix is an identity matrix.
#'
#' @description \code{is.identity.matrix} takes a matrix and tests if it is an identity matrix.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix is an identity matrix, FALSE if the matrix is not.



is.identity.matrix <- function(x,tol=.Machine$double.eps){

  if(is.null(dim(x))){
    stop("is.identity.matrix: x is not a matrix.")
  }

  x <- as.matrix(x)
  x[abs(x) < tol] <- 0

  if(is.diagonal.matrix(x)){
    if( all(diag(x)==1) ){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }

}


#' Objective and subjective wine data
#'
#' Data used for illustrative purposes in Abdi H., Guillemot, V., Eslami, A., & Beaton, D. (in Press, 2018). Canonical correlation analysis (CCA). In R. Alhajj and J. Rokne (Eds.), Encyclopedia of Social Networks and Mining (2nd Edition). New York: Springer Verlag.
#'
#' @format A list called 'wine' with three matrices:
#' \describe{
#' \item{objective}{a 36 wine x 5 attribute matrix of objective variables for wines.}
#' \item{subjective}{a 36 wine x 9 attribute matrix of subjective variables for wines.}
#' \item{supplemental}{a 36 wine x 3 attribute matrix of supplemental information on the wines (country, color, grape).}
#' }
#'
"wine"

#' Authors
#'
#' Data used for illustrative purposes in Abdi, H. & Williams, L.J. (2010). Correspondence analysis. In N.J. Salkind, D.M., Dougherty, & B. Frey (Eds.): Encyclopedia of Research Design. Thousand Oaks (CA): Sage. pp. 267-278.
#'
#' @format A matrix that contains the frequency of punctuation usage of 6 authors x 3 punctuations.
"authors"


#' SNPs and drug Use
#'
#' Data used for illustrative purposes in Beaton, D., Filbey, F., & Abdi H. (2013). Integrating partial least squares correlation and correspondence analysis for nominal data. In Abdi, H., Chin, W., Esposito Vinzi, V., Russolillo, G., & Trinchera, L. (Eds.), New Perspectives in Partial Least Squares and Related Methods. New York: Springer Verlag. pp.81-94.
#'
#'@format A list called 'snps.druguse' with two matrices:
#' \describe{
#' \item{DATA1}{a 50 x 3 matrix of two-level categorical ("yes" or "no") variables.}
#' \item{DATA2}{a 50 x 2 matrix of three-level categorical (genotypes) variables. Note: Some entries contains NAs.}
#' }
"snps.druguse"


