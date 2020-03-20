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


#' @export
#'
#' @title \code{matrix_exponent}: raise matrix to a power and rebuild lower rank version
#'
#' @description \code{matrix_exponent} takes in a matrix and will compute raise that matrix to some arbitrary power via the singular value decomposition.
#'  Additionally, the matrix can be computed for a lower rank estimate of the matrix.
#'
#' @param x data matrix
#' @param power the power to raise \code{x} by (e.g., 2 is squared)
#' @param ... parameters to pass through to \code{\link{tolerance_svd}}
#'
#' @return The (possibly lower rank) raised to an arbitrary \code{power} version of \code{x}
#'
#' @seealso \code{\link{tolerance_svd}}
#'
#' @examples
#'  data(wine)
#'  X <- as.matrix(wine$objective)
#'  X.power_1 <- matrix.exponent(X)
#'  X / X.power_1
#'
#'  ## other examples.
#'  X.power_2 <- matrix.exponent(X,power=2)
#'  X.power_negative.1.div.2 <- matrix.exponent(X,power=-1/2)
#'
#'  X.power_negative.1 <- matrix.exponent(X,power=-1)
#'  X.power_negative.1 / (X %^% -1)
#'
#' @author Derek Beaton
#'
#' @keywords multivariate, diagonalization, eigen

matrix_exponent <- function(x, power = 1, ...){

  ##stolen from MASS::ginv()
  if (length(dim(x)) > 2 || !(is.numeric(x) || is.complex(x))){
    stop("matrix_exponent: 'x' must be a numeric or complex matrix")
  }
  if (!is.matrix(x)){
    x <- as.matrix(x)
  }

  ## should be tested for speed.

  #res <- tolerance_svd(x,...)
  #comp.ret <- 1:min(length(res$d),k)
  #return( (res$u[,comp.ret] * matrix(res$d[comp.ret]^power,nrow(res$u[,comp.ret]),ncol(res$u[,comp.ret]),byrow=T)) %*% t(res$v[,comp.ret]) )


  ## the special cases:
  ## power = 0
  if(power==0){
    x <- diag(1,nrow(x),ncol(x))
    attributes(x)$message.to.user = "https://www.youtube.com/watch?v=9w1y-kMPNcM"
    return( x )
  }
  ## is diagonal
  if(is.diagonal.matrix(x)){
    return( diag( diag(x)^power ) )

  }
  ## is vector
  if( any(dim(x)==1) ){
    return( x^power )
  }

  res <- tolerance_svd(x, ...)
  # need to replace this sweep
  return( sweep(res$u,2,res$d^power,"*") %*% t(res$v) )

}



#' @export
#'
#' @title Matrix exponentiation
#'
#' @description takes in a matrix and will compute raise that matrix to some arbitrary power via the singular value decomposition.
#'  Additionally, the matrix can be computed for a lower rank estimate of the matrix.
#'
#' @param x data matrix
#' @param power the power to raise \code{x} by (e.g., 2 is squared)
#'
#' @return \code{x} raised to an arbitrary \code{power}
#'
#' @seealso \code{\link{matrix_exponent}}
#'
#' @examples
#'  data(wine)
#'  X <- as.matrix(wine$objective)
#'  X %^% 2 # power of 2
#'  X %^% -1 # (generalized) inverse
#'
#' @author Derek Beaton
#'
#' @keywords multivariate, diagonalization, eigen
#'
`%^%` <- function(x,power){
  matrix_exponent(x,power=power)
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

#' Personal beer tasting notes
#'
#' Tasting notes, preferences, breweries and styles of 38 different craft beers from various breweries, across various styles.
#'
#'@format A list called 'beer.tasting.notes' with four matrices:
#' \describe{
#'  \item{data}{Data matrix. Tasting notes (ratings) of 38 different beers (rows) described by 16 different flavor profiles (columns).}
#'  \item{brewery.design}{Design matrix. Source brewery of 38 different beers (rows) across 26 breweries (columns).}
#'  \item{style.design}{Design matrix. Style of 38 different beers (rows) across 20 styles (columns) (styles as listed from Beer Advocate website).}
#'  \item{sup.data}{Supplementary data matrix. ABV and overall preference ratings of 38 beers described by two features (ABV & overall) in original value and rounded value.}
#' }
#'
#' @author Jenny Rieck and Derek Beaton
#'
#' @source
#' Jenny Rieck and Derek Beaton laboriously ``collected'' these data for ``experimental purposes''.
#'
#' @references
#' http://www.beeradvocate.com
#'
"beer.tasting.notes"


