# imports from other packages
#' @importFrom MASS ginv
NULL


#' @export
#'
#' @title \code{is_diagonal_matrix}: test if a matrix is a diagonal matrix.
#'
#' @description \code{is_diagonal_matrix} takes a matrix and tests if it is a diagonal matrix.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix is a diagonal matrix, FALSE if the matrix is not.


## enhanced version of this: https://stackoverflow.com/questions/11057639/identifying-if-only-the-diagonal-is-non-zero
is_diagonal_matrix <- function(x,tol=.Machine$double.eps){

  if( length(dim(x)) != 2 ){
    stop("is_diagonal_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("is_diagonal_matrix: x is not numeric.")
  }
  if( dim(x)[1] != dim(x)[2] ){
    stop("is_diagonal_matrix: x is not a square matrix.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }

  x[ x^2 < tol ] <- 0
  return(all(x[lower.tri(x)] == 0, x[upper.tri(x)] == 0))
}


#' @export
#'
#' @title \code{is_empty_matrix}: test if a matrix contains all 0s.
#'
#' @description \code{is_empty_matrix} takes a matrix and tests if it contains all 0s.
#'
#' @param x A matrix to test.
#' @param tol Tolerance precision to eliminate all abs(x) values below \code{tol}. Default is \code{.Machine$double.eps}.
#'
#' @return A boolean. TRUE if the matrix contains all 0s, FALSE if the matrix does not.


is_empty_matrix <- function(x,tol=.Machine$double.eps){

  if( length(dim(x)) != 2 ){
    stop("is_empty_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("is_empty_matrix: x is not numeric.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }

  x[abs(x) < tol] <- 0

  if(sum(abs(x))==0){
    return(TRUE)
  }else{
    return(FALSE)
  }

}

#' @export
#'
#' @title \code{are_all_values_positive}: check for strictly positive values in object
#'
#' @description \code{are_all_values_positive} test presumably a vector, but any numeric object, for strictly positive values
#'
#' @param x A numeric object
#'
#' @return A boolean. TRUE if all values are positive, FALSE otherwise

are_all_values_positive <- function(x){

  !(any(is.null(x)) |
  any(is.infinite(x)) |
  any(is.na(x)) |
  any(is.nan(x)) |
  any(x < 0))

}

#' @export
#'
#' @title \code{sqrt_psd_matrix}: square root of a positive semi-definite matrix
#'
#' @description \code{sqrt_psd_matrix} takes a square, positive semi-definite matrix and returns the square root of that matrix (via the eigenvalue decomposition by way of \code{tolerance_eigen}).
#'
#' @details Note that \code{sqrt_psd_matrix} uses \code{tolerance_eigen(...,tol=1e-13)} with the tolerance fixed at 1e-13. This enforces an "acceptably zero" value for eigenvalues, and also stops the computation if any eigenvalues are not zero or positive (within the tolerance parameter). If the computation is halted, that means the matrix was not positive semi-definite.
#'
#' @param x A square matrix (presumably positive semi-definite)
#'
#' @return A matrix. The square root of the \code{x} matrix
#' @seealso \code{\link{tolerance_eigen}}

sqrt_psd_matrix <- function(x){

  ## checks: just that they are a matrix & square & numeric
  if( length(dim(x)) != 2 ){
    stop("sqrt_psd_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("sqrt_psd_matrix: x is not numeric.")
  }
  if( dim(x)[1] != dim(x)[2] ){
    stop("sqrt_psd_matrix: x is not a square matrix.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }

  ## tolerance_eigen
  res <- tolerance_eigen(x, tol = 1e-13)

  ## rebuild
  return(t(t(res$vectors) * sqrt(res$values) ) %*% t(res$vectors))

}


#' @export
#'
#' @title \code{invsqrt_psd_matrix}: inverse square root of a positive semi-definite matrix
#'
#' @description \code{invsqrt_psd_matrix} takes a square, positive semi-definite matrix and returns the inverse square root of that matrix (via the eigenvalue decomposition by way of \code{tolerance_eigen}).
#'
#' @details Note that \code{invsqrt_psd_matrix} uses \code{tolerance_eigen(...,tol=1e-13)} with the tolerance fixed at 1e-13. This enforces an "acceptably zero" value for eigenvalues, and also stops the computation if any eigenvalues are not zero or positive (within the tolerance parameter). If the computation is halted, that means the matrix was not positive semi-definite.
#'
#' @param x A square matrix (presumably positive semi-definite)
#'
#' @return A matrix. The inverse square root of the \code{x} matrix
#' @seealso \code{\link{tolerance_eigen}}

invsqrt_psd_matrix <- function(x){

  ## checks: just that they are a matrix & square & numeric
  if( length(dim(x)) != 2 ){
    stop("invsqrt_psd_matrix: x is not a matrix.")
  }
  if( !is.numeric(x) ){
    stop("invsqrt_psd_matrix: x is not numeric.")
  }
  if( dim(x)[1] != dim(x)[2] ){
    stop("invsqrt_psd_matrix: x is not a square matrix.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }

  ## tolerance_eigen
  res <- tolerance_eigen(x, tol = 1e-13)

  ## rebuild
  return(t(t(res$vectors) * (1/sqrt(res$values)) ) %*% t(res$vectors))

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
#'@format A list called 'beer.tasting.notes' with 10 matrices:
#'\describe{
#'  \item{data}{Data matrix. Tasting notes (ratings) of 38 different beers (rows) described by 16 different flavor profiles (columns).}
#'  \item{sup.data}{Supplementary data matrix. ABV and overall preference ratings of 38 beers described by two features (ABV & overall) in original value and rounded value.}
#'  \item{brewery.design}{Design matrix. Source brewery of 38 different beers (rows) across 26 breweries (columns).}
#'  \item{style.design}{Design matrix. Style of 38 different beers (rows) across 20 styles (columns) (styles as listed from Beer Advocate website).}
#'  \item{superordinate.styles}{Supplementary data matrix. Style of 38 different beers (rows) across 8 "superordinate" styles.}
#'  \item{city.design}{Supplementary data matrix. Style of 38 different beers (rows) across 18 cities of origin.}
#'  \item{region.design}{Supplementary data matrix. Style of 38 different beers (rows) across 5 regions of origin.}
#'  \item{inspired.style}{Supplementary data matrix. Style of 38 different beers (rows) across 6 styles that inspired these beers.}
#'  \item{pale.sour.style}{Supplementary data matrix. Style of 38 different beers (rows) across 4 styles to capture if it is pale, sour, a saison, or "miscellaneous".}
#'  \item{physical.distances}{Data matrix. Square symmetric distance matrix of 29 different beers (rows) from \code{data} that reflect the physical distance---as the Nazgul flies---between each beer (by brewery).}
#' }
#' @author Jenny Rieck and Derek Beaton
#'
#' @source
#' Jenny Rieck and Derek Beaton laboriously ``collected'' these data for ``experimental purposes''.
#'
#' @references
#' http://www.beeradvocate.com
#'
"beer.tasting.notes"



#' synthetic_ADNI
#'
#' A synthetic data set derived from the Alzheimer's Disease Neuroimaging Initiative (ADNI). Data were generated from a specific set of variables in the \code{ADNIMERGE} \code{R} package. Synthetic data were produced with the \code{synthpop} package with the "cart" method.
#'
#' @details see http://adni.loni.usc.edu/ for more details on the ADNI data and the \code{ADNIMERGE} package
#'
#' @format A matrix that contains 623 observations (rows) and 17 variables (columns) of various data types (i.e., a mixture of continuous, categorical, and ordinal)
"synthetic_ADNI"



#' synthetic_ONDRI
#'
#' A synthetic data set derived from the Ontario Neurodegenerative Disease Research Initiative (ONDRI). Data were generated from a specific set of variables in the \code{ADNIMERGE} \code{R} package. Synthetic data were produced with the \code{synthpop} package with the "cart" method.
#'
#' @details see http://ondri.ca for more details on the study, and see https://github.com/ondri-nibs/toy_data for details on the synthetic data. To note, the \code{data.frame} available here is in a more "prepared" format than available as raw data files.
#'
#' @format A matrix that contains 138 observations (rows) and 17 variables (columns) of various data types (i.e., a mixture of continuous, categorical, and ordinal), including some identification variables (e.g., IDs)
"synthetic_ONDRI"
