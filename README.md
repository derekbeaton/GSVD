The Generalized Singular Value Decomposition
================

# GSVD

The generalized singular value decomposition (GSVD) is a [name shared by
two different SVD
techniques](https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition).
This package is for the “weighted” or “vector-constrained” GSVD. For one
of the most straight forward introductions to the SVD and GSVD see the
Appendices of Greenacre (1984).

## Background

The GSVD generalizes in two ways:

1.  it generalizes the standard SVD to allow for constraints or weights
    to be applied to the left and right singular vectors, and

2.  it generalizes many multivariate techniques, e.g., principal
    components analysis (PCA), correspondence analysis (CA), canonical
    correlation analysis (CCA), partial least squares (PLS),
    multidimensional scaling (MDS), linear discriminant analysis (LDA),
    and many many more.

The GSVD is an extraordinarily powerful and flexible tool for
multivariate analyses and it is the core technique in the French school
of data science/analyses (Holmes & Josse, 2017).

## Package overview

This `GSVD` package is the first and most important package in the
family of ExPosition2 packages. For an exposition of ExPosition see
Beaton et al., (2014). The `GSVD` package is a focused package with one
goal: give `R` users better and simpler access to the GSVD. `GSVD`’s
companion packages will allow users more direct access to specific
methods.

`gsvd()` is an efficient pure `R` implementation of the GSVD with no
(current) dependencies. However, in the not-so-distant future, we plan
to make `GSVD` better, faster, and more efficient with the use of `Rcpp`
and `Matrix`.

## Installation

The `GSVD` package is not yet on CRAN but should be soon. For now, the
simplest approach to installation is through the `devtools` package:

``` r
# devtools install via github
devtools::install_github("derekbeaton/GSVD")
```

## Functions

### Core

  - `gsvd()` is the generalized SVD function. It passes through to
    `tolerance.svd()` and makes use of `matrix.exponent()` with `%^%`.

  - `tolerance.svd()` is an alternative to the SVD function to only
    return vectors and values *above* some precision threshold (e.g.,
    `.Machine$double.eps`)

  - `geigen()` is the generalized eigen function. It passes through to
    `tolerance.eigen()` and makes use of `matrix.exponent()` with `%^%`.

  - `tolerance.eigen()` is an alternative to the eigenvalue
    decomposition function to only return vectors and values *above*
    some precision threshold (e.g., `.Machine$double.eps`)

  - `gplssvd()` is the generalized partial least squares-singular value
    decomposition function. It passes through to `tolerance.svd()` and
    makes use of `matrix.exponent()` with `%^%`.

### Bells-and-whistles

  - `matrix.exponent()` is matrix exponentiation through the SVD (i.e.,
    raise the singular values to a power and then rebuild the matrix).

  - `%^%` is the same as `matrix.exponent()` but much more convenient
    and follows usual matrix operations in `R` (e.g., `%*%`)

  - `matrix.generalized.inverse()` is a specific implementation of
    `matrix.exponent()` that strictly performs the generalized inverse
    (i.e., singular values raised to -1).

  - `matrix.low.rank.rebuild()` is a heavy-duty function to rebuild
    lower rank versions of matrices with the SVD. It currently allows
    for a set of continuous components, arbitrary components, or to
    rebuild by percent of explained variance.

## Usage

### One table analyses

Here we provide three examples of “one table” analyses: principal
components analysis, multidimensional scaling (distances), and
correspondence analysis (with a smidgen of multidimensional scaling).
These make use of `gsvd()` and `geigen()`

``` r
library(GSVD)

# several examples of principal component analysis
 data(wine)
 wine.objective <- wine$objective
 ## "covariance" PCA
 cov.pca.data <- scale(wine.objective,scale=FALSE)
 cov.pca.res <- gsvd(cov.pca.data)
 ## "correlation" PCA
 cor.pca.data <- scale(wine.objective,scale=TRUE)
 cor.pca.res <- gsvd(cor.pca.data)
 ## an alternative approach to "correlation" PCA with GSVD constraints
 cor.pca.res2 <- gsvd(cov.pca.data,RW=1/apply(wine.objective,2,var))

# an example of multidimensional scaling
  squared.wine.attributes.distances <- as.matrix(dist(t(scale(wine.objective))))^2
  double.centered.wine <- (.Call(stats:::C_DoubleCentre,squared.wine.attributes.distances)/2)*-1
  mds.res <- geigen(double.centered.wine)

 
# an example of correspondence analysis.
 data(authors)
 Observed <- authors/sum(authors)
 row.w <- rowSums(Observed)
   row.W <- diag(1/row.w)
 col.w <- colSums(Observed)
   col.W <- diag(1/col.w)
 Expected <- row.w %o% col.w
 Deviations <- Observed - Expected
 ca.res <- gsvd(Deviations,row.W,col.W)
 
 
 # an alternate example of correspondence analysis by way of multidimensional scaling of Chi-squared distances
 Chi2DistanceMatrix <- t(Deviations) %*% diag(1/row.w) %*% Deviations
 ca.res_geigen <- geigen(Chi2DistanceMatrix, col.W)
 
```

### Two table analyses

Here we provide four examples of “two table” analyses all through
`gplssvd()`: partial least squares correlation, canonical correlation
analysis, reduced rank regression/redundancy analysis, and partial least
squares-correspondence analysis. This also requires data from the
`ExPosition` package. Each of these techniques can be expressed as
optimization of latent vectors.

``` r
library(GSVD)

  data(wine)
  X <- scale(wine$objective)
  Y <- scale(wine$subjective)
  
  ## an example of partial least squares (correlation)
  pls.res <- gplssvd(X, Y)
  
  pls.res$d
  diag( t(pls.res$lx) %*% pls.res$ly )
  
  ## an example of canonical correlation analysis
  cca.res <- gplssvd( X = X %^% (-1), Y = Y %^% (-1), XRW=crossprod(X), YRW=crossprod(Y))
    ## to note: cca.res$p and cca.res$q are the "canonical vectors" (see ?cancor and $xcoef and $ycoef)
  cca.res$d
  diag( t(cca.res$lx) %*% cca.res$ly )
  
  ## an alternate example of canonical correlation analysis
  cca.res_alt <- gplssvd(  X = X, Y = Y, XRW=crossprod(X %^% (-1)), YRW=crossprod(Y %^% (-1)))
    ## to note: cca.res_alt$fi %*% diag(1/cca.res_alt$d) and cca.res_alt$fj %*% diag(1/cca.res_alt$d) are the "canonical vectors" (see ?cancor and $xcoef and $ycoef)
  cca.res_alt$d
  diag( t(cca.res_alt$lx) %*% cca.res_alt$ly )
  
  
  ## an example of reduced rank regression/redundancy analysis
  rrr.res <- gplssvd(X, Y, XRW=crossprod(X %^% (-1))) 
    ## to note: rrr.res$fi is "beta" and rrr.res$v is alpha (see rrr.nonmiss: http://ftp.uni-bayreuth.de/math/statlib/S/rrr.s)
  rrr.res$d
  diag( t(rrr.res$lx) %*% rrr.res$ly )

  
  ## an example of pls-correspondence analysis (see https://utd.edu/~herve/abdi-bdAa2015_PLSCA.pdf)
  library(ExPosition)
  data("snps.druguse")
  X_nom <- makeNominalData(snps.druguse$DATA1)
    Ox <- X_nom / sum(X_nom)
    rx <- rowSums(Ox)
    cx <- colSums(Ox)
    Ex <- rx %o% cx
    Zx <- Ox - Ex
  Y_nom <- makeNominalData(snps.druguse$DATA2)
    Oy <- Y_nom / sum(Y_nom)
    ry <- rowSums(Oy)
    cy <- colSums(Oy)
    Ey <- ry %o% cy
    Zy <- Oy - Ey
  
  plsca.res <- gplssvd(Zx, Zy, 
                       XLW = 1/rx, YLW = 1/ry,
                       XRW = 1/cx, YRW = 1/cy)
  plsca.res$d
  diag(t(plsca.res$lx) %*% plsca.res$ly)
  
 
```

Et voila\! We have a unified generalized framework for many standard
multivariate analyses, all through the `g*()` family of functions here
in the `GSVD` package. (To note: the above two table analyses could also
have been done through the `gsvd()` but it is more convenient to do so
with `gplssvd()`.)

## And beyond\!

`GSVD` is the first package part of the larger ExPosition2 family. More
to come soon\!

## References

1.  Greenacre, M. (1984). Theory and applications of correspondence
    analysis. Academic Press.

2.  Holmes, S., & Josse, J. (2017). Discussion of “50 Years of Data
    Science”. Journal of Computational and Graphical Statistics, 26(4),
    768-769.

3.  Beaton, D., Fatt, C. R. C., & Abdi, H. (2014). An ExPosition of
    multivariate analysis with the singular value decomposition in R.
    Computational Statistics & Data Analysis, 72, 176-189.
