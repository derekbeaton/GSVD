The Generalized Singular Value Decomposition
================

GSVD
====

The generalized singular value decomposition (GSVD) is a [name shared by two different SVD techniques](https://en.wikipedia.org/wiki/Generalized_singular_value_decomposition). This package is for the "weighted" or "vector-constrained" GSVD. For one of the most straight forward introductions to the SVD and GSVD see the Appendices of Greenacre (1984).

Background
----------

The GSVD generalizes in two ways:

1.  it generalizes the standard SVD to allow for constraints or weights to be applied to the left and right singular vectors, and

2.  it generalizes many multivariate techniques, e.g., principal components analysis (PCA), correspondence analysis (CA), canonical correlation analysis (CCA), partial least squares (PLS), multidimensional scaling (MDS), linear discriminant analysis (LDA), and many many more.

The GSVD is an extraordinarily powerful and flexible tool for multivariate analyses and it is the core technique in the French school of data science/analyses (Holmes & Josse, 2017).

Package overview
----------------

This `GSVD` package is the first and most important package in the family of ExPosition2 packages. For an exposition of ExPosition see Beaton et al., (2014). The `GSVD` package is a focused package with one goal: give `R` users better and simpler access to the GSVD. `GSVD`'s companion packages will allow users more direct access to specific methods.

`gsvd()` is an efficient pure `R` implementation of the GSVD with no (current) dependencies. However, in the not-so-distant future, we plan to make `GSVD` better, faster, and more efficient with the use of `Rcpp` and `Matrix`.

Installation
------------

The `GSVD` package is not yet on CRAN but should be soon. For now, the simplest approach to installation is through the `devtools` package:

``` r
# devtools install via github
devtools::install_github("derekbeaton/GSVD")
```

Usage
-----

``` r
library(GSVD)

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
 
```

### Core

-   `gsvd()` is the generalized SVD function. It passes through to `tolerance.svd()` and makes use of `matrix.exponent()` with `%^%`.

-   `tolerance.svd()` is an alternative to the SVD function to only return vectors and values *above* some precision threshold (e.g., `.Machine$double.eps`)

### Bells-and-whistles

-   `matrix.exponent()` is matrix exponentiation through the SVD (i.e., raise the singular values to a power and then rebuild the matrix).

-   `%^%` is the same as `matrix.exponent()` but much more convenient and follows usual matrix operations in `R` (e.g., `%*%`)

-   `matrix.generalized.inverse()` is a specific implementation of `matrix.exponent()` that strictly performs the generalized inverse (i.e., singular values raised to -1).

-   `matrix.low.rank.rebuild()` is a heavy-duty function to rebuild lower rank versions of matrices with the SVD. It currently allows for a set of continuous components, arbitrary components, or to rebuild by percent of explained variance.

And beyond!
-----------

`GSVD` is the first package part of the larger ExPosition2 family. More to come soon!

References
----------

1.  Greenacre, M. (1984). Theory and applications of correspondence analysis. Academic Press.

2.  Holmes, S., & Josse, J. (2017). Discussion of “50 Years of Data Science”. Journal of Computational and Graphical Statistics, 26(4), 768-769.

3.  Beaton, D., Fatt, C. R. C., & Abdi, H. (2014). An ExPosition of multivariate analysis with the singular value decomposition in R. Computational Statistics & Data Analysis, 72, 176-189.
