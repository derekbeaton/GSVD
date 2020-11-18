## ----setup, include = FALSE---------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_data, echo = F------------------------------------------------------
library(GSVD)
data("synthetic_ONDRI", package="GSVD")

## ----pca_1--------------------------------------------------------------------

continuous_data <- synthetic_ONDRI[,c("TMT_A_sec", "TMT_B_sec",
                                      "Stroop_color_sec", "Stroop_word_sec", 
                                      "Stroop_inhibit_sec", "Stroop_switch_sec")]

cov_pca_geigen <- geigen( cov(continuous_data) )
cor_pca_geigen <- geigen( cor(continuous_data) )


## ----pca_1_print--------------------------------------------------------------

cov_pca_geigen

cor_pca_geigen


## ----mds_setup----------------------------------------------------------------

data_for_distances <- synthetic_ONDRI[,
                                      c("TMT_A_sec", "TMT_B_sec", 
                                        "Stroop_color_sec", "Stroop_word_sec",
                                        "Stroop_inhibit_sec","Stroop_switch_sec")]
scaled_data <- scale(data_for_distances, center = T, scale = T)

distance_matrix <- as.matrix( dist( scaled_data ) )

row_weights <- rep(1/nrow(distance_matrix), nrow(distance_matrix))

centering_matrix <- diag(nrow(distance_matrix)) - 
  ( rep(1,nrow(distance_matrix)) %o% row_weights )

matrix_to_decompose <- centering_matrix %*% 
  (-(distance_matrix^2)/2) %*% 
  t(centering_matrix)

mds_geigen <- geigen(matrix_to_decompose)


## ----mds_setup2, warning=FALSE------------------------------------------------

row_weights <- 1/synthetic_ONDRI$AGE

centering_matrix <- 
  diag(nrow(distance_matrix)) - (
  rep(1,nrow(distance_matrix)) %o% ( row_weights/sum(row_weights) )
)


matrix_to_decompose <- 
  -(centering_matrix %*% (distance_matrix^2) %*% t(centering_matrix))/2

mds_weighted_geigen <- geigen( matrix_to_decompose , W = row_weights, tol = NA )



## ----ics_prep-----------------------------------------------------------------


continuous_data <- synthetic_ONDRI[,c("TMT_A_sec", "TMT_B_sec",
                                      "Stroop_color_sec", "Stroop_word_sec", 
                                      "Stroop_inhibit_sec", "Stroop_switch_sec")]
scaled_data <- scale(continuous_data, center = T, scale = T)
cov_data <- cov(scaled_data)

### the following three lines help us compute the fourth moments covariance matrix.
sigma.data.sqrt <- GSVD::sqrt_psd_matrix(cov_data)
radius <- sqrt(rowSums((scaled_data %*% solve(sigma.data.sqrt))^2))
cov4_data <- crossprod(radius * scaled_data) / nrow(scaled_data)


## ----ics----------------------------------------------------------------------

ics_geigen <- GSVD::geigen(cov4_data, solve(cov_data))

ics_unmixing_matrix <- crossprod(ics_geigen$v, GSVD::invsqrt_psd_matrix(cov_data))
ics_generalzied_kurtosis <- ics_geigen$l_full / prod(ics_geigen$l_full)^(1/ncol(scaled_data))
ics_row_scores <- tcrossprod(scaled_data, ics_unmixing_matrix)


## ----pca_row_scores-----------------------------------------------------------

pca_row_scores <- tcrossprod(scaled_data, t(cor_pca_geigen$v))


## ----orthogonawesome, eval = F------------------------------------------------
#  
#  cor(pca_row_scores)
#  cor(ics_row_scores)
#  

## ----ics_vs_pca, echo = F-----------------------------------------------------

cor_pca_ics <- cor(pca_row_scores, ics_row_scores)
rownames(cor_pca_ics) <- paste0("PCA Component", 1:nrow(cor_pca_ics))
colnames(cor_pca_ics) <- paste0("ICS Component", 1:ncol(cor_pca_ics))

kable(cor_pca_ics, digits = 3)



