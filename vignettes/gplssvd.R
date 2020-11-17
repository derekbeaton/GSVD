## ----setup, include = FALSE---------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_data, echo = F------------------------------------------------------
library(GSVD)
data("synthetic_ONDRI", package="GSVD")

## ----two_tables---------------------------------------------------------------

X <- synthetic_ONDRI[,
                     c("TMT_A_sec", "TMT_B_sec",
                       "Stroop_color_sec","Stroop_word_sec",
                       "Stroop_inhibit_sec","Stroop_switch_sec")]
Y <- synthetic_ONDRI[,c("AGE","NAGM_PERCENT","NAWM_PERCENT")]

scaled_X <- scale(X, center = T, scale = T)
scaled_Y <- scale(Y, center = T, scale = T)


pls_gplssvd <- gplssvd(scaled_X, scaled_Y)

rrr_gplssvd <- gplssvd(scaled_X, scaled_Y, 
                       XRW = MASS::ginv(crossprod(scaled_X)))

cca_gplssvd <- gplssvd(scaled_X, scaled_Y, 
                       XRW = MASS::ginv(crossprod(scaled_X)), 
                       YRW = MASS::ginv(crossprod(scaled_Y)))


## ----cca_gplsvd_vis-----------------------------------------------------------

cca_gplssvd


## ----plsca_data, echo = T-----------------------------------------------------

pls_table_1 <- data.frame(
  MAPT = as.factor(synthetic_ONDRI$MAPT_DIPLOTYPE),
  APOE = as.factor(synthetic_ONDRI$APOE_GENOTYPE)
)
rownames(pls_table_1) <- synthetic_ONDRI$SUBJECT


pls_table_2 <- data.frame(
  SEX = as.factor(synthetic_ONDRI$SEX),
  NIHSS = as.factor(synthetic_ONDRI$NIHSS)
)
rownames(pls_table_2) <- synthetic_ONDRI$SUBJECT

## from: https://community.rstudio.com/t/how-to-get-a-full-set-of-dummy-variables/21682/2
disjunctive_data_X <- 
  model.matrix( ~. , 
    data=pls_table_1, 
    contrasts.arg = 
      lapply(
        data.frame(pls_table_1[,sapply(data.frame(pls_table_1), is.factor)]),
        contrasts, 
        contrasts = FALSE)
    )
disjunctive_data_X <- disjunctive_data_X[,-1]

disjunctive_data_Y <- 
  model.matrix( ~. , 
    data=pls_table_2, 
    contrasts.arg = 
      lapply(
        data.frame(pls_table_2[,sapply(data.frame(pls_table_2), is.factor)]),
        contrasts, 
        contrasts = FALSE)
    )
disjunctive_data_Y <- disjunctive_data_Y[,-1]


observed_matrix_X <- disjunctive_data_X / sum(disjunctive_data_X)
row_probabilities_X <- rowSums(observed_matrix_X)
col_probabilities_X <- colSums(observed_matrix_X)
expected_matrix_X <- row_probabilities_X %o% col_probabilities_X
deviations_matrix_X <- observed_matrix_X - expected_matrix_X

observed_matrix_Y <- disjunctive_data_Y / sum(disjunctive_data_Y)
row_probabilities_Y <- rowSums(observed_matrix_Y)
col_probabilities_Y <- colSums(observed_matrix_Y)
expected_matrix_Y <- row_probabilities_Y %o% col_probabilities_Y
deviations_matrix_Y <- observed_matrix_Y - expected_matrix_Y



plsca_gplssvd <- gplssvd( 
                  X = deviations_matrix_X,
                  Y = deviations_matrix_Y,
                  XLW = 1/row_probabilities_X,
                  YLW = 1/row_probabilities_Y,
                  XRW = 1/col_probabilities_X,
                  YRW = 1/col_probabilities_Y
                  )


