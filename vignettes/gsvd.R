## ----setup, include = FALSE---------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_data, echo = F------------------------------------------------------
library(GSVD)
data("synthetic_ONDRI", package="GSVD")

## ----pca_2--------------------------------------------------------------------

continuous_data <- synthetic_ONDRI[,c("TMT_A_sec", "TMT_B_sec",
                                      "Stroop_color_sec", "Stroop_word_sec", 
                                      "Stroop_inhibit_sec", "Stroop_switch_sec")]

scaled_data <- scale(continuous_data, center = T, scale = T)
degrees_of_freedom <- nrow(continuous_data)-1
row_constraints <- rep(1/(degrees_of_freedom), nrow(scaled_data))

cor_pca_gsvd <- gsvd( scaled_data, LW = row_constraints )

cor_pca_gsvd


## ----ca_setup, echo = T-------------------------------------------------------

mapt_by_apoe_table <- table(synthetic_ONDRI$MAPT_DIPLOTYPE, synthetic_ONDRI$APOE_GENOTYPE)

rownames(mapt_by_apoe_table) <- paste0("MAPT_",rownames(mapt_by_apoe_table))
colnames(mapt_by_apoe_table) <- paste0("APOE_",colnames(mapt_by_apoe_table))


## ----ca-----------------------------------------------------------------------

observed_matrix <- mapt_by_apoe_table / sum(mapt_by_apoe_table)
row_frequencies <- rowSums(observed_matrix)
col_frequencies <- colSums(observed_matrix)
expected_matrix <- row_frequencies %o% col_frequencies
deviations_matrix <- observed_matrix - expected_matrix

ca_gsvd <- gsvd( deviations_matrix, 
                 LW = 1/row_frequencies, 
                 RW = 1/col_frequencies )


## ----mca_data_prep, echo = T--------------------------------------------------

mapt_apoe_sex <- data.frame(
  MAPT = as.factor(synthetic_ONDRI$MAPT_DIPLOTYPE),
  APOE = as.factor(synthetic_ONDRI$APOE_GENOTYPE),
  SEX = as.factor(synthetic_ONDRI$SEX),
  NIHSS = as.factor(synthetic_ONDRI$NIHSS)
)
rownames(mapt_apoe_sex) <- synthetic_ONDRI$SUBJECT

## from: https://community.rstudio.com/t/how-to-get-a-full-set-of-dummy-variables/21682/2
disjunctive_data <- 
  model.matrix( ~. , 
    data=mapt_apoe_sex, 
    contrasts.arg = 
      lapply(
        data.frame(mapt_apoe_sex[,sapply(data.frame(mapt_apoe_sex), is.factor)]),
        contrasts, 
        contrasts = FALSE)
    )


disjunctive_data <- disjunctive_data[,-1]

kable(disjunctive_data[c(1,10,20,40,80),c(1:5)], booktabs = T, caption = "Illustration of complete disjunctive coding (a.k.a. dummy coding, one-hot encoding) where each level of a categorical variable is represented. A value of '1' indicates the presence of that level for that row ('0' otherwise).")


## ----mca----------------------------------------------------------------------

observed_matrix <- disjunctive_data / sum(disjunctive_data)
row_frequencies <- rowSums(observed_matrix)
col_frequencies <- colSums(observed_matrix)
expected_matrix <- row_frequencies %o% col_frequencies
deviations_matrix <- observed_matrix - expected_matrix

mca_gsvd <- gsvd( deviations_matrix, 
                  LW = 1/row_frequencies,
                  RW = 1/col_frequencies)


## ----ics_prep-----------------------------------------------------------------


continuous_data <- synthetic_ONDRI[,c("TMT_A_sec", "TMT_B_sec",
                                      "Stroop_color_sec", "Stroop_word_sec", 
                                      "Stroop_inhibit_sec", "Stroop_switch_sec")]
scaled_data <- scale(continuous_data, center = T, scale = T)
cov_data <- cov(scaled_data)

sigma.data.sqrt <- GSVD::sqrt_psd_matrix(cov_data)
radius_squared <- rowSums((scaled_data %*% solve(sigma.data.sqrt))^2)

ics_gsvd <- GSVD::gsvd(scaled_data, LW = radius_squared, RW = solve(cov_data))

ics_gsvd_unmix <- crossprod(ics_gsvd$q, solve(cov_data))
ics_generalzied_kurtosis <- ics_gsvd$l_full / prod(ics_gsvd$l_full)^(1/ncol(scaled_data))
ics_row_scores <- tcrossprod(scaled_data, ics_gsvd_unmix)


## ----rmca---------------------------------------------------------------------

omegas <- c(0, 1, 2, 3, 4, 5, 10, 25, 50)
rmca_results <- vector("list",length(omegas))

centered_data <- scale(disjunctive_data,scale=F)
projection_matrix <- t(centered_data) %*% 
  MASS::ginv( tcrossprod(centered_data) ) %*% 
  centered_data
  
for(i in 1:length(omegas)){
  
  LW <- diag(nrow(centered_data)) + 
    (omegas[i] * MASS::ginv(tcrossprod(centered_data)))
  RW <- diag(colSums(disjunctive_data)) + (omegas[i] * projection_matrix)
  invRW <- t(MASS::ginv(RW))
  
  rownames(LW) <- colnames(LW) <- rownames(centered_data)
  rownames(invRW) <- rownames(RW)
  colnames(invRW) <- colnames(RW)
  
  rmca_results[[i]] <- gsvd(centered_data, LW = LW, RW = invRW, k = 2)
  
}


