## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(knitr)
library(kableExtra)

## ---- echo = F----------------------------------------------------------------

test_tab <- rbind(c("","$\\bf{X}$","$\\bf{Y}$","${\\bf{M}_{\\bf X}}$","${\\bf{W}_{\\bf X}}$","${\\bf{M}_{\\bf Y}}$","${\\bf{W}_{\\bf Y}}$","$C$"),
c("`geigen()`","`X`","","","`W`","","","`k`"),
c("`gsvd()`","`X`",""," `LW`","`RW`","","","`k`"),
c("`gplssvd()`","`X`"," `Y`","`XRW`","`XLW`","`YRW`","`YLW`"," `k`"))

kable(test_tab, escape = FALSE, booktabs = TRUE, caption = "\\label{tab:arguments}Mapping between arguments (input) to functions (rows) and notation for the analysis tuples (columns).")


## ----echo=FALSE---------------------------------------------------------------
test_tab2 <- rbind(
  c("","What it is","Notation","`geigen`","`gsvd`","`gplssvd`"),
c("`d`","`k` singular values","${\\bf \\Delta}$","$\\checkmark$","$\\checkmark$", "$\\checkmark$"),
c("`d_full`","all singular values","${\\bf \\Delta}$","$\\checkmark$","$\\checkmark$","$\\checkmark$"),
c("`l`","`k` eigenvalues","${\\bf \\Lambda}$","$\\checkmark$", "$\\checkmark$","$\\checkmark$"),
c("`l_full`","all eigenvalues", "${\\bf \\Lambda}$", "$\\checkmark$", "$\\checkmark$", "$\\checkmark$"),
c("`u`", "`k` Left singular/eigen vectors", "${\\bf U}$", "",  "$\\checkmark$" , "$\\checkmark$"),
c("`v`","`k` Right singular/eigen vectors", "${\\bf V}$", "$\\checkmark$", "$\\checkmark$","$\\checkmark$"),
c("`p`", "`k` Left generalized singular/eigen vectors", "${\\bf P}$","","$\\checkmark$", "$\\checkmark$"),
c("`q`", "`k` Right generalized singular/eigen vectors", "${\\bf Q}$", "$\\checkmark$", "$\\checkmark$", "$\\checkmark$"),
c("`fi`","`k` Left component scores", "${\\bf F}_{I}$" ," ", "$\\checkmark$", "$\\checkmark$"),
c("`fj`", "`k` Right component scores", "${\\bf F}_{J}$", "$\\checkmark$", "$\\checkmark$", "$\\checkmark$"),
c("`lx`", "`k` Latent variable scores for `X`","${\\bf L}_{\\bf X}$", "","","$\\checkmark$"),
c("`ly`", "`k` Latent variable scores for `Y`", "${\\bf L}_{\\bf Y}$","","","$\\checkmark$"))

kable(test_tab2, escape = FALSE, caption = "\\label{tab:values}Mapping of values (output from functions; rows) to their conceptual meanings, notation used here, and which `GSVD` functions have these values.", booktabs = T)


