## some benchmarks on tolerance.svd() and t()

library(microbenchmark)
library(GSVD)
library(ggplot2)
times <- 100

n <- 100
p <- 100

X <- matrix(rnorm(n*p),n,p)
mbm.1 <- microbenchmark(
  just.svd = { tolerance.svd(X) },
  t.svd = {tolerance.svd(t(X))},
  times=times
)

X <- matrix(rnorm(n*10*p),n*10,p)
mbm.2 <- microbenchmark(
  just.svd = { tolerance.svd(X) },
  t.svd = {tolerance.svd(t(X))},
  times=times
)

X <- matrix(rnorm(n*100*p),n*100,p)
mbm.3 <- microbenchmark(
  just.svd = { tolerance.svd(X) },
  t.svd = {tolerance.svd(t(X))},
  times=times
)


# X <- matrix(rnorm(n*1000*p),n*1000,p)
# mbm.4 <- microbenchmark(
#   just.svd = { tolerance.svd(X) },
#   t.svd = {tolerance.svd(t(X))},
#   times=times
# )


X <- matrix(rnorm(n*10*p),n,p*10)
mbm.5 <- microbenchmark(
  just.svd = { tolerance.svd(X) },
  t.svd = {tolerance.svd(t(X))},
  times=times
)

X <- matrix(rnorm(n*100*p),n,p*100)
mbm.6 <- microbenchmark(
  just.svd = { tolerance.svd(X) },
  t.svd = {tolerance.svd(t(X))},
  times=times
)


# X <- matrix(rnorm(n*1000*p),n,p*1000)
# mbm.7 <- microbenchmark(
#   just.svd = { tolerance.svd(X) },
#   t.svd = {tolerance.svd(t(X))},
#   times=times
# )



