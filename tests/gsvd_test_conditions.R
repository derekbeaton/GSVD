## correspondence analysis (CA)

#source('./R/gsvd.R')
#source('./R/tolerance.svd.R')
#source('./R/power.rebuild_matrix.R')
#source('./R/invert.rebuild_matrix.R')
#source('./R/isDiagonal.matrix.R')

  #authors data for CA
#load('./data/authors.rda')


## let's just do this directly through the GSVD and do all the preprocessing here.
sum.data <- sum(authors$ca$data)
wi <- rowSums(authors$ca$data)/sum.data
wj <- colSums(authors$ca$data)/sum.data

ca.res <- gsvd( (authors$ca$data/sum.data) - (wi %o% wj), 1/wi, 1/wj)

ca.res2 <- gsvd( (authors$ca$data/sum.data) - (wi %o% wj), diag(1/wi), diag(1/wj))

ca.res3 <- gsvd( (authors$ca$data/sum.data) - (wi %o% wj), 1/wi, diag(1/wj))

ca.res4 <- gsvd( (authors$ca$data/sum.data) - (wi %o% wj), diag(1/wi), 1/wj)

epca.res <- epCA(authors$ca$data,graphs=F)



### next I need CCA tests...

obj.wine <- scale(wine$objective,scale=T)
subj.wine <- scale(wine$subjective,scale=T)
base.cca.res <- cancor(obj.wine,subj.wine,xcenter = F,ycenter = F)


cca.res <- gsvd(
  matrix.generalized.inverse(crossprod(obj.wine)) %*% t(obj.wine) %*% subj.wine %*% matrix.generalized.inverse(crossprod(subj.wine)),
  crossprod(obj.wine),
  crossprod(subj.wine)
)

cca.res$lx <- (obj.wine %*% cca.res$p)
cca.res$ly <- (subj.wine %*% cca.res$q)





pca.res <- gsvd(subj.wine)
eppca.res <- epPCA(subj.wine,graphs=F,center=F,scale=F)

