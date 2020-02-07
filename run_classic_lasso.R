# the code in this file can be used to run a LASSO fit on the processed mu transposition data, as in the accompanying paper
# the model is fitted as specified in STAR methods, and fitted coefficients will be written to the file lasso_cv_cust_symm.txt
# this code has been tested and runs successfully under both R 3.4.4 and 3.6.2; we have observed that minor numerical differences in the results
# may be observed using different R versions, but these differences do not appear to discernibly change the fits

library(glmnet)


# read and preprocess the data
dat.counts = read.csv("all_interaction_counts_with_norm.csv")
dat.counts = subset(dat.counts, is.finite(dat.counts$count_raw))
dat.counts$from_bin = as.factor(paste0( "Bin", dat.counts$Bin_ID) )
dat.counts$to_bin = as.factor(paste0( "Bin", dat.counts$Other_bin) )
distances.plus = abs( dat.counts$Other_bin - dat.counts$Bin_ID )
distances.minus = abs( (pmin( dat.counts$Other_bin, dat.counts$Bin_ID) + 100 ) - pmax( dat.counts$Other_bin, dat.counts$Bin_ID) )
dat.counts$dist = pmin(distances.plus, distances.minus)
dat.counts$logdist = log(1+dat.counts$dist)
dat.counts$binpair = as.factor(paste0("int", pmin(dat.counts$Bin_ID, dat.counts$Other_bin), "_", pmax(dat.counts$Bin_ID, dat.counts$Other_bin), sep=""))

dat.for.mat=dat.counts
dat.for.mat$X = NULL
dat.for.mat$Bin_ID = NULL
dat.for.mat$Other_bin = NULL
dat.for.mat$dist = NULL
dat.for.mat$norm_val = NULL
dat.for.mat$count_normed = NULL
dat.for.mat$count_raw = NULL
options(contrasts = c("contr.sum", "contr.sum"))

mm = model.matrix(~0+from_bin + to_bin + binpair+ logdist, dat.for.mat)

# change the model matrix to the type of contrasts that we want
mm.from = model.matrix(~0+from_bin, dat.for.mat)
mm_to = model.matrix(~0+to_bin, dat.for.mat)
mm_int = model.matrix(~0+binpair, dat.for.mat)
mm_dist = model.matrix(~0+logdist, dat.for.mat)
mm.alt = cbind( mm.from, mm_to, mm_int, mm_dist)


require(doMC)
registerDoMC(cores=10)
# run an initial fit
my.model.symm.cv = cv.glmnet( mm.alt, dat.counts$count_raw, family="poisson", offset=log(1.0/dat.counts$norm_val), standardize=F, intercept=T, parallel=T, nfolds=100)

# use a custom set of lambda values since the fit doesn't saturate using the defaults
cust.lambda = 0.8912509**seq(150)
my.mod.cv.symm.cust = cv.glmnet( mm.alt, dat.counts$count_raw, family="poisson", offset=log(1.0/dat.counts$norm_val), standardize=F, intercept=T, nfolds=100, lambda=cust.lambda, parallel=T)


# write the outputs
sink('lasso_cv_cust_symm.txt')
coef(my.mod.cv.symm.cust, s="lambda.1se")
sink()

sink('lasso_cv_cust_symm_min.txt')
coef(my.mod.cv.symm.cust, s="lambda.min")
sink()
