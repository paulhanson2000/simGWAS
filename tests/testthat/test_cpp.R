library(simGWAS)
library(testthat)
data(simGWAS_test_data) # loads a matrix of haplotypes
nsnps = ncol(haps)
snps <- colnames(haps) <- paste0("s",1:nsnps)
freq <- as.data.frame(haps+1)
freq$Probability <- 1/nrow(freq)
  CV=snps[1]
logOR=log(1.5)
nreps=100
set.seed(42)
  args=list(N0=10000, # number of controls
            N1=10000, # number of cases
            snps=snps, # column names in freq of SNPs for which Z scores should be generated
            W=CV, # causal variants, subset of snps
            gamma.W=logOR, # odds ratios
            freq=freq, # reference haplotypes
            nrep=nreps)
zsim <- do.call(simulated_z_score,args)
vbetasim=do.call(simulated_vbeta,args)
tmp=args; tmp$nrep=NULL
zexp <- do.call(expected_z_score,tmp)


test_that("mean.cpp", {
  x=rnorm(100)
  expect_equal(mean(x),simGWAS:::meanC(x))
})

test_that("wsumsq", {
  x=rnorm(100)
  y=rnorm(100)
  w=runif(100)
  expect_equal(sum( (x-y)^2 * w)/sum(w), simGWAS:::wsumsq(x,y,w))
})

test_that("vcf2haps", {
  x=matrix(c("0|0","0|1","1|0","1|1"),2,2)
  y=matrix(c(0,0,1,1,0,1,0,1),2,4)
  expect_equal(simGWAS:::vcf2haps(x),y)
})

test_that("psum", {
  x=rnorm(100)
  y=rnorm(100)
  expect_equal(sum(x*y), simGWAS:::psum(x,y))
})
