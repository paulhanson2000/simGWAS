## load some data
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


test_that("errors as expected", {
  tmp=args; tmp$snps=c(snps,"notasnp") # snp not in freq
  expect_error(do.call(simulated_z_score, tmp))
  tmp=args; tmp$W=NULL # W missing
  expect_error(do.call(simulated_z_score, tmp))
  tmp=args; tmp$gamma.W=c(logOR,1) # gamma.W wrong length
  expect_error(do.call(simulated_z_score, tmp))
  tmp=args; tmp$freq=matrix() # freq not a data.frame
  expect_error(do.call(simulated_z_score, tmp))
})

test_that("correct types returned", {
  expect_type(zsim, "double")
  expect_length(zsim, nreps * ncol(haps))
  expect_type(vbetasim, "double")
  expect_length(vbetasim, nreps * ncol(haps))
})

test_that("simulated beta and z scores give expected values", {
  betasim <- zsim * sqrt(vbetasim)
  mu=mean(betasim[,1])
  expect_equal(abs(mu - logOR), 0.0060940571)
  expect_equal(mean(zsim[,1]), 16.103446)
  expect_equal(mean(vbetasim[,1]), 0.00061506872)
})

test_that("make_GenoProbList/fastextractsnps", {
  expect_equal(make_GenoProbList(snps[1:2], CV, freq[,c(snps[1:2],"Probability")]),
               list(list(c(`0` = 0.688900000000001, 0, 0), c(0, `0` = 0.2822,
0), c(0, 0, `1` = 0.0289)), list(c(0.688900000000001, 0.1826,
                                   0.0121), c(0, 0.0996000000000001, 0.0132), c(0, 0, 0.0036))))
})
