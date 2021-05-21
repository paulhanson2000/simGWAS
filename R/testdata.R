#' Simulated data to use in testing and vignettes in the coloc package
#'
#' @title simGWAS_test_data
#'
#' @docType data
#'
#' @usage data(simGWAS_test_data)
#'
#' @format A matrix of haplotypes, one row per haplotype, one column per SNP
#'
#' @keywords datasets
#' @examples
#' data(simGWAS_test_data) # loads object called haps
#' str(haps)
"haps"

if(FALSE) {

  library(mvtnorm)
  library(bindata)
  simx <- function(nsamples,S,maf) {
    rmvbin(n=nsamples, margprob=maf, sigma=S)
    ## mu <- rep(0,nsnps)
    ## rawvars <- rmvnorm(n=nsamples, mean=mu, sigma=S)
    ## pvars <- pnorm(rawvars)
    ## x <- qbinom(1-pvars, 2, maf)
  }

  nsnps=100;nhaps=1000
  cat("Generate a small set of data\n")
  ntotal <- nsnps * nsamples
  ## S <- toeplitz(sample(1:10,nsnps,replace=TRUE)/10)
  S <- toeplitz((nsnps:1)/nsnps)
  R=rWishart(1,2*nsnps,S)[,,1]
  ## (1 - (abs(outer(1:nsnps,1:nsnps,`-`))/nsnps))
  maf=runif(nsnps,0.2,0.8)^2
  ## causals=c(1:nsnps)[order(abs(maf-0.5) * abs(((1:nsnps) - 12.5)/nsnps))][1:ncausals]
  set.seed(42)
  haps=simx(nsamples,S,maf)
  colnames(haps) <- paste0("s",1:nsnps)
  save(haps, file="data/simGWAS_test_data.rda", version=2)
}
