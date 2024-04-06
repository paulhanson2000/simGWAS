This fork of simGWAS makes it faster, and ARM-compatible. It can be used the same way as the original simGWAS, so see the original vignette. However, when using the `simulated_z_score()` function, do the following for maximum speed gains:

```{r}
# Pre-calculate GenoProbList
gpl <- make_GenoProbList(..., n_thread=parallel::detectCores())

# Pre-calculate weighted LD
LD <- simGWAS:::wcor2(as.matrix( freq[!colnames(freq) %in% "Probability"] ), freq$Probability)

# Trick to make the LD matrix positive definite, so Cholesky decomposition can be used
diag(LD) <- 1.0001

# Pass in pre-calculated GenoProbList
# Pass in pre-calculated LD
# Choose Cholesky decomposition instead of the default Eigenvalue decomposition
simulated_z_score(..., GenoProbList=gpl, LD=LD, rmvnorm_method="chol")
```

Compiling R with BLAS is highly recommended and will speed everything up by another 50 times. This is tricky but if you are using a compute cluster it has likely already been done for you.

# simGWAS
Simulating GWAS output with a given causal model.

## Installation

This is an R package, and can be install directly from github using:

```{r}
install.packages("devtools") # if not already installed
devtools::install_github("chr1swallace/simGWAS")
```

## Documentation

Please see the vignette in package, or at https://chr1swallace.github.io/simGWAS/articles/intro.html for example use

See [FAQ](FAQ.md) also.

## What is this about?

From the abstract:

> Methods for analysis of GWAS summary statistics have encouraged data sharing and democratised the analysis of different diseases. Ideal validation for such methods is application to simulated data, where some "truth" is known. As GWAS increase in size, so does the computational complexity of such evaluations; standard practice repeatedly simulates and analyses genotype data for all individuals in an example study. We have developed a novel method based on an alternative approach, directly simulating GWAS summary data, without individual data as an intermediate step. We mathematically derive the expected statistics for any set of causal variants and their effect sizes, conditional upon control haplotype frequencies (available from public reference datasets). Simulation of GWAS summary output can be conducted independently of sample size by simulating random variates about these expected values. Across a range of scenarios, our method, available as an open source R package, produces very similar output to that from simulating individual genotypes with a substantial gain in speed even for modest sample sizes. Fast simulation of GWAS summary statistics will enable more complete and rapid evaluation of summary statistic methods as well as opening new potential avenues of research in fine mapping and gene set enrichment analysis.

See more at the full manuscript (open access)\ 
Fortune and Wallace (2019) Bioinformatics https://doi.org/10.1093/bioinformatics/bty898

## Changelog

see NEWS.md

