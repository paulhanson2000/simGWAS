# simGWAS
Simulating GWAS output with a given causal model.

## Installation

This is an R package, and can be install directly from github using:

```{r}
install.packages("devtools") # if not already installed
devtools::install_github("chr1swallace/simGWAS")
```

## What is this about?

From the abstract:

> Methods for analysis of GWAS summary statistics have encouraged data sharing and democratised the analysis of different diseases. Ideal validation for such methods is application to simulated data, where some "truth" is known. As GWAS increase in size, so does the computational complexity of such evaluations; standard practice repeatedly simulates and analyses genotype data for all individuals in an example study. We have developed a novel method based on an alternative approach, directly simulating GWAS summary data, without individual data as an intermediate step. We mathematically derive the expected statistics for any set of causal variants and their effect sizes, conditional upon control haplotype frequencies (available from public reference datasets). Simulation of GWAS summary output can be conducted independently of sample size by simulating random variates about these expected values. Across a range of scenarios, our method, available as an open source R package, produces very similar output to that from simulating individual genotypes with a substantial gain in speed even for modest sample sizes. Fast simulation of GWAS summary statistics will enable more complete and rapid evaluation of summary statistic methods as well as opening new potential avenues of research in fine mapping and gene set enrichment analysis.

See more at the full manuscript https://www.biorxiv.org/content/early/2018/05/02/313023
