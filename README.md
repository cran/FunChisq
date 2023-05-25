---
title: The 'FunChisq' R package
bibliography: inst/REFERENCES.bib
---
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/FunChisq)](https://cran.r-project.org/package=FunChisq)
[![CRAN_latest_release_date](https://www.r-pkg.org/badges/last-release/FunChisq)](https://cran.r-project.org/package=FunChisq)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/FunChisq)](https://cran.r-project.org/package=FunChisq)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/FunChisq)](https://cran.r-project.org/package=FunChisq)



### Overview

The package provides statistical hypothesis testing methods for inferring model-free functional dependency. Functional test statistics are asymmetric and functionally optimal, unique from other related statistics. The test significance is based on either asymptotic chi-squared or exact distributions.

The tests include an asymptotic *functional chi-squared test* [@zhang2013deciphering], *an adapted functional chi-squared test* [@Kumar2022AFT], and an *exact functional test* [@zhong2019eft;@Nguyen2020EFT]. The *normalized* functional chi-squared test was used by Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges (Hill et al., 2016) [<10.1038/nmeth.3773>](https://doi.org/10.1038/nmeth.3773). 

To measure the effect size, one can use the asymmetric *function index* [@Zhong2019FANTOM5;@KumarZSLS18]. Its value is minimized to 0 by perfectly independent patterns and maximized to 1 by perfect non-constant functions.

A simulator [@sharma2017simulating] can generate functional, non-functional, and independent patterns as contingency tables. The simulator provides options to control row and column marginal distributions and the noise level.

### When to use the package

Tests in this package can be used to reveal evidence for causality based on the causality-by-functionality principle. They target model-free inference without assuming a parametric model. For continuous data, these tests offer an advantage over regression analysis when a parametric functional form cannot be assumed. Data can be first discretized, e.g., by R packages ['Ckmeans.1d.dp'](https://cran.r-project.org/package=Ckmeans.1d.dp) or ['GridOnClusters'](https://cran.r-project.org/package=GridOnClusters). For categorical data, they provide a novel means to assess directional dependency not possible with symmetrical Pearson's chi-squared or Fisher's exact tests. They are a better alternative to conditional entropy in many aspects.

### To download and install the package

```{r}
install.packages("FunChisq")
```

### Citing the package

