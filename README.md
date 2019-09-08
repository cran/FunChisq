README
======

## What?

The R package **'FunChisq'** provides statistical hypothesis testing methods for inferring model-free functional dependency. Functional test statistics are asymmetric and functionally optimal, unique from other related statistics. The test significance is based on either asymptotic chi-squared or exact distributions.

The tests include asymptotic *functional chi-squared tests* (Zhang & Song, 2013) [<arXiv:1311.2707>](https://arxiv.org/pdf/1311.2707v3.pdf) and an *exact functional test* (Zhong & Song, 2019) [<10.1109/TCBB.2018.2809743>](https://doi.org/10.1109/TCBB.2018.2809743). The *normalized* functional chi-squared test was used by Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges (Hill et al., 2016) [<10.1038/nmeth.3773>](https://doi.org/10.1038/nmeth.3773). 

A *function index* (Zhong & Song, in press) (Kumar et al., 2018) [<10.1109/BIBM.2018.8621502>](https://doi.org/10.1109/BIBM.2018.8621502) derived from the functional test statistic offers a new effect size measure of the strength of functional dependency

## Why?

Tests in this package reveal evidence for causality based on the causality-by-functionality principle. For continuous data, these tests offer an advantage over regression analysis when a parametric functional form cannot be assumed. Data can be first discretized, e.g., by the R package [**'Ckmeans.1d.dp'**](https://cran.r-project.org/package=Ckmeans.1d.dp). For categorical data, they provide a novel means to assess
directional dependency not possible with symmetrical Pearson's chi-squared or Fisher's exact tests. They are a better alternative to conditional entropy in many aspects.
