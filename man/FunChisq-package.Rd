\name{FunChisq-package}
\alias{FunChisq-package}
\docType{package}
\title{
Model-Free Functional Chi-Squared and Exact Tests
}
\description{
Statistical hypothesis testing methods for model-free functional dependency using asymptotic chi-squared or exact distributions. Functional chi-squared test statistics (Zhang and Song, 2013; Zhang, 2014; Nguyen, 2018; Zhong and Song, 2019) are asymmetric, functionally optimal, and model-free, unique from other related statistical measures.

Tests in this package reveal evidence for causality based on the causality-by-functionality principle (Simon and Rescher, 1966). The tests require data from two or more variables be formatted as a contingency table. Continuous variables need to be discretized first, for example, using the R package \pkg{Ckmeans.1d.dp}.

The package implements asymptotic functional chi-squared tests (Zhang and Song, 2013; Zhang, 2014), an exact functional test (Zhong and Song, 2019; Nguyen 2018), a comparative functional chi-squared test (Zhang, 2014), and also a comparative chi-squared test (Song et al., 2014; Zhang et al., 2015). The normalized functional chi-squared test was used by Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges (Hill et al., 2016).

A function index derived from the functional chi-squared offers a new effect size measure for the strength of function dependency. It is asymmetrically functionally optimal, different from the symmetric Cramer's V, also a better alternative to conditional entropy in many aspects.

A simulator is provided to generate functional, dependent non-functional, and independent patterns (Sharma et al., 2017).

For continuous data, these tests offer an advantage over regression analysis when a parametric form cannot be reliably assumed for the underlying function. For categorical data, they provide a novel means to assess directional dependency not possible with symmetrical Pearson's chi-squared test, G-test, or Fisher's exact test.
}
\details{
\tabular{ll}{
Package: \tab FunChisq\cr
Type: \tab Package\cr
Current version: \tab 2.4.9.1\cr
Initial release version: \tab 1.0\cr
Initial release date: \tab 2014-03-08\cr
License: \tab LGPL (>= 3) \cr
}

}
\seealso{
  For data discretization, an option is optimal univariate clustering. See package \pkg{Ckmeans.1d.dp}.

  For symmetric dependency tests on discrete data, see Pearson's chi-squared test (\code{\link[stats]{chisq.test}}), Fisher's exact test (\code{\link[stats]{fisher.test}}), mutual information (package \pkg{entropy}), and G-test, implemented in packages \pkg{DescTools} and \pkg{RVAideMemoire}.
}

\author{
Yang Zhang, Hua Zhong, Hien Nguyen, Ruby Sharma, Sajal Kumar and Joe Song
}

\references{
Hill, S. M., Heiser, L. M., Cokelaer, T., Unger, M., Nesser, N. K., Carlin, D. E., Zhang, Y., Sokolov, A., Paull, E. O., Wong, C. K., Graim, K., Bivol, A., Wang, H., Zhu, F., Afsari, B., Danilova, L. V., Favorov, A. V., Lee, W. S., Taylor, D., Hu, C. W., Long, B. L., Noren, D. P., Bisberg, A. J., HPN-DREAM Consortium, Mills, G. B., Gray, J. W., Kellen, M., Norman, T., Friend, S., Qutub, A. A., Fertig, E. J., Guan, Y., Song, M., Stuart, J. M., Spellman, P. T., Koeppl, H., Stolovitzky, G., Saez-Rodriguez, J. and Mukherjee, S. (2016) Inferring causal molecular networks: empirical assessment through a community-based effort. \emph{Nature Methods} \bold{13}(4), 310--318.

Nguyen, H. H. (2018) \emph{Inference of Functional Dependency via Asymmetric, Optimal, and Model-free Statistics.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State University, Las Cruces, USA.

Sharma, R., Kumar, S., Zhong, H. and Song, M. (2017) Simulating noisy, nonparametric, and multivariate discrete patterns. \emph{The R Journal} \bold{9}(2), 366--377. Retrieved from \url{https://doi.org/10.32614/RJ-2017-053}

Simon, H. A. and Rescher, N. (1966) Cause and counterfactual. \emph{Philosophy of Science} \bold{33}(4), 323--340.

Song M., Zhang, Y., Katzaroff, A. J., Edgar, B. A. and Buttitta, L. (2014) Hunting complex differential gene interaction patterns across molecular contexts. \emph{Nucleic Acids Research} \bold{42}(7), e57. Retrieved from \url{https://doi.org/10.1093/nar/gku086}

Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State University, Las Cruces, USA.

Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393--4407. Retrieved from \url{https://doi.org/10.1093/nar/gkv358}

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707. Retrieved from
\url{https://arxiv.org/abs/1311.2707}

Zhong, H. and Song, M. (2019) A fast exact functional test for directional association and cancer biology applications. \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics} \bold{16}(3), 818--826.  Retrieved from \url{https://doi.org/10.1109/TCBB.2018.2809743}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{package}
\keyword{htest}
\keyword{datagen}
\keyword{nonparametric}
