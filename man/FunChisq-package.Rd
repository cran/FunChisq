\name{FunChisq-package}
\alias{FunChisq-package}
\docType{package}
\title{
Chi-Square and Exact Tests for Non-Parametric Functional Dependencies
}
\description{
Statistical hypothesis testing methods for non-parametric functional dependencies using asymptotic chi-square or exact statistics. These tests were motivated to reveal evidence for causality based on functional dependencies (Simon and Rescher, 1966). The package implements asymptotic functional chi-square tests (Zhang and Song, 2013; Zhang, 2014), an exact functional test (Zhong, 2014), a comparative functional chi-square test (Zhang, 2014), and also a comparative chi-square test (Song et al., 2014; Zhang et al., 2015). The tests require data from two or more variables be formatted as a contingency table. Continuous variables need to be discretized first, for example, using the R package Ckmeans.1d.dp. The normalized functional chi-square test was used by Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges (Hill et al., 2016). For continuous data, these tests offer an advantage over regression analysis when a parametric form cannot be reliably assumed for the underlying function. For categorical data, they provide a novel means to assess directional dependencies not possible with symmetrical Pearson's chi-square or Fisher's exact tests.
}
\details{
\tabular{ll}{
Package: \tab FunChisq\cr
Type: \tab Package\cr
Version: \tab 2.3.3\cr
Initial version: \tab 1.0\cr
Initial date: \tab 2014-03-08\cr
License: \tab LGPL (>= 3) \cr
}

}
\seealso{
  The \code{Ckmeans.1d.dp} package.
}

\author{
Yang Zhang, Hua Zhong and Joe Song
}

\references{
Hill, S. M., Heiser, L. M., Cokelaer, T., Unger, M., Nesser, N. K., Carlin, D. E., Zhang, Y., Sokolov, A., Paull, E. O., Wong, C. K., Graim, K., Bivol, A., Wang, H., Zhu, F., Afsari, B., Danilova, L. V., Favorov, A. V., Lee, W. S., Taylor, D., Hu, C. W., Long, B. L., Noren, D. P., Bisberg, A. J., HPN-DREAM Consortium, Mills, G. B., Gray, J. W., Kellen, M., Norman, T., Friend, S., Qutub, A. A., Fertig, E. J., Guan, Y., Song, M., Stuart, J. M., Spellman, P. T., Koeppl, H., Stolovitzky, G., Saez-Rodriguez, J. and Mukherjee, S. (2016) Inferring causal molecular networks: empirical assessment through a community-based effort. \emph{Nature Methods} \bold{13}(4), 310-318.

Simon, H. A. and Rescher, N. (1966) Cause and counterfactual. \emph{Philosophy of Science} \bold{33}(4), 323-340.

Song M., Zhang, Y., Katzaroff, A. J., Edgar, B. A. and Buttitta, L. (2014). Hunting complex differential gene interaction patterns across molecular contexts. \emph{Nucleic Acids Research} \bold{42}(7), e57. Retrieved from \url{http://nar.oxfordjournals.org/content/42/7/e57.long}

Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State University, Las Cruces, USA.

Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393-4407. Retrieved from \url{http://nar.oxfordjournals.org/content/43/9/4393.long}

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707,
\url{http://arxiv.org/abs/1311.2707}

Zhong, H. (2014) \emph{An Exact and Fast Statistical Test for Nonparametric Functional Dependencies.} Unpublished M.S. thesis, Department of Computer Science, New Mexico State University, Las Cruces, USA.
}

\keyword{ package }
