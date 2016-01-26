\name{FunChisq-package}
\alias{FunChisq-package}
\docType{package}
\title{
Chi-Square and Exact Tests for Non-Parametric Functional Dependencies
}
\description{
Statistical hypothesis testing methods for non-parametric functional dependencies based on asymptotic chi-square and exact statistics. These tests were designed to reveal evidence for causality based on functional dependencies. It implements asymptotic functional chi-square tests (Zhang and Song, 2013; Zhang, 2014), an exact functional test (Zhong, 2014), a comparative functional chi-square test (Zhang, 2014), and also a comparative chi-square test (Song et al., 2014; Zhang et al., 2015). It requires data from two or more variables be formatted as a contingency table. Continuous variables must be discretized first, for example, using the R package Ckmeans.1d.dp. The normalized functional chi-square test was used by the Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges. These tests offer an advantage over regression analysis when a parametric form cannot be reliably assumed for the underlying function.
}
\details{
\tabular{ll}{
Package: \tab FunChisq\cr
Type: \tab Package\cr
Version: \tab 2.2.1\cr
Date: \tab 2016-01-25\cr
Initial version: \tab 1.0\cr
Initial date: \tab 2014-03-08\cr
License: \tab LGPL (>= 3) \cr
}

}
\seealso{
  The \code{Ckmeans.1d.dp} package.
}

\author{
Yang Zhang, Hua Zhong and Joe Song.
}

\references{
Song M., Zhang Y., Katzaroff A. J., Edgar B. A., and Buttitta L. (2014). Hunting complex differential gene interaction patterns across molecular contexts. \emph{Nucleic Acids Research} \bold{42}(7), e57. Retrieved from \url{http://nar.oxfordjournals.org/content/42/7/e57.long}

Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State Unviersity, Las Cruces, USA.

Zhang, Y., Liu, Z. L., and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393-4407. Retrieved from \url{http://nar.oxfordjournals.org/content/43/9/4393.long}

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707,
\url{http://arxiv.org/abs/1311.2707}

Zhong, H. (2014) \emph{An Exact and Fast Statistical Test for Nonparametric Functional Dependencies.} Unpublished M.S. thesis, Department of Computer Science, New Mexico State University, Las Cruces, USA.
}

\keyword{ package }
