\name{FunChisq-package}
\alias{FunChisq-package}
\docType{package}
\title{
Chi-Square Tests for Nonparametric Functional Dependencies
}
\description{
Statistical hypothesis testing methods for nonparametric functional dependencies. These tests were designed to reveal causality based on evidence from nonparametric functional dependencies based on the functional chi-square statistic. It implements asymptotic non-constant functional chi-square tests (Zhang and Song, 2013; Zhang, 2014), an exact functional test (Zhong, 2014), and a comparative functional chi-square test to compare contingency tables across multiple conditions (Zhang, 2014). Such tests offer an advantage over regression analysis when a parametric form cannot be reliably assumed for the underlying function. It requires data from two or more variables be formatted as a contingency table. Continuous variables must be discretized first using the Ckmeans.1d.dp function. The normalized functional chi-square test was used by the Best Performer NMSUSongLab in DREAM8 Breast Cancer Network Inference Challenges.
}
\details{
\tabular{ll}{
Package: \tab FunChisq\cr
Type: \tab Package\cr
Version: \tab 2.0.1\cr
Date: \tab 2015-02-23\cr
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
Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State Unviersity, Las Cruces, USA

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707, 
\url{http://arxiv.org/abs/1311.2707}

Zhong, H. (2014) \emph{An exact and fast statistical test for nonparametric functional dependencies.} Unpublished M.S. thesis, Department of Computer Science, New Mexico State University, Las Cruces, USA.

}

\keyword{ package }
