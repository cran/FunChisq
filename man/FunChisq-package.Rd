\name{FunChisq-package}
\alias{FunChisq-package}
\docType{package}
\title{
Chi-Square and Exact Tests for Non-Parametric Functional Dependencies
}
\description{
Statistical hypothesis testing methods for non-parametric functional dependencies based on asymptotic chi-square and exact statistics. These tests were designed to reveal causality based on evidence from functional dependencies. It implements an asymptotic non-constant functional chi-square test (Zhang and Song, 2013; Zhang, 2014), an exact functional test (Zhong, 2014), and a comparative functional chi-square test (Zhang, 2014). It requires data from two or more variables be formatted as a contingency table. Continuous variables must be discretized first, for example, using the R package Ckmeans.1d.dp. The normalized functional chi-square test was used by the Best Performer NMSUSongLab in DREAM8 Breast Cancer Network Inference Challenges. These tests offer an advantage over regression analysis when a parametric form cannot be reliably assumed for the underlying function.
}
\details{
\tabular{ll}{
Package: \tab FunChisq\cr
Type: \tab Package\cr
Version: \tab 2.1.0\cr
Date: \tab 2015-06-29\cr
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
Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State Unviersity, Las Cruces, USA.

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707, 
\url{http://arxiv.org/abs/1311.2707}

Zhong, H. (2014) \emph{An Exact and Fast Statistical Test for Nonparametric Functional Dependencies.} Unpublished M.S. thesis, Department of Computer Science, New Mexico State University, Las Cruces, USA.
}

\keyword{ package }
