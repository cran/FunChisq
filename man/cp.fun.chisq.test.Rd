\name{cp.fun.chisq.test}
\alias{cp.fun.chisq.test}
\title{
Comparative Chi-Square Test for Non-Parametric Functional Heterogeneity
}
\description{
Comparative functional chi-square tests on two or more contingency tables.
}

\usage{
cp.fun.chisq.test(x, method = "fchisq", log.p = FALSE)
}

\arguments{
  \item{x}{
  a list of at least two matrices representing contingency tables of the same dimensionality.
}

  \item{method}{
  a character string to specify the method to compute the functional chi-square statistic and its p-value. The default is \code{"fchisq"} (equivalent to \code{"default"}). See Details.
}

  \item{log.p}{
  logical; if \code{TRUE}, the p-value is given as \code{log(p)}. Taking the log improves the accuracy when p-value is close to zero. The default is \code{FALSE}.
  }
}

\details{
The comparative functional chi-square test determines whether the patterns underlying the contingency tables are heterogeneous in a functional way. Specifically, it evaluates whether the column variable is a changed function of the row variable across the contingency tables.

Two methods are provided to compute the functional chi-square statistic and its p-value. When \code{method = "fchisq"} (or \code{"default"}), the p-value is computed using the chi-square distribution; when \code{method =} \code{"nfchisq"} (or \code{"normalized"}) a normalized functional chi-square is obtained by shifting and scaling the original chi-square and a p-value is computed using the standard normal distribution (Box et al., 2005). The normalized test is more conservative on the degrees of freedom.
}

\value{
A list with class "\code{htest}" containing the following components:

\item{statistic}{functional heterogeneity chi-square if \code{method = "fchisq"} (equivalent to \code{"default"}), or normalized functional chi-square if \code{method = "nfchisq"} (equivalent to \code{"normalized"}).}
\item{parameter}{degrees of freedom.}
\item{p.value}{p-value of the comparative functional chi-square test. By default, it is computed by the chi-square distribution. If \code{method = "normalized"}, it is the p-value of the normalized functional chi-square computed by the standard normal distribution.}
}

\references{
Box, G. E., Hunter, J. S., and Hunter, W. G. (2005). \emph{Statistics for Experimenters: Design, Innovation and Discovery}, 2nd Edition. Wiley-Interscience, New York.

Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State University, Las Cruces, USA.

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707. \url{http://arxiv.org/abs/1311.2707}
}

\author{
Yang Zhang and Joe Song
}

\seealso{
For comparative chi-square test that does not consider functional dependencies, \code{\link{cp.chisq.test}}.
}

\examples{
x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
y <- t(x)
z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
data <- list(x,y,z)
cp.fun.chisq.test(data)
cp.fun.chisq.test(data, method="nfchisq")
}
