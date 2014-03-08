\name{fun.chisq.test}
\alias{fun.chisq.test}
\title{

Functional chi-square tests on contingency tables

}
\description{

fun.chisq.test performs non-constant functional chi-square tests on contingency tables.

}
\usage{
fun.chisq.test(x, method="default")
}

\arguments{
  \item{x}{
  a matrix representing a contingency table. The row variable is the independent variable and the column variable is the dependent variable
}

  \item{method}{
  method to compute the functional chi-square statistic and its p-value. See Details.
}
}
\details{

The functional chi-square test (Zhang and Song, 2013) is performed to test whether the column variable is a function of the row variable in contingency table x.

The methods to compute the functional chi-square value include a functional chi-square (\code{method="default"}) in which the p-value is computed using the chi-square distribution, and a normalized functional chi-square (\code{method="normalized"}) obtained by shifting and scaling the original chi-square and a p-value computed using the standard normal distribution (Box et al., 2005). The normalized test is more conservative on the degrees of freedom.
}
\value{

A list with class "\code{htest}" containing the following components:

\item{statistic}{functional chi-square (by default), or normalized functional chi-square if method="\code{normalized}"}
\item{parameter}{degrees of freedom for the functional chi-square statistic}
\item{p.value}{p-value of the functional chi-square test. By default, it is computed by the chi-square distribution. If \code{method="normalized"}, it is the p-value of the normalized functional chi-square computed by the standard normal distribution}
}
\references{

Zhang, Y. and Song, M. (2013) \emph{Deciphering interactions in causal networks without parametric assumptions.} 
\url{http://arxiv.org/abs/1311.2707}

Box, G. E., Hunter, J. S., and Hunter, W. G. (2005). \emph{Statistics for Experimenters: Design, Innovation and Discovery, 2nd Edition.} Wiley-Interscience, New York.
}
\author{
Yang Zhang and Joe Song
}

\examples{

x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
fun.chisq.test(x) #functional relationship gives significant p-value
fun.chisq.test(t(x)) #non-functional relationship gives relatively weak p-value

}
