\name{cp.chisq.test}
\alias{cp.chisq.test}
\title{
Comparative Chi-Square Test for Association Heterogeneity
}
\description{
Comparative chi-square tests on two or more contingency tables. This test does not consider functional dependencies.
}

\usage{
cp.chisq.test(x, method="default")
}

\arguments{
  \item{x}{
  a list of at least two matrices representing contingency tables of the same dimensionality.
}

  \item{method}{
  a character string to specify the method to compute the chi-square statistic and its p-value. See Details.
}
}

\details{
The comparative chi-square test determines whether the patterns underlying the contingency tables are heterogeneous. This null test statistic is proved to asymptotically follow the chi-square distribution (Song et al., 2014), different from the widely used heterogeneity chi-square test (Zar, 2010).

Two methods are provided to compute the chi-square statistic and its p-value. When \code{method = "default"}, the p-value is computed using the chi-square distribution; when \code{method =} \code{"normalized"} a normalized chi-square is obtained by shifting and scaling the original chi-square and a p-value is computed using the standard normal distribution (Box et al., 2005). The normalized test is more conservative on the degrees of freedom.
}

\value{
A list with class "\code{htest}" containing the following components:

\item{statistic}{heterogeneity chi-square (by default), or normalized chi-square if \code{method = "normalized"}.}
\item{parameter}{degrees of freedom.}
\item{p.value}{p-value of the comparative chi-square test. By default, it is computed by the chi-square distribution. If \code{method = "normalized"}, it is the p-value of the normalized chi-square computed by the standard normal distribution.}
}

\references{
Box, G. E., Hunter, J. S., and Hunter, W. G. (2005). \emph{Statistics for Experimenters: Design, Innovation and Discovery}, 2nd Edition. Wiley-Interscience, New York.

Song M., Zhang Y., Katzaroff A. J., Edgar B. A., and Buttitta L. (2014). Hunting complex differential gene interaction patterns across molecular contexts. \emph{Nucleic Acids Research} \bold{42}(7), e57. Retrieved from \url{http://nar.oxfordjournals.org/content/42/7/e57.long}

Zar, J. H. (2010). \emph{Biostatistical Analysis}, 5th Edition. Prentice Hall, New Jersey.
}

\author{
Joe Song
}

\seealso{
For comparative functional chi-square test, \code{\link{cp.fun.chisq.test}}.
}

\examples{
x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
y <- t(x)
z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
data <- list(x,y,z)
cp.chisq.test(data)
cp.chisq.test(data, method="normalized")

}
