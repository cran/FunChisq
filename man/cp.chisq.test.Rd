\name{cp.chisq.test}
\alias{cp.chisq.test}
\title{
Comparative Chi-Squared Test for Association Heterogeneity
}
\description{
Comparative chi-squared tests on two or more contingency tables. This test does not consider functional dependencies.
}

\usage{
cp.chisq.test(
  x, method=c("chisq", "nchisq", "default", "normalized"),
  log.p = FALSE
)
}

\arguments{
  \item{x}{
  a list of at least two matrices representing contingency tables of the same dimensionality.
}

  \item{method}{
  a character string to specify the method to compute the chi-squared statistic and its p-value. The default is \code{"chisq"}. See Details.

    Note: \code{"default"} and \code{"normalized"} are deprecated.
}
  \item{log.p}{
  logical; if \code{TRUE}, the p-value is given as \code{log(p)}. Taking the log improves the accuracy when p-value is close to zero. The default is \code{FALSE}.
  }
}

\details{
The comparative chi-squared heterogeneity test determines whether the patterns underlying multiple contingency tables are heterogeneous. Its null test statistic is proved to asymptotically follow the chi-squared distribution (Song et al., 2014; Zhang et al., 2015), different from the widely used chi-squared heterogeneity test (Zar, 2010).

Two methods are provided to compute the chi-squared statistic and its p-value. When \code{method = "chisq"} (or \code{"default"}), the p-value is computed using the chi-squared distribution; when \code{method =} \code{"nchisq"} (or \code{"normalized"}) a normalized statistic is obtained by shifting and scaling the original chi-squared and a p-value is computed using the standard normal distribution (Box et al., 2005). The normalized test is more conservative on the degrees of freedom.
}

\value{
A list with class "\code{htest}" containing the following components:

\item{statistic}{heterogeneity statistic if \code{method = "chisq"} (equivalent to \code{"default"}), or normalized statistic if \code{method = "nchisq"} (equivalent to \code{"normalized"}).}
\item{parameter}{degrees of freedom of the chi-squared statistic.}
\item{p.value}{p-value of the comparative chi-squared test. By default, it is computed by the chi-squared distribution (\code{method = "chisq"} or \code{"default"}). If \code{method = "nchisq"} (or \code{"normalized"}), it is the p-value of the normalized chi-squared statistic using the standard normal distribution.}
}

\references{
Box, G. E., Hunter, J. S. and Hunter, W. G. (2005) \emph{Statistics for Experimenters: Design, Innovation and Discovery}, 2nd Ed., New York: Wiley-Interscience.

Song M., Zhang Y., Katzaroff A. J., Edgar B. A. and Buttitta L. (2014) Hunting complex differential gene interaction patterns across molecular contexts. \emph{Nucleic Acids Research} \bold{42}(7), e57. Retrieved from \url{https://nar.oxfordjournals.org/content/42/7/e57.long}

Zar, J. H. (2010) \emph{Biostatistical Analysis}, 5th Ed., New Jersey: Prentice Hall.

Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393--4407. Retrieved from \url{https://nar.oxfordjournals.org/content/43/9/4393.long}
}

\author{
Joe Song
}

\seealso{
For comparative \emph{functional} chi-squared test, \code{\link[FunChisq]{cp.fun.chisq.test}}.
}

\examples{
\dontrun{
x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
y <- t(x)
z <- matrix(c(1,0,1,4,0,4,0,4,0), 3)
data <- list(x,y,z)
cp.chisq.test(data)
cp.chisq.test(data, method="nchisq")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{htest}
\keyword{nonparametric}
