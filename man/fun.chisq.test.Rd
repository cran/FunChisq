\name{fun.chisq.test}
\alias{fun.chisq.test}
\title{
Chi-Square and Exact Tests for Model-Free Functional Dependency
}

\description{
Asymptotic chi-square, normalized chi-square or exact tests on contingency tables to determine model-free functional dependency of the column variable on the row variable.
}

\usage{
fun.chisq.test(
  x,
  method = c("fchisq", "nfchisq", "exact", "default",
             "normalized", "simulate.p.value"),
  alternative = c("non-constant", "all"), log.p=FALSE,
  index.kind = c("unconditional", "conditional"),
  simulate.nruns = 2000
)
}

\arguments{
  \item{x}{
  a matrix representing a contingency table. The row variable represents the independent variable or all unique combinations of multiple independent variables. The column variable is the dependent variable.
}

  \item{method}{
  a character string to specify the method to compute the functional chi-square test statistic and its p-value. The options are \code{"fchisq"} (equivalent to \code{"default"}, the default), \code{"nfchisq"} (equivalent to \code{"normalized"}), \code{"exact"} or \code{"simulate.p.value"}. See Details.

  Note: \code{"default"} and \code{"normalized"} are deprecated.
}
  \item{alternative}{
  a character string to specify the alternative hypothesis. The options are \code{"non-constant"} (default, non-constant functions) and \code{"all"} (all types of functions including constant ones).
  }
  \item{log.p}{
  logical; if \code{TRUE}, the p-value is given as \code{log(p)}. Taking the log improves the accuracy when p-value is close to zero. The default is \code{FALSE}.
  }
  \item{index.kind}{
  a character string to specify the kind of function index xi.f to be estimated. The options are \code{"unconditional"} (default) and \code{"conditional"}. See Details.
  }
  \item{simulate.nruns}{
   A number to specify the number of tables generated to simulate the null distribution. Default is \code{2000}. Only used when \code{method="simulate.p.value"}.
  }
}

\details{

The functional chi-square test determines whether the column variable is a function of the row variable in contingency table \code{x} (Zhang and Song, 2013; Zhang, 2014). This function supports three hypothesis testing methods:

\code{index.kind} specifies the kind of function index to be computed. If the experimental design controls neither the row nor column marginal sums, \code{index.kind = "unconditional"} (default) is recommended; If the column marginal sums are controlled, \code{index.kind = "conditional"} is recommended. The choice of \code{index.kind} affects only the function index xi.f value, but not the test statistic or p-value.

When \code{method="fchisq"} (equivalent to \code{"default"}, the default), the test statistic is computed as described in (Zhang and Song, 2013; Zhang, 2014) and the p-value is computed using the chi-square distribution.

When \code{method="nfchisq"} (equivalent to \code{"normalized"}), the test statistic is a normalized functional chi-square obtained by shifting and scaling the original chi-square (Zhang and Song, 2013; Zhang, 2014); and the p-value is computed using the standard normal distribution (Box et al., 2005). The normalized chi-square, more conservative on the degrees of freedom, was used by the Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges.

When \code{method="exact"}, an exact functional test (Zhong and Song, 2018) is performed. It computes an exact p-value and is fast when both the sample and table sizes are small. If the sample size is greater than 200 or the table size is larger than 5 by 5,  the exact test may not complete within a reasonable amount of time and the asymptotic functional chi-square test (\code{method="fchisq"}) is used instead. %When the sample size of an input contingency table is large, the exact functional test and the functional chi-square test will return similar p-values. When the expected values for all entries in a contingency table are greater than 5, the asymptotic tests will perform similarly with the exact test.
% On a 2-by-2 contingency tables, Fisher's exact test (Fisher, 1922) will be applied.
For 2-by-2 contingency tables, the asymptotic test options (\code{method="fchisq"} or \code{"nfchisq"}) are recommended to test functional dependency.

When \code{method="simulate.p.value"}, a simulated null distribution is used to calculate \code{p-value}. The null distribution is a multinomial distribution that is the product of two marginal distributions. Like other Monte Carlo based methods, this method is slower but may be more accurate than other methods based on asymptotic distributions.

%The following is an example suitable for the exact functional test. An epigenetic indicator called CpG island methylator phenotype (CIMP) is strongly associated with liver cancers. Specimens (Shen et al, 2002) to study CIMP were collected and divided into three groups with different CIMP status: negative (no methylation gene), intermediate (1~2 methylated genes) and positive (>2 methylated genes). The following table represents the frequencies of observed tumor protein p53 mutations and CIMP status:
%\tabular{lcccc}{
%  \tab \tab  CpG Island Methylator\cr
%  \tab \tab Phenotype (CIMP)\cr
%  \tab  Negative  \tab  Intermediate	\tab	Positive\cr
%%	\bold{Cirrhosis} \tab  \tab  \tab  \tab  0.038 (Pearson`s chisq)\cr
%%	  Negative  \tab	12	\tab	16	\tab	10\cr
%%	  Positive	\tab	5	\tab	18	\tab	21\cr
%%    \cr
%%	\bold{Hepatitis} \tab \tab  \tab  \tab 0.010 (Pearson`s chisq$)\cr
%%	  Negative	\tab	12	\tab	12	\tab	8\cr
%%	  Positive	\tab	 5	\tab	22	\tab	22\cr
%%    \cr
%%	\bold{Country risk} \tab  \tab  \tab \tab  0.021 (Fisher`s exact)\cr
%%	  Low risk	\tab	14	\tab	17	\tab	14\cr
%%	  High risk	\tab	 3	\tab	19	\tab	18\cr
%%    \cr
%	\bold{p53 mutation}\cr
%    No	\tab	12	\tab	26	\tab	18\cr
%	  Yes	\tab	 0	\tab	 8	\tab	12
%}
%Example 4 below performs the exact functional test on this table.

%\tabular{lccc}{
%\tab  Exact functional test\cr
%Cirrhosis->CIMP \tab  0.0706  \tab  0.0702  \tab  \bold{0.038*(P)}\cr
%CIMP->Cirrhosis \tab  \bold{0.0424*} \tab  \bold{0.0395*} \tab  \bold{0.038*(P)}\cr
%\cr
%Hepatitis->CIMP \tab  \bold{0.0301*} \tab  \bold{0.0311*} \tab  \bold{0.010*(P)}\cr
%CIMP->Hepatitis \tab  \bold{0.0103*} \tab  \bold{0.0123*} \tab  \bold{0.010*(P)}\cr
%\cr
%Country risk->CIMP  \tab  0.0706  \tab  0.0683  \tab  \bold{0.021*(F)}\cr
%CIMP->Country risk  \tab  \bold{0.0243*} \tab  \bold{0.0243*} \tab  \bold{0.021*(F)}\cr
%\cr
%\bold{Interactions}\tab  p-value\cr

%p53 -> CIMP \tab  0.0426\cr% \tab  0.0595  \tab  \bold{0.017*(F)}\cr
%CIMP -> p53 \tab  0.0273% \tab  0.0585  \tab  \bold{0.017*(F)}
%}

}

\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the functional chi-square statistic if \code{method = "fchisq"}, \code{"default"}, or \code{"exact"}; or the normalized functional chi-square statistic if \code{method = "nfchisq"} or \code{"normalized"}.}
\item{parameter}{degrees of freedom for the functional chi-square statistic.}
\item{p.value}{p-value of the functional test. If \code{method = "fchisq"} (or \code{"default"}), it is computed by an asymptotic chi-square distribution; if \code{method = "nfchisq"} (or \code{"normalized"}), it is computed by the standard normal distribution; if \code{method} \code{= "exact"}, it is computed by an exact hypergeometric distribution.}
\item{estimate}{an estimate of function index between 0 and 1. The value of 1 indicates a strictly mathematical function. It is asymmetrical with respect to transpose of the input contingency table, different from the symmetrical Cramer's V for Pearson's chi-square.}
}

\references{
Box, G. E., Hunter, J. S. and Hunter, W. G. (2005) \emph{Statistics for Experimenters: Design, Innovation and Discovery}, 2nd ed., New York: Wiley-Interscience.

Fisher, R. A. (1922) On the interpretation of chi-square from contingency tables, and the calculation of P. \emph{Journal of the Royal Statistical Society} \bold{85}(1), 87--94.

Pearson, K. (1900) On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. \emph{Philosophical Magazine Series 5} \bold{50}(302), 157--175.

%Shen, L., Ahuja, N., Shen, Y., Habib, N. A., Toyota, M., Rashid, A., and Issa, J.-P. J. (2002). DNA methylation and environmental exposures in human hepatocellular carcinoma. \emph{Journal of the National Cancer Institute}, 94(10), 755-761.

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707,
\url{https://arxiv.org/abs/1311.2707}

Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State University, Las Cruces, USA.

Zhong, H. and Song, M. (2018) A fast exact functional test for directional association and cancer biology applications. \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics.} In press.
}

\author{
Yang Zhang, Hua Zhong and Joe Song
}

\seealso{
  For data discretization by optimal univariate \var{k}-means clustering, see \pkg{Ckmeans.1d.dp}.

  For symmetrical dependency tests on discrete data, see Pearson's chi-square test \code{\link[stats]{chisq.test}}, Fisher's exact test \code{\link[stats]{fisher.test}}, and mutual information \pkg{entropy}.
}

\examples{
\dontrun{
# Example 1. Asymptotic functional chi-square test
x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
fun.chisq.test(x) # strong functional dependency
fun.chisq.test(t(x)) # weak functional dependency

# Example 2. Normalized functional chi-square test
x <- matrix(c(8,0,8,0,8,0,2,0,2), 3)
fun.chisq.test(x, method="nfchisq") # strong functional dependency
fun.chisq.test(t(x), method="nfchisq") # weak functional dependency

# Example 3. Exact functional chi-square test
x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
fun.chisq.test(x, method="exact") # strong functional dependency
fun.chisq.test(t(x), method="exact") # weak functional dependency

# Example 4. Exact functional chi-square test on a real data set
#            (Shen et al., 2002)
# x is a contingency table with row variable for p53 mutation and
#   column variable for CIMP
x <- matrix(c(12,26,18,0,8,12), nrow=2, ncol=3, byrow=TRUE)

# Test the functional dependency: p53 mutation -> CIMP
fun.chisq.test(x, method="exact")

# Test the functional dependency CIMP -> p53 mutation
fun.chisq.test(t(x), method="exact")

# Example 5. Asymptotic functional chi-square test with simulated distribution
x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
fun.chisq.test(x, method="simulate.p.value")
fun.chisq.test(x, method="simulate.p.value", simulate.n = 1000)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\keyword{htest}
\keyword{nonparametric}

