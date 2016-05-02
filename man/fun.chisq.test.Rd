\name{fun.chisq.test}
\alias{fun.chisq.test}
\title{
Chi-Square and Exact Tests for Non-Parametric Functional Dependencies
}

\description{
Asymptotic chi-square, normalized chi-square or exact tests on contingency tables to determine non-parametric functional dependencies of the column variable on the row variable.
}

\usage{
fun.chisq.test(x, method = "fchisq", alternative = "non-constant",
               log.p = FALSE, index.kind = "unconditional")
}

\arguments{
  \item{x}{
  a matrix representing a contingency table. The row variable represents the independent variable or all unique combinations of multiple independent variables. The column variable is the dependent variable.
}

  \item{method}{
  a character string to specify the method to compute the functional chi-square test statistic and its p-value. The options are \code{"fchisq"} (equivalent to \code{"default"}, the default), \code{"nfchisq"} (equivalent to \code{"normalized"}), or \code{"exact"}. See Details.
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
}

\details{

The functional chi-square test determines whether the column variable is a function of the row variable in contingency table \code{x} (Zhang and Song, 2013; Zhang, 2014). This function supports three hypothesis testing methods:

\code{index.kind} specifies the kind of estimated function index to be computed. If the experimental design controls neither the row nor column marginal sample size, \code{index.kind = "unconditional"} (default) is recommended; otherwise, \code{index.kind = "conditional"} is recommended. The choice of \code{index.kind} affects only the estimated function index xi.f value, but not the test statistic or p-value.

When \code{method="fchisq"} (equivalent to \code{"default"}, the default), the test statistic is computed as described in (Zhang and Song, 2013; Zhang, 2014) and the p-value is computed using the chi-square distribution.

When \code{method="nfchisq"} (equivalent to \code{"normalized"}), the test statistic is a normalized functional chi-square obtained by shifting and scaling the original chi-square (Zhang and Song, 2013; Zhang, 2014); and the p-value is computed using the standard normal distribution (Box et al., 2005). The normalized chi-square, more conservative on the degrees of freedom, was used by the Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges.

When \code{method="exact"}, an exact functional test (Zhong, 2014) is performed. It computes an exact p-value and is fast when both the sample and table sizes are small. If the sample size is greater than 200 or the table size is larger than 5 by 5,  the exact test may not complete within a reasonable amount of time and the asymptotic functional chi-square test (\code{method="fchisq"}) is used instead. %When the sample size of an input contingency table is large, the exact functional test and the functional chi-square test will return similar p-values. When the expected values for all entries in a contingency table are greater than 5, the asymptotic tests will perform similarly with the exact test.
% On a 2-by-2 contingency tables, Fisher's exact test (Fisher, 1922) will be applied.
For 2-by-2 contingency tables, the asymptotic test options (\code{method="fchisq"} or \code{"nfchisq"}) are recommended to test functional dependency.

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
Box, G. E., Hunter, J. S., and Hunter, W. G. (2005). \emph{Statistics for Experimenters: Design, Innovation and Discovery}, 2nd Edition. Wiley-Interscience, New York.

Fisher, R. A. (1922). On the interpretation of chi-square from contingency tables, and the calculation of P. \emph{Journal of the Royal Statistical Society}, 85(1), 87-94.

Pearson, K. (1990). On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. \emph{Philosophical Magazine Series 5}, 50(302), 157-175.

%Shen, L., Ahuja, N., Shen, Y., Habib, N. A., Toyota, M., Rashid, A., and Issa, J.-P. J. (2002). DNA methylation and environmental exposures in human hepatocellular carcinoma. \emph{Journal of the National Cancer Institute}, 94(10), 755-761.

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707,
\url{http://arxiv.org/abs/1311.2707}

Zhang, Y. (2014) \emph{Nonparametric Statistical Methods for Biological Network Inference.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State University, Las Cruces, USA.

Zhong, H. (2014) \emph{An Exact and Fast Statistical Test for Nonparametric Functional Dependencies.} Unpublished M.S. thesis, Department of Computer Science, New Mexico State University, Las Cruces, USA.

}


\author{
Yang Zhang, Hua Zhong and Joe Song
}

\examples{
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
}
