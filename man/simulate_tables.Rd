\name{simulate_tables}
\alias{simulate_tables}
\title{
Simulate Noisy, Nonparametric, and Discrete-Valued Contingency Tables
}

\description{
Generate random contingency tables representing various functional, non-functional, dependent, or independent patterns.
}

\usage{
simulate_tables(
  n=100, nrow=3, ncol=3,
  type = c("functional", "many.to.one",
           "discontinuous", "independent",
           "dependent.non.functional"),
  noise.model = c("house", "candle"), noise=0.0,
  n.tables=1,
  row.marginal=rep(1/nrow, nrow),
  col.marginal=rep(1/ncol, ncol)
)
}

\arguments{

   \item{n}{
   an integer specifying the sample size to be distributed in the table. For \code{"functional"},  \code{"many.to.one"}, and \code{"discontinuous"} tables, \code{n} must be no less than \code{nrow}. For \code{"independent"} and \code{"dependent.non.functional"} tables, \code{n} must be no less than \code{nrow}*\code{ncol}.
  }

  \item{nrow}{
  an integer specifying the number of rows in output tables. The value must be no less than 2. For \code{"many.to.one"} tables, \code{nrow} must be no less than 3.
  }

  \item{ncol}{
  an integer specifying the number of columns in output table. \code{ncol} must be no less than 2.
  }

  \item{type}{
  a character string to specify the type of pattern underlying the table. The options are \code{"functional"} (default),  \code{"many.to.one"}, \code{"discontinuous"}, \code{"independent"}, and \code{"dependent.non.functional"}. See Details.
%The options are \code{"functional"} (default, where y=f(x) but x may or may not be function of y),  \code{"many.to.one"} (where y=f(x) but x != f(y)), \code{"independent"} (x and y are statistically independent),  \code{"dependent.non.functional"} (x and y are statistically dependent but y!=f(x) and x!=f(y)).
  }
  \item{noise.model}{
  a character string indicating the noise model of either \code{"house"} for ordinal variables (Zhang et al., 2015) or \code{"candle"} for categorical variables. See \code{\link{add.noise}} for details.
  }

   \item{noise}{
  a numeric value between 0 and 1 specifying the noise level to be added to the table using function \code{\link{add.noise}}. The noise is applied along the rows of the table. See \code{\link{add.noise}} for details.
  }
  \item{n.tables}{
  an integer value specifying the number of tables to be generated.
  }
  \item{row.marginal}{
  a numeric vector of length \code{nrow} specifying row marginal probabilities. The default is a uniform distribution.
  %For \code{"functional table"}, \code{"dependent.non.functional"} and \code{"independent"} tables, \code{row.marginal} should be greater than or equal to 2.
  For \code{"many.to.one"} tables, the length of \code{row.marginal} vector must be no less than 3.
  }
  \item{col.marginal}{
  a numeric vector of length \code{ncol} specifying column marginal probabilities. It is only applicable in generating independent tables. The default is a uniform distribution.
  }

}

\details{
This function can generate five types of table representing different interaction patterns between row and column discrete random variables \eqn{X}{X} and \eqn{Y}{Y}. Three of the five types are non-constant functional patterns (\eqn{Y}{Y} is a non-constant function of \eqn{X}{X}):

\code{type="functional"}: \eqn{Y}{Y} is a function of \eqn{X}{X} but \eqn{X}{X} may or may not be a function of \eqn{Y}{Y}. The samples are distributed using the given row marginal probabilities.

\code{type="many.to.one"}: \eqn{Y}{Y} is a many-to-one function of \eqn{X}{X} but \eqn{X}{X} is not a function of \eqn{Y}{Y}. The samples are distributed on the basis of row probabilities.

\code{type="discontinuous"}: \eqn{Y}{Y} is a function of \eqn{X}{X}, where the function value of X must differ from its neighbors. \eqn{X}{X} may or may not be a function of \eqn{Y}{Y}. The samples are distributed using the given row marginal probabilities. A discontinuous function forms a contrast with those that are close to constant functions.

The fourth \code{type="dependent.non.functional"} is non-functional patterns where \eqn{X}{X} and \eqn{Y}{Y} are statistically dependent but not function of each other.

The fifth \code{type="independent"} is a pattern where \eqn{X}{X} and \eqn{Y}{Y} are statistically independent whose joint probability mass function is the product of their marginal probability mass functions.

Random noise can be optionally applied to the tables using either the house or the candle noise model. See \code{\link{add.noise}} for details.
%(Zhang et al., 2015) by specifying a non-zero value (\eqn{\le 1}{<= 1}) for argument \code{noise}.

Sharma et al. (2017) give full mathematical and statistical details of the simulation strategies for the above table types except the "discontinuous" type.
}

\value{
A list containing the following components:

\item{pattern.list}{a list of tables containing binary patterns in 0's and 1's. Each table is created by setting all non-zero entries in the corresponding sampled contingency table from \code{sample.list} to 1. Each table strictly satisfies the functional relationship for a given pattern \code{type} requested. This table does not meet the statistical requirements. As each table represents the truth regarding the mathematical relationship between the row and column variables, they can be used as the ground truth or gold standard for benchmarking.}

\item{sample.list}{a list of tables satisfying both the functional and statistical requirements. These tables are noise free.}

\item{noise.list}{a list of tables after applying noise to the corresponding tables in \code{sample.list}. Each table is the noisy version of the sampled contingency table. Due to the added noise, each table may no longer strictly satisfy the required functional or statistical relationships. These tables are the main output to be used for the evaluation of a discrete pattern discovery algorithm.}

\item{pvalue.list}{a list of p-values reporting the statistical significance of the generated tables for the required type. When the pattern type specifies a functional relationship, the p-values are computed by the  functional chi-square test (Zhang and Song, 2013); otherwise, the Pearson's chi-square test of independence is used to calculate the p-value.}
}

\references{
%% ~put references to the literature/web site here ~
Sharma, R., Kumar, S., Zhong, H. and Song, M. (2017) Simulating noisy, nonparametric, and multivariate discrete patterns. \emph{The R Journal} \bold{9}(2), 366--377. Retrieved from \url{https://journal.r-project.org/archive/2017/RJ-2017-053/index.html}

Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393--4407. Retrieved from \url{https://nar.oxfordjournals.org/content/43/9/4393.long}

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707,
\url{https://arxiv.org/abs/1311.2707}
}

\author{
Ruby Sharma, Sajal Kumar, Hua Zhong and Joe Song
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{add.noise}} for details of the noise model.
}

\examples{
\dontrun{
# In all examples, x is the row variable and y is the column
#    variable of a table.

# Example 1. Simulating a noise free function where y=f(x),
#            x may or may not be g(y)

simulate_tables(n=100, nrow=4, ncol=5, type="functional",
                noise=0.0, n.tables = 1,
                row.marginal = c(0.3,0.2,0.3,0.2))

# Example 2. Simulating a noisy functional pattern where
#            y=f(x), x may or may not be g(y)

simulate_tables(n=100, nrow=4, ncol=5, type="functional",
                noise=0.1, n.tables = 1,
                row.marginal = c(0.3,0.2,0.3,0.2))

# Example 3. Simulating a noise free many.to.one function where
#            y=f(x), x!=f(y).

simulate_tables(n=100, nrow=4, ncol=5, type="many.to.one",
                noise=0.0, n.tables = 1,
                row.marginal = c(0.4,0.3,0.1,0.2))

# Example 4. Simulating a pattern where x and y are
#            statistically independent.

simulate_tables(n=100, nrow=4, ncol=5, type="independent",
                noise=0.0, n.tables = 1,
                row.marginal = c(0.4,0.3,0.1,0.2),
                col.marginal = c(0.1,0.2,0.4,0.2,0.1))

# Example 5. Simulating noise-free dependent.non.functional
#            pattern where y!=f(x) and x and y are statistically
#            dependent.

simulate_tables(n=100, nrow=4, ncol=5,
                type="dependent.non.functional", noise=0.0,
                n.tables = 1, row.marginal = c(0.2,0.4,0.2,0.2))

# Example 6. Simulating noise-free discontinuous
#            pattern where  y=f(x), x may or may not be g(y)

simulate_tables(n=100, nrow=4, ncol=5,
                type="discontinuous", noise=0.0,
                n.tables = 1, row.marginal = c(0.2,0.4,0.2,0.2))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{ datagen }