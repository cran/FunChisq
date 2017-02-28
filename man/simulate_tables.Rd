\name{simulate_tables}
\alias{simulate_tables}
\title{
Simulate Noisy, Nonparametric, and Discrete-Valued Tables
}

\description{
Generate random contingency tables representing functional, non-functional, dependent, or independent interaction patterns.
}

\usage{
simulate_tables(n=100, nrow=3, ncol=3,
                type = c("functional",
                         "nonmonotonic",
                         "independent",
                         "dependent.non.functional"),
                noise=0.0, n.tables=1,
                row.marginal=rep(1/nrow, nrow),
                col.marginal= rep(1/ncol, ncol))

}

\arguments{

   \item{n}{
   an integer specifying the sample size to be distributed in the table. For \code{"functional"} and \code{"nonmonotonic"} tables, \code{n} must be no less than \code{nrow}. For \code{"independent"} tables, \code{n} must be no less than \code{nrow}*\code{ncol}. For \code{"dependent.non.functional"} tables, \code{n} must be greater than \code{nrow}.
  }

  \item{nrow}{
  an integer specifying the number of rows of output tables. The value must be no less than 2. For \code{"nonmonotonic"} tables, \code{nrow} must be no less than 3.
  }

  \item{ncol}{
  a integer value specifying the number of columns in desired table. \code{ncol} must be no less than 2.
  }

  \item{type}{
  a character string to specify the type of table. The options are \code{"functional"} (default),  \code{"nonmonotonic"}, \code{"independent"}, and \code{"dependent.non.functional"}. See Details.
%The options are \code{"functional"} (default, where y=f(x) but x may or may not be function of y),  \code{"nonmonotonic"} (where y=f(x) but x != f(y)), \code{"independent"} (x and y are statistically independent),  \code{"dependent.non.functional"} (x and y are statistically dependent but y!=f(x) and x!=f(y)).
  }

   \item{noise}{
  a numeric value between 0 and 1 specifying the factor of noise to be added to the table using the house noise model (Zhang et al., 2015). The house noise is applied along the rows of the table. See \code{\link{add.house.noise}} for details.
  }
  \item{n.tables}{
  a integer value specifying the number of tables to be generated.
  }
  \item{row.marginal}{
  a numeric vector of length at least 2 specifying row marginal probabilites.
  %For \code{"functional table"}, \code{"dependent.non.functional"} and \code{"independent"} tables, \code{row.marginal} should be greater than or equal to 2.
  For \code{"nonmonotonic"} tables, the length of \code{row.marginal} vector must be no less than 3.
  }
  \item{col.marginal}{
  a numeric vector specifying column marginal probabilites. It is only applicable in generating independent tables.
  }

}

\details{
This function can generate four types of table representing different interaction patterns between row variable \eqn{x}{x} and column varible \eqn{y}{y}:

\code{type="functional"}. \eqn{y}{y} is a function of \eqn{x}{x} (\eqn{y=f(x)}{y=f(x)}) but \eqn{x}{x} may or may not be a function of \eqn{y}{y}. The samples are distributed using the given row marginal probabilites.

\code{type="nonmonotonic"}.  This table type generates random nonmonotonic functional tables where \eqn{y=f(x)}{y=f(x)} but \eqn{x\neq f(y)}{x!=f(y)}. The samples are distributed on the basis of row probabilites.

\code{type="independent"}.  This table type generates independent tables where \eqn{y\neq f(x)}{y!=f(x)} and \eqn{x\neq f(y)}{x!=f(y)} and \eqn{x}{x} and \eqn{y}{y} are statistically independent.


\code{type="dependent.non.functional"}.  This table type generates dependent.non.functional tables  where \eqn{y\neq f(x)}{y!=f(x)} and \eqn{x\neq f(y)}{x!=f(y)} and \eqn{x}{x} and \eqn{y}{y} are statistically dependent to each other.

Random perturbations can be optionally applied to the tables using the house noise model (Zhang et al., 2015).

}

\value{
A list containing the following components:
\item{pattern.list}{a list of tables containing patterns in the form of 0's and 1's.}
\item{sample.list}{a list of tables after distributing the sample in pattern.table.}
\item{noise.list}{a list of tables after applying noise to the sample.table.}
\item{pvalue.list}{a list of P-values reporting the statitical significance of the generated tables for the required type. Depending on the table type, the P-values are from either Pearson's chi-square tests of independence or functional chi-square tests (Zhang and Song, 2013).}
}

\references{
%% ~put references to the literature/web site here ~
Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393-4407. Retrieved from \url{https://nar.oxfordjournals.org/content/43/9/4393.long}

Zhang, Y. and Song, M. (2013) Deciphering interactions in causal networks without parametric assumptions. \emph{arXiv Molecular Networks}, arXiv:1311.2707,
\url{https://arxiv.org/abs/1311.2707}
}

\author{
Ruby Sharma, Sajal Kumar, Hua Zhong and Joe Song
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{add.house.noise}} for details of the house noise model.
}

\examples{
# Example 1. Simulating noise free functional table where y=f(x), x may or may
#            not be f(y)

simulate_tables(n=100, nrow=4, ncol=5, type="functional", noise=0.0,
                n.tables = 1, row.marginal = c(0.3,0.2,0.3,0.2))

# Example 2. Simulating functional table with noise where y=f(x), x may or may
#            not be f(y)

simulate_tables(n=100, nrow=4, ncol=5, type="functional", noise=0.1,
                n.tables = 1, row.marginal = c(0.3,0.2,0.3,0.2))

# Example 3. Simulating noise free nonmonotonic table where y=f(x), x!=f(y).

simulate_tables(n=100, nrow=4, ncol=5, type="nonmonotonic", noise=0.0,
                n.tables = 1, row.marginal = c(0.4,0.3,0.1,0.2))

# Example 4. Simulating noise free Independent tables where y!=f(x), x!=f(y) but
#            x and y are statistically dependent.

simulate_tables(n=100, nrow=4, ncol=5, type="independent", noise=0.0,
                n.tables = 1, row.marginal = c(0.4,0.3,0.1,0.2),
                col.marginal = c(0.1,0.2,0.4,0.2,0.1))

# Example 5. Simulating noise free dependent.non.functional tables where y!=f(x)
#            and x and y are statistically dependent.

simulate_tables(n=100, nrow=4, ncol=5, type="dependent.non.functional", noise=0.0,
                n.tables = 1, row.marginal = c(0.2,0.4,0.2,0.2))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{ datagen }
