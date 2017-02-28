\name{add.house.noise}
\alias{add.house.noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Apply House Noise to Discrete-Valued Tables
}
\description{
    The house noise model is designed for contingency tables of discrete values (Zhang et al., 2015). The model defines the condition probability for a given integer value to jump to all possible integer levels. The probability decreases as the new value deviates further away from the original value. Noise is applied indepdently for each data point in each table.
}
\usage{
add.house.noise(tables, u, margin=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tables}{
    A list of tables or one table. A table can be either a matrix or a data frame of integer values.
  }
  \item{u}{
    A numeric value between 0 and 1 to specify the noise level to be applied to the input tables. See Details.
  }
  \item{margin}{
    A value of either 0, 1, or 2. Default is 0.

    0: noise is applied along both rows and columns in a table. The sum of values in the table is the same before and after noise application.

    1: noise is applied along each row. The sum of each row is the same before and after noise application.

    2: noise is applied along each column. The sum of each column is the same before and after noise application;

  }
}
\details{
  At a minimum level of 0, no noise is applied on the input table(s). A maximum level of 1 indicates that the original sample will be changed to some other values with a probability of 1. For a discrete random variable of two possible values, a noise level of 1 will flip the values; a noise level of 0.5 is considered more appropriate if the interaction pattern is to be studied.
}

\value{
    If \code{tables} is a list, the function returns a list of tables with noised applied. If \code{tables} is a numeric matrix or a data frame, the function returns one table with noise applied.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
%% ~put references to the literature/web site here ~
Zhang, Y., Liu, Z. L. and Song, M. (2015) ChiNet uncovers rewired transcription subnetworks in tolerant yeast for advanced biofuels conversion. \emph{Nucleic Acids Research} \bold{43}(9), 4393-4407. Retrieved from \url{https://nar.oxfordjournals.org/content/43/9/4393.long}
}

\author{
    Hua Zhong, Yang Zhang and Joe Song.
}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{simulate_tables}}.
}
\examples{
  # Example 1. Add noise to a single table

  # Create a 4x4 table
  t <- matrix(c(3,0,0,0,
                0,2,2,0,
                0,0,0,4,
                3,3,2,0),
              nrow=4, ncol=4, byrow=TRUE)
  # Apply noise at level 0.1 along both rows and columns of the table.
  add.house.noise(t, 0.1, 0)

  # Example 2. Add noise to a list of tables

  # Create a list of tables
  t.list <- list(t+5, t*10, t*2)
  # Apply noise at level 0.2 along the rows of the table.
  add.house.noise(t.list, 0.2, 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{ datagen }
\concept{ house noise model }
