\name{Exact Functional Test}
\alias{Exact Functional Test}
\alias{EFTDP}
\alias{EFTDQP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Exact Functional Test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Perform the exact functional test on a contingency table
}
\usage{
EFTDP(nm)
EFTDQP(nm)
}
\arguments{
  \item{nm}{
%%     ~~Describe \code{x} here~~
A matrix of nonnegative integers representing a contingency table.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The functional exact test is performed using branch-and-bound with two algorithms (DP and DQP) to avoid re-calculation of bounds.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The exact P-value of the test.
}
\references{
%% ~put references to the literature/web site here ~
Nguyen, H. H. (2018) \emph{Inference of Functional Dependency via Asymmetric, Optimal, and Model-free Statistics.} Unpublished doctoral dissertation, Department of Computer Science, New Mexico State University, Las Cruces, USA.
}
\author{
%%  ~~who you are~~
Hien Nguyen, Hua Zhong, Joe Song
}
\note{
%%  ~~further notes~~
The functions provide a direct entry into the C++ implementations of the exact functional test.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{fun.chisq.test}}
}
\examples{
x = matrix(c(0, 6, 3, 0, 10, 5, 4, 4, 1), nrow=3)
EFTDQP(x)
EFTDQP(t(x))

EFTDP(x)
EFTDP(t(x))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of RShowDoc("KEYWORDS")
\keyword{ nonparametric }% __ONLY ONE__ keyword per line