\name{plot_table}
\alias{plot_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a contingency table
}
\description{
The input contingency table is visualized by shades of a given color for better understanding of the underlying pattern.
}
\usage{
plot_table(table, xlab = "Column", ylab = "Row", col = "green3",
           xaxt = "n", yaxt = "n", main = NA, show.value = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{A data frame or a matrix.
}
  \item{xlab}{
  The lable of the horizontal axis.
}
  \item{ylab}{
  The lable of the vertical axis.
}
  \item{col}{
  The color corresponding to the maximum value in the table.
}
  \item{xaxt}{

  The style of the horizontal axis. See \code{par()}.
}
  \item{yaxt}{
  The style of the vertical axis. See \code{par()}.
}
  \item{main}{
  The title of the plot.
}
  \item{show.value}{
  logical. Show the value of each cell in the table on the plot.
}
  \item{\dots}{
  Parameters acceptable to \code{image()} function in the \pkg{graphics} package.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Joe Song
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
plot_table(matrix(1:6, nrow=2), col="mediumspringgreen")
plot_table(matrix(rnorm(25), nrow=5), col="orange", show.value=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }% use one of  RShowDoc("KEYWORDS")
