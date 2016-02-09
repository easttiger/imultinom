\name{imultinom}
\alias{imultinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  imultinom S4 object constructor
}
\description{
Constructs an imultinom S4 object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imultinom(listInput, matrixInput, a, b, De)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listInput}{
%%     ~~Describe \code{listInput} here~~
  This is a list in the shape of list(a = numeric vector, b = numeric vector, De = 0-1 matrix). If supplied, all subsequent inputs will be ignored.
}
  \item{matrixInput}{
%%     ~~Describe \code{matrixInput} here~~
  Input in adjacency matrix format. Precedence is after list input. If supplied, all subsequent inputs will be ignored.
}
  \item{a}{
%%     ~~Describe \code{a} here~~
  Individual input of the complete data exponents as a numeric vector. When using this input, make sure both listInput and matrixInput are not supplied.
}
  \item{b}{
%%     ~~Describe \code{b} here~~
  Individual input of the composite data exponents as a numeric vector. When using this input, make sure both listInput and matrixInput are not supplied.
}
  \item{De}{
%%     ~~Describe \code{De} here~~
  Individual input of the subset pattern indicator matrix as a 0-1 matrix (class matrix, mode numeric). Requires: nrow = length(a), ncol = length(b). When using this input, make sure both listInput and matrixInput are not supplied.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  Returns an S4 object having 3 slots: @a, @b, @De
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
  Fanghu Dong <easttiger@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
obj1_imultinom <- imultinom(a=c(3,3,3),b=c(2,9),De = matrix(c(1,0,1,1,1,0),ncol=2))
mle(obj1_imultinom)

cat ("Press [enter] to continue"); line <- readline()
obj2_imultinom <- imultinom(Weaver.Examples.Hankin2010Volleyball.ConvertInput())
mle(obj2_imultinom)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
