\name{imultinom}
\alias{imultinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imultinom(listInput, matrixInput, a, b, De)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listInput}{
%%     ~~Describe \code{listInput} here~~
}
  \item{matrixInput}{
%%     ~~Describe \code{matrixInput} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{De}{
%%     ~~Describe \code{De} here~~
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
x = imultinom(a=c(3,3,3),b=c(2,9),De = matrix(c(1,0,1,1,1,0),ncol=2))
mle(x)
y = imultinom(Weaver.Examples.Hankin2010Volleyball.ConvertInput())
mle(y)


## The function is currently defined as
function (listInput, matrixInput, a, b, De) 
{
    if (!missing(listInput)) {
        a = listInput$a
        b = listInput$b
        De = listInput$De
    }
    else if (!missing(matrixInput)) {
        .matrix.input.validate(matrixInput)
        listInput = .input.matrix(matrixInput)
        a = listInput$a
        b = listInput$b
        De = t(listInput$tDe)
    }
    stopifnot("matrix" == class(De))
    stopifnot("numeric" == class(a))
    stopifnot("numeric" == class(b))
    stopifnot(length(a) > 1L)
    stopifnot(length(b) >= 1L)
    stopifnot((nrow(De) == length(a)) && (ncol(De) == length(b)))
    stopifnot(all(a >= 0))
    new("imultinom", a = a, b = b, De = De)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
