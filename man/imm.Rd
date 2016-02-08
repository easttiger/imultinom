\name{imm}
\alias{imm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imm(p, a, b, De, listInput, logscale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
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
  \item{listInput}{
%%     ~~Describe \code{listInput} here~~
}
  \item{logscale}{
%%     ~~Describe \code{logscale} here~~
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

## The function is currently defined as
function (p, a, b, De, listInput, logscale = TRUE) 
{
    {
        if (.Revo.tb("E:/Dropbox/GitHub/imultinom/R/weaver.R|390|0")) 
            browser()
    }
    if (!missing(listInput)) {
        {
            if (.Revo.tb("E:/Dropbox/GitHub/imultinom/R/weaver.R|391|0")) 
                browser()
        }
        a = listInput$a
        {
            if (.Revo.tb("E:/Dropbox/GitHub/imultinom/R/weaver.R|392|0")) 
                browser()
        }
        b = listInput$b
        {
            if (.Revo.tb("E:/Dropbox/GitHub/imultinom/R/weaver.R|393|0")) 
                browser()
        }
        De = listInput$De
    }
    {
        if (.Revo.tb("E:/Dropbox/GitHub/imultinom/R/weaver.R|395|0")) 
            browser()
    }
    p = p/sum(p)
    {
        if (.Revo.tb("E:/Dropbox/GitHub/imultinom/R/weaver.R|396|0")) 
            browser()
    }
    loglik = sum(a * log(p)) + sum(b * log(drop(t(p) \%*\% De)))
    {
        if (.Revo.tb("E:/Dropbox/GitHub/imultinom/R/weaver.R|397|0")) 
            browser()
    }
    ifelse(logscale, loglik, exp(loglik))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
