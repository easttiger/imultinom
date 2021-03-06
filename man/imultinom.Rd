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
options(iterrec=T)
cat("\n=================\n")
subject = "MLE of Ford (1957) adjacency matrix format"
cat(c(subject,"\n"))
Ford1957_adjmat = matrix(
	c( 0,15,15, 0,
	  11, 0,10,20,
	  11,10, 0,20,
	   0, 1, 1, 0),
	nrow = 4,
	byrow=TRUE);

Ford1957 <- imultinom(matrixInput = Ford1957_adjmat)
Ford1957_out <- mle(Ford1957)
if(!is.null(getOption("iterrec")) && getOption("iterrec")){
  with(Ford1957_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  plot(Ford1957_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(Ford1957_out)
}

cat ("Press [enter] to continue"); line <- readline()
cat("\n=================\n")
subject = "MLE of Hunter (2004) NASCAR2002 Plackett-Luce Modelling"
cat(c(subject,"\n"))
nascar2002PL <- imultinom(Weaver.Examples.Hunter2004PL.ConvertInput())
nascar2002PL_out = mle(nascar2002PL)
if(!is.null(getOption("iterrec")) && getOption("iterrec")){
  with(nascar2002PL_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  plot(nascar2002PL_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(nascar2002PL_out)
}

cat ("Press [enter] to continue"); line <- readline()
cat("\n=================\n")
subject = "MLE of Hunter (2004) NASCAR2002 Bradley-Terry Modelling"
cat(c(subject,"\n"))
data(nascar2002BT)
nascar2002BT_out = mle(nascar2002BT)
if(!is.null(getOption("iterrec")) && getOption("iterrec")){
  with(nascar2002BT_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  plot(nascar2002BT_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(nascar2002BT_out)
}

cat ("Press [enter] to continue"); line <- readline()
cat("\n=================\n")
subject = "MLE of Hankin (2010) Volleyball dataset";
cat(c(subject,"\n"))
Hankin2010 = imultinom(Weaver.Examples.Hankin2010Volleyball.ConvertInput())
Hankin2010_out = mle(Hankin2010)
if(!is.null(getOption("iterrec")) && getOption("iterrec")){
  with(Hankin2010_out, print(list(x=x,iterCount=iter$count,lnLik=iter$lnLik[iter$count])))
  plot(Hankin2010_out$iter$lnLik, type="b",pch=20,main=subject)
}else{
  print(Hankin2010_out)
}
options(iterrec=F)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
