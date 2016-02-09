\name{mle}
\alias{mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Maximum Likelihood Estimation of Incomplete Multinomial Model
}
\description{
  Find the Maximum Likelihood Estimate and its Asymptotic Covariance of an imultinom S4 object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  mle(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{listInput} here~~
  This is an imultinom S4 object. The function also reads the following global options for iteration related settings. \cr
  getOption("tolerance"), can be set by options(tolerance=1e-15), if you want a finer tolerance other than the default 1e-10 \cr
  getOption("iterrec"), can be set by options(iterrec = TRUE) if you want to turn on iteration recording. \cr
  getOption("maxit"), can be set by options(maxit=1000) if the default 500 maximum iteration steps is insufficient. \cr
  getOption("PriorThickness"), can be set by options(PriorThickness=1000). This parameter is only for weaver.bayes, and it requires some skills to set. \cr
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
  Returns $x as the MLE, $obs.info.mat as the asymptotic covariance matrix of the MLE, along with other iteration details.
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
