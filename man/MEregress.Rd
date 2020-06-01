\name{MEregress}
\alias{MEregress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized Maximum Entropy for estimating linear regression.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MEregress(y,x,number,Z,V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{dependent variable}
  \item{x}{One dimension of dependent variable}
  \item{number}{number of supports}
  \item{Z}{bound of coefficient}
  \item{V}{bound of error}
}
\details{
This funciton is used to estimated the linear regression
}
\value{
\item{beta}{intercept,beta}
\item{Maxent}{Maximum entropy}
}
\references{
Golan, A., Judge, G. G., & Miller, D. (1996). Maximum entropy econometrics. Iowa State University, Department of Economics.
Jaynes, E. T. (1957). Information theory and statistical mechanics. Physical review, 106(4), 620.
Maneejuk, P. and Yamaka, W. (2020). Entropy Inference in Smooth Transition Kink Regression
}
\author{
Dr.Woraphon Yamaka
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("Rsolnp")
set.seed(1)
n=100
e=rnorm(n)
x0=rnorm(n)
x1=rnorm(n)
y=1+2*x0+3*x1+e
x=cbind(x0,x1)
MEregress(y,x,number="3",Z=10,V=5)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
