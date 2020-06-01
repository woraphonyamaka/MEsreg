\name{MEskink}
\alias{MEskink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized Maximum Entropy for estimating the smooth transition kink regression
}
\description{
GME inference method for the smooth transition kink regression model with under kink point. The advantage of GME method is that it is robust even when we have ill-posed or ill-conditioned problems, and thus, it has higher estimation accuracy and robustness, especially when the probability distribution of errors is unknown.
}
\usage{
MEskink(y,x,number,Z,V)
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
Entropy refers to the amount of uncertainty represented by a discrete probability distribution. The maximum entropy method was proposed by Jaynes (1957) and developed in the early 1990s by Golan, Judge, and Miller (1996) for estimating the unknown probabilities of a discrete probability distribution. This estimator uses the entropy-information measure of Shannon (1948) to recover those unknown probability distributions of underdetermined problems. This function is a simple estimation function for one covariate.
}
\value{

\item{beta }{intercept,beta_regime1,beta_regime2 }
\item{threshold }{kink point}
\item{smooth}{kink point}
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
set.seed(1)
n=100
thres=3
gam=1.2
e=rnorm(n)
x=rnorm(n,thres,5)
alpha=c(0.5,1,-1)

y=alpha[1]+(alpha[2]*(x*(1-logis(gam,x,thres))))+(alpha[3]*(x*(logis(gam,x,thres))))+e

MEskink(y,x,number="5",Z=10,V=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
