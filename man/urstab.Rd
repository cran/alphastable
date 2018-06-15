\name{urstab}
\alias{urstab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{urstab}
\description{simulates \code{iid} realizations from univariate stable distribution based on formulas given by Chambers et al. (1976) <DOI: 10.1080/01621459.1976.10480344> and Weron (1996) <doi.org/10.1016/0167-7152(95)00113-1>.}
\usage{urstab(n,alpha,beta,sigma,mu,param)}
\arguments{
  \item{n}{sample size}
  \item{alpha}{tail index parameter}
  \item{beta}{skewness parameter}
  \item{sigma}{scale parameter}
  \item{mu}{location parameter}
  \item{param}{kind of parameterization must; be 0 or 1 for \code{S_0} and \code{S_1} parameterizations, respectively}
}
%\details{}
\value{a vector of \code{n} numeric values}
\references{
Chambers, J. M., Mallows, C. L., and Stuck, B. W. (1976). A method for simulating stable random variables, Journal of the american statistical association, 71(354), 340-344.

Weron, R. (1996). On the Chambers-Mallows-Stuck method for simulating skewed stable random variables, Statistics & probability letters, 28(2), 165-171.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
%\note{}
%\seealso{}
\examples{
# By the following example, we simulate n=200 iid realizations from univariate stable
# distribution with parameters alpha=1.2, beta=0.5, sigma=2, and mu=0 in S_0 parameterization.
x <- urstab(200, 1.2, 0.5, 2, 0, 0)
}
%hist(x, freq=FALSE, breaks=20, col='green')}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{  }% use one of  RShowDoc("KEYWORDS")
%\keyword{  }% __ONLY ONE__ keyword per line
