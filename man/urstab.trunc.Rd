\name{urstab.trunc}
\alias{urstab.trunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{urstab.trunc}
\description{using the methodology given by Soltan and Shirvani (2010), Shirvani and Soltani (2013) for simulating \code{iid} truncated stable random variable, it simulates truncated stable realizations.}
\usage{urstab.trunc(n, alpha, beta, sigma, mu, a, b, param)}
\arguments{
  \item{n}{sample size}
  \item{alpha}{tail index parameter}
  \item{beta}{skewness parameter}
  \item{sigma}{scale parameter}
  \item{mu}{location parameter}
  \item{a}{lower bound of truncation}
  \item{b}{upper bound of truncation}
  \item{param}{kind of parameterization; must be 0 or 1 for \code{S_0} and \code{S_1} parameterizations, respectively}
}
%\details{}
\value{a vector of \code{n} numeric values}
\references{
Shirvani, A. and Soltani, A. R. (2013). A characterization for truncated Cauchy random variables with nonzero skewness parameter, Computational Statistics, 28(3), 1011-1016.

Soltani, A. R. and Shirvani, A. (2010). Truncated stable random variables: characterization and simulation, Computational Statistics, 25(1), 155-161.

Teimouri, M. and Nadarajah, S. (2013). On simulating truncated stable random variables, Computational Statistics, 28(5), 2367-2377.

Teimouri, M. and Nadarajah, S. (2017). On simulating truncated skewed Cauchy random variables, Communications in Statistics-Simulation and Computation, 46(2), 1318-1321.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
%\note{}
%\seealso{}
\examples{
# We simulate n=200 iid realizations from truncated stable distribution with parameters
# alpha=1.3, beta=0.5, sigma=2, and mu=0 which is truncated over (-5,5) in S_0 parameterization.
urstab.trunc(200,1.3,0.5,2,0,-5,5,0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
