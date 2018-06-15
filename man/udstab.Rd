\name{udstab}
\alias{udstab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{udstab}
\description{computes the probability density function (pdf) of the univariate stable distribution based on formulas given by Nolan (1997) <doi.org/10.1080/15326349708807450> and asymptotic series, see Teimouri and Amindavar (2008).}
\usage{udstab(x, alpha, beta, sigma, mu, param)}
\arguments{
  \item{x}{point at which the pdf is computed}
  \item{alpha}{tail index parameter}
  \item{beta}{skewness parameter}
  \item{sigma}{scale parameter}
  \item{mu}{location parameter}
  \item{param}{kind of parameterization; must be 0 or 1 for \code{S_0} and \code{S_1} parameterizations, respectively}
}
%\details{}
\value{a numeric value}
\references{
Nolan, J. P. (1997). Numerical calculation of stable densities and distribution functions, Communications in statistics-Stochastic models, 13(4), 759-774.

Teimouri, M. and Amindavar, H. (2008). A novel approach to calculate stable densities, Proceedings of the World Congress on Engineering, 1, 710-714.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{\code{udstab()} computes the pdf of univariate stable distribution using asymptotic series within their convergence regions. For points outside of convergence regions, the pdf is computed using \code{stabledist} package based on formulas given by Nolan (1997). So, to compute the pdf using the \code{upstab()} we may need to install \code{stabledist} package.}
%\seealso{}
\examples{
# In the following, we compute the pdf of a univariate stable distribution at point 2
# with parameters alpha=1.2, beta=0.9, sigma=1, and mu=0 in S_{0} parameterization.
library("stabledist")
udstab(2,1.2,0.9,1,0,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
