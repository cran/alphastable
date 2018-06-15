\name{mdstab.elliptical}
\alias{mdstab.elliptical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mdstab.elliptical}
\description{computes the probability density function of a \code{d}-dimensional elliptically contoured stable distribution at a given point in \code{R^{d}}, see Teimouri et al. (2018).}
\usage{mdstab.elliptical(x, alpha, Sigma, Mu)}
\arguments{
  \item{x}{ vector of real values in \code{R^d}}
  \item{alpha}{tail index parameter}
  \item{Sigma}{ \code{d} by \code{d} positive definite dispersion matrix}
  \item{Mu}{location vector in \code{R^d}}
}
%\details{}
\value{a numeric value}
\references{
Teimouri, M., Rezakhah, S., and Mohammadpour, A. (2018). Parameter estimation using the EM algorithm for symmetric stable random variables and sub-Gaussian random vectors, Journal of Statistical Theory and Applications, 17(3),1-20.}

\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{
\code{mdstab.elliptical()} computes the probability density function of an \code{d}-dimensional elliptically contoured stable distribution using either asymptotic series or Monte Carlo approximation.}
%\seealso{}
\examples{
# In the following example, we compute the pdf of a two-dimensional elliptically contoured
# stable distribution with parameters alpha=1.3, Sigma=matrix(c(1,.5,.5,1),2,2), and mu=(0,0)^T.
library("stabledist")
mdstab.elliptical(c(5,5),1.2,matrix(c(1,0.5,0.5,1),2,2),c(0,0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
