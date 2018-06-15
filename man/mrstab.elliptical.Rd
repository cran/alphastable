\name{mrstab.elliptical}
\alias{mrstab.elliptical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mrstab.elliptical}
\description{generates \code{iid} realizations from \code{d}-dimensional elliptically contoured stable distribution, see Nolan (2013) <doi.org/10.1007/s00180-013-0396-7>.}
\usage{mrstab.elliptical(n, alpha, Sigma, Mu)}
\arguments{
  \item{n}{sample size}
  \item{alpha}{tail index parameter}
  \item{Sigma}{\code{d} by \code{d} positive definite dispersion matrix}
  \item{Mu}{location vector in \code{R^d}}
}
\details{\code{mrstab.elliptical()} needs to install the \code{mvtnorm} package}
\value{an \code{n} by \code{d} matrix of numeric values}
\references{
Nolan. J. P. (2013). Multivariate elliptically contoured stable distributions: theory and estimation, Computational Statistics, 28(5), 2067-2089.

Samorodnitsky, G. and Taqqu, M. S. (1994). Stable Non-Gaussian Random Processes: Stochastic Models and Infinite Variance, Chapman and Hall, London.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{\code{mrstab.elliptical()} generates \code{iid} realizations from \code{d}-dimensional elliptically contoured stable distribution based on definitions given by Nolan (2013) and Samorodnitsky and Taqqu (1994)}

%\seealso{}
\examples{
# In the following example, we simulate n=200 iid vectors of a two-dimensional elliptically
# contoured stable distribution with parameters alpha=1.3, Sigma=matrix(c(1,.5,.5,1),2,2),
# and mu=(0,0)^T.
library("mvtnorm")
library("stabledist")
mrstab.elliptical(200,1.3,matrix(c(1,.5,.5,1),ncol=2,nrow=2),c(0,0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
