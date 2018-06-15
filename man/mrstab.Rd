\name{mrstab}
\alias{mrstab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mrstab}
\description{generates \code{iid} realizations from bivariate stable vectors using the methodology proposed by Modarres and Nolan (1994).}
\usage{mrstab(n, m, alpha, Gamma, Mu)}
\arguments{
  \item{n}{sample size}
  \item{m}{number of masses}
  \item{alpha}{tail index parameter}
  \item{Gamma}{vector of masses}
  \item{Mu}{location vector}
}
%\details{}
\value{a vector of \code{n} numeric values}
\references{
Modarres, R. and Nolan, J. P. (1994). A method for simulating stable random vectors, Computational Statistics, 9(1), 11-19.}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{\code{mrstab()} assumes that masses are located at unit sphere with addresses \code{s_j=(cos(2*pi(j-1)/m), sin(2*pi(j-1)/m))}; for \code{j=1,...,4}.
}
%\seealso{}
\examples{
# We use the following command to simulate n=200 iid vectors from a two-dimensional stable
# distribution with alpha=1.3, with a vector of 4 masses as gamma=(0.1,0.5,0.5,0.1)^T,
# and mu=(0,0)^T.
library("stabledist")
mrstab(200,4,1.3,c(0.1,0.5,0.5,0.1),c(0,0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
