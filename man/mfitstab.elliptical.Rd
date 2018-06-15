\name{mfitstab.elliptical}
\alias{mfitstab.elliptical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mfitstab.elliptical}
\description{estimates the parameters of a \code{d}-dimensional elliptically contoured stable distribution, see Teimouri et al. (2018).}
\usage{mfitstab.elliptical(yy, alpha0, Sigma0, Mu0)}
\arguments{
  \item{yy}{vector of d-dimensional realizations}
  \item{alpha0}{initial value of the tail index parameter to start the EM algorithm}
  \item{Sigma0}{initial value of the dispersion matrix to start the EM algorithm}
  \item{Mu0}{initial value of the location vector to start the EM algorithm}
}
%\details{}
\value{
 \item{alpha}{estimated value of the tail index parameter}
 \item{Sigma}{estimated value of the dispersion matrix}
 \item{Mu}{estimated value of the location vector}
}
\references{
Teimouri, M., Rezakhah, S., and Mohammadpour, A. (2018). Parameter estimation using the EM algorithm for symmetric stable random variables and sub-Gaussian random vectors, Journal of Statistical Theory and Applications, 17(3),1-20,
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
%\note{The \code{S&P} and \code{IPC} data are among 9 years of daily returns consists of 2535 observations %from 1/4/2000 to 9/22/2009. This data set is available at the website of Yahoo finance. We note that the EM %algorithm is robust with respect to the initial values.}
%\seealso{}
\examples{
# Here we follow for applying the EM algorithm to Z=(x1, x2)^T using the
# initial values alpha0=1, Sigma0=matrix(c(0.75,0.25,0.25,0.75),2,2), and
# Mu0=(0.5,0.5)^T.
library("stabledist")
x1<-urstab(100,1.2,0,1,2,0)
x2<-urstab(100,1.2,0,0.5,2,0)
z<-cbind(x1,x2)
\donttest{mfitstab.elliptical(z,1,matrix(c(0.75,0.25,0.25,0.75),2,2),c(0.5,0.5))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
