\name{ufitstab.sym}
\alias{ufitstab.sym}
\title{ufitstab.sym}
\description{estimates the parameters of a symmetric stable distribution through the EM algorithm, see Teimouri et al. (2018).}
\usage{ufitstab.sym(yy, alpha0, sigma0, mu0)}
\arguments{
  \item{yy}{a vector of observations}
  \item{alpha0}{initial value for the tail index parameter}
  \item{sigma0}{initial values for the scale parameter}
  \item{mu0}{initial values for the location parameter}
}
%\details{}
\value{
 \item{alpha  }{  estimated value of the tail index parameter}
 \item{sigma  }{  estimated value of the scale parameter}
 \item{mu  }{   estimated value of the location parameter}
}
\references{
Teimouri, M., Rezakhah, S., and Mohammadpour, A. (2018). Parameter estimation using the EM algorithm for symmetric stable random variables and sub-Gaussian random vectors, Journal of Statistical Theory and Applications, 17(3),1-20.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
%\note{
%The \code{CAC40} data are among 9 years of daily returns consists of 2535 observations from 1/4/2000 to 9/22%/2009. This set of data is available at the website of Yahoo finance. We note that the EM algorithm is robust %with respect to the initial values.}
%\seealso{}
%# By the following example, we apply the EM algorithm to CAC40 (Bourse de Paris) data.
%# The initial values are alpha_0=1.2, sigma_0=1, and mu_0=1.
%# CAC40<-c(-0.020807,0.003035,-0.036817,...,-0.054140,-14.87798)
\examples{
# By the following example, we apply the EM algorithm to n=50 iid realization of symmetric
# stable distribution with parameters alpha=1.2, sigma=1, and mu=1. The initial values
# are alpha_0=1.2, sigma_0=1, and mu_0=1.
library("stabledist")
y<-urstab(50,1.2,0,1,1,0)
\donttest{ufitstab.sym(y,1.2,1,1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
