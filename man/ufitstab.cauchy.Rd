\name{ufitstab.cauchy}
\alias{ufitstab.cauchy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ufitstab.cauchy}
\description{estimates the parameters of the Cauchy distribution. Given the initial values of the skewness, scale, and location parameters, it uses the EM algorithm to estimate the parameters of the Cauchy distribution.}
\usage{ufitstab.cauchy(y, beta0, sigma0, mu0, param)}
\arguments{
  \item{y}{vector of observations}
  \item{beta0}{initial value of skewness parameter to start the EM algorithm}
  \item{sigma0}{initial value of scale parameter to start the EM algorithm}
  \item{mu0}{initial value of location parameter to start the EM algorithm}
  \item{param}{kind of parameterization; must be 0 or 1 for \code{S_0} and \code{S_1} parameterizations, respectively}
}
\details{
Generally the EM algorithm seeks for the ML estimations when the log-likelihood function is not tractable mathematically. This is done by considering an extra missing (or latent) variable when the conditional expectation of the complete data log-likelihood given observed data and a guess of unknown parameter(s) is maximized. So, first we look for a stochastic representation. The representation given by the following proposition is valid for Cauchy distribution. Suppose \code{Y~S_0(1,beta,sigma,mu)} and \code{T~S_{1}(1,1,1,0)} (here \code{S_0} and \code{S_1} refer to parameterizations \code{S_0} and \code{S_1}, respectively). Then \code{Y=sigma*(1-|beta|)*N/Z+sigma*beta*T+mu} where \code{N~Z~N(0,1)}. The random variables \code{N}, \code{Z}, and \code{T} are mutually independent.}
\value{
 \item{beta}{estimated value of the skewness parameter}
 \item{sigma}{estimated value of the scale parameter}
 \item{mu}{estimated value of the location parameter}
}
\references{
Davidian, M. and Giltinan, D.M. (1995). Nonlinear Mixed Effects Models for Repeated Measurement Data, Chapman and Hall.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{
The set of data considered here is large recorded intensities (in Richter scale) of the earthquake at seismometer locations in western North America between 1940 and 1980, see Davidian and Giltinan (1995). Among the features, we focus on the 182 distances from the seismological measuring station to the epicenter of the earthquake (in km) as the variable of interest. This set of data can be found in package \code{nlme}. We note that \code{ufitstab.cauchy()} is robust with respect to the initial values.}
%\seealso{}
\examples{
# In the following example, using the initial values beta_0=0.5, sigma_0=5, and mu_0=10,
# we apply the EM algorithm to estimate the parameters of Cauchy distribution fitted to
# the earthquake data given by the vector y.
y<-c(7.5,  8.8,   8.9,   9.4,   9.7,   9.7,   10.5,  10.5,  12.0,  12.2,  12.8,  14.6,
     14.9,  17.6,  23.9,  25.0,  2.9,   3.2 ,  7.6,   17.0,  8.0,   10.0,  10.0,  8.0,
     19.0,  21.0,  13.0,  22.0,  29.0,  31.0,  5.8,   12.0,  12.1,  20.5,  20.5,  25.3,
     35.9,  36.1,  36.3,  38.5,  41.4,  43.6,  44.4,  46.1,  47.1,  47.7,  49.2,  53.1,
     4.0,   10.1,  11.1,  17.7,  22.5,  26.5,  29.0,  30.9,  37.8,  48.3,  62.0,  50.0,
     16.0,  62.0,  1.2,   1.6,   9.1,   3.7,   5.3,   7.4,   17.9,  19.2,  23.4,  30.0,
     38.9,  10.8,  15.7,  16.7,  20.8,  28.5,  33.1,  40.3,  8.0,   32.0,  30.0,  31.0,
     16.1,  63.6,  6.6,   9.3,   13.0,  17.3,  105.0, 112.0, 123.0, 5.0,   23.5,  26.0,
     0.5,   0.6,   1.3,   1.4,   2.6,   3.8,   4.0,   5.1,   6.2,   6.8,   7.5,   7.6,
     8.4,   8.5,   8.5,   10.6,  12.6,  12.7,  12.9,  14.0,  15.0,  16.0,  17.7,  18.0,
     22.0,  22.0,  23.0,  23.2,  29.0,  32.0,  32.7,  36.0,  43.5,  49.0,  60.0,  64.0,
     105.0, 122.0, 141.0, 200.0, 45.0,  130.0, 147.0, 187.0, 197.0, 203.0, 211.0, 17.0,
     19.6,  20.2,  21.1,  88.0,  91.0,  12.0,  148.0, 42.0,  85.0,  21.9,  24.2,  66.0,
     87.0,  23.4,  24.6,  25.7,  28.6,  37.4,  46.7,  56.9,  60.7,  61.4,  62.0,  64.0,
     82.0,  107.0, 109.0, 156.0, 224.0, 293.0, 359.0, 370.0, 25.4,  32.9,  92.2,  45.0,
     145.0, 300.0)
library("stabledist")
\donttest{ufitstab.cauchy(y,0.5,5,10,0)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
