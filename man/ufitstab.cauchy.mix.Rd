\name{ufitstab.cauchy.mix}
\alias{ufitstab.cauchy.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ufitstab.cauchy.mix}
\description{estimates the parameters of a \code{k}-component mixture of Cauchy distributions. Assuming that \code{k} is known, given the vector of initial values of entire parameter space, it uses the EM algorithm to estimate the parameters of the \code{k}-component mixture of Cauchy distributions.}
\usage{ufitstab.cauchy.mix(y, k, omega0, beta0, sigma0, mu0)}
\arguments{
  \item{y}{vector of observations}
  \item{k}{number of components}
  \item{omega0}{initial value for weight vector to start the EM algorithm}
  \item{beta0}{initial value for skewness vector to start the EM algorithm}
  \item{sigma0}{initial value for scale vector to start the EM algorithm}
  \item{mu0}{initial value for location vector to start the EM algorithm}
}
\details{
Generally the EM algorithm seeks for the ML estimations when the log-likelihood function is not tractable mathematically. This is done by considering an extra missing (or latent) variable when the conditional expectation of the complete data log-likelihood given observed data and a guess of unknown parameter(s) is maximized. So, first we look for a stochastic representation. The representation given by the following proposition is valid for Cauchy distribution. Suppose \code{Y~ S_{0}(1,beta,sigma,mu)} and \code{T~S_{1}(1,1,1,0)} (here \code{S_0} and \code{S_1} refer to parameterizations \code{S_0} and \code{S_1}, respectively). Then \code{Y=sigma*(1-|beta|)*N/Z+sigma*beta*T+mu} where \code{N~Z~N(0,1)}. The random variables \code{N}, \code{Z}, and \code{P} are mutually independent.}
\value{
 \item{omega-bar}{a \code{k}-component vector of estimated values for the weight vector}
 \item{beta-bar}{a \code{k}-component vector of estimated values for the skewness vector}
 \item{sigma-bar}{a \code{k}-component vector of estimated values for the scale vector}
 \item{mu-bar}{a \code{k}-component vector of estimated values for the location vector}
}
\references{
Bjerkedal, T. (1960) Acquisition of resistance in guinea pigs infected with different doses of virulent tubercle bacilli, American Journal of Epidemiology, 72, 130-148.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{
We use the survival times in days of 72 guinea pigs infected with different doses of tubercle bacilli, see Bjerkedal (1960). We note that the EM algorithm is robust with respect to the initial values.}
%\seealso{}
\examples{
# In the following, we give an example that fits a two-component mixture of Cauchy distributions
# to the survival times (in days) of 72 guinea pigs through the EM algorithm. For this, the initial
# values are: omega_0=(0.65,0.35), sigma_0=(20,50), beta_0=(0.20,0.05), and mu_0=(95,210).
library("stabledist")
y<-c(10,33,44,56,59,72,74,77,92,93,96,100,100,102,105,107,107,108,108,108,
     109,112,121,122,122,124,130,134,136,139,144,146,153,159,160,163,163,
     168,171,172,176,113,115,116,120,183,195,196,197,202,213,215,216,222,
     230,231,240,245,251,253,254,255,278,293,327,342,347,361,402,432,458,
     555)
\donttest{ufitstab.cauchy.mix(y,2,c(0.65,0.35),c(0.20,0.05),c(20,50),c(95,210))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
