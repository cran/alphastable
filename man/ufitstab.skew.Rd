\name{ufitstab.skew}
\alias{ufitstab.skew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ufitstab.skew}
\description{using the EM algorithm, it estimates the parameters of skew stable distribution.}
\usage{ufitstab.skew(y, alpha0, beta0, sigma0, mu0, param)}
\arguments{
 \item{y}{vector of observations}
 \item{alpha0}{initial value of tail index parameter to start the EM algorithm}
 \item{beta0}{initial value of skewness parameter to start the EM algorithm}
 \item{sigma0}{initial value of scale parameter to start the EM algorithm}
 \item{mu0}{initial value of location parameter to start the EM algorithm}
 \item{param}{kind of parameterization; must be 0 or 1 for \code{S_0} and \code{S_1} parameterizations, respectively}
 }
\details{
For any skew stable distribution we give a new representation by the following. Suppose \code{Y~ S_{0}(alpha, beta, sigma, mu), P~ S_{1}(alpha/2,1,(cos(pi*alpha/4))^(2/alpha),0)}, and \code{V~ S_{1}(alpha,1,1,0)}. Then, \code{Y=eta*(2P)^(1/2)*N+theta*V+ mu-lambda}, where \code{eta=sigma*(1-|beta|)^(1/alpha)}, \code{theta=sigma*sign(beta)*|beta|^(1/alpha)}, \code{lambda=sigma*beta*tan(pi*alpha/2)}, and \code{N~N(0,1)} follows a skew stable distribution. All random variables \code{N}, \code{P}, and \code{V} are mutually independent.}
\value{
 \item{alpha}{  estimated value of the tail index parameter}
 \item{beta}{   estimated value of the skewness parameter}
 \item{sigma}{  estimated value of the scale parameter}
 \item{mu}{estimated value of the location parameter}
}
\references{
Buckle, D. J. (1995). Bayesian inference for stable distributions, Journal of the American Statistical Association, 90(430), 605-613.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{
Daily price returns of Abbey National shares between 31/7/91 and 8/10/91 (including \code{n=50} business days). By assuming that \code{p_{t}} denotes the price at \code{t}-th day, the price return at \code{t}-th day is defined as \code{(p_{t-1}-p_{t})/p_{t-1}}; for \code{t=2,...,n}, see Buckle (1995). We note that the EM algorithm is robust with respect to the initial values.}
%\seealso{}
\examples{
# For example, We use the daily price returns of Abbey National shares. Using the initial
# values as alpha_{0}=0.8, beta_{0}=0, sigma_{0}=0.25, and mu_{0}=0.25.
price<-c(296,296,300,302,300,304,303,299,293,294,294,293,295,287,288,297,
         305,307,304,303,304,304,309,309,309,307,306,304,300,296,301,298,
         295,295,293,292,307,297,294,293,306,303,301,303,308,305,302,301,
         297,299)
x<-c()
n<-length(price)
for(i in 2:n){x[i]<-(price[i-1]-price[i])/price[i-1]}
library("stabledist")
\donttest{ufitstab.skew(x[2:n],0.8,0,0.25,0.25,1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{skew stable distribution }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
