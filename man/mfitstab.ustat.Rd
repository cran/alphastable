\name{mfitstab.ustat}
\alias{mfitstab.ustat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mfitstab.ustat}
\description{estimates the parameters of a strictly bivariate stable distribution using approaches
proposed by Mohammadi et al. (2015)<doi.org/10.1007/s00184-014-0515-7> and Teimouri et al. (2017)<doi.org/10.1155/2017/3483827>. The estimated parameters are tail index and discretized spectral measure on \code{m} equidistant points located on unit sphere in \code{R^2}.}
\usage{mfitstab.ustat(u,m,method=method)}
\arguments{
  \item{u}{an \code{n} by 2 vector of observations}
  \item{m}{number of masses located on unit circle in \code{R^2}}
  \item{method}{integer values 1 or 2, respectively, corresponds to the method given by Teimouri et al. (2017) and Mohammadi et al. (2015)}
}
%\details{}
\value{
\item{alpha }{estimated value of tail index}
\item{mass}{estimated value of discrete spectral measure}
}
\references{
Mohammadi, M., Mohammadpour, A., and Ogata, H. (2015). On estimating the tail index and the spectral measure of multivariate alpha-stable distributions, Metrika, 78(5), 549-561.

Nolan. J. P. (2013). Multivariate elliptically contoured stable distributions: theory and estimation, Computational Statistics, 28(5), 2067-2089.

Teimouri, M., Rezakhah, S., and Mohammadpour, A. (2017). U-Statistic for multivariate stable distributions, Journal of Probability and Statistics, https://doi.org/10.1155/2017/3483827.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
%\note{
%For daily log-return (in percent) of 1247 closing prices of \code{AXP} (American Express Company) and %\code{MRK (Merck & Co. Inc.)} we refer readers to Nolan (2013). This set of data are between January 3, 2000 %and December 31, 2004.}
%\seealso{}
\examples{
# Here, for example, we are interested to estimate the parameters of a bivariate
# stable distribution. For this, two sets of n=400 iid realizations which are
# assumed to distributed jointly as a strictly bivariate stable distribution with
# tail index alpha=1.2 are simulated. Considering m=4, masses of the discrete spectral
# measure are addressed by s_j=(cos(2*pi(j-1)/m), sin (2*pi(j-1)/m)); for j=1,...,4.
library("nnls")
x1<-urstab(400,1.2,-0.50,1,0,0)
x2<-urstab(400,1.2,0.50,0.5,0,0)
z<-cbind(x1,x2)
mfitstab.ustat(z,4,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
