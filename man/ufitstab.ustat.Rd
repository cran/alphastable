\name{ufitstab.ustat}
\alias{ufitstab.ustat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ufitstab.ustat}
\description{estimates the tail index and scale parameters of a symmetric and zero-location stable distribution using U-statistic proposed by Fan (2006) <DOI: 10.1080/03610920500439992>.}
\usage{ufitstab.ustat(x)}
\arguments{
   \item{x}{vector of observations}
}
%\details{}
\value{
 \item{alpha }{estimated value of the tail index parameter}
 \item{sigma }{estimated value of the scale parameter}
}
\references{
Fan, Z. (2006). Parameter estimation of stable distributions, Communications in Statistics-Theory and Methods, 35(2), 245-255.}

\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{The \code{ufitstab.ustat()} must be applied to a symmetric and zero-location stable distribution.}
%\seealso{}
\examples{
# We are estimating the parameters of a symmetric stable distribution. For this, firstly,
# we simulate a sample of n=100 iid realizations from stable distribution in S_1 parameterization
# with parameters alpha=1.2, beta=0, sigma=1, and mu=0.
x<-urstab(100,1.2,0,1,0,1)
ufitstab.ustat(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
