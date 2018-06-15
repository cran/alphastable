\name{ufitstab.sym.mix}
\alias{ufitstab.sym.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ufitstab.sym.mix}
\description{estimates the parameters of a \code{k}-component mixture of symmetric stable distributions, Teimouri et al. (2018) <doi.org/10.1080/03610918.2017.1288244>. Having \code{k} and given the vector of initial values of entire parameter space, it uses some type of the EM algorithm (ECME) to estimate the parameters of mixture of symmetric stable distributions.}
\usage{ufitstab.sym.mix(yy, k, omega0, alpha0, sigma0, mu0)}
\arguments{
  \item{yy}{vector of observations}
  \item{k}{ number of components}
  \item{omega0}{vector of  initial values for weights}
  \item{alpha0}{vector of  initial values for tail indices}
  \item{sigma0}{vector of  initial values for scale parameters}
  \item{mu0}{vector of initial values for location parameters}
}
%\details{}
\value{
 \item{omega-bar }{a \code{k}-component vector of estimated values for the weight vector}
 \item{alpha-bar }{a \code{k}-component vector of estimated values for the tail index vector}
 \item{sigma-bar }{a \code{k}-component vector of estimated values for the scale vector}
 \item{mu-bar    }{a \code{k}-component vector of estimated values for the location vector}
 \item{membership}{a \code{k} by \code{n} matrix whose entries are 1 and 0; for \code{i}-th row, if \code{j}-th column is one; then \code{i}-th observation belongs to the \code{j}-th component}
}
\references{
Teimouri, M., Rezakhah, S., and Mohammadpour, A. (2018). EM algorithm for symmetric stable mixture model, Communications in Statistics-Simulation and Computation, 47(2), 582-604.
}
\author{Mahdi Teimouri, Adel Mohammadpour, and Saralees Nadarajah}
\note{
The set of data, used here, is the velocities of 82 distant galaxies, diverging from our own galaxy, called here galaxy. These data are available at \code{https://people.maths.bris.ac.uk/mapjg/mixdata}. We note that \code{ufitstab.sym.mix()} is robust with respect to the initial values.}
%\seealso{}
\examples{
# In what follows, we apply the EM algorithm to estimate the parameters of the
# mixture of symmetric stable distributions. For this, the initial values for
# fitting a three-component mixture of symmetric stable distribution to the
# galaxy data are: (0.1,0.35,0.55) for weight vector, (1.2,1.2,1.2) for tail
# index vector, (1,1,1) for scale vector, and (8,20,22) for the location vector.
galaxy<-c(9.172,9.350,9.483,9.558,9.775,10.227,10.406,16.084,16.170,18.419,
          18.552,18.600,18.927,19.052,19.070,19.330,19.343,19.349,19.440,
          19.473,19.529,19.541,19.547,19.663,19.846,19.856,19.863,19.914,
          19.918,19.973,19.989,20.166,20.175,20.179,20.196,20.215,20.221,
          20.415,20.629,20.795,20.821,20.846,20.875,20.986,21.137,21.492,
          21.701,21.814,21.921,21.960,22.185,22.209,22.242,22.249,22.314,
          22.374,22.495,22.746,22.747,22.888,22.914,23.206,23.241,23.263,
          23.484,23.538,23.542,23.666,23.706,23.711,24.129,24.285,24.289,
          24.366,24.717,24.990,25.633,26.960,26.995,32.065,32.789,34.279)
library("stabledist")
\donttest{ufitstab.sym.mix(galaxy,3,c(0.1,0.35,0.55),c(1.2,1.2,1.2),c(1,1,1),c(8,20,22))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
