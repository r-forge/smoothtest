\name{orth.poly}
\alias{orth.poly}
\alias{orthpoly}
\title{Orthonormal polynomials}
\description{This function generates polynomials that are orthonormal to the density function of a specified distribution. Some distributions are indexed by nuisance parameters which may be estimated from sample observations. Both maximum likelihood estimation and method of moments estimation are implemented.}
\usage{orth.poly(degree, distr = "unif", pars = c(0, 1), sample = NA, method = "NONE", f = NA, moments = NA, typedistr = "cont", ntrials = NA)}
\arguments{\item{degree}{an integer that specifies the degree of the orthonormal polynomial.}
\item{distr}{a character string that indicates the distribution to which the orthonormal polynomials have to be generated. It has to be chosen from a list; see *Details*.}
\item{pars}{a vector containing the nuissance parameters of the concerning distribution. See *Details*.}
\item{sample}{a vector that contains sample observations from which the nuissance parameters can be estimated.}
\item{method}{a character string indicating the method of estimation of the nuisance parameters. See *Details*.}
\item{f}{a function. Here you can provide a density function directly. This only works in combination with \code{distr}="otherwise".}
\item{moments}{a vector of length 2*\code{degree} that specifies the first 2*\code{degree} non-central moments of the density function.}
\item{typedistr}{either "cont" (continuous) or "disc" (discrete), indicating the type of the distribution. This is only necessary when the density function is given through \code{f}.}
\item{ntrials}{a positive natural number representing the known and fixed number of trials in the beta binomial distribution. This is only used when \code{distr}="betab".}}
\details{The character string \code{distr} must be chosen from "logis" (logistic), "exp" (exponential), "norm" (normal), "unif" (uniform), "pois" (Poisson), "ZIP" (zero-inflated Poisson), "negb" (negative binomial), "logar" (logarithmic), "laplace" (Laplace), "extrval" (extreme value), "geom" (geometric), "ZTP" (zero-truncated Poisson), "genpareto" (generalised Pareto), "betab" (beta-binomial), "gamma" (gamma) or "otherwise". The latter is used when the density function is specified through \code{f} or \code{moments}.\cr
The nuisance parameters must be specified through \code{pars} when they are known, or through \code{sample} and \code{method} when they have to be estimated. If \code{f} or \code{moments} are used, and thus \code{distr}="otherwise", then \code{pars}, \code{sample} and \code{method} are not used.\cr

The meaning of \code{pars} depends on the distribution:\cr
If \code{distr}="logis", then \code{pars} = c('location', 'scale') in the density function \code{dlogis}. \cr
If \code{distr}="exp", then \code{pars} = 'rate' in the density function \code{dexp}.  \cr
If \code{distr}="norm", then \code{pars} = c('mean', 'sd') in the density function \code{dnorm}. \cr
If \code{distr}="unif", then \code{pars} = c('min', 'max') in the density function \code{dunif}. \cr
If \code{distr}="pois", then \code{pars} = 'lambda' in the density function \code{dpois}. \cr
If \code{distr}="negb", then \code{pars} = c('prob', 'size') in the density function \code{dnbinom}. \cr
If \code{distr}="ZIP", then \code{pars} = c('phi', 'lambda') in the density function \code{dzipois}. \cr
If \code{distr}="logar", then \code{pars} = 'prob' in the density function \code{dlog}. \cr
If \code{distr}="laplace", then \code{pars} = c('location', 'scale') in the density function \code{dlaplace}. \cr
If \code{distr}="extrval", then \code{pars} = c('location', 'scale') in the density function \code{dgumbel}. \cr
If \code{distr}="geom", then \code{pars} = 'prob' in the density function \code{dgeom}. \cr
If \code{distr}="ZTP", then \code{pars} = '\eqn{\lambda}{lambda}' in the density function of the zero-truncated Poisson distribution:\cr
\deqn{f(x) = \lambda^x/(x!(\exp(\lambda)-1))}{f(x) = lambda^x/(x!(exp(lambda)-1))} 
for \eqn{x > 0}, and \eqn{f(0)=0}.\cr
If \code{distr}="genpareto", then \code{pars} = c('scale', -'shape') in the density function \code{dgpd}.\cr
If \code{distr}="betab", then \code{pars} = c('shape1', 'shape2') in the density function \code{dbetabin.ab}.\cr
If \code{distr}="gamma", then \code{pars} = c('shape', 'scale') in the density function \code{dgamma}.\cr

When the nuisance parameters are known and specified through \code{pars}, you should set \code{method}="NONE". If the nuisance parameters have to be estimated from a sample (which is given as the vector \code{sample}), the estimation method must be chosen. For \code{method}="MLE" maximum likelihood estimation is used, for \code{method}="MME" method of moments estimation is used.\cr

The vector \code{moments} is only relevant when \code{distr}="otherwise" and thus none of the implemented distributions is used. It specifies the first 2*\code{degree} non-central moments of the density function.\cr

The function \code{f} and the character string \code{typedistr} are only relevant when \code{distr}="otherwise" and \code{moments}=NA and thus none of the implemented distributions is used. In this case not only the density function \code{f} has to be specified, but also the type of the density function; i.e. \code{typedistr}="cont" for a continuous distribution (default), and \code{typedistr}="disc" for a discrete distribution.}
\value{A function: the orthonormol polynomial of order \code{degree}.}
\references{Rayner, J.C.W., Thas, O. and De Boeck B. (2007). A Generalised Emerson Recurrence Relation. To be published.}
\author{Olivier Thas and Bert De Boeck}
\section{Warning}{For discrete distributions the support is limited to 0:500.}
\examples{
# generates the second order Legendre polynomial
h2<-orth.poly(degree=2,distr="unif",pars=c(0,1))
# plot this function
u<-seq(0,1,0.01)
plot(u,h2(u),type="l")

# the third order Hermite polynomial based on estimates of the mean and the variance of the normal distribution
x<-rnorm(100)
h3<-orth.poly(degree=3,distr="norm",method="MLE",sample=x)
# plot this function
plot(sort(x),h3(sort(x)),type="l")

# the fourth order polynomial which is orthonormal to the function given by
d<-function(x) ifelse((x<=1)&(x>=0),(1+2*x)/2,0)
h4<-orth.poly(degree=4,distr="otherwise",f=d)
# plot this function
plot(x,h4(x),type="l")
}
\seealso{\code{\link{make.poly}}}
\keyword{math}
\keyword{symbolmath}
