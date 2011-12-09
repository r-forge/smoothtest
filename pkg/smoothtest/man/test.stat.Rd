\name{test.stat}
\alias{test.stat}
\title{Smooth test statistics}
\description{Function that calculates the smooth test statistics asked for when using \code{smooth.test}}
\usage{test.stat(degree, distr = "unif", pars = c(0, 1), sample, method = "NONE", f = NA, moments = NA, typedistr = "cont", chol = FALSE, ntrials = NA)}
\arguments{\item{degree}{the order of the calculated smooth test statistic}
\item{distr}{specifies the distribution specified in the null hypothesis and must be one of "unif" (default), "pois", "logar", "exp", "norm", "negb", "ZIP", "laplace", "extrval", "geom", "ZTP", "genpareto", "betab", "gamma" or "logis"}
\item{pars}{if \code{method}="NONE", this is a numeric vector which contains the nuisance parameters as speficied under the null hypothesis.\cr
if \code{method}="MLE" or \code{method}="MME", these are the estimated nuisance parameters (when \code{pars} is NA they are estimated by the corresponding method).\cr
Default is \code{pars}=c(0,1) (the lower and upper bound of the interval over which the uniform distribution is defined).}
\item{sample}{a numeric vector of sample values. The sample of which the smooth test statistic is calculated.}
\item{method}{indicated the method for parameter estimation and must be one of "NONE" (default), "MLE", or "MME".}
\item{f}{a density function. This must only be provided when \code{distr} is none of the distributions listed above and when \code{moments}=NA. One then has to change the default value of \code{distr} ("unif") into "otherwise".}
\item{moments}{a vector specifying the the first 2*\code{degree} non-central moments of the null density function. This must only be provided when \code{distr} is none of the distributions listed above. One then has to change the default value of \code{distr} ("unif") into "otherwise".}
\item{typedistr}{must only be specified if \code{f} is given. It specifies whether \code{f} refers to a continuous ("cont", default) or a discrete ("disc") distribution.}
\item{chol}{logical; if TRUE a Cholesky decomposition is used.}
\item{ntrials}{a positive natural number representing the known and fixed number of trials in the beta binomial distribution, only used when \code{distr}="betab".}}
\details{The smooth test statistic is calculated for the random sample of observations given through \code{sample}.\cr
For *Details* how the null hypothesis is specified through \code{distr}, \code{method} and \code{pars} (and possibly \code{moments} or \code{f} and \code{typedistr}) see \code{\link{orth.poly}}.}
\value{A list with the following smooth test statistics.\cr
\item{Tstat}{a numeric vector containing the theoretically rescaled components. The last position in the vector contains the smooth test statistic of the requested order.}
\item{comp}{a numeric vector containing the unscaled components.}
\item{Sigma}{the theoretical variance covariance matrix.}
\item{Evar}{a numeric vector containing the empirical variances of the components.}}
\references{Rayner, J.C.W. and Best, D.J. (1989). Smooth Tests of Goodness-of-fit. Oxford University Press, New York.}
\author{Olivier Thas and Bert De Boeck}
\seealso{\code{\link{smooth.test}},\code{\link{orth.poly}}}	
\keyword{internal}
