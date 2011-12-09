\name{smooth.test}
\alias{smooth.test}
\alias{smoothtest}
\title{Smooth Tests}
\description{Performs the smooth test for the one-sample goodness-of-fit problem as described by Rayner et al.(2009). Both simple and composite null hypotheses can be tested. The maximum likelihood (MLE) and the method of moments (MME) methods for nuisance parameter estimation are implemented.}
\usage{smooth.test(sample, order, distr = "unif", method = "NONE", pars = c(0,1), B = 200, rescale=FALSE, f = NA, moments = NA, typedistr = "cont", chol=FALSE, output=TRUE, ntrials=NA)}
\arguments{\item{sample}{a numeric vector of sample observations}
\item{order}{the order of the test}
\item{distr}{specifies the distribution to be tested, and must be one of "unif" (default), "pois", "logar", "exp", "norm", "negb", "ZIP", "laplace", "extrval", "geom", "ZTP", "logis", "genpareto", "betab" or "gamma"}
\item{method}{indicates the method for parameter estimation and must be one of "NONE" (default), "MLE", or "MME"}
\item{pars}{a numeric vector which contains the nuisance parameters as speficied under the null hypothesis. This is only used if \code{method}="NONE". Default is \code{pars}=c(0,1) (the lower and upper bound of the interval over which the uniform distribution is defined).}
\item{B}{the number of bootstrap runs for p-value calculation. When \code{B}=NULL the asymptotic chi-squared distribution is used for p-value calculation.}
\item{rescale}{if \code{rescale}=TRUE the empirical variance is used to rescale the components; if \code{rescale}=FALSE (default) the theoretical variance under the null hypothesis is used.}
\item{f}{a density function. This must only be provided when \code{distr} is none of the distributions listed above and when \code{moments}=NA.}
\item{moments}{a vector specifying the first 2*\code{degree} non-central moments of the null density function. This must only be provided when \code{distr} is none of the distributions listed above and when \code{f}=NA.}
\item{typedistr}{must only be specified if \code{f} is given. It specifies whether \code{f} refers to a continuous ("cont", default) or a discrete ("disc") distribution.}
\item{chol}{logical; if TRUE the components result from a Cholesky decomposition. Default is FALSE.}
\item{output}{logical; if TRUE (default) an extensive output of the smooth test is given, otherwise, no output is given.}
\item{ntrials}{a positive natural number representing the fixed number of trials in the beta binomial distribution, only used when \code{distr}="betab".}}
\details{This function performs a smooth test of goodness of fit for both simple and composite null hypotheses. The nuisance parameters can be estimated by means of maximum likelihood (\code{method}="MLE) or method of moments (\code{method}="MME). The function also performs the individual component tests. These components may also be obtained from a Cholesky decomposition (\code{chol}=TRUE). The p-values are computed using the asymptotic chi-squared approximation (\code{B}=NULL) or using a parametric bootstrap approximation (\code{B}>0). The \code{smooth.test} function only works with a fixed order (\code{order}). See \code{\link{ddsmooth.test}} for the data-driven version.\cr

See the *Details* section of \code{\link{orth.poly}} for details on the use of \code{distr}, \code{pars}, \code{f} and \code{moments}.}
\value{A list with components:\cr
\item{statistics}{a numeric vector containing the (rescaled) components. The last position in the vector contains
the smooth test statistic.}
\item{p.value}{a numeric vector containing the p-values corresponding to \code{statistics}}
\item{par.est}{a numeric vector containing the nuisance parameter estimates}
\item{comp}{a numeric vector containing the (unscaled) components (this is exactly equal to \code{statistics} without the last entry).}
\item{Sigma}{the theoretical variance covariance matrix of the (unscaled) components, calculated under the null hypothesis}
\item{EVar}{a numeric vector containing the empirical variances of the (unscaled) components}}
\references{Rayner, J.C.W., Thas, O. and Best, D.J. (2009). Smooth Tests of Goodness-of-fit. Wiley, Singapore.}
\author{Olivier Thas and Bert De Boeck}
\seealso{\code{\link{orth.poly}}}
\examples{
x<-rnorm(100,mean=2,sd=2)
## does x come from a normal distribution? Perform a smooth test of order 4 and use MLE for the 
## estimation of the mean and the variance of the normal distribution.
smooth.test(sample=x,distr="norm",order=4,method="MLE")
## Smooth goodness-of-fit test
## Null hypothesis:  norm  against  4 th order alternative
## Nuisance parameter estimation:  MLE 
## Parameter estimates:  2.256117 4.308846  ( MEAN VAR )

## Smooth test statistic S_k =  1.857947  p-value =  0.26 
##     3 th compenent V_k =  0.5342814  p-value =  0.51 
##     4 th compenent V_k =  1.323665  p-value =  0.17 

## All p-values are obtained by the boostrap with  200  runs

## does x come from a logistic distribution? Perform a smooth test of order 4 and use MME for the 
## estimation of the locatin and the scale parameters of the logistic distribution.
smooth.test(sample=x,distr="logis",order=4,method="MME")
## Smooth goodness-of-fit test
## Null hypothesis:  logis  against  4 th order alternative
## Nuisance parameter estimation:  MME 
## Parameter estimates:  2.256117 1.144435  ( MU SIGMA )

## Smooth test statistic S_k =  1.190191  p-value =  0.19 
##     3 th compenent V_k =  0.1355062  p-value =  0.62 
##     4 th compenent V_k =  1.054685  p-value =  0.04 

## All p-values are obtained by the boostrap with  200  runs
}
\keyword{smooth}
