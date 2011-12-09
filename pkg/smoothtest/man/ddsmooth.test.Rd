\name{ddsmooth.test}
\alias{ddsmooth.test}
\alias{ddsmoothtest}
\title{Data-driven smooth tests}
\description{Performs the data-driven smooth test for the one-sample goodness-of-fit problem. See Rayner et al.(2009) for more details. Both simple and composite null hypotheses can be tested. The maximum likelihood (MLE) and the method of moments (MME) methods for nuisance parameter estimation are implemented. \cr
Several model selection criteria can be used: AIC, BIC and MISE.}
\usage{ddsmooth.test(sample, max.order, horizon = "order", criterion = "AIC", distr = "norm", method = "NONE", pars = c(0, 1), B = 200, f = NA, typedistr = "cont", plot = FALSE, plot.range = NULL, ntrials = NA, ...)}
\arguments{\item{sample}{a numeric vector of sample observations}
\item{max.order}{the maximal order of the test that can be selected by the model selection criterion}
\item{horizon}{\code{horizon}="order" (default) refers to order-selection, and \code{horizon}="subset" refers to subset selection}
\item{criterion}{a character specifying the model selection criterion: "AIC" (default), "BIC" or "MISE"}
\item{distr}{specifies the distribution to be tested, and must be one of "unif" (default), "pois", "logar", "exp", "norm", "negb", "ZIP", "laplace", "extrval", "geom", "ZTP", "logis", "genpareto", "betab" or "gamma"}
\item{method}{indicates the method for parameter estimation and must be one of "NONE" (default), "MLE", or "MME"}
\item{pars}{a numeric vector which contains the nuisance parameters as speficied under the null hypothesis. This is only used if \code{method}="NONE". Default is \code{pars}=c(0,1) (the lower and upper bound of the interval over which the uniform distribution is defined).}
\item{B}{the number of bootstrap runs for p-value calculation. When \code{B}=NULL the asymptotic chi-squared distribution is used for p-value calculation.}
\item{f}{a density function. This must only be provided when \code{distr} is none of the distributions listed
 above and when \code{moments}=NA. }
\item{typedistr}{must only be specified if \code{f} is given. It specifies whether \code{f} refers to a continuous ("cont", default) or a discrete ("disc") distribution.}
\item{plot}{logical; when \code{plot}=TRUE a graph is the improved density estimate is plotted. Default is \code{plot}=FALSE.}
\item{plot.range}{a vector of length 2, specifying the lower and upper bound of the range over which the improved density estimated has to be plotted.} 
\item{ntrials}{a positive natural number representing the fixed number of trials in the beta binomial distribution, only used when \code{distr}="betab". }
\item{...}{additional arguments forwarded to the plot function.}}
\details{This function performs a data-driven smooth test of goodness of fit for both simple and composite null hypotheses. The nuisance parameters can be estimated by means of maximum likelihood (\code{method}="MLE) or method of moments (\code{method}="MME). The function also performs the individual component tests. The p-values are computed using the asymptotic chi-squared approximation (\code{B}=NULL) or using a parametric bootstrap approximation (\code{B}>0).\cr

The components are selected from a horizon that is specified by a maximal order (\code{max.order}), a model selection criterion (\code{criterion}) and a type of selection (\code{horizon}="order" for order selection, and \code{horizon}="subset" for subset selection).\cr

See \code{\link{smooth.test}} for the smooth test with a fixed order. See the *Details* section of \code{\link{orth.poly}} for details on the use of \code{distr}, \code{pars}, \code{f} and \code{moments}.\cr
}
\value{A list with components: \cr
  \item{stat}{the test statistic}
  \item{model}{a vector with the selected orders}
  \item{p.value}{the p-value}
  \item{par.est}{a numeric vector containing the nuisance parameter estimates}
  \item{dnull}{a vector with the densities of the hypothesised distribution, evaluated on a set of points within \code{plot.range}; this is only non-empty when \code{plot}=TRUE.}
  \item{dimp}{a vector with the densities of the improved density estimate, evaluated on a set of points within \code{plot.range}; this is only non-empty when \code{plot}=TRUE.}
}
\references{Rayner, J.C.W., Thas, O. and Best, D.J. (2009). Smooth Tests of Goodness-of-fit. Wiley, Singapore.}
\author{Olivier Thas and Bert De Boeck }
\seealso{\code{\link{orth.poly}}, \code{link{smooth.test}}}
\examples{
# Example 10.4.1: Industrial data series (4)
data(industrial4)
ddsmooth.test(industrial4,distr="norm",method="MME",horizon="subset",criterion="AIC",max.order=6,B=1000,plot=TRUE)
# the output
#
#Data-Driven Smooth goodness-of-fit test
#Null hypothesis: norm against 6 th order alternative
#Nuisance parameter estimation: MME 
#Parameter estimates: 0.63364 0.3388172  ( mean sd )
#Horizon: subset 
#Selection criterion: AIC 
#
#Data-Driven Smooth test statistic S_k = 3.999744  p-value = 0.0522 
#Selected model: 5 
#
#All p-values are obtained by the bootstrap with 1000 runs

data(aluminium)
ddsmooth.test(aluminium,distr="pois",method="MME",horizon="subset",max.order=5,criterion="MISE",B=1000,plot=TRUE)
# the output
#
#Data-Driven Smooth goodness-of-fit test
#Null hypothesis: pois against 5 th order alternative
#Nuisance parameter estimation: MME 
#Parameter estimates: 1.999  ( lambda )
#Horizon: subset 
#Selection criterion: MISE 
#
#Data-Driven Smooth test statistic S_k = 4.709159  p-value = 0.142 
#Selected model: 2 
#
#All p-values are obtained by the bootstrap with 1000 runs
}
\keyword{smooth}
\keyword{htest}

