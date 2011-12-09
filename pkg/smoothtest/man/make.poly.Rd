\name{make.poly}
\alias{make.poly}
\alias{make.poly.disc}
\title{Generalised Emerson Recurrence Relation}
\description{Constructs polynomials that are orthonormal to a given density function (continuous or discrete) by using the generalised Emerson recurrence relations of Rayner et al. (2008). These functions are basically internally used functions, used by \code{orth.poly}.}
\usage{make.poly(f, k, moments = NA, coef=F, norm=T)
make.poly.disc(f, k, moments = NA, coef=F, norm=T)}
\arguments{\item{f}{a function which is typically a density function against to which the orthonormal polynomials have to be generated. Only used when \code{moments}=NA (default).}
\item{k}{an integer which specifies up to which degree the orthonormal polynomials are generated.}
\item{moments}{The default is \code{moments}=NA. When \code{moments} is a vector of length 2*\code{k}, it specifies the first 2*\code{k} non-central moments of a density function against to which the orthonormal polynomials have to be generated.}
\item{coef}{logical; if FALSE (default) the function returns a function-object; otherwise a matrix with the coefficients of the polynomials are returned.}
\item{norm}{logical; if TRUE (default) normalised polynomials are computed. \code{norm}=FALSE can only be used with \code{coef}=TRUE.}}
\details{If the density function \code{f} is continuous \code{make.poly} should be used; if it is discrete \code{make.poly.disc} should be used.}
\value{\code{make.poly} and \code{make.poly.disc} with \code{coef}=FALSE give a list with \code{k}+1 elements that contain the orthonormal polynomials as functions, starting form degree 0 up to degree \code{k}.\cr
When \code{coef}=FALSE the functions return a matrix with columns corresponding to the coefficients of the orthonormal (orthogonal) polynomials, starting form degree 0 up to degree \code{k}.\cr
The rows of the matrix correpond with the coefficients of the terms of 0-th order up to \code{k}-th order of those polynomials.}
\references{Rayner, J.C.W., Thas, O. and De Boeck B. (2008). A Generalized Emerson Recurrence Relation. Aust.N.Z.J.Stat.50(3), 235-240}
\author{Olivier Thas and Bert De Boeck}
\seealso{\code{\link{orth.poly}}}
\keyword{math}
\keyword{symbolmath}
