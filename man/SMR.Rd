\name{SMR}
\alias{dSMR}
\alias{pSMR}
\alias{qSMR}
\alias{rSMR}
\docType{package}
\title{The externally studentized normal midrange distribution}

\description{
Computes the probability density, the cumulative distribution function and the quantile function and generates random samples for the  externally studentized normal midrange distribution with the numbers means equal to \code{n}, the degrees of freedom equal to \code{nu} and the number of 
points of the Gauss-Legendre quadrature equal to \code{np}. 
}

\usage{
dSMR(x, n, nu, np=32)
pSMR(x, n, nu, np=32)
qSMR(p, n, nu, np=32)
rSMR(N, n,  nu=Inf)
}

\arguments{
  \item{x}{
vector of quantiles \eqn{x \in R}.
}
  \item{p}{
vector of probabilities \eqn{(0, 1)}.
}
  \item{n}{
number of means. Only for  \eqn{n > 1}.
}
 \item{N}{
vector size to be simulated \eqn{N > 1}.
}
 \item{nu}{
degrees of freedom \eqn{\nu > 0}.
}
 \item{np}{
number of points of the gaussian quadrature \eqn{np > 2}. 
}     
}

\details{
Assumes np = 32 as default value for dSMR,  pSMR and qSMR. If \code{nu} is not specified, it assumes the default value \code{Inf} in \code{rSMR}. When \code{nu=1}, the convergence of the routines requires \code{np>250} to obtain the desired result  accurately.
The Midrange distribution has density
\deqn{f(\overline{q};n,\nu) =\int^{\infty}_{0}
\int^{x\overline{q}}_{-\infty}2n(n-1)x\phi(y) \phi(2x\overline{q}-y)[\Phi(2x\overline{q}-y)-\Phi(y)]^{n-2}f(x;\nu)dydx,}
where, \eqn{q} is the quantile of  externally studentized midrange distribution, \eqn{n} is the numbers of means and \eqn{\nu} is the degrees of freedon.

The  externally studentized midrange distribution function is given by
\deqn{
F(\overline{q};n,\nu)=\int^{\overline{q}}_{\infty}
\int^{\infty}_{0}\int^{x\overline{q}}_{-\infty}2n(n-1)x\phi(y) \phi(2xz-y)[\Phi(2xz-y)-\Phi(y)]^{n-2}f(x;\nu)dydxdz.
}
where, \eqn{q} is the quantile of  externally studentized  midrange distribution, \eqn{n} is the numbers of means and \eqn{\nu} is the degrees of freedon.
}

\value{\code{dSMR} gives the density, \code{pSMR} gives the cumulative distribution function, \code{qSMR} gives the quantile function, and \code{rSMR}
generates random deviates.}

\references{
Batista, BDO; Ferreira, DF. The  externally studentized normal midrange distribution. Submmited for publications. 2012.
}

\author{
Ben Deivide de Oliveira Batista and Daniel Furtado Ferreira.

Maintainer: Daniel Furtado Ferreira \code{<danielff@dex.ufla.br>}
}

\seealso{
Package homepage: \code{<www.dex.ufla.br/~danielff/r_resources.html>}
}

\examples{
library(SMR)
x  <- 2
p  <- 0.9
N  <- 100
n  <- 5
nu <- 3
np <- 32
dSMR(x, n, nu, np)
pSMR(x, n, nu, np)
qSMR(p, n, nu, np)
rSMR(N, n, nu)
}

\keyword{midrange}