\name{SMR}
\alias{dSMR}
\alias{pSMR}
\alias{qSMR}
\alias{rSMR}
\docType{package}
\title{The externally studentized normal midrange distribution}

\description{
Computes the probability density, the cumulative distribution function and the quantile function and generates random samples for the  externally studentized normal midrange distribution with the numbers means equal to \code{size}, the degrees of freedom equal to \code{df} and the number of
points of the Gauss-Legendre quadrature equal to \code{np}.}

\usage{
dSMR(x, size, df, np=32, log = FALSE)
pSMR(q, size, df, np=32, lower.tail = TRUE, log.p = FALSE)
qSMR(p, size, df, np=32, eps = 1e-13, maxit = 5000, lower.tail = TRUE, log.p = FALSE)
rSMR(n, size,  df = Inf)}

\arguments{
  \item{x, q}{
vector of quantiles \eqn{x \in R} and \eqn{q \in R}.}

  \item{p}{
vector of probabilities \eqn{(0, 1)}.
}
  \item{size}{
sample size. Only for  \eqn{size > 1}.
}
 \item{n}{
vector size to be simulated \eqn{n > 1}.
}
 \item{df}{
degrees of freedom \eqn{df > 0}.
}
 \item{np}{
number of points of the gaussian quadrature \eqn{np > 2}.
}
 \item{log, log.p}{
logical argument; if \code{TRUE}, the probabilities \eqn{p} are given as \eqn{log(p)}.
}
 \item{lower.tail}{
logical argument; if \code{TRUE}, the probabilities are \eqn{P[X \leq x]} otherside, \eqn{P[X \geq x].}
}
 \item{eps}{
 stopping criterion for Newton-Raphson's iteraction method.
 }
 \item{maxit}{
 maximum number of interaction in the Newton-Raphson method.
 }
}

\details{
Assumes np = 32 as default value for dSMR,  pSMR and qSMR. If \code{df} is not specified, it assumes the default value \code{Inf} in \code{rSMR}. When \code{df=1}, the convergence of the routines requires \code{np>250} to obtain the desired result  accurately.
The Midrange distribution has density

\deqn{f(\overline{q};n,\nu) =\int^{\infty}_{0}
\int^{x\overline{q}}_{-\infty}2n(n-1)x\phi(y) \phi(2x\overline{q}-y)[\Phi(2x\overline{q}-y)-\Phi(y)]^{n-2}f(x;\nu)dydx,}{f(q;n,\nu) = int^{\infty}_{0} int^{xq}_{- \infty} 2n(n-1)x\phi(y)\phi(2xq-y)[\Phi(2xq-y)-\Phi(y)]^{n-2}f(x;\nu)dydx}

where, \eqn{q} is the quantile of  externally studentized midrange distribution, \eqn{n} (\code{size}) is the sample size and \eqn{\nu} is the degrees of freedon.

The  externally studentized midrange distribution function is given by
\deqn{
F(\overline{q};n,\nu)=\int^{\overline{q}}_{-\infty}
\int^{\infty}_{0}\int^{x\overline{q}}_{-\infty}2n(n-1)x\phi(y) \phi(2xz-y)[\Phi(2xz-y)-\Phi(y)]^{n-2}f(x;\nu)dydxdz.
}{F(q;n,\nu) = int^{q}_{- \infty} int^{\infty}_{0} int^{xq}_{- \infty} 2n(n-1)x\phi(y)\phi(2xz-y)[\Phi(2xz-y)-\Phi(y)]^{n-2}f(x;\nu)dydxdz.}
where, \eqn{q} is the quantile of  externally studentized  midrange distribution, \eqn{n} (\code{size}) is the sample size and \eqn{\nu} is the degrees of freedon.
}

\value{\code{dSMR} gives the density, \code{pSMR} gives the cumulative distribution function, \code{qSMR} gives the quantile function, and \code{rSMR}
generates random deviates.}

\references{
Batista, BDO; Ferreira, DF. The  externally studentized normal midrange distribution. Submmited for publications. 2012.
}



\examples{
library(SMR)

#example 1:
x  <- 2
q  <- 1
p  <- 0.9
n  <- 30
size  <- 5
df <- 3
np <- 32
dSMR(x, size, df, np)
pSMR(q, size, df, np)
qSMR(p, size, df, np)
rSMR(n, size, df)

#example 2:
x  <- c(-1, 2, 1.1)
q  <- c(1, 0, -1.5)
p  <- c(0.9, 1, 0.8)
n  <- 10
size  <- 5
df <- 3
np <- 32
dSMR(x, size, df, np)
pSMR(q, size, df, np)
qSMR(p, size, df, np)
rSMR(n, size, df)
}

\keyword{midrange}
