% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{mymcmc}
\alias{mymcmc}
\title{Adaptive Multivariate MCMC sampler}
\usage{
mymcmc(
  x,
  nlogf,
  m1,
  m2 = m1,
  m3,
  scl1 = 0.5,
  scl2 = 2,
  skip = 1,
  covm = 0,
  nfcn = 0,
  plot = FALSE,
  plot.range = 0
)
}
\arguments{
\item{x}{a list with components 'label' (of mode character), 'est' (the
parameter vector with the initial guess), 'low' (vector with lower bounds),
and 'upp' (vector with upper bounds)}

\item{nlogf}{negative log of the density function (not necessarily
normalized) for which samples are to be obtained}

\item{m1}{length of first pilot run (not used when covm supplied)}

\item{m2}{length of second pilot run (not used when covm supplied )}

\item{m3}{length of final run}

\item{scl1}{scale for covariance of mv normal proposal (second pilot run)}

\item{scl2}{scale for covariance of mv normal proposal (final run)}

\item{skip}{number of cycles skipped for graphical output}

\item{covm}{covariance matrix for multivariate normal proposal
distribution. If supplied, all pilot runs will be skipped and a run of
length m3 will be produced. Useful to continue a simulation from a given
point with specified covm}

\item{nfcn}{number of function calls}

\item{plot}{logical variable. If TRUE the chain and the negative log
density (nlogf) is plotted. The first m1+m2 cycles are shown in green,
other cycles in red}

\item{plot.range}{[Not documented. Leave as default]}
}
\value{
list with the following components: \item{f }{ values of nlogf for
the samples obtained } \item{mcmc }{ the chain (samples obtained) }
\item{covm }{ current covariance matrix for mv normal proposal
distribution}
}
\description{
This function generates MCMC-based samples from a (posterior) density f
(not necessarily normalized). It uses a Metropolis algorithm in conjunction
with a multivariate normal proposal distribution which is updated
adaptively by monitoring the correlations of succesive increments of at
least 2 pilot chains. The method is described in De Gunst, Dewanji and
Luebeck (submitted). The adaptive method is similar to the one proposed in
Gelfand and Sahu (JCGS 3:261--276, 1994).
}
\details{
standard output reports a summary of the acceptance fraction, the current
values of nlogf and the parameters for every (100*skip) th cycle. Plotted
chains show values only for every (skip) th cycle.
}
\note{
This function is part of the Bhat exploration tool
}
\examples{

# generate some Poisson counts on the fly
  dose <- c(rep(0,50),rep(1,50),rep(5,50),rep(10,50))
  data <- cbind(dose,rpois(200,20*(1+dose*.5*(1-dose*0.05))))

# neg. log-likelihood of Poisson model with 'linear-quadratic' mean: 
  nlogf <- function (x) { 
  ds <- data[, 1]
  y  <- data[, 2]
  g <- x[1] * (1 + ds * x[2] * (1 - x[3] * ds)) 
  return(sum(g - y * log(g)))
  }

# start MCMC near mle
  x <- list(label=c("a","b","c"), est=c(20, 0.5, 0.05), low=c(0,0,0), upp=c(100,10,.1))
# samples from posterior density (~exp(-nlogf))) with non-informative
# (random uniform) priors for "a", "b" and "c".
out <- mymcmc(x, nlogf, m1=2000, m2=2000, m3=10000, scl1=0.5, scl2=2, skip=10, plot=TRUE)
# start MCMC from some other point: e.g. try x$est <- c(16,.2,.02)

}
\references{
too numerous to be listed here
}
\seealso{
\code{\link{dfp}}, \code{\link{newton}},
\code{\link{logit.hessian}}
}
\author{
E. Georg Luebeck (FHCRC)
}
\keyword{iteration}
\keyword{methods}
\keyword{optimize}
