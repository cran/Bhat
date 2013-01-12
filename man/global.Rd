\name{global}
\alias{global}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Random search for a global function minimum}
\description{ This function generates MCMC samples from a (posterior)
  density function f (not necessarily normalized) in search of a global minimum of f. 
  It uses a simple Metropolis algorithm to generate the samples. 
  Global monitors the mcmc samples and returns the minimum 
  value of f, as well as a sample covariance (covm) that can be used as input
  for the Bhat function mymcmc.  
}
\usage{
global(x, nlogf, beta = 1., mc, scl=2, skip=1, nfcn = 0, plot=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a list with components 'label' (of mode character), 'est'
    (the parameter vector with the initial guess), 'low' (vector with
    lower bounds), and 'upp' (vector with upper bounds) }
  \item{nlogf}{ negative log of the density function (not necessarily
    normalized) }
  \item{beta}{ 'inverse temperature' parameter } 
  \item{mc}{ length of MCMC search run }
  \item{scl}{ not used } 
  \item{skip}{ number of cycles skipped for graphical output }
  \item{nfcn}{ number of function calls }
  \item{plot}{ logical variable. If TRUE the chain and the negative
    log density (nlogf) is plotted }
}
\details{ 
  standard output reports a summary of the acceptance fraction, the
  current values of nlogf and the parameters for every (100*skip) th
  cycle. Plotted chains show values only for every (skip) th cycle.
}
\value{
  list with the following components:
  \item{fmin }{ minimum value of nlogf for the samples obtained }
  \item{xmin }{ parameter values at fmin }
  \item{covm }{ covariance matrix of differences between consecutive samples in chain } 
}
\references{ too numerous to be listed here }
\author{ E. Georg Luebeck (FHCRC) }
\note{ This function is part of the Bhat package
}
%  ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{dfp}}, \code{\link{newton}}, \code{\link{logit.hessian}} \code{\link{mymcmc}} }

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

# initialize global search
  x <- list(label=c("a","b","c"), est=c(10, 0.25, 0.05), low=c(0,0,0), upp=c(100,10,.1))
# samples from posterior density (~exp(-nlogf))) with non-informative
# (random uniform) priors for "a", "b" and "c".
out <- global(x, nlogf, beta = 1., mc=1000, scl=2, skip=1, nfcn = 0, plot=TRUE)
# start MCMC from some other point: e.g. try x$est <- c(16,.2,.02)
}
\keyword{ iteration }% __ONLY ONE__ keyword per line
\keyword{ methods }% __ONLY ONE__ keyword per line
\keyword{ optimize }% __ONLY ONE__ keyword per line

