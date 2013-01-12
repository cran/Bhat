\name{dfp} 
\alias{dfp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Function minimization with box-constraints
}
\description{
  This Davidon-Fletcher-Powell optimization algorithm has been `hand-tuned'
  for minimal setup configuration and for efficency. It uses an internal
  logit-type transformation based on the pre-specified box-constraints.
  Therefore, it usually does not require rescaling (see help for the R
  optim function). \code{dfp} automatically computes step sizes for each
  parameter to operate with sufficient sensitivity in the functional
  output. Performance is comparable to the BFGS algorithm in the R
  function \code{optim}. \code{dfp} interfaces with \code{newton} to
  ascertain convergence, compute the eigenvalues of the Hessian, and 
  provide 95\% confidence intervals when the function to be minimized is a
  negative log-likelihood.
 
  The dfp function minimizes a function \code{f} over the parameters
  specified in the input list \code{x}. The algorithm is based on Fletcher's
  "Switching Method" (Comp.J. 13,317 (1970))
}
\usage{
dfp(x, f, tol=1e-05, nfcn=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a list with components 'label' (of mode character), 'est'
    (the parameter vector with the initial guess), 'low' (vector with
    lower bounds), and 'upp' (vector with upper bounds) }
  \item{f}{ the function that is to be minimized over the parameter 
    vector defined by the list \code{x} }
  \item{tol}{ a tolerance used to determine when convergence
	     should be indicated }
  \item{nfcn}{ number of function calls }
}
\details{
  the code has been 'transcribed' from Fortran source code into R
}
\value{
  list with the following components:
  \item{fmin }{ the function value f at the minimum }
  \item{label }{ the labels taken from list \code{x} }
  \item{est }{ a vector of the estimates at the minimum. dfp
    does not overwrite \code{x} }
  \item{status }{ 0 indicates convergence, 1 indicates non-convergence }
  \item{nfcn }{ no. of function calls }
}
\references{ Fletcher's Switching Method (Comp.J. 13,317, 1970)
}
\author{ E. Georg Luebeck (FHCRC)
}
\note{ This function is part of the Bhat exploration tool
}
% ~Make other sections like WARNING with \section{WARNING }{....} ~
\seealso{ optim, \code{\link{newton}}, \code{\link{ftrf}}, \code{\link{btrf}}, \code{\link{logit.hessian}}
}
\examples{
        # generate some Poisson counts on the fly
          dose <- c(rep(0,50),rep(1,50),rep(5,50),rep(10,50))
          data <- cbind(dose,rpois(200,20*(1+dose*.5*(1-dose*0.05))))

        # neg. log-likelihood of Poisson model with 'linear-quadratic' mean: 
          lkh <- function (x) { 
          ds <- data[, 1]
          y  <- data[, 2]
          g <- x[1] * (1 + ds * x[2] * (1 - x[3] * ds)) 
          return(sum(g - y * log(g)))
          }

	# for example define
          x <- list(label=c("a","b","c"),est=c(10.,10.,.01),low=c(0,0,0),upp=c(100,20,.1))

	# call:
	  results <- dfp(x,f=lkh)
}
\keyword{ optimize }% at least one, from doc/KEYWORDS
\keyword{ methods }% __ONLY ONE__ keyword per line

