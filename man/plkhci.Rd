\name{plkhci}
\alias{plkhci}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Profile-likelihood based confidence intervals }
\description{
  function to find \code{prob}*100\% confidence intervals using
  profile-likelihood. Numerical solutions are obtained via a  
  modified  Newton-Raphson algorithm. The method is described in
  Venzon and Moolgavkar, Journal of the Royal Statistical Society,
  Series C vol 37, no.1, 1988, pp. 87-94. 
}
\usage{
plkhci(x, nlogf, label, prob=0.95, eps=.001, nmax=10, nfcn=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a list with components 'label' (of mode character), 'est'
    (the parameter vector with the initial guess), 'low' (vector with
    lower bounds), and 'upp' (vector with upper bounds) }
  \item{nlogf}{ the negative log of the density function (not necessarily
    normalized) for which samples are to be obtained }
  \item{label}{ parameter for which confidence bounds are computed }
  \item{prob}{ probability associated with the confidence interval }
  \item{eps}{ a numerical value. Convergence results when all
    (logit-transformed) derivatives are smaller \code{eps} }  
  \item{nmax}{ maximum number of Newton-Raphson iterations in each direction }
  \item{nfcn}{ number of function calls }
}
\value{
  2 component vector giving lower and upper p\% confidence bounds 
}
\author{ E. Georg Luebeck (FHCRC) }
\note{ At this point, only a single parameter label can be passed to
  plkhci.  This function is part of the Bhat exploration tool}

%  ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{dfp}}, \code{\link{newton}}, \code{\link{logit.hessian}} }

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

	# for example define
          x <- list(label=c("a","b","c"),est=c(10.,10.,.01),low=c(0,0,0),upp=c(100,20,.1))

	# get MLEs using dfp:
	  r <- dfp(x,f=nlogf)
          x$est <- r$est
          plkhci(x,nlogf,"a")
          plkhci(x,nlogf,"b")
          plkhci(x,nlogf,"c")
	# e.g. 90\% confidence bounds for "c" 
          plkhci(x,nlogf,"c",prob=0.9)
        
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line




