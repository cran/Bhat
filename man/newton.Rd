\name{newton}
\alias{newton}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Function minimization with box-constraints }
\description{
  Newton-Raphson algorithm for minimizing a function \code{f} over the
  parameters specified in the input list \code{x}. Note, a
  Newton-Raphson search is very efficient in the 'quadratic region'
  near the optimum. In higher dimensions it tends to be rather
  unstable and may behave chaotically. Therefore, a (local or global)
  minimum should be available to begin with. Use the \code{optim} or
  \code{dfp} functions to search for optima.
}
\usage{
newton(x, f, eps=0.1, itmax=10, relax=0, nfcn=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a list with components 'label' (of mode character), 'est'
    (the parameter vector with the initial guess), 'low' (vector with
    lower bounds), and 'upp' (vector with upper bounds) }
  \item{f}{ the function that is to be minimized over the parameter 
             vector defined by the list \code{x} }
  \item{eps}{ converges when all (logit-transformed) derivatives are
    smaller \code{eps} }  
  \item{itmax}{ maximum number of Newton-Raphson iterations }
  \item{relax}{ numeric. If 0, take full Newton step, otherwise
    'relax' step incrementally until a better value is found }
  \item{nfcn}{ number of function calls }
}
\value{
  list with the following components:
  \item{fmin }{ the function value f at the minimum }
  \item{label }{ the labels }
  \item{est }{ a vector of the parameter estimates at the minimum. newton
    does not overwrite \code{x} }
  \item{low }{ lower 95\% (Wald) confidence bound }
  \item{upp }{ upper 95\% (Wald) confidence bound }
  The confidence bounds assume that the function \code{f} is a negative
  log-likelihood 
}
\author{ E. Georg Luebeck (FHCRC) }
\note{
  \code{newton} computes the (logit-transformed) Hessian of \code{f}
  (using logit.hessian). This function is part of the Bhat exploration
  tool
}
%  ~Make other sections like WARNING with \section{WARNING }{....} ~
\seealso{ \code{\link{dfp}}, \code{\link{ftrf}}, \code{\link{btrf}}, \code{\link{logit.hessian}}, \code{\link{plkhci}} }
\examples{
        # generate some Poisson counts on the fly
          dose <- c(rep(0,100),rep(1,100),rep(5,100),rep(10,100))
          data <- cbind(dose,rpois(400,20*(1+dose*.5*(1-dose*0.05))))

        # neg. log-likelihood of Poisson model with 'linear-quadratic' mean: 
          lkh <- function (x) { 
          ds <- data[, 1]
          y  <- data[, 2]
          g <- x[1] * (1 + ds * x[2] * (1 - x[3] * ds)) 
          return(sum(g - y * log(g)))
          }

	# for example define
          x <- list(label=c("a","b","c"),est=c(10.,10.,.01),low=c(0,0,0),upp=c(100,20,.1))

	# calls:
	  r <- dfp(x,f=lkh)
          x$est <- r$est
          results <- newton(x,lkh)
}
\keyword{ optimize }% at least one, from doc/KEYWORDS
\keyword{ methods }% __ONLY ONE__ keyword per line



