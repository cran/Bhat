\name{dqstep}
\alias{dqstep}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ step size generator }
\description{
  \code{dqstep} determines the smallest steps ds from s so that
  abs(f(s+ds)-f(s)) equals a pre-specified sensitivity 
}
\usage{
dqstep(x, f, sens)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a list with components 'label' (of mode character), 'est'
    (the parameter vector with the initial guess), 'low' (vector with
    lower bounds), and 'upp' (vector with upper bounds) }
  \item{f}{ the function that is to be minimized over the parameter 
             vector defined by the list \code{x} }
  \item{sens}{ target sensitivity (i.e. the value of f(s+ds)-f(s)) }
}
\details{ 
  uses simple quadratic interpolation
}
\value{
  returns a vector with the desired step sizes
}
\author{ E. Georg Luebeck (FHCRC) }
\note{ This function is part of the Bhat exploration tool
}
% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{dfp}}, \code{\link{newton}}, \code{\link{logit.hessian}} }

\examples{
  ## Rosenbrock Banana function
   fr <- function(x) {
         x1 <- x[1]
         x2 <- x[2]
         100 * (x2 - x1 * x1)^2 + (1 - x1)^2
    }
  ## define
   x <- list(label=c("a","b"),est=c(1,1),low=c(0,0),upp=c(100,100))
   dqstep(x,fr,sens=1)
}
\keyword{ optimize }% at least one, from doc/KEYWORDS
\keyword{ iteration }% __ONLY ONE__ keyword per line
