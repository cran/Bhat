\name{logit.hessian}
\alias{logit.hessian}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Hessian (curvature matrix) }
\description{
  Numerical evaluation of the Hessian of a real function f: \eqn{R^n
    \rightarrow R}{R^n -> R}
  on a generalized logit scale, i.e. using transformed parameters
  according to x'=log((x-xl)/(xu-x))), with xl < x < xu.
}
\usage{
logit.hessian(x=x, f=f, del=rep(0.002, length(x$est)), dapprox=FALSE, nfcn=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a list with components 'label' (of mode character), 'est'
    (the parameter vector with the initial guess), 'low' (vector with
    lower bounds), and 'upp' (vector with upper bounds) }
  \item{f}{ the function for which the Hessian is to be computed at
    point x }
  \item{del}{ step size on logit scale (numeric) }
  \item{dapprox}{ logical variable. If TRUE the off-diagonal elements
    are set to zero. If FALSE (default) the full Hessian is computed }
  \item{nfcn}{ number of function calls }
}
\details{
  This version uses a symmetric grid for the numerical evaluation
  computation of first and second derivatives.
}
\value{
  returns list with
  \item{df }{ first derivatives (logit scale) }
  \item{ddf }{ Hessian (logit scale) }
  \item{nfcn }{ number of function calls }
  \item{eigen }{ eigen values (logit scale) }
}
\author{ E. Georg Luebeck (FHCRC) }
\note{ This function is part of the Bhat exploration tool }
\seealso{ \code{\link{dfp}}, \code{\link{newton}}, \code{\link{ftrf}}, \code{\link{btrf}}, \code{\link{dqstep}} } 
\examples{
  ## Rosenbrock Banana function
   fr <- function(x) {
         x1 <- x[1]
         x2 <- x[2]
         100 * (x2 - x1 * x1)^2 + (1 - x1)^2
    }
  ## define
   x <- list(label=c("a","b"),est=c(1,1),low=c(-100,-100),upp=c(100,100))
   logit.hessian(x,f=fr,del=dqstep(x,f=fr,sens=0.01))
  ## shows the differences in curvature at the minimum of the Banana
  ## function along principal axis (in a logit-transformed coordinate system)
}
\keyword{ array }% at least one, from doc/KEYWORDS
\keyword{ iteration }% __ONLY ONE__ keyword per line
\keyword{ methods }% __ONLY ONE__ keyword per line
\keyword{ optimize }% __ONLY ONE__ keyword per line


