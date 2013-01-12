\name{ftrf}
\alias{ftrf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generalized logit transform }
\description{
  maps a bounded parameter x onto the real line according to
  y=log((x-xl)/(xu-x))), with xl < x < xu. If this constraint is
  violated, an error occurs. x may be vector  
}
\usage{
ftrf(x, xl, xu)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector }
  \item{xl}{ a numeric vector of same length as x with x > xl }
  \item{xu}{ a numeric vector of same length as x with x < xu }
}
\value{
  returns numerical vector of transforms
}
\author{ E. Georg Luebeck (FHCRC) }
\seealso{ \code{\link{btrf}} }
\keyword{ optimize }% at least one, from doc/KEYWORDS
\keyword{ misc }% __ONLY ONE__ keyword per line
