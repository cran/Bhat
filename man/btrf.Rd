\name{btrf}
\alias{btrf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generalized inverse-logit transform }
\description{ maps real line onto open interval (xl, xu) using the
  transform y = (exp(xt) * xu + xl)/(1.+exp(xt)) where xt is a numeric
  vector with -Inf < xt < Inf
}
\usage{
btrf(xt, xl, xu)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{xt}{ a numeric vector }
  \item{xl}{ a numeric vector of same length as x}
  \item{xu}{ a numeric vector of same length as x, and xu > xl}
}
\value{
  returns the inverse-logit transform (numeric) of xt
}
\author{ E. Georg Luebeck (FHCRC) }
\seealso{ \code{\link{ftrf}} }
\keyword{ optimize }% at least one, from doc/KEYWORDS
\keyword{ misc }% __ONLY ONE__ keyword per line
