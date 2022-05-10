##' Generalized logit transform
##' 
##' maps a bounded parameter x onto the real line according to
##' y=log((x-xl)/(xu-x))), with xl < x < xu. If this constraint is violated, an
##' error occurs. x may be vector
##' 
##' 
##' @param x a numeric vector
##' @param xl a numeric vector of same length as x with x > xl
##' @param xu a numeric vector of same length as x with x < xu
##' @return returns numerical vector of transforms
##' @author E. Georg Luebeck (FHCRC)
##' @seealso \code{\link{btrf}}
##' @keywords optimize misc
##' @export
##' 
"ftrf" <-
function(x,xl,xu) {

####  forward transformation
####  this assumes logit transformations of the parameters
####  bounded from below by xl and from above by xu
  if(any((x-xl) <= 0)) stop('ftrf requires x > xl')
  if(any((xu-x) <= 0)) stop('ftrf requires x < xu')
  return(log((x-xl)/(xu-x)))
}
