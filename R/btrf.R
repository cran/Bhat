##' Generalized inverse-logit transform
##' 
##' maps real line onto open interval (xl, xu) using the transform y = (exp(xt)
##' * xu + xl)/(1.+exp(xt)) where xt is a numeric vector with -Inf < xt < Inf
##' 
##' 
##' @param xt a numeric vector
##' @param xl a numeric vector of same length as x
##' @param xu a numeric vector of same length as x, and xu > xl
##' @return returns the inverse-logit transform (numeric) of xt
##' @author E. Georg Luebeck (FHCRC)
##' @seealso \code{\link{ftrf}}
##' @keywords optimize misc
##' @export
##' 
"btrf" <-
function(xt,xl,xu) {

####  back transformation
####  this assumes logit transformations of the parameters
####  bounded from below by xl and from above by xu

  rho <- exp(xt)
  return((rho * xu + xl)/(1.+rho))
}
