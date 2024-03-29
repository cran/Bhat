##' step size generator
##' 
##' \code{dqstep} determines the smallest steps ds from s so that
##' abs(f(s+ds)-f(s)) equals a pre-specified sensitivity
##' 
##' uses simple quadratic interpolation
##' 
##' @param x a list with components 'label' (of mode character), 'est' (the
##' parameter vector with the initial guess), 'low' (vector with lower bounds),
##' and 'upp' (vector with upper bounds)
##' @param f the function that is to be minimized over the parameter vector
##' defined by the list \code{x}
##' @param sens target sensitivity (i.e. the value of f(s+ds)-f(s))
##' @return returns a vector with the desired step sizes
##' @note This function is part of the Bhat exploration tool
##' @author E. Georg Luebeck (FHCRC)
##' @seealso \code{\link{dfp}}, \code{\link{newton}},
##' \code{\link{logit.hessian}}
##' @keywords optimize iteration
##' @examples
##' 
##'   ## Rosenbrock Banana function
##'    fr <- function(x) {
##'          x1 <- x[1]
##'          x2 <- x[2]
##'          100 * (x2 - x1 * x1)^2 + (1 - x1)^2
##'     }
##'   ## define
##'    x <- list(label=c("a","b"),est=c(1,1),low=c(0,0),upp=c(100,100))
##'    dqstep(x,fr,sens=1)
##' 
##' @export
##' 
"dqstep" <-
function(x,f,sens) {
  # fix nfcn counter later
  npar  <- length(x$est)
  step <- .001; dsteps <- rep(step,npar)

  # OBTAIN REFERENCE POINT VALUE
  xt <- ftrf(x$est, x$low, x$upp)
  f0 <- f(x$est)
  
  for(i in 1:npar) {
    stepi <- step
    flag <- 0
    xt.new <- xt
    while(flag==0) {
      flag <- 1
      xt.new[i] <- xt[i]-stepi; x1 <- -stepi
      f1 <- f(btrf(xt.new, x$low, x$upp))
      xt.new[i] <- xt[i]+stepi; x2 <- stepi
      f2 <- f(btrf(xt.new, x$low, x$upp))

      # handle exceptions
      if (is.na(f2)) f2 <- Inf
      if(f2==Inf | f2==-Inf) {
        warning('Infs - reducing step size')
        stepi <- stepi/10; flag <- 0}
      if(f2==f0 & f1==f0) {
        cat('increasing step size','\n')
        stepi <- stepi*10; flag <- 0}
      if(abs(f2-f0) > (.5 * sens) | abs(f1-f0) > (.5 * sens)) {
        cat('reducing step size','\n')
        stepi <- stepi/10; flag <- 0
        # cat(stepi,f1,f2,'\n')
      }
    }
    
    b <- ((f1-f0)*x2-(f2-f0)*x1)/(x1*x1*x2-x2*x2*x1)
    a <- (f1-f0)/x1-b*x1

    # ***roots
    r <- a*a+4*b*sens
    if(r < 0 | is.na(r) | b ==0) {
      warning('oops: unable to find stepsize, use default')
      cat('problem with ',x$label[i],'\n')
      break
    }
      
    xs1 <- 0.5*(-a-sqrt(a*a+4*b*sens))/b
    xs2 <- 0.5*(-a+sqrt(a*a+4*b*sens))/b

    if(abs(xs1) <= abs(xs2)) {xs <- xs1} else { xs <- xs2}

    # *** see where we end up
    xt.new[i] <- xt[i]+xs
    f2 <- f(btrf(xt.new, x$low, x$upp))
    if (is.na(f2)) f2 <- Inf
    if(f2==Inf | f2==-Inf) {
      warning('oops: unable to find stepsize, use default')
      break
    }
    dsteps[i] <- xs
    # cat('DSTEP:',x$label[i],signif(xs,6),(f2-f0),'\n')
  }
  return(dsteps)
}

