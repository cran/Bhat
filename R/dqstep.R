"dqstep" <-
function(x,f,sens) {
  # fix nfcn counter later
  npar  <- length(x$est)
  step <- .0001; dsteps <- rep(.001,npar)

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

      if(f2==Inf | f2==-Inf) {
        warning('Infs - reducing step size')
        stepi <- stepi/10; flag <- 0}
      if(abs(f2-f0) > (.5 * sens) | abs(f1-f0) > (.5 * sens)) {
        # warning('reducing step size')
        stepi <- stepi/10; flag <- 0
        # cat(stepi,f1,f2,'\n')
      }
    }
    
    b <- ((f1-f0)*x2-(f2-f0)*x1)/(x1*x1*x2-x2*x2*x1)
    a <- (f1-f0)/x1-b*x1

    # ***roots
    xs1 <- 0.5*(-a-sqrt(a*a+4*b*sens))/b
    xs2 <- 0.5*(-a+sqrt(a*a+4*b*sens))/b

    if(abs(xs1) <= abs(xs2)) {xs <- xs1} else { xs <- xs2}

    # *** see where we end up
    xt.new[i] <- xt[i]+xs
    f2 <- f(btrf(xt.new, x$low, x$upp))
    if(f2==Inf | f2==-Inf) {
      warning('oops: unable to find stepsize, use default')
      break
    }
    dsteps[i] <- xs
    # cat('DSTEP:',x$label[i],signif(xs,6),(f2-f0),'\n')
  }
  return(dsteps)
}

