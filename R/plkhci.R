#       Program:     plkhci.R
#                    modified newton-raphson optimizer
#                    to obtain profile-likelihood-based confidence intervals
#
#       Method:      D.J. VENZON and S.H. MOOLGAVKAR
#                    " A method for computing profile-likelihood-based
#                      confidence intervals "
#                    Journal of the Royal Statistical Society, Series C
#                       volume 37, no.1, 1988, pp. 87-94
#                    not yet fully implemented, but may converge
#                    
#       Notes:       p:  identifies which parameter (among npar) is to be used for profile
#                    likelihood calc. 
#                    npar:	total number of parameters in model

plkhci <- function(x, nlogf, label, prob=0.95, eps=.001, nmax=10, nfcn = 0) {
  
  if (!is.list(x)) {
  	    cat("x is not a list! see help file", "\n")
  	    return()
  }
  names(x)[1] <- "label"
  names(x)[2] <- "est"
  names(x)[3] <- "low"
  names(x)[4] <- "upp"
  npar <- length(x$est)

  if (npar <= 0) {
  	    warning("no. of parameters < 1")
  	    exit
  }

  small <- 1.e-8
  nd1 <- npar-1    
  sens <- .5
  rhs <- numeric(npar)
  disc <- numeric(npar)
  qchi <- qchisq(prob,1)

### check selected parameter
  if(!is.character(label)) stop('label needs to be of type character')
  p <- match(label,x$label)
  if(is.na(p)) stop('label not found')
  
  xt <- ftrf(x$est, x$low, x$upp)
  fmle <- nlogf(x$est); nfcn <- nfcn + 1

  cat(' neg. log. likelihood: ',fmle,'\n','\n')
  cat(' will atempt to compute both bounds (+/- direction)','\n')

        del <- dqstep(x,nlogf,sens)

        # *** shoot along tangent?
        # *** initialize first step

        h <- logit.hessian(x,nlogf,del,dapprox=FALSE,nfcn)  # returns list: $df,$ddf,$nfcn
        nfcn <- h$nfcn
        ddf <- h$ddf; df <- h$df

        # ****  construct ddl
	ddl <- ddf[-p,-p]

        # ****  construct dbo
        dbo <- ddf[p,-p]

        # ****  construct dbb
	dbb0 <- ddf[p,p]

        # ***   matrix inversion: ddl inverse
        b <- diag(1,nd1)
        # *** if inversion fails need to return or reinitialize NR (not implemented)
        ggl <- solve(ddl,b,tol=1.e-16)
        dbl0 <- ggl %*% dbo
        tmp0 <- t(dbo) %*% dbl0 
        cat('\n')
  
      for(idir in c(1,-1)) {

        xt0 <- xt
        xt0[p] <- xt0[p]  +idir *       sqrt(qchi/(dbb0-tmp0))/2
        xt0[-p] <- xt0[-p]-idir * as.vector(dbl0) * sqrt(qchi/(dbb0-tmp0))/2

        f0 <- nlogf(btrf(xt0, x$low, x$upp)); nfcn <- nfcn + 1
        
        if(idir ==1)  cat('trying lower bound ------------------------', '\n')
        if(idir ==-1) cat('trying upper bound ------------------------', '\n')
        cat('starting at:   ',2*(f0-fmle),'\n')
        cat('initial guess: ',btrf(xt0, x$low, x$upp),'\n')

# *** search via newton-raphson

      cat('\n')
      cat('begin Newton-Raphson search for profile lkh conf. bounds:','\n')
      cat('eps value for stop criterium:',eps,'\n')
      cat('nmax                        :',nmax,'\n')  

        for (n in 1:nmax) {

        x0 <- x; x0$est <- btrf(xt0, x$low, x$upp)
        f0 <- nlogf(x0$est); nfcn <- nfcn + 1
        h <- logit.hessian(x0,nlogf,del,dapprox=FALSE,nfcn)  # returns list: $df,$ddf,$nfcn
        nfcn <- h$nfcn
        ddf <- h$ddf; df <- h$df

        # ****  construct ddl
	ddl <- ddf[-p,-p]

        # ****  construct dbo
        dbo <- ddf[p,-p]

        # ***  construct matrix fdd
        fdd <- ddf 

        fdd[1,1] <- -df[p]
        rhs[1]  <- (fmle - f0 + 0.5*qchi)
        fdd[1,2:npar] <- -df[-p]
        rhs[2:npar] <- df[-p]
        fdd[-1,-1] <- ddl #lower right block   
        fdd[2:npar,1] <- dbo
      
        b <- diag(1,npar)
        ggl <- solve(fdd,b,tol=1.e-16)

        # *** shoot ahead

        if(n <=2) rneps <- .2
        if(n>2 && n <=4) rneps <- .8
        if(n > 4) rneps <- 1.
        
        nepsc <- 0
        xt00 <- xt0
      
        disc0 <- disc <- as.vector(ggl %*% rhs)

        rhs1 <- 1000.
        #while(nepsc < 8 &&  rhs1 > qchi) {
        disc <- rneps * disc0

        nz <- 0
        xt0[p]  <- xt00[p]+disc[1]
        xt0[-p] <- xt00[-p]+disc[2:npar]
                tmp <- numeric(npar)
                tmp[abs(disc) <= eps] <- 1
                nz <- sum(tmp)

        # *** update f0
        f0 <- nlogf(btrf(xt0, x$low, x$upp)); nfcn <- nfcn + 1
        rhs1 <- abs(fmle-f0+0.5*qchi)

        #if(rhs1 > qchi) {
        #   rneps <- .75 * rneps
        #   cat('reducing step size:',rneps,rhs1,'\n') 
        #   nepsc <- nepsc + 1
        # }

        #if(nepsc ==8) cat('problem closing in, will continue.','\n')

        estimate <- btrf(xt0, x$low, x$upp)
        # ***   stop-criterion for iteration
        if(nz ==  npar && nepsc  ==  0) {
        cat('\n','CONVERGENCE: ',trunc(n),' iterations','\n')
        status <- 'converged'
        cat('\n')
        cat('chisquare value is: ',2*(f0-fmle),'\n')
        cat('confidence bound of ',x$label[p],' is ',estimate[p],'\n')
        cat('log derivatives:    ',df[-p],'\n')
        m.out <- cbind(x$label, signif(estimate,6), signif(df,6), signif(diag(ddf),6))
        dimnames(m.out) <- list(1:npar, c("label", "estimate", "log deriv", "log curv"))
        print(m.out, quote = FALSE)
        cat('\n');break
      }
        # cat('\n')
        # cat('chisquare value is: ',2*(f0-fmle),'\n')
        # cat('confidence bound of ',x$label[p],' is ',estimate[p],'\n')
        # cat('log derivatives:    ',df[-p],'\n')
      }
      if(idir==-1) {xlow <- estimate[p]} else {xupp <- estimate[p]}
      }
        return(c(xlow,xupp))
}








