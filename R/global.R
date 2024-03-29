##' Random search for a global function minimum
##' 
##' This function generates MCMC samples from a (posterior) density function f
##' (not necessarily normalized) in search of a global minimum of f.  It uses a
##' simple Metropolis algorithm to generate the samples.  Global monitors the
##' mcmc samples and returns the minimum value of f, as well as a sample
##' covariance (covm) that can be used as input for the Bhat function mymcmc.
##' 
##' standard output reports a summary of the acceptance fraction, the current
##' values of nlogf and the parameters for every (100*skip) th cycle. Plotted
##' chains show values only for every (skip) th cycle.
##' 
##' @param x a list with components 'label' (of mode character), 'est' (the
##' parameter vector with the initial guess), 'low' (vector with lower bounds),
##' and 'upp' (vector with upper bounds)
##' @param nlogf negative log of the density function (not necessarily
##' normalized)
##' @param beta 'inverse temperature' parameter
##' @param mc length of MCMC search run
##' @param scl not used
##' @param skip number of cycles skipped for graphical output
##' @param nfcn number of function calls
##' @param plot logical variable. If TRUE the chain and the negative log
##' density (nlogf) is plotted
##' @return list with the following components: \item{fmin }{ minimum value of
##' nlogf for the samples obtained } \item{xmin }{ parameter values at fmin }
##' \item{covm }{ covariance matrix of differences between consecutive samples
##' in chain }
##' @note This function is part of the Bhat package
##' @author E. Georg Luebeck (FHCRC)
##' @seealso \code{\link{dfp}}, \code{\link{newton}},
##' \code{\link{logit.hessian}} \code{\link{mymcmc}}
##' @references too numerous to be listed here
##' @keywords iteration methods optimize
##' @examples
##' 
##' # generate some Poisson counts on the fly
##'   dose <- c(rep(0,50),rep(1,50),rep(5,50),rep(10,50))
##'   data <- cbind(dose,rpois(200,20*(1+dose*.5*(1-dose*0.05))))
##' 
##' # neg. log-likelihood of Poisson model with 'linear-quadratic' mean: 
##'   nlogf <- function (x) { 
##'   ds <- data[, 1]
##'   y  <- data[, 2]
##'   g <- x[1] * (1 + ds * x[2] * (1 - x[3] * ds)) 
##'   return(sum(g - y * log(g)))
##'   }
##' 
##' # initialize global search
##'   x <- list(label=c("a","b","c"), est=c(10, 0.25, 0.05), low=c(0,0,0), upp=c(100,10,.1))
##' # samples from posterior density (~exp(-nlogf))) with non-informative
##' # (random uniform) priors for "a", "b" and "c".
##' out <- global(x, nlogf, beta = 1., mc=1000, scl=2, skip=1, nfcn = 0, plot=TRUE)
##' # start MCMC from some other point: e.g. try x$est <- c(16,.2,.02)
##' 
##' @export
##' 
"global" <-
function (x, nlogf, beta = 1., mc=1000, scl=2, skip=1, nfcn = 0, plot=FALSE)
{
	    #     simple MCMC/MH sampler for R
	    #     This module is part of the Bhat likelihood exploration tool.

	    #     This program is free software; you can redistribute it and/or modify
	    #     it under the terms of the GNU General Public License as published by
	    #     the Free Software Foundation; either version 2 of the License, or
	    #     (at your option) any later version.
	    #     This program is distributed in the hope that it will be useful,
	    #     but WITHOUT ANY WARRANTY; without even the implied warranty of
	    #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	    #     GNU General Public License for more details.

            #     E. Georg Luebeck (gluebeck@fhcrc.org)

	    if (!is.list(x)) {
	    	    cat("x is not a list! see help file", "\n")
	    	    return()
	    }
	    names(x)[1] <- "label"
	    names(x)[2] <- "est"
	    names(x)[3] <- "low"
	    names(x)[4] <- "upp"
	    npar <- length(x$est)
            x0 <- x
            
	    ####  objects:
	    if (npar <= 0) stop('no. of parameters must be >= 1')
            eps <- .05

            # *** initialize graphical output
            if(plot==TRUE) {
            graphics::par(mfrow=c(npar+1,1),bg="grey")}

            # *** initialize f.old
	    cat(date(), "\n")
	    fp <- nlogf(x$est); nfcn <- nfcn + 1

            # initialize monitors
            f.mon <- numeric(mc); x.mon <- matrix(0,ncol=npar,nrow=mc)
            f.mon[1] <- fp; x.mon[1,] <- x$est

            # *** scale dsc: used in defining proposal
            dsc <- rep(.05,npar)
            fmin <- fp

            # --- start M chain ------------------------------------------------------------

            sum1 <- sum2 <- sum3 <- n.accept2 <- sm <- c.sum <- 0
            n.tried <- n.accepted <- numeric(npar)
            
            acc <- 1+numeric(npar) #??
            for(n in 1:mc) {                            
              # ** parameters in turn
              for(i in 1:npar) {
              acc.p1 <- acc.p2 <- 0
              # *** create suitable random uniform proposal kernel (width)
              di <- dsc[i] * min(abs(x$est[i]-x$low[i]),abs(x$est[i]-x$upp[i]))
              dih <- di/2
              x0$est[i] = x$est[i] - dih + di*stats::runif(1)

              if(x0$est[i] > x$low[i] & x0$est[i] < x$upp[i]) {
                f0 <- nlogf(x0$est); nfcn <- nfcn + 1
                # IF(FP.EQ.nan .or. FP.EQ.-nan) GOTO 4
                # IF(FP.EQ.inf .or. FP.EQ.-inf) GOTO 4
                acc[i] = exp(beta*(fp-f0))
                acc.p1 <- max(eps,min(1,acc[i])) # max(eps,...)  avoids freezing at initialization
                n.tried[i] <- n.tried[i]+1
              }

              if(stats::runif(1) <= acc.p1) {
                x$est[i] <- x0$est[i]; fp <- f0; n.accept2 <- n.accept2+1
                n.accepted[i] <- n.accepted[i]+1
              }
              if(fmin>f0) {fmin <- f0; xmin <- x0$est}
            }
              
              f.mon[n] <- fp; x.mon[n,] <- x$est

              # *** adjust di to obtain acceptance near 50%

              if(n%%100==0) {
                # *** update proposal width array DI every 100 updates for acceptance near 50 %
                cat("Parameter","no. tried","no. accepted","ratio",'\n')
                freq <- n.accepted/n.tried
                for (i in 1:npar) {
                  cat(i,' ',n.tried[i],' ',n.accepted[i],' ',freq[i],'\n')
                  dsc[i] <- max(1.e-10,dsc[i]*(freq[i]+.5)^2)
                }

                n.accepted <- n.tried <- numeric(npar); n.accept2 <- 0
                n.skip <- seq(skip,n,skip)

                if(plot==TRUE) {
                graphics::par(mar=c(0, 5, 0.1, 4) + 0.1)
                plot(f.mon[n.skip], type='l', xlab = " ", ylab = "-log-density",col=2)
                for (i in 1:(npar-1)) {
                  graphics::par(mar=c(0, 5, 0, 4) + 0.1)
                  plot(x.mon[n.skip,i], type='l', xlab = " ", ylab = x$label[i], col=2)
                }
                graphics::par(mar=c(0.1, 5, 0, 4) + 0.1)
                plot(x.mon[n.skip,npar], type='l', xlab = " ", xaxt='n', ylab = x$label[npar], col=2)
              }
              }
            }

            covm <- stats::cov(x.mon[2:mc,]-x.mon[1:(mc-1),]) 
            cat('\n','\n')
            # print(covm,quote=FALSE)
            # cat('\n','\n')

            # *** pilot run 2 *** includes m2 cycles for incremental adjustment of covm
            cat('second pilot run:','\n','\n')
            eig <- eigen(covm)
            cat('eigen values:',eig$values,'\n','\n')

            return(list(fmin=fmin,xmin=xmin,covm=covm))
            }
