##' Adaptive Multivariate MCMC sampler
##' 
##' This function generates MCMC-based samples from a (posterior) density f
##' (not necessarily normalized). It uses a Metropolis algorithm in conjunction
##' with a multivariate normal proposal distribution which is updated
##' adaptively by monitoring the correlations of succesive increments of at
##' least 2 pilot chains. The method is described in De Gunst, Dewanji and
##' Luebeck (submitted). The adaptive method is similar to the one proposed in
##' Gelfand and Sahu (JCGS 3:261--276, 1994).
##' 
##' standard output reports a summary of the acceptance fraction, the current
##' values of nlogf and the parameters for every (100*skip) th cycle. Plotted
##' chains show values only for every (skip) th cycle.
##' 
##' @param x a list with components 'label' (of mode character), 'est' (the
##' parameter vector with the initial guess), 'low' (vector with lower bounds),
##' and 'upp' (vector with upper bounds)
##' @param nlogf negative log of the density function (not necessarily
##' normalized) for which samples are to be obtained
##' @param m1 length of first pilot run (not used when covm supplied)
##' @param m2 length of second pilot run (not used when covm supplied )
##' @param m3 length of final run
##' @param scl1 scale for covariance of mv normal proposal (second pilot run)
##' @param scl2 scale for covariance of mv normal proposal (final run)
##' @param skip number of cycles skipped for graphical output
##' @param covm covariance matrix for multivariate normal proposal
##' distribution. If supplied, all pilot runs will be skipped and a run of
##' length m3 will be produced. Useful to continue a simulation from a given
##' point with specified covm
##' @param nfcn number of function calls
##' @param plot logical variable. If TRUE the chain and the negative log
##' density (nlogf) is plotted. The first m1+m2 cycles are shown in green,
##' other cycles in red
##' @param plot.range [Not documented. Leave as default]
##' @return list with the following components: \item{f }{ values of nlogf for
##' the samples obtained } \item{mcmc }{ the chain (samples obtained) }
##' \item{covm }{ current covariance matrix for mv normal proposal
##' distribution}
##' @note This function is part of the Bhat exploration tool
##' @author E. Georg Luebeck (FHCRC)
##' @seealso \code{\link{dfp}}, \code{\link{newton}},
##' \code{\link{logit.hessian}}
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
##' # start MCMC near mle
##'   x <- list(label=c("a","b","c"), est=c(20, 0.5, 0.05), low=c(0,0,0), upp=c(100,10,.1))
##' # samples from posterior density (~exp(-nlogf))) with non-informative
##' # (random uniform) priors for "a", "b" and "c".
##' out <- mymcmc(x, nlogf, m1=2000, m2=2000, m3=10000, scl1=0.5, scl2=2, skip=10, plot=TRUE)
##' # start MCMC from some other point: e.g. try x$est <- c(16,.2,.02)
##' 
##' @export
##' 
"mymcmc" <-
function (x, nlogf, m1, m2=m1, m3, scl1=0.5, scl2=2, skip=1, covm=0, nfcn = 0, plot=FALSE, plot.range=0)
{
            #     MCMC/MH sampler for R
            #     This module is part of the Bhat likelihood function exploration tool.

            #     This program is free software; you can redistribute it and/or modify
            #     it under the terms of the GNU General Public License as published by
            #     the Free Software Foundation; either version 2 of the License, or
            #     (at your option) any later version.
            #     This program is distributed in the hope that it will be useful,
            #     but WITHOUT ANY WARRANTY; without even the implied warranty of
            #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            #     GNU General Public License for more details.

            #     E. Georg Luebeck (gluebeck@fhcrc.org)
            #     last updated: 07-27-2015

            if (!is.list(x)) {
                    cat("x is not a list! see help file", "\n")
                    return()
            }
            names(x)[1] <- "label"
            names(x)[2] <- "est"
            names(x)[3] <- "low"
            names(x)[4] <- "upp"
            npar <- length(x$est)

            ####  objects:
            if (npar <= 0) stop('no. of parameters must be >= 1')

            small <- 1.e-8
            sens <- 1.
            xinf <- 20.

            # *** initialize graphical output
            if(plot==TRUE) {
              if(plot.range==0) {
                graphics::par(mfrow=c(npar+1,1),bg="grey")
              }
            }

            # *** initialize f.old
            cat(date(), "\n")
            xt <- ftrf(x$est, x$low, x$upp)
            f.old <- nlogf(x$est); nfcn <- nfcn + 1
            f.first <- f.old

            acc.count <- 0

            # *** if COVM is not provided
            if(!is.matrix(covm)) {

              # *** initialize monitors
              f.mon <- numeric(m1+m3); x.mon <- matrix(0,ncol=npar,nrow=m1+m3)

              cat('trying a proposal COVM using the Hessian','\n')
              del <- dqstep(x,nlogf,sens)
              h <- logit.hessian(x,nlogf,del,dapprox=FALSE,nfcn)  # returns list: $df,$ddf,$nfcn
              nfcn <- h$nfcn
              ddf <- h$ddf

              # ***   EIGEN VALUES and MATRIX INVERSION
              eig <- eigen(ddf)

              if(any(eig$values < small)) {
                  cat('Hessian may not be positive definite','\n')
                  cat('trying a proposal COVM using dqstep with sensitivity sens','\n')
                  del <- dqstep(x,nlogf,sens)
                  # eig <- .1/del/del; ddf <- diag(eig); eig <- eigen(ddf)
                  inv.ddf <- diag(del*del)
              } else {
                  # we've just eigendecomposed the matrix, let's use it for inverting
                  inv.ddf  <- 0.5 * eig$vectors %*% diag(1/eig$values, nrow=length(eig$values)) %*% t(eig$vectors)
              }

               # ggf <- 0.5 * solve(ddf,diag(1,npar),tol=1.e-10)

               # cat('unity test:','\n'); print(format(ggf %*% ddf),quote=FALSE)
              cat('first pilot chain:','\n','\n')

              nc <- 0
              for (n in 1:m1) {

                accept <- 1

                # *** compute proposal x' (=yt). If unacceptable, reject
                # dx <-  t(eig$vectors) %*% rnorm(npar,0,1/(0.5*sqrt(eig$values)))
                dx <-  MASS::mvrnorm(n = 1, mu = rep(0,npar), Sigma = inv.ddf, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
                yt <- xt + scl1 * dx
                Xt <- btrf(xt, x$low, x$upp)
                Yt <- btrf(yt, x$low, x$upp)

                # *** get log-likelihood
                f.new <- nlogf(Yt); nfcn <- nfcn + 1

                # *** Jacobian of logit-transform
                djac <- sum(log((Yt-x$low)*(x$upp-Yt)) - log((Xt-x$low)*(x$upp-Xt)))

                # *** boundary checks from within func ...

                # *** R ratio
                if(any(abs(yt) > xinf)) {
                  cat('cycle',n,': mymcmc close to boundary, move rejected!','\n'); accept <- 0} else {
                    accept <- min(accept,exp(-f.new+f.old+djac))
                }

                if(accept == 1) {xt <- yt; f.old <- f.new; acc.count <- acc.count+1
                             xt.maxl <- yt; f.maxl <- f.new # approximate/search max-likelihood
                             } else {
                               if(stats::runif(1) <= accept) {
                                         xt <- yt; f.old <- f.new; acc.count <- acc.count+1}
                             }

                # *** record to monitor
                f.mon[n] <- f.old; x.mon[n,] <- btrf(xt, x$low, x$upp)

                # *** regular output and graphical monitor
                if(n%%(100*skip) == 0) {
                  m.out <- c("n:",signif(n,3),"acceptance:",signif(acc.count/100/skip,3),"-log lkh:",signif(f.new,8),signif(x.mon[n,],6))
                  cat(m.out,'\n')
                  acc.count <- 0
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
                    nc <- nc + 1
                  }
                }
              }

            # *** obtain empirical covariance of increments
            covm <- stats::cov(x.mon[2:m1,]-x.mon[1:(m1-1),]) # now redundant
            cat('\n','\n')
            # print(covm,quote=FALSE)
            # cat('\n','\n')

            # *** pilot run 2 *** includes m2 cycles for incremental adjustment of covm
            cat('second pilot run:','\n','\n')
            eig <- eigen(covm)
            cat('eigen values:',eig$values,'\n','\n')

            } else {
            # covm for mvn proposal distribution given as input
            if(nrow(covm)!=length(x$est)) {stop('dimension of covm not specified correctly')}
            nc <- 0
            # re-initialize monitors
            m1 <- 1
            f.mon <- numeric(m1+m3); x.mon <- matrix(0,ncol=npar,nrow=m1+m3)
            # eig <- eigen(covm);
            x.mon[1,] <- x$est; f.mon <- f.first
          }


              # *** continue 'production' run on original scale
              acc.count <- 0
              for (n in (m1+1):(m1+m3)) {

                x.mon[n,] <- x.mon[n-1,]; f.mon[n] <- f.mon[n-1]
                accept <- 1

                # *** compute proposal x' (=yt). If unacceptable, reject
                # dx <-  t(eig$vectors) %*% rnorm(npar,0,sqrt(eig$values))
                dx <-  MASS::mvrnorm(n = 1, mu = rep(0,npar), Sigma = covm, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
                y <- x.mon[n-1,] + scl2 * dx
                if(any((y-x$low) < small)) {cat('move rejected (lower bound)','\n'); accept <- 0}
                if(any((x$upp-y) < small)) {cat('move rejected (upper bound)','\n'); accept <- 0}

                # *** boundary checks from within func ...

                if(accept > 0) {
                # *** get log-likelihood
                  f.new <- nlogf(y); nfcn <- nfcn + 1

                  # *** acceptance ratio and updates
                  accept <- min(accept,exp(-f.new+f.mon[n]))
                  if(accept == 1) {x.mon[n,] <- y; f.mon[n] <- f.new; acc.count <- acc.count+1
                             x.maxl <- y; f.maxl <- f.new # approximate/search max-likelihood
                             } else {
                               if(stats::runif(1) <= accept) {
                                 x.mon[n,] <- y; f.mon[n] <- f.new; acc.count <- acc.count+1
                               }
                             }
                }

                # *** regular output and graphical monitor
                if(n%%(100*skip) == 0) {
                  m.out <- c("n:",signif(n,3),"acceptance:",signif(acc.count/100/skip,3),"-log lkh:",signif(f.new,8),signif(x.mon[n,],6))
                  cat(m.out,'\n')
                  acc.count <- 0

                  n.skip.1 <- seq(skip,n,skip)
                  n.skip.2 <- seq(skip,min(n,m1+m2),skip)

                  if(plot==TRUE) {
                    if(m1 > 1) {brncol <- 3} else {brncol <- 2}

                    graphics::par(mar=c(0, 5, 0.1, 4) + 0.1)
                    plot(f.mon[n.skip.1], type='l', xlab = " ", ylab = "-log-density",col=2)
                    graphics::lines(f.mon[n.skip.2],col=brncol) #pilot cycles
                    for (i in 1:(npar-1)) {
                      graphics::par(mar=c(0, 5, 0, 4) + 0.1)
                      plot(x.mon[n.skip.1,i], type='l', xlab = " ", ylab = x$label[i], col=2)
                      graphics::lines(x.mon[n.skip.2,i],col=brncol) #pilot cycles
                    }
                    graphics::par(mar=c(0.1, 5, 0, 4) + 0.1)
                    plot(x.mon[n.skip.1,npar], type='l', xlab = " ", xaxt='n', ylab = x$label[npar], col=2)
                    graphics::lines(x.mon[n.skip.2,npar],col=brncol) #pilot cycles

                    nc <- nc + 1
                  }


                  if(m1 > 1) { #note: when covm is passed to mymcmc, m1 is set to 1
                  # update covariance using sampled increments (m1+1):n
                    if(n <= m1+m2) {
                      covm <- stats::cov(x.mon[2:n,]-x.mon[1:(n-1),])
                      eig <- eigen(covm)
                      if(any(eig$values < small)) {warning('covm nearly singular')}
                    }
                  }
                }
            }
            return(list(f=f.mon,mcmc=x.mon,covm=covm))
            }
