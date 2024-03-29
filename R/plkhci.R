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

#                    modified 04/26/2004 by Russell Millar Dept of Stats U. Akld
#                    to handle case of npar=1



##' Profile-likelihood based confidence intervals
##' 
##' function to find \code{prob}*100\% confidence intervals using
##' profile-likelihood. Numerical solutions are obtained via a modified
##' Newton-Raphson algorithm. The method is described in Venzon and Moolgavkar,
##' Journal of the Royal Statistical Society, Series C vol 37, no.1, 1988, pp.
##' 87-94.
##' 
##' 
##' @param x a list with components 'label' (of mode character), 'est' (the
##' parameter vector with the initial guess), 'low' (vector with lower bounds),
##' and 'upp' (vector with upper bounds)
##' @param nlogf the negative log of the density function (not necessarily
##' normalized) for which samples are to be obtained
##' @param label parameter for which confidence bounds are computed
##' @param prob probability associated with the confidence interval
##' @param eps a numerical value. Convergence results when all
##' (logit-transformed) derivatives are smaller \code{eps}
##' @param nmax maximum number of Newton-Raphson iterations in each direction
##' @param nfcn number of function calls
##' @return 2 component vector giving lower and upper p\% confidence bounds
##' @note At this point, only a single parameter label can be passed to plkhci.
##' This function is part of the Bhat exploration tool
##' @author E. Georg Luebeck (FHCRC)
##' @seealso \code{\link{dfp}}, \code{\link{newton}},
##' \code{\link{logit.hessian}}
##' @keywords distribution multivariate
##' @examples
##' 
##'         # generate some Poisson counts on the fly
##'           dose <- c(rep(0,50),rep(1,50),rep(5,50),rep(10,50))
##'           data <- cbind(dose,rpois(200,20*(1+dose*.5*(1-dose*0.05))))
##' 
##'         # neg. log-likelihood of Poisson model with 'linear-quadratic' mean: 
##'           nlogf <- function (x) { 
##'           ds <- data[, 1]
##'           y  <- data[, 2]
##'           g <- x[1] * (1 + ds * x[2] * (1 - x[3] * ds)) 
##'           return(sum(g - y * log(g)))
##'           }
##' 
##' 	# for example define
##'           x <- list(label=c("a","b","c"),est=c(10.,10.,.01),low=c(0,0,0),upp=c(100,20,.1))
##' 
##' 	# get MLEs using dfp:
##' 	  r <- dfp(x,f=nlogf)
##'           x$est <- r$est
##'           plkhci(x,nlogf,"a")
##'           plkhci(x,nlogf,"b")
##'           plkhci(x,nlogf,"c")
##' 	# e.g. 90% confidence bounds for "c" 
##'           plkhci(x,nlogf,"c",prob=0.9)
##'         
##' 
##' @export
##' 
plkhci <- function (x, nlogf, label, prob = 0.95, eps = 0.001, nmax = 10, nfcn = 0)
{
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
        stop() 
    }
    small <- 1e-08
    nd1 <- npar - 1
    sens <- 0.5
    rhs <- numeric(npar)
    disc <- numeric(npar)
    qchi <- stats::qchisq(prob, 1)
    if (!is.character(label)) 
        stop("label needs to be of type character")
    p <- match(label, x$label)
    if (is.na(p)) 
        stop("label not found")
    xt <- ftrf(x$est, x$low, x$upp)
    fmle <- nlogf(x$est)
    nfcn <- nfcn + 1
    cat(" neg. log. likelihood: ", fmle, "\n", "\n")
    cat(" will attempt to compute both bounds (+/- direction)", 
        "\n")
    del <- dqstep(x, nlogf, sens)
    h <- logit.hessian(x, nlogf, del, dapprox = FALSE, nfcn)
    nfcn <- h$nfcn
    ddf <- as.matrix(h$ddf)
    df <- h$df
    dbb0 <- ddf[p, p]
    tmp0 <- 0
    if(npar > 1) {
        ddl <- ddf[-p, -p]
        dbo <- ddf[p, -p]
        b <- diag(1, nd1)
        ggl <- solve(ddl, b, tol = 1e-16)
        dbl0 <- ggl %*% dbo
        tmp0 <- t(dbo) %*% dbl0
    }
    cat("\n")
    for (idir in c(1, -1)) {
        xt0 <- xt
        xt0[p] <- xt0[p] + idir * sqrt(qchi/(dbb0 - tmp0))/2
        if(npar > 1)
            xt0[-p] <- xt0[-p] - idir * as.vector(dbl0) * sqrt(qchi/(dbb0 - tmp0))/2
        f0 <- nlogf(btrf(xt0, x$low, x$upp))
        nfcn <- nfcn + 1
        if (idir == 1) 
            cat("trying lower bound ------------------------", 
                "\n")
        if (idir == -1) 
            cat("trying upper bound ------------------------", 
                "\n")
        cat("starting at:   ", 2 * (f0 - fmle), "\n")
        cat("initial guess: ", btrf(xt0, x$low, x$upp), "\n")
        cat("\n")
        cat("begin Newton-Raphson search for profile lkh conf. bounds:", 
            "\n")
        cat("eps value for stop criterium:", eps, "\n")
        cat("nmax                        :", nmax, "\n")
        for (n in 1:nmax) {
            x0 <- x
            x0$est <- btrf(xt0, x$low, x$upp)
            f0 <- nlogf(x0$est)
            nfcn <- nfcn + 1
            h <- logit.hessian(x0, nlogf, del, dapprox = FALSE, 
                nfcn)
            nfcn <- h$nfcn
            ddf <- as.matrix(h$ddf)
            df <- h$df
            fdd <- ddf   
            fdd[1, 1] <- -df[p]
            rhs[1] <- (fmle - f0 + 0.5 * qchi)
            if(npar > 1) {
                ddl <- ddf[-p, -p]
                dbo <- ddf[p, -p]
                fdd[1, 2:npar] <- -df[-p]
                rhs[2:npar] <- df[-p]
                fdd[-1, -1] <- ddl
                fdd[2:npar, 1] <- dbo
            }
            ggl <- solve(fdd, tol = 1e-16)
            if (n <= 2) 
                rneps <- 0.2
            if (n > 2 && n <= 4) 
                rneps <- 0.8
            if (n > 4) 
                rneps <- 1
            nepsc <- 0
            xt00 <- xt0
            disc0 <- disc <- as.vector(ggl %*% rhs)
            rhs1 <- 1000
            disc <- rneps * disc0
            nz <- 0
            xt0[p] <- xt00[p] + disc[1]
            if(npar > 1)
                xt0[-p] <- xt00[-p] + disc[2:npar]
            tmp <- numeric(npar)
            tmp[abs(disc) <= eps] <- 1
            nz <- sum(tmp)
            f0 <- nlogf(btrf(xt0, x$low, x$upp))
            nfcn <- nfcn + 1
            rhs1 <- abs(fmle - f0 + 0.5 * qchi)
            estimate <- btrf(xt0, x$low, x$upp)
            if (nz == npar && nepsc == 0) {
                cat("\n", "CONVERGENCE: ", trunc(n), " iterations", 
                  "\n")
                status <- "converged"
                cat("\n")
                cat("chisquare value is: ", 2 * (f0 - fmle), 
                  "\n")
                cat("confidence bound of ", x$label[p], " is ", 
                  estimate[p], "\n")
                if(npar > 1)
                    cat("log derivatives:    ", df[-p], "\n")
                m.out <- cbind(x$label, signif(estimate, 6), 
                  signif(df, 6), signif(diag(ddf), 6))
                dimnames(m.out) <- list(1:npar, c("label", "estimate", 
                  "log deriv", "log curv"))
                print(m.out, quote = FALSE)
                cat("\n")
                break
            }
        }
        if (idir == -1) {
            xlow <- estimate[p]
        }
        else {
            xupp <- estimate[p]
        }
    }
    return(c(xlow, xupp))
}

