parabolic_well = function(x, a,b,c) { 3*(x[1]-a)^2 + 0.5*(x[2]-b)^2 + c }

rosenbrock_banana_ev = function(x, a,b) {
    n = length(x)
    if (n%%2 != 0) stop("Vector length must be even")
    x_odd = x[seq.int(1,n,2)]
    x_even = x[seq.int(2,n,2)]
    return (sum(a*(x_odd^2-x_even)^2 + (x_odd-b)^2))
}


##' Double-well polynomial function
##'
##' Has the following critical points and the diagonal Hessian elements:
##'  (x1,x2)  [H_{11}, H_{22}]
##'  (0,b1) [12*a1*a2, 6*(b1-b2)]
##'  (0,b2) [12*a1*a2, 6*(b2-b1)]
##'  (a1,b1) [12*a1*(a1-a2), 6*(b1-b2)]
##'  (a2,b1) [12*a2*(a2-a1), 6*(b1-b2)]
##'  (a1,b2) [12*a1*(a1-a2), 6*(b2-b1)]
##'  (a2,b2) [12*a2*(a2-a1), 6*(b2-b1)]
##' 
double_well = function(x, a1,a2,b1,b2) {
    (3*x[1]^4 - 4*(a1+a2)*x[1]^3 + 6*a1*a2*x[1]^2
        + 2*x[2]^3 - 3*(b1+b2)*x[2]^2 + 6*b1*b2*x[2])
}




test_that("dfp works", {
    fn = function(x) parabolic_well(x, a=0.25, b=-1.75, c=5)

    x=list(label = c("x1", "x2"),
           est = c(0.1, 0.2),
           low = c(-10, -10),
           upp = c(10, 10))

    res = dfp(x, f=fn)

    ##str(res)

    expect_equal(res$status, 0)
    expect_equal(res$fmin, 5, tolerance=1e-6)
    expect_equal(res$label, c("x1", "x2"))
    expect_equal(res$est, c(0.25, -1.75), tolerance=1e-6)

    expect_equal(res$low, c(-0.55, -3.6186), tolerance=1e-3)
    expect_equal(res$high, c(1.047362, 0.254630), tolerance=1e-3)
    
    expect_lte(res$nfcn, 30)
})


test_that("dfp on Rosebrock works", {
    fn = function(x) rosenbrock_banana_ev(x, a=100, b=1)

    x=list(label = c("x1", "x2", "x3", "x4"),
           est = c(0.1, 0.2, 0.3, 0.4),
           low = c(-5, -5, -5, -5),
           upp = c(10, 10, 10, 10))

    res = dfp(x, f=fn)

    ##str(res)

    expect_equal(res$status, 0)
    expect_equal(res$fmin, 0, tolerance=1e-5)
    expect_equal(res$est, c(1, 1, 1, 1), tolerance=1e-2)

    expect_equal(res$low, c(-0.32, -1.464, -0.319, -1.463), tolerance=1e-2)
    expect_equal(res$high, c(2.42, 3.85, 2.42, 3.85), tolerance=1e-2)
    
    expect_lte(res$nfcn, 405)
})


test_that("dfp on double-well works", {
    fn = function(x) double_well(x, a1=-0.5, a2=1, b1=0, b2=1)

    x=list(label = c("x", "y"),
           est = c(0, 0),
           low = c(-5, -5),
           upp = c(10, 10))

    res = dfp(x, f=fn)

    str(res)

    expect_equal(res$status, 0)
    expect_equal(res$fmin, -3, tolerance=1e-5)
    expect_equal(res$est, c(1,1), tolerance=1e-5)

    expect_equal(res$low, c(0.545, 0.220), tolerance=5e-3)
    expect_equal(res$high, c(1.47, 1.81), tolerance=5e-3)
    
    expect_lte(res$nfcn, 60)
})


test_that("confined dfp on double-well finds the other minimum", {
    ## 100x horizontally stretched double-well 
    fn = function(x) double_well(0.01*x, a1=-0.5, a2=1, b1=0, b2=1)

    x=list(label = c("x", "y"),
           est = c(-10, -10),
           low = c(-150, -150),
           upp = c(0, 300)) # do not let x[1] go above 0

    res = dfp(x, f=fn)

    str(res)

    expect_equal(res$status, 0)
    expect_equal(res$fmin, -1.3125, tolerance=5e-5)
    expect_equal(res$est, c(-50,100), tolerance=1e-2) # should find another minimum

    expect_equal(res$low, c(-116.9, 20.2), tolerance=5e-3)
    expect_equal(res$high, c(-9.92, 173.92), tolerance=5e-3)
    
    expect_lte(res$nfcn, 65)
})


test_that("dfp can pass additional arguments", {

    x=list(label = c("x", "y"),
           est = c(0, 0),
           low = c(-5, -5),
           upp = c(10, 10))

    res = dfp(x, f=double_well, a1=-0.5, a2=1, b1=0, b2=1)

    str(res)

    expect_equal(res$status, 0)
    expect_equal(res$fmin, -3, tolerance=1e-5)
    expect_equal(res$est, c(1,1), tolerance=1e-5)

    expect_equal(res$low, c(0.545, 0.220), tolerance=5e-3)
    expect_equal(res$high, c(1.47, 1.81), tolerance=5e-3)
    
    expect_lte(res$nfcn, 60)
})

