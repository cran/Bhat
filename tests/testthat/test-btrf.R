test_that("btrf works", {
  expect_equal(btrf(c(2,-5), c(1,10), c(5,50)), c(4.5231883,10.267714))
})

test_that("btrf clamps 0 to mid", {
  expect_equal(btrf(c(0,0,0), c(0,1,-5), c(2,3,11)), c(1,2,3))
})

