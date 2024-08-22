test_that(
  "Critical values. Expected = 1.444563 ",
  {
    expected <- 1.421182
    actual <- cv.kappa(100,3.5,0.05)
    expect_equal(round(actual,5), round(expected,5))
  }
)
