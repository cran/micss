test_that(
  "p-values. Expected = 0.05 ",
  {
    expected <- 0.04592267
    actual <- p.val.kappa(1.444563,100,3.5)
    expect_equal(round(actual,5), round(expected,5))
  }
)
