test_that(
  "kappa.test with data example. Returns 3.548",
  {
    expected <- 2.047063
    data(logReturnsRandDollar)
    e <- whitening(data$rand.dollar)$e # whitening
    actual <- kappa_test(e)$kappa
    expect_equal(round(actual,5), round(expected,5))
  }
)
