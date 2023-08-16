test_that(
  "icss with data example. Returns number of breaks = 2 ",
  {
    expected <- 2
    data(logReturnsRandDollar)
    e <- whitening(data$rand.dollar)$e # whitening
    actual <- icss(e)$nb
    expect_equal(actual,expected)
  }
)
