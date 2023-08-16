test_that(
  "micss with data example. Returns number of breaks = 2 ",
  {
    expected <- 2
    data(logReturnsRandDollar)
    e <- whitening(data$rand.dollar)$e # whitening
    actual <- micss(e)$icss$nb
    expect_equal(actual,expected)
  }
)
