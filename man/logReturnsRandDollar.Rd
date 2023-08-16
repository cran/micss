\name{logReturnsRandDollar}
\alias{logReturnsRandDollar}
\docType{data}
\title{Data used in the examples}
\description{Log returns of the exchange rate South African Rand versus United States Dollar.}
\usage{data(logReturnsRandDollar)}
\value{Time series with 7705 observations.}
\source{Paulo Rodrigues}
\references{J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.}
\author{J.L. Carrion-i-Silvestre and A. Sansó}
\examples{
data(logReturnsRandDollar)
names(data)
# The following example replicates some of the results of Table 6 in
# Carrion-i-Silvestres & Sanso (2023)
\donttest{
data(logReturnsRandDollar)
e <- whitening(data$rand.dollar)$e  # pre-whitening
m <- micss(e)
print.micss(m)}
}
\keyword{Datasets}
