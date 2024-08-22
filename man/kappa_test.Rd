\name{kappa_test}
\alias{kappa_test}
\title{CUMSUMQ test to test for changes in the unconditional variance}
\description{Computes the CUMSUMQ test to test for changes in the unconditional variance and reports the p-value adapted to the tail index and sample size}
\usage{
kappa_test(e,sig.lev=0.05,alpha=NULL,kmax=NULL)
}
\arguments{
  \item{e}{A numeric vector. Stationary variable on which the constancy of unconditional variance is tested.}
  \item{sig.lev}{Significance level. The default value is 0.05.}
  \item{alpha}{Tail index. Must be a number between 2 and 4. The default value is 4.}
  \item{kmax}{Maximum lag to be used for the estimation of the long-run fourth order moment. If not reported, an automatic procedure computes it depending on the sample size.}
}
\details{It is only computed if the sample size is greater than 25 observations.
}
\value{
  \item{kappa}{CUMSUMQ test.}
  \item{tb}{Possible time of the break (with maximum value of the statistic).}
  \item{cv}{critical value at the specified significance level.}
  \item{p.val}{p-value.}
}
\references{J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.}
\author{J.L. Carrion-i-Silvestre and A. Sanso.}
\seealso{
 \code{\link{micss}}
}
\examples{
data(logReturnsRandDollar)
e <- whitening(data$rand.dollar)$e # whitening
kappa_test(e)
}
