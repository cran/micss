\name{micss}
\alias{micss}
\title{Modiffied Iterative Cumulative Sum of Squares Algorithm}
\description{Implements the MICSS algorithm of Carrion-i-Silvestre & Sansó (2023)}
\usage{
micss(e,sig.lev=0.05,kmax=NULL,alpha=NULL,tail.est="NR",k=0.1)
}
\arguments{
  \item{e}{A numeric vector. Stationary variable on which the constancy of unconditional variance is tested.}
  \item{sig.lev}{Significance level. The default value is 0.05.}
  \item{kmax}{Maximum lag to be used for the estimation of the long-run fourth order moment. If not reported, an automatic procedure computes it depending on the sample size.}
  \item{alpha}{Tail index. If not reported, it is estimated automatically.}
  \item{tail.est}{Estimator of the tail index. The default value is "NR", which uses Nicolau & Rodrigues (2019) estimator. "Hill" uses the Hill's (1975) estimator.}
  \item{k}{Fraction of the upper tail to be used to estimate of the tail index.}
}
\details{The tail index is estimated using the absolute values.
}
\value{
  \item{icss}{An object with the output of the \link{icss} algorithm.}
  \item{alpha}{An object with the output of the \link{estimate.alpha}.}
}
\references{J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.

  B. Hill (1975): A Simple General Approach to Inference About the Tail of a Distribution. The Annals of Mathematical Statistics 3, 1163-1174.

  J. Nicolau & P.M.M. Rodrigues (2019): A new regression-based tail index estimator. The Review of Economics and Statistics 101, 667-680.}
\author{J.L. Carrion-i-Silvestre and A. Sansó}

\seealso{
 \code{\link{icss}}
 \code{\link{estimate.alpha}}
 \code{\link{print.micss}}
 \code{\link{plot.icss}}
}
\examples{
set.seed(2)
e <- c(stats::rnorm(200),3*stats::rnorm(200))
o <- micss(e)
print.micss(o)
#
# The following example replicates some of the results of Table 6 in
# Carrion-i-Silvestres & Sanso (2023)
\donttest{
data(logReturnsRandDollar)
e <- whitening(data$rand.dollar)$e  # pre-whitening
m <- micss(e)
print.micss(m)}
}
