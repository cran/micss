\name{icss}
\alias{icss}
\title{Iterative Cumulative Sum of Squares Algorithm}
\description{Implements the ICSS algorithm of Inclan and Tiao (1994) using the CUMSUMQ test detailed in Carrion-i-Silvestre & Sansó (2023)}
\usage{
icss(e,sig.lev=0.05,kmax=NULL,alpha=NULL)
}
\arguments{
  \item{e}{A numeric vector. Stationary variable on which the constancy of unconditional variance is tested.}
  \item{sig.lev}{Significance level. The default value is 0.05}
  \item{kmax}{Maximum lag to be used for the estimation of the long-run fourth order moment. If not reported, an automatic procedure computes it depending on the sample size.}
  \item{alpha}{Tail index. If not reported, it is set at 4, which was the implicit assumption of Inclan & Tiao (1994). Values between 2 and 4 are allowed because this function is used by \link{micss}}
}
\details{Used internally by \link{micss}.
}
\value{
  \item{nb}{Number of breaks found by the algorithm.}
  \item{tb}{Vector with the time breaks.}
}
\references{J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.

  C. Inclan & G.C. Tiao (1994): Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance. Journal of the American Statistical Association 89, 913-923.}
\author{J.L. Carrion-i-Silvestre and A. Sanso}
\seealso{
 \code{\link{micss}}
 \code{\link{print.icss}}
 \code{\link{plot.icss}}
}
\examples{
set.seed(2)
e <- c(stats::rnorm(200),3*stats::rnorm(200))
o <- icss(e)
print.icss(o)
}
