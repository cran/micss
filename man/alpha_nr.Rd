\name{alpha_nr}
\alias{alpha_nr}
\title{alpha_nr}
\usage{
alpha_nr(y, k)
}
\arguments{
\item{y}{A numeric vector.}

\item{k}{Fraction of the upper tail to be used to estimate of the tail index.}
}
\value{
\itemize{
  \item \code{alpha}: Estimated tail index.
  \item \code{sd.alpha}: Standard error.
}
}
\description{
Computes the estimator of the tail index proposed by Nicolau & Rodrigues (2019).
}
\examples{
alpha_nr(rnorm(500),k=0.1)
}
\references{
J. Nicolau and P.M.M. Rodrigues (2019): A new regression-based tail index estimator. The Review of Economics and Statistics 101, 667-680.
}
\seealso{
\link{alpha_hill}
}
