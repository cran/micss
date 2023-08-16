\name{alpha_hill}
\alias{alpha_hill}
\title{alpha_hill}
\usage{
alpha_hill(x, k)
}
\arguments{
\item{x}{A numeric vector.}
\item{k}{Fraction of the upper tail to be used to estimate of the tail index.}
}
\value{
\itemize{
  \item \code{alpha}: Estimated tail index.
  \item \code{sd.alpha}: Standard error.
  \item \code{s}: Number of observations used in the estimation.
}
}
\description{
Computes the estimator of the tail index proposed by Hill (1975).
}
\examples{
alpha_hill(rnorm(500),k=0.1)
}
\references{
B. Hill (1975): A Simple General Approach to Inference About the Tail of a Distribution. The Annals of Mathematical Statistics 3, 1163-1174.
}
\seealso{
\link{alpha_nr}
}
