\name{whitening}
\alias{whitening}
\title{Whitening}
\usage{
whitening(y, kmax = NULL)
}
\arguments{
\item{y}{A numeric vector. Variable to be whiten.}

\item{kmax}{Maximum lag to be used for the long-run estimation of the variance.
If not specified uses [12*(t/100)^(1/4)].}
}
\value{
\itemize{
  \item \code{e}: Whiten variable.
  \item \code{rho}: Vector of autoregressive parameters.
  \item \code{lag}: number of lags used.
}
}
\description{
Eliminates the autocorrelation of a variable using an AR model.
}
\details{
Selects the model using the Bayes Information Criterium.
}
\examples{
whitening(rnorm(100))
}
