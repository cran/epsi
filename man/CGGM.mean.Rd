\name{CGGM.mean}
\alias{CGGM.mean}
\alias{CGGM.lts}
\alias{CGGM.autoscale}
\title{Corner-Preserving Smoothing}
\description{
  Computation of two dimensional redescending M-kernel estimators.
}

\usage{
  CGGM.mean(data, h, g=NULL, silent=FALSE)

  CGGM.lts(data, h, g=NULL, trim=0, silent=FALSE)

  CGGM.autoscale(data, h, silent=FALSE)
}

\arguments{
  \item{data}{numerical matrix of observed data.}
  \item{h}{positive number. Bandwidth for the kernel.}
  \item{g}{optional positive number. Scale parameter. If \code{g} is \code{NULL}
    (default), the scale parameter is determined automatically by the
    function \code{CGGM.autoscale}.}
  \item{trim}{optional number within [0,1). Part of observations trimmed
    by \code{CGGM.lts}}
  \item{silent}{optional boolean. If true, \code{CGGM.autoscale}
    produces no output.}
}
\details{
  \code{CGGM.mean} implements a corner-preserving smoothing method
  introduced by Chu et al. (1998) which is based on a redescending
  M-kernel estimator. As kernel and score function the density of the
  standard normal distribution is used. A robust version of this
  estimator is introduced by Hillebrand (2002) and implemented in
  \code{CGGM.lts}.
  
  \code{CGGM.autoscale} calculates the median of the interquartile
  ranges within the 'windows' used in \code{CGGM.mean} and
  \code{CGGM.lts}. This can be used as scale parameter.
}

\value{
  Return value is a numerical matrix containing the smoothed data.
}
\references{
  Chu, C.K., Glad, I.K., Godtliebsen, F., Marron, J.S. (1998)
  Edge-Preserving Smoothers for Image Processing,
  \emph{J. Amer. Statis. Assoc.} 93, 526-541.

  Hillebrand, M. (2002) On Robust Corner-Preserving Smoothing in Image
  Processing, \emph{Carl-von-Ossietzky-Universität Oldenburg, Dissertation}
  (http://docserver.bis.uni-oldenburg.de/publikationen/dissertation/2003/hilonr03/hilonr03.html).
}

\author{
  Tim Garlipp (garlipp@mathematik.uni-oldenburg.de)
}

\examples{

  y <- matrix(rep(0,60*60),nrow=60)
  y[21:40,21:40]<-1
  y <- y + matrix(rnorm(60*60,0,0.1),nrow=60)
  image(y,col=gray(seq(0,1,1/255)))

  ymean <- CGGM.mean(y,0.04)
  image(ymean,col=gray(seq(0,1,1/255)))
}

\keyword{ cluster }
\keyword{ regression }
\keyword{ robust }
\keyword{ nonparametric }
\keyword{ smooth }

