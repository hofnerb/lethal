\name{geoderm.uv}
\alias{geoderm.uv}
\docType{data}

\title{UV tolerance of \emph{Geodermatophilus}}

\description{
  A novel actinobacterium \emph{Geodermatophilus poikilotrophi}
  (G18T = DSM 44209 = CCUG 63018) was isolated from dolomitic marble
  collected from outcrops in Samara, near of Namib Desert (Namibia).
  This bacterium was compared to \emph{Geodermatophilus obscurus} (DSM
  43160T). The two strains of \emph{Geodermatophilus} were exposed to
  ultra-violet radiation with a dose of 5-10 \eqn{J/s/m^2} for
  various time periods and the number of colony-forming units per ml
  (a proxy for the number of surviving bacteria) was counted.
}

\usage{
data(geoderm.uv)
}

\format{
  A data frame with 72 observations and the following 4 variables:
  \describe{
    \item{ \code{time}}{time period of UV exposure in minutes.}
    \item{ \code{strain}}{type of strain: 44209 (\emph{Geodermatophilus
	poikilotrophi}) or 43160 (\emph{Geodermatophilus obscurus}).}
    \item{ \code{replicate}}{biological replicate (A or B).}
    \item{ \code{value}}{number of colony-forming units per ml
      (c.f.u. ml^{-1}).}
  }
}

\details{
  For usage and results see other manual pages.
}

\references{
  M. C. Montero-Calasanz, B. Hofner, M. Göker,M. Rohde, C. Spröer, K.
  Hezbri, M. Gtari,  P. Schumann, H.-P. Klenk, \emph{Geodermatophilus
  poikilotrophi sp. nov., a multi-tolerant actinomycete isolated from
  dolomitic marble}, submitted.

}
\examples{
data(geoderm.uv)
summary(geoderm.uv)
}
\keyword{datasets}
