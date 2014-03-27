\name{methods}
\alias{predict.LD}
\alias{print.LD}
\alias{print.LDconfint}
\alias{summary.LD}
\alias{anova.LD}

\alias{get_x}

\title{
  Methods for \code{"LD"} models
}
\description{
  Methods for models fitted by \code{\link{LD}}
}
\usage{
## predict function
\method{predict}{LD}(object, group = NULL,
        type = c("response", "link", "terms", "lpmatrix"),
        newdata = NULL, ...)

## print function
\method{print}{LD}(x, ...)
## summary function
\method{summary}{LD}(object, ...)
## anova function
\method{anova}{LD}(object, ...)

## print function for confidence intervals
\method{print}{LDconfint}(x, ...)

## get dose values from model
get_x(object, ...)
}
\arguments{
  \item{object, x}{
    object fitted via \code{\link{LD}}.
  }
  \item{group}{
    (optional) specify a factor level of \code{groups} to extract only
    LDxx values for one of the groups.
  }
  \item{type}{
    type of prediction. See \code{\link{predict.gam}} in package
    \code{mgcv} for details.
  }
  \item{newdata}{
    (optional) new data used for predictions.
  }
  \item{\dots}{
    additional arguments passed to callies.
  }
}
\details{

}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  Benjamin Hofner
}
\seealso{
  \code{\link{LD}}
}

\examples{
## model survival fractions for UV radiation experiment
data("geoderm.uv")
mod.uv <- LD(value ~ time, groups = "strain", experiment = "replicate",
             dose_trafo = "sqrt", data = geoderm.uv,
             family = negbin(theta = c(0.5, 10)))

## extract model and LDs
mod.uv
## a richer representation of the model (with LDs)
summary(mod.uv)

## extract LDs only
LD(mod.uv)

## plot the results
plot(mod.uv)
## with better labels:
plot(mod.uv, xlab = "Time (min)", ylab = expression(c.f.u.ml^-1))

\dontrun{## automatically as this takes a while.
## Users can run this code at any time
ci.uv <- confint(mod.uv)
## extract confidence intervals (and LDs)
ci.uv

## graphic with survival fractions and confidence intervals:
plot(ci.uv, xlab = "Time (min)", upper_ylab = expression(c.f.u.ml^-1),
     mar = c(4, 9.3, 2, 2.5))
## add labels
mtext(rep(c("LD10", "LD50"), 3), side = 4,
      at = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2),
      cex = 0.75, las = 2)
}
}
\keyword{methods}
