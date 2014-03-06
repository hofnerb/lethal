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
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
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
## currently none available
}
\keyword{methods}
