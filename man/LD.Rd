\name{LD}
\alias{LD}
\alias{LD.formula}
\alias{LD.LD}
\alias{confint}
\alias{confint.LD}
\alias{lethal.dose}
\alias{lethal.dose.formula}
\alias{lethal.dose.LD}

\title{
  Lethal dose
}
\description{
  Compute generalized additive model for lethal dose finding and
  corresponding confidence intervals using parametric bootstrap
  approaches
}
\usage{
\method{LD}{formula}(formula, groups = NULL, experiment = NULL,
   lethal.dose = c(50, 10),
   dose_trafo = c("sqrt", "log", "none"),
   data, family = nb(), ...)

\method{LD}{LD}(object, lethal.dose = NULL, group = NULL, ...)

\method{confint}{LD}(object, level = 0.95, lethal.dose = NULL,
       B1 = 20, B2 = 100, newdata = NULL,
       myapply = mclapply, ...)

### just wrappers to LD functions as defined above
### (i.e., with same arguments)
\method{lethal.dose}{formula}(...)
\method{lethal.dose}{LD}(...)
}

\arguments{
  \item{formula}{
    Specify the main dose-response model in the form
    \code{outcome ~ dose}. Only used to extract \code{outcome} and
    \code{dose} varibles. Internally a more complex model is fitted. See
    details.
  }
  \item{groups}{
    (optional) define the name of the grouping factor (if available in
    the data set). This could be for example the species, the strain, or
    any other factor (with two or more groups).
  }
  \item{experiment}{
    (optional) define the name of the experiment variable if multiple
    replicates (biological or technical) are used. This variable is
    included as randondom effect in the lethal dose model.
  }
  \item{lethal.dose}{
    vector of LDXX values. Per default, \code{c(50, 10)} for LD50 and
    LD10. If no value is specified in functions for class \code{"LD"},
    the LD values of the fitted \code{"LD"} object are used.
  }
  \item{dose_trafo}{
    define a character string for the transformation that is used to
    stabilize the estimation of the dose effect. Note that this does not
    influence the form of the dose-response relation ship as we use a
    smooth function of dose (or a smooth function of the transformation
    of dose) in the model. This is really more a technical option.
  }
  \item{data}{
    data set with an outcome variable (i.e., number of surviving
    cells/bacteria/...), a dose variable and optionally a grouping
    variable (with two or more groups) that is used to comparte the
    groups. Further more, the data set can optionally contain an
    indicator for the experiment, which is used to model possible
    experiment effects.
  }
  \item{family}{
    distributional assumption. Currently only \code{negbin} is
  }
  \item{object}{
    Lethal dose model fitted using function \code{LD}.
  }
  \item{group}{
    (optional) specify a factor level of \code{groups} to extract only
    LDxx values for one of the groups.
  }
  \item{level}{
    coverage probability of confidence interval
  }
  \item{newdata}{
    (optional) list of data frames, one per strain, with a fine grid for dose.
    Only for expert users. Usually this data set is built internally.
  }
  \item{B1}{
    number of new sampled data sets.
  }
  \item{B2}{
    number of draws from the distribution of the fitted coefficients.
  }
  \item{myapply}{
    either \code{\link{lapply}} or \code{\link{mclapply}} (for parallel
    computing on Linux machines). In the latter case, use the additional
    argument \code{mc.cores} to specifiy the number of processes to use.
  }
  \item{\dots}{
    further arguments.
  }
}
\details{
  Explain the model that is used to derive the lethal doses here...
}
\author{
  Benjamin Hofner <benjamin.hofner@fau.de>
}
\seealso{
  methods for fitted \code{LD} models:
  \code{\link{predict.LD}}, \code{\link{summary.LD}}, etc.
}
\examples{
## model survival fractions for UV radiation experiment
data("geoderm.uv")
mod.uv <- LD(value ~ time, groups = "strain", experiment = "replicate",
             dose_trafo = "sqrt", data = geoderm.uv)

## for more examples see ?summary.LD
}

\keyword{models}
\keyword{nonlinear}
\keyword{nonparametric}
