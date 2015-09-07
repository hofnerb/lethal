lethal
======

[![Build Status](https://travis-ci.org/hofnerb/lethal.svg)](https://travis-ci.org/hofnerb/lethal)
[![Build status](https://ci.appveyor.com/api/projects/status/uoe5e6o43932a2u8?svg=true)](https://ci.appveyor.com/project/hofnerb/lethal)

`lethal`  computes lethal doses for count data based on generalized additive models (GAMs) together with parametric bootstrap confidence intervals for the lethal dose. `lethal` is designed for experiments with counts as outcome, which need a separate preparation for each measurment. Examples for such experiments are survival experiments where the survival is measured as the number of colony forming units (c.f.u.). In this case, one cannot measure one prepartation multiple times with various doses but one needs one experiment (with one or more biological replicates) for each dose.

## Installation:

- Latest development version from GitHub:

        library("devtools")
        install_github("hofnerb/lethal")

  To be able to use the `install_github()` command, one needs to install `devtools` first:
  
        install.packages("devtools")

## Using lethal

Afer installation (see above) we first need to load the package:

    library("lethal")

To illustrate the usage of the package `lethal` we will use a data set on the UV tolerance of _Geodermatophilus_ ([Montero-Calasanz, Hofner _et al._, 2014](https://github.com/hofnerb/lethal/blob/master/README.md#references)):

    data(geoderm.uv, package = "lethal")
    summary(geoderm.uv)
    
The data set has 72 observations and the following 4 variables: 
- `time` time period of UV exposure in minutes.
- `strain` type of strain: 44209 (_Geodermatophilus poikilotrophi_) or 43160 (_Geodermatophilus obscurus_).
- `replicate` biological replicate (A or B).
- `value` number of colony-forming units (c.f.u.) per ml.

Now, we want to model this data

    ## model survival fractions for UV radiation experiment
    mod.uv <- LD(value ~ time, groups = "strain", experiment = "replicate",
                 dose_trafo = "sqrt", data = geoderm.uv,
                 family = negbin(theta = c(0.5, 10)))

To display the results, one can use

    ## extract model and LDs
    mod.uv
    ## a richer representation of the model (with LDs)
    summary(mod.uv)
    
or we can extract lethal doses only    
    
    ## extract LDs only
    LD(mod.uv)

To plot the results is very simple:

    ## plot the results
    plot(mod.uv)
    ## with better labels:
    plot(mod.uv, xlab = "Time (min)", ylab = expression(c.f.u.ml^-1))

### Confidence intervals

To obtain confidence intervals, we need to rely on parametric bootstrap approaches. To simplify things, this is already implemented

    ## compute confidence intervals
    ci.uv <- confint(mod.uv)
    ## extract confidence intervals (and LDs)
    ci.uv

As bevore, a simple way to plot the results is available

    ## graphic with survival fractions and confidence intervals:
    plot(ci.uv, xlab = "Time (min)", upper_ylab = expression(c.f.u.ml^-1),
    mar = c(4, 9.3, 2, 2.5))
    ## add labels
    mtext(rep(c("LD10", "LD50"), 3), side = 4,
          at = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2),
          cex = 0.75, las = 2)

### References 

M. d. C. Montero-Calasanz, B. Hofner, M. Göker, M. Rohde, C. Spröer, K. Hezbri, M. Gtari, P. Schumann, H.-P. Klenk (2014). "Geodermatophilus poikilotrophi sp. nov., a multi-tolerant actinomycete isolated from dolomitic marble." BioMed Research International. 2014(914767):1-11. [[PDF]](http://downloads.hindawi.com/journals/bmri/2014/914767.pdf)
