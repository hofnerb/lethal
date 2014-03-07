################################################################################
## compute lethal dose model

LD <- lethal.dose <- function(...)
    UseMethod("LD")

LD.formula <- function(formula, groups = NULL, experiment = NULL,
               lethal.dose = c(50, 10),
               dose_trafo = c("sqrt", "log", "none"),
               data, family = negbin(theta = c(0.01, 1000)), ...) {

    if (!all(grepl("Negative Binomial", family$family)))
        stop("family must currently be negbin")

    if (!is.null(groups) && nlevels(data[, groups]) != 2)
        stop("currently only two groups are supported")

    dose_trafo <- match.arg(dose_trafo)

    dose <- formula[[3]]
    outcome <- formula[[2]]

    if (length(dose) != 1)
        stop("Specify a formula of the type ", sQuote("outcome ~ dose"))

    if (!all(cc <- complete.cases(data))) {
        warning("Any cases with missing values are dropped")
        data <- data[cc, ]
    }

    fm_dose <- dose
    if (dose_trafo != "none") {
        fm_dose <- paste0("I(", dose_trafo, "(", fm_dose, "))")
    }

    fm_groups <- fm_experiment <- ""
    if (!is.null(groups)) {
        fm_groups <- paste("+", groups)
        fm_dose <- paste0(fm_dose, ", by = ", groups)
    }
    if (!is.null(experiment))
        fm_experiment <- paste0("+ s(", experiment, ", bs = 're')")

    fm_dose <- paste0("s(", fm_dose, ", k = 10, bs = 'ps')")
    formula <- paste(outcome, "~", fm_dose, fm_groups, fm_experiment)

    model <- gam(as.formula(formula, env = .GlobalEnv), data = data, family = family, ...)

    if (dose_trafo != "none") {
        if (dose_trafo == "log") {
            fct <- exp
        } else {
            fct <- function(x) x^2
        }
    } else {
        fct <- function(x) x
    }

    RET <- list(model = model,
                data = data,
                trafo = dose_trafo,
                family = family,
                additional_arguments = list(...))
    RET$variables <- list(outcome = as.character(outcome),
                          dose = as.character(dose),
                          groups = groups,
                          experiment = experiment)
    class(RET) <- c("LD")
    RET$lethal.dose <- LD(RET, lethal.dose)
    return(RET)
}

LD.LD <- function(object, lethal.dose = NULL, group = NULL, ...) {
    if (is.null(lethal.dose) ||
        (!is.null(attr(object$lethal.dose, "values")) &&
         all(attr(object$lethal.dose, "values") == lethal.dose))) {

        RET <- object$lethal.dose
        if (!is.null(group)) {
            if (is.list(RET)) {
                RET <- lapply(RET, function(x, group) x[group, , drop = FALSE],
                              group = group)
            } else {
                RET <- RET[group, , drop = FALSE]
            }
        }
        return(RET)
    }

    pr <- predict(object, group)
    RET <- LD(pr, dose = get_x(pr), lethal.dose = lethal.dose, ...)
    return(RET)
}

LD.LDpred <- function(pred, dose = NULL, lethal.dose = c(50, 10), ...) {
    ## pred can be a vector or a matrix of predictions
    if (is.null(dose)) {
        warning("get_x needs to be fixed")
        dose <- get_x(pred)
    }

    ## make sure pred is a matrix
    if (!is.matrix(pred))
        pred <- matrix(pred, ncol = 1)

    LD <- matrix(NA, nrow = ncol(pred), ncol = 2)
    colnames(LD) <- c("dose", "predicted count")
    if (class(pred) == "LDpred")
        rownames(LD) <- paste0("group:", colnames(pred))
    if (!is.null(colnames(pred)))
        rownames(LD) <- colnames(pred)
    if (length(lethal.dose) > 1) {
        LD <- lapply(1:length(lethal.dose), function(i) return(LD))
        names(LD) <- paste("LD", lethal.dose, sep = "")
    }

    if (is.matrix(LD)) {
        for (i in 1:ncol(pred)){
            idx <- pred[,i] < (pred[1, i] * lethal.dose / 100)
            LD[i, ] <- c(dose[idx][1], pred[idx, i][1])
        }
    } else {
        for (j in 1:length(lethal.dose)) {
            for (i in 1:ncol(pred)){
                idx <- pred[,i] < (pred[1, i] * lethal.dose[j] / 100)
                LD[[j]][i, ] <- c(dose[idx][1], pred[idx, i][1])
            }
        }
    }
    attr(LD, "values") <- lethal.dose
    class(LD) <- "LDmatrix"
    return(LD)
}

LD.matrix <- function(pred, dose, lethal.dose = c(50, 10), ...) {
    LD <- LD.LDpred(pred = pred, dose = dose, lethal.dose = lethal.dose, ...)
    return(LD)
}

lethal.dose.formula <- lethal.dose.matrix <- lethal.dose.LD <- lethal.dose.LDpred <- function(...)
    LD(...)


#    if (is.null(cll[["lethal.dose"]]))
#        lethal.dose <- attr(object$lethal.dose, "values")

################################################################################

## upr <- p$fit + get_tval(mgcv) * (p$se.fit)
## lwr <- p$fit - get_tval(mgcv) * (p$se.fit)

predict.LD <- function(object, group = NULL,
                       type = c("response", "link", "terms", "lpmatrix"),
                       newdata = NULL, ...) {

    type <- match.arg(type)
    variables <- object$variables

    if (!is.null(newdata)) {
        pr <- predict(object$model, newdata = newdata, type = type, ...)
        attr(pr, "newdata") <- newdata
        class(pr) <- "LDpred"
        return(pr)
    }

    ## if no new data is given:

    if (is.null(group) && !is.null(variables$groups)) {
        ## make new data set for each group:
        g <- nlevels(object$data[, variables$groups])
        newdata <- lapply(1:g, make_newdata, object = object)
    } else {
        newdata <- make_newdata(object, group = group)
    }

    ## use predict function of mgcv
    if (is.data.frame(newdata)) {
        pr <- predict(object$model, newdata = newdata, type = type, ...)
        ## average over experiments (if available)
        if (!is.null(variables$experiment)) {
            pr <- tapply(pr, newdata[, variables$dose], mean)
            newdata <- newdata[, !(names(newdata) %in% variables$experiment)]
            newdata <- unique(newdata)
        }
        pr <- matrix(pr, ncol = 1)
        if (!is.null(group))
            colnames(pr) <- levels(object$data[, variables$groups])[group]
    } else {
        pr <- sapply(newdata, predict, object = object$model, type = type, ...)
        ## average over experiments (if available)
        if (!is.null(variables$experiment)) {
            pr <- sapply(1:g, function(i)
                         tapply(pr[, i], newdata[[i]][, variables$dose], mean))
            newdata <- lapply(newdata, function(ND) {
                              ND <- ND[, !(names(ND) %in% variables$experiment)]
                              unique(ND)
                          })
        }
        colnames(pr) <- levels(object$data[, variables$groups])
    }

    ## fixme se.fit = TRUE will kill this!

    RET <- pr
    attr(RET, "newdata") <- newdata
    attr(RET, "variables") <- variables
    class(RET) <- "LDpred"
    return(RET)
}

################################################################################
## Inference based on Wood p.261

## parametric bootstrap that takes the uncerteinty of smoothing-parameter
## selection into account

confint <- function(object, level, ...)
    UseMethod("confint")

confint.LD <- function(object, level = 0.95, lethal.dose = NULL, B1 = 20, B2 = 100,
                       newdata1 = NULL, newdata2 = NULL,
                       myapply = mclapply, ...) {

    data <- object$data
    fm <- object$model$formula
    responsefct <- object$model$family$linkinv
    dose <- object$variables$dose
    outcome <- object$variables$outcome

    if (is.null(lethal.dose))
        lethal.dose <- attr(object$lethal.dose, "values")

    ## <FIXME> What about experiments with more than two groups?
    if (!is.null(object$variables$groups)) {
        if (is.null(newdata1) && is.null(newdata2)) {
            newdata1 <- make_newdata(object, group = 1)
            LD1.pe <- LD(object, lethal.dose = lethal.dose, group = 1)
            newdata2 <- make_newdata(object, group = 2)
            LD2.pe <- LD(object, lethal.dose = lethal.dose, group = 2)
        } else {
            if (!is.null(newdata1)) {
                pred <- predict(object, newdata = newdata1)
                LD1.pe <- LD(pred, lethal.dose = lethal.dose,
                             dose = newdata1[, dose])
            }
            if (!is.null(newdata2)) {
                pred <- predict(object, newdata = newdata2)
                LD2.pe <- LD(pred, lethal.dose = lethal.dose,
                             dose = newdata2[, dose])
            }
        }
        if (!is.null(LD2.pe) && !is.null(LD2.pe)) {
            LD_diff.pe <- ld.diff(LD1.pe, LD2.pe)
        }
    } else {
        if (is.null(newdata1)) {
            newdata1 <- make_newdata(object, group = 1)
            LD1.pe <- LD(object, group = 1)
        } else {
            pred <- predict(object, newdata = newdata1)
            LD1.pe <- LD(pred, lethal.dose = lethal.dose,
                         dose = newdata1[, dose])
        }
    }

    probs <- c((1 - level)/2, 1 - (1 - level)/2)

    ## Inner parametric bootstrap (draw coefficients from the multivariate
    ## normal distribution they are assumed to follow as they are REML
    ## estimates)
    if (B2 > 0)
        rbeta <- rbind(mvrnorm(n = B2, mu = coef(object),
                               Sigma = vcov(object$model)))
    else
        rbeta <- coef(object)

    refit_model <- function(counter, model, data, fm) {
        ## Outer parametric bootstrap (draw new observations from fitted model):
        ## - We assume that the model holds, i.e., the data follows a negative
        ##   binomial distribution and the coefficients are consistent.
        ## - We simpulate data from the negative binomial distribution with the
        ##   conditional means as fitted bevor.

        ## <FIXME>: Add alternative sample functions for different families
        res <- rnegbin(rep(1, length(fitted(model))), mu = fitted(model),
                       theta = model$family$getTheta())
        ## make shure that new outcome is really a count variable > 0
        y <- round(res)
        y[y < 0] <- 0
        data.BS <- data
        data.BS[, outcome] <- y


        ## <FIXME>: Add alternative families:

        ## estimate smoothing parameters and overdispersion parameter theta from
        ## simulated data
        arguments <- list(formula = fm, data = data.BS, family = object$family)
        if (!is.null(object$additional_arguments))
            arguments <- c(arguments, object$additional_arguments)

        mod.BS <- do.call("gam", arguments)
        ## now estimate coefficients from original data but with the above
        ## smoothing parameters
        arguments$data <- data
        arguments$sp <- mod.BS$sp
        arguments$family <- negbin(theta = mod.BS$family$getTheta())
        mod <- do.call("gam", arguments)

        ## Inner parametric bootstrap for model on new data
        ## (draw coefficients from the multivariate normal distribution they are
        ## assumed to follow as they are REML estimates)
        if (B2 > 0)
            return( mvrnorm(n = B2, mu = coef(mod), Sigma = vcov(mod)) )
        else
            return(coef(mod))
    }
    ## use B1 - 1 loops as we already used the original model
    if (B1 > 1) {
        res <- myapply(1:(B1 - 1), refit_model, model = object$model,
                       data = data, fm = fm, ...)
        rbeta <- rbind(rbeta, do.call(rbind, res))
    }

    compute_CIs <- function(pred, probs, newdata, lethal.dose) {

        LD <- LD(pred, dose = newdata[, dose], lethal.dose = lethal.dose)
        if (any(sapply(LD, function(x) any(is.na(x))))) {
            if (is.matrix(LD)) {
                NAs <- sum(is.na(LD[,1]))
                if (NAs > probs[1] * nrow(LD))
                    warning("Upper limit of confidence interval might be to small")
                LD[is.na(LD[, 1]), 1] <- max(newdata[, dose])
            } else {
                NAs <- sapply(LD, function(x) sum(is.na(x[,1])))
                if (any(NAs > probs[1] * nrow(LD[[1]])))
                    warning("Upper limit of confidence interval might be to small")
                LD <- lapply(LD, function(x) {
                    x[is.na(x[, 1]), 1] <- max(newdata[, dose])
                    return(x)
                })


            }

        }
        if (is.list(LD)) {
            CI_LD <- lapply(LD, function(x) {
                            mat <- matrix(quantile(x[,1], probs = probs),
                                          nrow = 1)
                            colnames(mat) <- names(quantile(x[,1], probs = probs))
                            return(mat)
                        })
        } else {
            CI_LD <- matrix(quantile(LD[,1], probs = probs), nrow = 1)
            colnames(CI_LD) <- names(quantile(LD[,1], probs = probs))
        }
        return(list(ci = CI_LD, bootsrap_LDs = LD))
    }

    ## now we need to obtain y = X'beta on a fine grid for all different
    ## coefficients and then obtain the maximum for each group and the LDXX
    ## values. From these values we then are able to obtain the CIs by using the
    ## empirical quantiles (either of LDXX directly or for the difference
    ## LDXX(group1) - LDXX(group2)
    X1 <- predict(object$model, newdata = newdata1, type = "lpmatrix")
    pred1 <- responsefct(tcrossprod(X1, rbeta))
    CI_LD1 <- compute_CIs(pred1, probs, newdata1, lethal.dose)
    LD1 <- CI_LD1$bootsrap_LDs
    CI_LD1 <- CI_LD1$ci

    ## merge LD1.pe with CI
    CI <- vector("list", 3)
    CI[[1]] <- append.CI(LD1.pe, CI_LD1)

    if (!is.null(newdata2)) {
        X2 <- predict(object$model, newdata = newdata2, type = "lpmatrix")
        pred2 <- responsefct(tcrossprod(X2, rbeta))
        CI_LD2 <- compute_CIs(pred2, probs, newdata2, lethal.dose)
        LD2 <- CI_LD2$bootsrap_LDs
        CI_LD2 <- CI_LD2$ci

        ## merge LD2.pe with CI
        CI[[2]] <- append.CI(LD2.pe, CI_LD2)

        LD_diff <- ld.diff(LD1, LD2)

        if (is.list(CI_LD1)) {
            CI_Diff <- lapply(LD_diff, function(x) {
                            mat <- matrix(quantile(x, probs = probs),
                                          nrow = 1)
                            colnames(mat) <- names(quantile(x, probs = probs))
                            return(mat)
                        })
            names(CI_Diff) <- names(CI_LD1)
        } else {
            CI_Diff <- matrix(quantile(LD_diff, probs = probs), nrow = 1)
            colnames(CI_Diff) <- names(quantile(LD_diff, probs = probs))
        }

        ## merge LD_Diff.pe with CI
        CI[[3]] <- append.CI(LD_diff.pe, CI_Diff)

    }

    CI <- combine.results(CI)
    attr(CI, "values") <- lethal.dose
    attr(CI, "model") <- object
    class(CI) <- "LDconfint"
    return(CI)
}

