################################################################################
## compute lethal dose model

LD <- lethal.dose <- function(...)
    UseMethod("LD")

LD.default <- function(formula, groups = NULL, experiment = NULL,
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
                trafo = dose_trafo)
    RET$variables <- list(outcome = as.character(outcome),
                          dose = as.character(dose),
                          groups = groups,
                          experiment = experiment)
    class(RET) <- c("LD")
    RET$lethal.dose <- LD(RET, lethal.dose)
    return(RET)
}

LD.LD <- function(object, lethal.dose = c(50, 10), group = NULL, ...) {
    cll <- match.call()
    if (is.null(cll[["lethal.dose"]]))
        return(object$lethal.dose)

    pr <- predict(object, group)
    RET <- LD(pr, dose = get_x(pr), lethal.dose = lethal.dose, ...)
    return(RET)
}

LD.LDpred <- function(pred, dose = NULL,
                                            lethal.dose = c(50, 10), ...) {
    ## pred can be a vector or a matrix of predictions
    if (is.null(dose))
        dose <- get_x(pred)

    LD <- matrix(NA, nrow = ncol(pred), ncol = 2)
    colnames(LD) <- c("dose", "predicted count")
    rownames(LD) <- paste0("group:", colnames(pred))
    if (length(lethal.dose) > 1) {
        LD <- lapply(1:length(lethal.dose), function(i) LD)
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

lethal.dose.default <- lethal.dose.LD <- lethal.dose.LDpred <- function(...)
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
        return(pr)
    }

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
        if (!is.null(variables$experiment))
            pr <- tapply(pr, newdata[, variables$dose], mean)
    } else {
        pr <- sapply(newdata, predict, object = object$model, type = type, ...)
        ## average over experiments (if available)
        if (!is.null(variables$experiment))
            pr <- sapply(1:g, function(i)
                         tapply(pr[, i], newdata[[i]][, variables$dose], mean))
        colnames(pr) <- levels(object$data[, variables$groups])
    }

    ## fixme se.fit = TRUE will kill this!

    RET <- pr
    attr(RET, "newdata") <- newdata
    class(RET) <- "LDpred"
    return(RET)
}

################################################################################
## Inference based on Wood p.261

## parametric bootstrap that takes the uncerteinty of smoothing-parameter
## selection into account

confint <- function(object, level, ...)
    UseMethod("confint")

confint.LD <- function(object, level = 0.95, lethal.dose = c(50, 10),
                       outcome = "value", dose = "time", B1 = 20, B2 = 100,
                       newdata1 = NULL, newdata2 = NULL,
                       myapply = mclapply, ...) {

    data <- object$data
    fm <- object$model$formula
    responsefct <- object$model$family$linkinv

    ## <FIXME> What about experiments with more than two groups?
    if (is.null(newdata1))
        newdata1 <- make_newdata(object, group = 1)
    if (is.null(newdata2))
        newdata2 <- make_newdata(object, group = 2)

    probs <- c((1 - level)/2, 1 - (1 - level)/2)

    ## Inner parametric bootstrap (draw coefficients from the multivariate
    ## normal distribution they are assumed to follow as they are REML
    ## estimates)
    if (B2 > 0)
        rbeta <- rbind(mvrnorm(n = B2, mu = coef(object$model),
                               Sigma = vcov(object$model)))
    else
        rbeta <- coef(object$model)

    refit_model <- function(counter, object, data, fm) {
        ## Outer parametric bootstrap (draw new observations from fitted model):
        ## - We assume that the model holds, i.e., the data follows a negative
        ##   binomial distribution and the coefficients are consistent.
        ## - We simpulate data from the negative binomial distribution with the
        ##   conditional means as fitted bevor.

        ## <FIXME>: Add alternative sample functions for different families
        res <- rnegbin(rep(1, length(fitted(object))), mu = fitted(object),
                       theta = object$family$getTheta())
        ## make shure that new outcome is really a count variable > 0
        y <- round(res)
        y[y < 0] <- 0
        data.BS <- data
        data.BS[, outcome] <- y


        ## <FIXME>: Add alternative families:

        ## estimate smoothing parameters and overdispersion parameter theta from
        ## simulated data
        mod.BS <- gam(fm, data = data.BS, family = negbin(theta = c(0.5, 10)))
        ## now estimate coefficients from original data but with the above
        ## smoothing parameters
        mod <- gam(fm, data = data, sp = mod.BS$sp,
                   family = negbin(theta = mod.BS$family$getTheta()))

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
        res <- myapply(1:(B1 - 1), refit_model, object = object$model,
                       data = data, fm = fm, ...)
        rbeta <- rbind(rbeta, do.call(rbind, res))
    }

    compute_CIs <- function(pred, probs, newdata) {


#####
###        LD <- matrix(NA, nrow = ncol(pred), ncol = 2)
###        if (length(lethal.dose) > 1) {
###            LD <- list(LD, LD)
###            names(LD) <- paste("LD", lethal.dose, sep = "")
###        }
###         for (val in lethal.dose) {
###           j <- paste("LD", val, sep = "")
###           for (i in 1:ncol(pred))
###               LD[[j]][i, ] <- compute_LD(pred[, i], ld_val = val,
###                                        dose = newdata[, dose])
###         }
###
#####

        if (is.list(LD)) {
            CI_LD <- sapply(LD, function(x)
                            quantile(x[,1], probs = probs))
        } else {
            CI_LD <- quantile(LD[,1], probs = probs)
        }
        return(list(LD = LD, CI_LD = CI_LD))
    }

    ## now we need to obtain y = X'beta on a fine grid for all different
    ## coefficients and then obtain the maximum for each group and the LDXX
    ## values. From these values we then are able to obtain the CIs by using the
    ## empirical quantiles (either of LDXX directly or for the difference
    ## LDXX(group1) - LDXX(group2)
    X1 <- predict(object$model, newdata = newdata1, type = "lpmatrix")
    pred1 <- responsefct(tcrossprod(X1, rbeta))
    LD1 <- compute_CIs(pred1, probs, newdata1)
    CI_LD1 <- LD1$CI_LD
    LD1 <- LD1$LD

    if (!is.null(newdata2)) {
        X2 <- predict(object$model, newdata = newdata2, type = "lpmatrix")
        pred2 <- responsefct(tcrossprod(X2, rbeta))
        LD2 <- compute_CIs(pred2, probs, newdata2)
        CI_LD2 <- LD2$CI_LD
        LD2 <- LD2$LD

        if (is.list(LD1)) {
            LD_diff <- lapply(1:length(LD1), function(i)
                              LD1[[i]][, 1] - LD2[[i]][, 1])
            CI_Diff <- sapply(LD_diff, function(x)
                              quantile(x, probs = probs))
            colnames(CI_Diff) <- colnames(CI_LD2)
        } else {
            CI_Diff <- quantile(LD1[, 1] - LD2[, 1], probs = probs)
        }
    } else {
        pred2 <- NA
        CI_LD2 <- NA
        LD2 <- NA
    }

    return(list(CI_LD1 = CI_LD1, CI_LD2 = CI_LD2, CI_Diff = CI_Diff,
                pred1 = pred1, pred2 = pred2, LD1 = LD1, LD2 = LD2))
}

