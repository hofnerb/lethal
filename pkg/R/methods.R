################################################################################

## Methods for LD objects

print.LD <- function(x, ...) {
    cat("\n\tFitted model:\n")
    print(x$model, ...)
    cat("\n\n")
    print(LD(x))
}

print_LDvals <- function(x, ...) {
    if (length(attr(x, "values")) == 1) {
        cat("LD", attr(x, "values"), "\n", sep = "")
        class(x) <- "matrix"
        attr(x, "values") <- NULL
        print(x)
        cat("\n")
    } else {
        for (i in 1:length(attr(x, "values"))) {
            cat(" LD", attr(x, "values")[i], "\n", sep = "")
            print(x[[i]])
            cat("\n")
        }
    }
}

print.LDmatrix <- function(x, ...) {
    cat("\tLethal dose(s):\n")
    print_LDvals(x, ...)
}

print.LDconfint <- function(x, ...) {
    cat("\tLethal dose(s) with confidence intervals:\n")
    print_LDvals(x)
}

summary.LD <- function(object, ...) {
    cat("\n\tFitted model:\n")
    print(summary(object$model, ...))
    cat("\n\n")
    print(LD(object))
}

coef.LD <- function(object, ...) {
    coef(object$model, ...)
}

anova.LD <- function(object, ...) {
    anova(object$model, ...)
}

get_x <- function(object, ...)
    UseMethod("get_x")

get_x.default <- function(object, ...)
    1:length(object)

get_x.LDpred <- function(object, ...) {
    newdata <- attr(object, "newdata")
    if (!is.data.frame(newdata))
        newdata <- newdata[[1]]
    newdata[, attr(object, "variables")$dose]
}
