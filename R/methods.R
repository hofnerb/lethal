################################################################################

## Methods for LD objects

print.LD <- function(x, ...) {
    cat("\n\tFitted model:\n")
    print(x$model, ...)
    cat("\n\n\tLethal dose(s):\n")
    print(LD(x))
}

print.LDmatrix <- function(x, ...) {
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

summary.LD <- function(object, ...) {
    cat("\n\tFitted model:\n")
    print(summary(object$model, ...))
    cat("\n\n\tLethal dose(s):\n")
    print(LD(object))
}

anova.LD <- function(object, ...) {
    anova(object$model, ...)
}

get_x <- function(object, ...)
    UseMethod("get_x")

get_x.default <- function(object, ...)
    1:length(object)

## <FIXME>
#get_x.LDpred <- function(object, ...)
#    #object$newdata[, DOSE]
