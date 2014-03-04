################################################################################

## (non-exported) helper functions

## for condidence intervals of the fitted splines
get_tval <- function(object, alpha = 0.05) {
    residual.df <- length(object$y) - sum(object$edf)
    tVal <- qt(1 - (alpha/2), residual.df)
    tVal
}


################################################################################

## helper function to create new data sets for prediction and compuation of LD
## values

make_newdata <- function(object, group = NULL) {
    if (!is.null(group) && length(group) != 1)
        stop("group must be an integer or a character string")

    ## extract info from LGmodel object
    variables <- object$variables

    ## extract data from fitted model
    newdata <- object$data

    ## get indices for variables
    idx <- !sapply(variables, is.null)
    terms <- lapply(variables[idx], function(x) grep(x, names(newdata)))
    terms <- lapply(terms, function(i) names(newdata)[i])

    ## get doses
    dose <- newdata[, terms$dose]

    ## set grid for dose
    dose <- seq(min(dose), max(dose), by = 1)

    ## remove response and dose from data
    newdata[, terms$outcome] <- NULL
    newdata[, terms$dose] <- NULL

    ## subset data to group
    if (!is.null(variables$groups)) {
        if (is.null(group)) {
            warning("No group selected, prediction might be misleading")
        } else {
            if (is.numeric(group)) {
                group <- levels(newdata[, terms$groups])[group]
            } # else {  ## make sure that the group exists
              #   group <- levels(newdata[, terms$groups])
              # }
            newdata[, terms$groups] <- group
        }
    }
    newdata <- unique(newdata)
    newdata <- newdata[rep(1:nrow(newdata), each = length(dose)),]
    newdata[, as.character(variables$dose)] <- dose
    return(newdata)
}
