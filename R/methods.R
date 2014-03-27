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



################################################################################

plot.LDconfint <- function(x, xlab = NULL,
                           upper_log = "y", upper_ylab = "total cell count",
                           col = c("black", "red"), lty = NULL,
                           legend = "topright", max.shift = 0.1,
                           mar  = c(4, 9.3, 2, 0.1),
                           ...) {

    LDmod <- attr(x, "model")

    outcome <- LDmod$variables$outcome
    dose <- LDmod$variables$dose
    groups <- LDmod$variables$groups
    data <- LDmod$data

    if (is.null(xlab))
        xlab <- dose

    ## get LD values
    lds <- LD(LDmod)
    n.lds <- length(attr(lds, "values"))

    layout(matrix(c(1, 2), ncol = 1), heights = c(2, 1))

    plot(LDmod, xlab = xlab, log = upper_log, ylab = upper_ylab,
         col = col, lty = lty, legend = legend, mar = mar, ...)

    op <- par(mar = mar, las  = 1)

    ny <- ifelse(n.lds > 1, nrow(x[[1]]), nrow(x))

    if (n.lds > 1) {
        main <- "95 % confidence intervals for lethal doses"
    } else {
        main <- "95 % confidence intervals for lethal dose"
    }

    plot(1:ny,
         xlim = range(x),
         ylim = c(0.75, ny + 0.25),
         type = "n", xlab = xlab, ylab = "", yaxt = "n",
         main = main, ...)

    if (n.lds > 1) {
        inc <- seq(-max.shift, max.shift, length = n.lds)
        for (i in 1:n.lds) {
            points(c(x[[i]][, 1]), (ny:1) + inc[i], pch = 20)
            for (j in 1:ny)
                lines(x[[i]][j, -1], rep(ny - j + 1, 2) + inc[i])
        }
    } else {
        points(c(x[, 1]), ny:1, pch = 20)
        for (j in 1:ny)
            lines(x[j, -1], rep(ny - j + 1, 2))
    }

    abline(v = 0, lty = "dashed")
    ## labels  y-axis
    if (n.lds > 1) {
        txt <- rownames(x[[1]])
    } else {
        txt <- rownames(x)
    }
    prefix <- ifelse(!is.null(groups), paste0(groups, ": "), "")
    axis(2, at = ny:1, labels = paste0(prefix, txt))
    par(op)
}


draw_lds <- function(object, col, lty, group = NULL) {
    if (is.null(group))
        group <- 1

    if (!is.matrix(object)) {
        for (i in 1:length(attr(object, "values"))) {
            lines(rep(object[[i]][group, 1], 2),
                  c(1e-10, object[[i]][group, 2]),
                  col = col, lty = lty[i])
        }
    } else {
        lines(rep(object[group, 1], 2), c(1e-10, object[group, 2]),
              col = col, lty = lty[1])
    }

}

plot.LD <- function(x, xlab = NULL,
                    log = "y", ylab = "total cell count",
                    col = c("black", "red"), lty = NULL,
                    legend = "topright", mar  = c(4, 9.3, 2, 0.1), ...) {

    outcome <- x$variables$outcome
    dose <- x$variables$dose
    groups <- x$variables$groups
    data <- x$data

    ## get LD values
    lds <- LD(x)
    n.lds <- length(attr(lds, "values"))

    if (is.null(xlab))
        xlab <- dose

    if (is.null(lty)) {
        lty <- 1:n.lds
    } else {
        if (length(lty) == 1)
            lty <- rep(lty, n.lds)
    }

    if (length(lty) != n.lds)
        stop("lty must be a vector of length 1 or",
             " equal to the number of LD values")

    if (!is.null(groups)) {
        cols <- col[data[, groups]]
    } else {
        cols <- col[1]
    }

    op <- par(mar = mar, las  = 1)
    fm <- as.formula(paste(outcome, dose, sep = "~"))
    plot(fm, data = data, pch = 20,
         xlab = xlab, ylab = "", log = log,
         col = cols, ...)
    mtext(ylab, side = 2, line = 5, las = 0)

    if (!is.null(groups)) {
        p1 <- predict(x, group = 1)
        p2 <- predict(x, group = 2)
        lines(get_x(p1), p1, lwd = 2, col = col[1])
        lines(get_x(p2), p2, lwd = 2, col = col[2])
        draw_lds(lds, col = col[1], lty = lty, group = 1)
        draw_lds(lds, col = col[2], lty = lty, group = 2)
    } else {
        p1 <- predict(x)
        lines(get_x(p1), p1, lwd = 2, col = col[1])
        draw_lds(lds, col = col[1], lty = lty)
    }

    if (legend != "none") {
        if (!is.null(groups)) {
            txt <- paste0(", ", groups, ": ",
                          rep(levels(data[, groups]), each = n.lds))
        } else {
            txt <- ""
        }

        txt <- paste0(rep(paste0("LD", attr(lds, "values")),
                          ifelse(is.null(groups), 1, 2)),
                      txt)
        legend(legend, legend = txt, bty = "n",
               title = "Lethal dose(s)",
               col = rep(col, each = n.lds),
               lty = rep(lty, ifelse(is.null(groups), 1, 2)))
    }
    par(op)
    #layout(matrix(c(1, 1), ncol = 1))
}
