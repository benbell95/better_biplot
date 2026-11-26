#' @rdname bb_biplot
#' @export
bb_biplot.prcomp <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                             varimax.args = list(), pc.biplot = FALSE, ...) {
    # Scaling code modified from base R biplot() [Copyright (C) 1995-2012 The R Core Team]
    if (is.null(x$x)) {
        stop("Error! prcomp object does not contain scores. Re-run prcomp() with retx = TRUE")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc
    n <- nrow(x$x)
    raw_sco <- x$x[, pc, drop = FALSE]
    raw_rot <- x$rotation[, pc, drop = FALSE]
    scaled <- .bb_scale_pca(raw_sco, raw_rot, x$sdev, n, pc, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = scaled$sdev, pc = pc), pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.princomp <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                               varimax.args = list(), pc.biplot = FALSE, ...) {
    # Scaling code modified from base R biplot() [Copyright (C) 1995-2012 The R Core Team]
    if (is.null(x$scores)) {
        stop("Error! princomp object does not contain scores. Re-run princomp() with scores = TRUE")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc
    n <- x$n.obs
    raw_sco <- x$scores[, pc, drop = FALSE]
    raw_rot <- x$loadings[, pc, drop = FALSE]
    scaled <- .bb_scale_pca(raw_sco, raw_rot, x$sdev, n, pc, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = scaled$sdev, pc = pc), pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.rda <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                          varimax.args = list(), pc.biplot = FALSE, ...) {
    if (!is.null(x$CCA)) {
        stop("Error! bb_biplot.rda() only supports unconstrained PCA models (no constraints)")
    }
    if (is.null(x$CA) || is.null(x$CA$u) || is.null(x$CA$v)) {
        stop("Error! rda object does not contain unconstrained PCA components")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc
    sdev <- sqrt(pmax(x$CA$eig, 0))
    n <- nrow(x$CA$u)
    raw_sco <- t(t(x$CA$u[, pc, drop = FALSE]) * (sdev[pc] * sqrt(n - 1)))
    raw_rot <- x$CA$v[, pc, drop = FALSE]
    scaled <- .bb_scale_pca(raw_sco, raw_rot, sdev, n, pc, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = scaled$sdev, pc = pc), pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.svyprcomp <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                                varimax.args = list(), pc.biplot = FALSE, group3, ...) {
    if (is.null(x$x)) {
        stop("Error! svyprcomp object does not contain scores. Re-run svyprcomp() with scores = TRUE")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc
    n <- nrow(x$x)
    raw_sco <- x$x[, pc, drop = FALSE]
    raw_rot <- x$rotation[, pc, drop = FALSE]
    scaled <- .bb_scale_pca(raw_sco, raw_rot, x$sdev, n, pc, scale,
                            pc.biplot, varimax.rotate, varimax.args)

    base_args <- list(x = scaled$sco, y = scaled$rot, xsd = scaled$sdev, pc = pc)
    if (!missing(group3) && !is.null(group3)) {
        if (identical(group3, "scaled")) {
            if (is.null(x$weights)) {
                stop("Error! svyprcomp object does not contain weights")
            }
            base_args$group3 <- structure(as.numeric(x$weights), class = "bb_svy_weights")
        } else {
            base_args$group3 <- group3
        }
    }
    do.call(bb_biplot.default, c(base_args, pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.PCA <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                          varimax.args = list(), pc.biplot = FALSE, ...) {
    if (is.null(x$ind$coord) || is.null(x$var$coord) || is.null(x$eig)) {
        stop("Error! PCA object does not contain expected FactoMineR::PCA components")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc
    sdev <- sqrt(pmax(x$eig[, 1], 0))
    sdev_safe <- ifelse(sdev[pc] > 0, sdev[pc], 1)
    n <- nrow(x$ind$coord)
    raw_sco <- as.matrix(x$ind$coord[, pc, drop = FALSE])
    raw_rot <- t(t(as.matrix(x$var$coord[, pc, drop = FALSE])) / sdev_safe)
    scaled <- .bb_scale_pca(raw_sco, raw_rot, sdev, n, pc, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = scaled$sdev, pc = pc), pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.pca <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                          varimax.args = list(), pc.biplot = FALSE, ...) {
    if (is.null(x$li) || (is.null(x$c1) && is.null(x$co)) || is.null(x$eig)) {
        stop("Error! pca object does not contain expected ade4::dudi.pca components")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc
    sdev <- sqrt(pmax(x$eig, 0))
    sdev_safe <- ifelse(sdev[pc] > 0, sdev[pc], 1)
    n <- nrow(x$li)
    raw_sco <- as.matrix(x$li[, pc, drop = FALSE])
    raw_rot <- if (!is.null(x$c1)) {
        as.matrix(x$c1[, pc, drop = FALSE])
    } else {
        t(t(as.matrix(x$co[, pc, drop = FALSE])) / sdev_safe)
    }
    scaled <- .bb_scale_pca(raw_sco, raw_rot, sdev, n, pc, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = scaled$sdev, pc = pc), pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.principal <- function(x, pc = c(1, 2), ...) {
    if (is.null(x$scores)) {
        stop("Error! principal object does not contain scores. Re-run principal() with scores = TRUE")
    }
    if (is.null(x$loadings)) {
        stop("Error! principal object does not contain loadings")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc
    sdev <- if (!is.null(x$Vaccounted) && is.matrix(x$Vaccounted) &&
                "Proportion Var" %in% rownames(x$Vaccounted)) {
        pv <- as.numeric(x$Vaccounted["Proportion Var", ])
        rem <- max(0, 1 - sum(pv))
        sqrt(pmax(c(pv, rem), 0))
    } else if (!is.null(x$values)) {
        sqrt(pmax(x$values, 0))
    } else {
        NULL
    }
    sco <- as.matrix(x$scores[, pc, drop = FALSE])
    rot <- unclass(x$loadings)[, pc, drop = FALSE]
    do.call(bb_biplot.default, c(list(x = sco, y = rot, xsd = sdev, pc = pc), pc_res$dots))
}
