#' Helper to apply base R biplot() scaling and optional varimax rotation
#'
#' @param sco Raw observation scores matrix (2 columns).
#' @param rot Raw variable loadings/rotation matrix (2 columns).
#' @param sdev Full vector of component standard deviations.
#' @param n Number of observations.
#' @param pc1,pc2 Indices of plotted principal components.
#' @param scale Numeric scaling exponent between 0 and 1.
#' @param pc.biplot Logical; additional scaling of variables.
#' @param varimax.rotate Logical; apply varimax rotation to scores.
#' @param varimax.args Named list of arguments passed to `stats::varimax()`.
#' @return A list with scaled `sco` and `rot` matrices.
#' @noRd
.bb_scale_pca <- function(sco, rot, sdev, n, pc1, pc2, scale, pc.biplot,
                          varimax.rotate, varimax.args = list()) {
    vm_opt <- .bb_merge_args(list(normalize = TRUE, eps = 1e-5), varimax.args, "varimax.args")
    lam <- sdev[c(pc1, pc2)] * sqrt(n)
    if (scale != 0) lam <- lam^scale else lam <- 1
    if (isTRUE(pc.biplot)) lam <- lam / sqrt(n)
    sco <- t(t(sco) / lam)
    rot <- t(t(rot) * lam)
    if (isTRUE(varimax.rotate)) {
        rawl <- rot %*% diag(sdev[c(pc1, pc2)])
        vm <- do.call(varimax, c(list(x = rawl), vm_opt))
        sco <- scale(sco) %*% vm$rotmat
    }
    list(sco = sco, rot = rot)
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
    pc1 <- pc_res$pc[1]
    pc2 <- pc_res$pc[2]
    sdev <- sqrt(pmax(x$CA$eig, 0))
    n <- nrow(x$CA$u)
    raw_sco <- t(t(x$CA$u[, c(pc1, pc2), drop = FALSE]) * (sdev[c(pc1, pc2)] * sqrt(n - 1)))
    raw_rot <- x$CA$v[, c(pc1, pc2), drop = FALSE]
    scaled <- .bb_scale_pca(raw_sco, raw_rot, sdev, n, pc1, pc2, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = sdev, pc = pc_res$pc), pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.svyprcomp <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                                varimax.args = list(), pc.biplot = FALSE, group3, ...) {
    if (is.null(x$x)) {
        stop("Error! svyprcomp object does not contain scores. Re-run svyprcomp() with scores = TRUE")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc1 <- pc_res$pc[1]
    pc2 <- pc_res$pc[2]
    n <- nrow(x$x)
    raw_sco <- x$x[, c(pc1, pc2), drop = FALSE]
    raw_rot <- x$rotation[, c(pc1, pc2), drop = FALSE]
    scaled <- .bb_scale_pca(raw_sco, raw_rot, x$sdev, n, pc1, pc2, scale,
                            pc.biplot, varimax.rotate, varimax.args)

    base_args <- list(x = scaled$sco, y = scaled$rot, xsd = x$sdev, pc = pc_res$pc)
    if (!missing(group3) && identical(group3, "scaled")) {
        if (is.null(x$weights)) {
            stop("Error! svyprcomp object does not contain weights")
        }
        base_args$group3 <- structure(as.numeric(x$weights), class = "bb_svy_weights")
    } else if (!missing(group3)) {
        base_args$group3 <- group3
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
    pc1 <- pc_res$pc[1]
    pc2 <- pc_res$pc[2]
    sdev <- sqrt(pmax(x$eig[, 1], 0))
    n <- nrow(x$ind$coord)
    raw_sco <- as.matrix(x$ind$coord[, c(pc1, pc2), drop = FALSE])
    raw_rot <- t(t(as.matrix(x$var$coord[, c(pc1, pc2), drop = FALSE])) / sdev[c(pc1, pc2)])
    scaled <- .bb_scale_pca(raw_sco, raw_rot, sdev, n, pc1, pc2, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = sdev, pc = pc_res$pc), pc_res$dots))
}

#' @rdname bb_biplot
#' @export
bb_biplot.pca <- function(x, pc = c(1, 2), scale = 1, varimax.rotate = FALSE,
                          varimax.args = list(), pc.biplot = FALSE, ...) {
    if (is.null(x$li) || (is.null(x$c1) && is.null(x$co)) || is.null(x$eig)) {
        stop("Error! pca object does not contain expected ade4::dudi.pca components")
    }
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc1 <- pc_res$pc[1]
    pc2 <- pc_res$pc[2]
    sdev <- sqrt(pmax(x$eig, 0))
    n <- nrow(x$li)
    raw_sco <- as.matrix(x$li[, c(pc1, pc2), drop = FALSE])
    raw_rot <- if (!is.null(x$c1)) {
        as.matrix(x$c1[, c(pc1, pc2), drop = FALSE])
    } else {
        t(t(as.matrix(x$co[, c(pc1, pc2), drop = FALSE])) / sdev[c(pc1, pc2)])
    }
    scaled <- .bb_scale_pca(raw_sco, raw_rot, sdev, n, pc1, pc2, scale,
                            pc.biplot, varimax.rotate, varimax.args)
    do.call(bb_biplot.default, c(list(x = scaled$sco, y = scaled$rot, xsd = sdev, pc = pc_res$pc), pc_res$dots))
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
    pc1 <- pc_res$pc[1]
    pc2 <- pc_res$pc[2]
    sdev <- if (!is.null(x$values)) sqrt(pmax(x$values, 0)) else NULL
    sco <- as.matrix(x$scores[, c(pc1, pc2), drop = FALSE])
    rot <- unclass(x$loadings)[, c(pc1, pc2), drop = FALSE]
    do.call(bb_biplot.default, c(list(x = sco, y = rot, xsd = sdev, pc = pc_res$pc), pc_res$dots))
}
