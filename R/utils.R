#' Helper to apply base R biplot() scaling and optional varimax rotation
#'
#' @param sco Raw observation scores matrix (2 columns).
#' @param rot Raw variable loadings/rotation matrix (2 columns).
#' @param sdev Full vector of component standard deviations.
#' @param n Number of observations.
#' @param pc Integer vector of length 2 giving the plotted principal components.
#' @param scale Numeric scaling exponent between 0 and 1.
#' @param pc.biplot Logical; additional scaling of variables.
#' @param varimax.rotate Logical; apply varimax rotation to scores.
#' @param varimax.args Named list of arguments passed to `stats::varimax()`.
#' @return A list with scaled `sco` and `rot` matrices, plus updated `sdev`.
#' @noRd
.bb_scale_pca <- function(sco, rot, sdev, n, pc, scale, pc.biplot,
                          varimax.rotate, varimax.args = list()) {
    vm_opt <- .bb_merge_args(list(normalize = TRUE, eps = 1e-5), varimax.args, "varimax.args")
    sdev_out <- sdev
    sdev_pc <- sdev[pc]

    if (isTRUE(varimax.rotate)) {
        rn_sco <- rownames(sco)
        cn_sco <- colnames(sco)
        rn_rot <- rownames(rot)
        cn_rot <- colnames(rot)
        rawl <- t(t(rot) * sdev_pc)
        vm <- do.call(varimax, c(list(x = rawl), vm_opt))
        if (!is.null(vm$rotmat)) {
            sdev_pc_safe <- ifelse(sdev_pc > 0, sdev_pc, 1)
            rot_l <- rawl %*% vm$rotmat
            s_rot <- sqrt(colSums(rot_l^2))
            s_rot_safe <- ifelse(s_rot > 0, s_rot, 1)
            sco <- t(t(t(t(sco) / sdev_pc_safe) %*% vm$rotmat) * s_rot)
            rot <- t(t(rot_l) / s_rot_safe)
            rownames(sco) <- rn_sco
            colnames(sco) <- cn_sco
            rownames(rot) <- rn_rot
            colnames(rot) <- cn_rot
            sdev_pc <- s_rot
            sdev_out[pc] <- s_rot
        }
    }

    lam <- sdev_pc * sqrt(n)
    if (scale != 0) lam <- lam^scale else lam <- rep(1, length(pc))
    if (isTRUE(pc.biplot)) lam <- lam / sqrt(n)
    lam_safe <- ifelse(lam > 0, lam, 1)
    sco <- t(t(sco) / lam_safe)
    rot <- t(t(rot) * lam)
    list(sco = sco, rot = rot, sdev = sdev_out)
}

#' Compute label rotation angles and text justification
#'
#' @param x A two-column numeric matrix of coordinates.
#' @param valign Vertical alignment (`0`, `0.5`, or `1`).
#' @param ran_adj Logical; whether to randomise justification.
#' @return A three-column numeric matrix (`angle`, `adj1`, `adj2`).
#' @noRd
.bb_lab_rotation <- function(x, valign, ran_adj) {
    x_an <- ifelse(x[, 1] == 0 & x[, 2] == 0, 0, atan(x[, 2] / x[, 1]) * 180 / pi)
    ad1 <- rep(0, times = nrow(x))
    ad2 <- rep(valign, times = nrow(x))

    left <- which(x[, 1] < 0)
    if (length(left) > 0) ad1[left] <- 1

    if (isTRUE(ran_adj)) {
        ad1 <- sample(c(0, 0.5, 1), length(ad1), replace = TRUE)
        ad2 <- sample(c(0, 0.5, 1), length(ad2), replace = TRUE)
    }

    cbind(x_an, ad1, ad2)
}

#' Compute covariance or Hotelling-t confidence ellipse coordinates for a group
#'
#' @param x A data frame with at least two numeric coordinate columns.
#' @param method Ellipse method: `"cov"` (bivariate normal / chi-square) or
#'   `"t"` (small-sample F / Hotelling's T-squared).
#' @param level Confidence level in `(0, 1)` (default `0.95`).
#' @return A data frame with columns `x` and `y`, or `NULL` if fewer than 3 rows.
#' @noRd
.bb_ellipse <- function(x, method = "cov", level = 0.95) {
    pts <- as.matrix(x[, 1:2, drop = FALSE])
    pts <- pts[is.finite(pts[, 1]) & is.finite(pts[, 2]), , drop = FALSE]
    n_pts <- nrow(pts)
    if (n_pts < 3) return(NULL)

    cen <- colMeans(pts)
    S <- cov(pts)
    if (anyNA(S)) return(NULL)

    eig <- eigen(S, symmetric = TRUE)
    vals <- pmax(eig$values, 0)
    r <- if (method == "t") {
        sqrt(2 * (n_pts - 1) / (n_pts - 2) * qf(level, df1 = 2, df2 = n_pts - 2))
    } else {
        sqrt(qchisq(level, df = 2))
    }

    t_seq <- seq(0, 2 * pi, length.out = 100)
    unit_circle <- rbind(cos(t_seq), sin(t_seq))
    ell <- t(eig$vectors %*% (sqrt(vals) * r * unit_circle))

    data.frame(x = ell[, 1] + cen[1], y = ell[, 2] + cen[2])
}

#' Validate and merge user styling options with defaults
#'
#' @param defaults Named list of default options.
#' @param user_args User-supplied list (or `NULL`).
#' @param arg_name Argument name for error messages.
#' @return Merged named list.
#' @noRd
.bb_merge_args <- function(defaults, user_args, arg_name) {
    if (is.null(user_args)) return(defaults)
    if (!is.list(user_args)) {
        stop("Error! ", arg_name, " must be a named list")
    }
    if (length(user_args) == 0) return(defaults)

    nms <- names(user_args)
    if (is.null(nms) || any(!nzchar(nms))) {
        stop("Error! All elements of ", arg_name, " must be named")
    }
    unknown <- setdiff(nms, names(defaults))
    if (length(unknown) > 0) {
        stop("Error! Unknown option(s) in ", arg_name, ": ", paste(unknown, collapse = ", "))
    }

    for (nm in nms) {
        defaults[nm] <- list(user_args[[nm]])
    }
    defaults
}

#' Resolve principal component indices with legacy pc1/pc2 deprecation shim
#'
#' @param pc Integer/numeric vector of length 2.
#' @param dots Captured `...` list.
#' @return List with `pc` (integer vector of length 2) and cleaned `dots`.
#' @noRd
.bb_resolve_pc <- function(pc, dots) {
    has_pc1 <- "pc1" %in% names(dots)
    has_pc2 <- "pc2" %in% names(dots)
    if (has_pc1 || has_pc2) {
        warning("Arguments 'pc1'/'pc2' are deprecated in betterbiplot 2.0.0; use 'pc = c(pc1, pc2)' instead.",
                call. = FALSE)
        if (has_pc1) {
            pc[1] <- dots$pc1
            dots$pc1 <- NULL
        }
        if (has_pc2) {
            pc[2] <- dots$pc2
            dots$pc2 <- NULL
        }
    }
    if (!is.numeric(pc) || length(pc) != 2 || anyNA(pc) ||
        any(pc != as.integer(pc)) || any(pc < 1)) {
        stop("Error! pc must be a positive integer vector of length 2")
    }
    list(pc = as.integer(pc), dots = dots)
}

#' Extract chull.args alias and deprecated v1.x arguments from dots
#'
#' @param dots Captured `...` list.
#' @param lab.args,var.args,ellipse.args,leg.args User-supplied option lists.
#' @return List with cleaned `dots`, updated option lists, and `f1_a`/`f2_a`.
#' @noRd
.bb_extract_legacy_args <- function(dots, lab.args, var.args, ellipse.args, leg.args) {
    if (is.null(lab.args)) lab.args <- list()
    if (is.null(var.args)) var.args <- list()
    if (is.null(ellipse.args)) ellipse.args <- list()
    if (is.null(leg.args)) leg.args <- list()

    # chull.args is a supported non-deprecated alias for ellipse.args
    if ("chull.args" %in% names(dots)) {
        hull_defaults <- list(
            method = "cov", level = 0.95, autolim = TRUE,
            nofill = FALSE, alpha = 0.25, border.alpha = 0.9, lwd = 2
        )
        hull_valid <- .bb_merge_args(hull_defaults, dots$chull.args, "chull.args")
        for (nm in names(dots$chull.args)) {
            if (!(nm %in% names(ellipse.args))) {
                ellipse.args[nm] <- list(hull_valid[[nm]])
            }
        }
        dots$chull.args <- NULL
    }

    lists <- list(
        lab.args = lab.args,
        var.args = var.args,
        ellipse.args = ellipse.args,
        leg.args = leg.args
    )

    legacy_map <- list(
        col.labd     = c("lab.args", "col"),
        cex.labd     = c("lab.args", "cex"),
        font.labd    = c("lab.args", "font"),
        lab_rotation = c("lab.args", "rotate"),
        ran_adj      = c("lab.args", "ran.adj"),
        valign       = c("lab.args", "valign"),
        vname        = c("var.args", "vname"),
        whichv       = c("var.args", "which"),
        expand       = c("var.args", "expand"),
        arrow.len    = c("var.args", "arrow.len"),
        lwd.v        = c("var.args", "lwd"),
        col.v        = c("var.args", "col"),
        col.labv     = c("var.args", "col.lab"),
        cex.labv     = c("var.args", "cex.lab"),
        font.labv    = c("var.args", "font.lab"),
        axes.v       = c("var.args", "axes"),
        autolim      = c("ellipse.args", "autolim"),
        nofill       = c("ellipse.args", "nofill"),
        lwd.e        = c("ellipse.args", "lwd"),
        title.leg    = c("leg.args", "title"),
        cex.leg      = c("leg.args", "cex"),
        horiz        = c("leg.args", "horiz"),
        lpx          = c("leg.args", "x"),
        lpy          = c("leg.args", "y")
    )

    for (old_name in intersect(names(dots), names(legacy_map))) {
        target <- legacy_map[[old_name]]
        list_name <- target[1]
        key <- target[2]
        warning("Argument '", old_name, "' is deprecated in betterbiplot 2.0.0; use '",
                list_name, " = list(", key, " = ...)' instead.", call. = FALSE)
        if (!(key %in% names(lists[[list_name]]))) {
            lists[[list_name]][key] <- list(dots[[old_name]])
        }
        dots[[old_name]] <- NULL
    }

    for (old_ell in intersect(names(dots), c("angle", "se"))) {
        warning("Argument '", old_ell, "' is deprecated in betterbiplot 2.0.0; ellipses now use 'ellipse.args = list(method = \"cov\", level = 0.95)' (or method = \"t\").",
                call. = FALSE)
        dots[[old_ell]] <- NULL
    }

    f1_a <- NULL
    if ("f1.a" %in% names(dots)) {
        warning("Argument 'f1.a' is deprecated in betterbiplot 2.0.0; pass a wrapper function such as 'f1 = function() fn(...)' instead.",
                call. = FALSE)
        f1_a <- dots$f1.a
        dots$f1.a <- NULL
    }

    f2_a <- NULL
    if ("f2.a" %in% names(dots)) {
        warning("Argument 'f2.a' is deprecated in betterbiplot 2.0.0; pass a wrapper function such as 'f2 = function() fn(...)' instead.",
                call. = FALSE)
        f2_a <- dots$f2.a
        dots$f2.a <- NULL
    }

    list(
        dots = dots,
        lab.args = lists$lab.args,
        var.args = lists$var.args,
        ellipse.args = lists$ellipse.args,
        leg.args = lists$leg.args,
        f1_a = f1_a,
        f2_a = f2_a
    )
}
