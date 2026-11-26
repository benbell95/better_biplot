### Ben's Better Biplot
### Plot function for creating better looking biplots in R (using only base R)
### Copyright 2025-2027, Benjamin Bell.
### Code and updates: https://github.com/benbell95/better_biplot

#' Better Biplot
#'
#' Creates customisable biplots for Principal Component Analysis (PCA) using
#' only base R graphics. Supports [stats::prcomp()], [stats::princomp()],
#' `vegan::rda()` (unconstrained PCA), `survey::svyprcomp()`, `FactoMineR::PCA()`,
#' `ade4::dudi.pca()`, and `psych::principal()` objects, as well as custom score
#' and loading matrices.
#'
#' @param x Observations matrix with at least 2 columns, or a supported PCA
#'   object (`"prcomp"`, `"princomp"`, `"rda"`, `"svyprcomp"`, `"PCA"`, `"pca"`,
#'   or `"principal"`).
#' @param y Variables (loadings/rotation) matrix with at least 2 columns.
#'   Automatically determined for supported PCA objects.
#' @param xsd Optional numeric vector of standard deviations for the variables,
#'   used for automatic axis variance labels. Automatically determined for
#'   supported PCA objects.
#' @param pc Integer vector of length 2 giving the principal components to plot
#'   (default `c(1, 2)`).
#' @param scale Numeric scaling exponent between `0` and `1` for PCA objects
#'   (default `1`), matching [stats::biplot()].
#' @param pc.biplot Logical; if `TRUE`, applies additional scaling of variables
#'   matching [stats::biplot()] (default `FALSE`).
#' @param varimax.rotate Logical; if `TRUE`, applies [stats::varimax()] rotation
#'   to the plotted components for supported PCA objects (default `FALSE`).
#' @param varimax.args Named list of arguments passed to [stats::varimax()]
#'   when `varimax.rotate = TRUE`:
#'   * `normalize`: Logical; perform Kaiser normalization (default `TRUE`).
#'   * `eps`: Numeric convergence tolerance (default `1e-5`).
#' @param limx Optional numeric multiplier to expand or shrink the automatic
#'   `xlim` and `ylim` values (e.g., `1.5`).
#' @param grid Logical; draw dashed reference lines through the origin (default `TRUE`).
#' @param col Optional vector of colours for observations or groups, or a
#'   single palette name from [grDevices::palette.pals()] (e.g. `"Okabe-Ito"`)
#'   or [grDevices::hcl.pals()] (e.g. `"Viridis"`, `"Dark 3"`). A palette is
#'   generated automatically if omitted.
#' @param pch Plotting symbol(s) for observations (default `21`).
#' @param cex.pt Numeric size multiplier for observation points (default `1`).
#' @param alpha Numeric point opacity in `[0, 1]` (default `0.85`).
#' @param xlab Optional x-axis label overriding the automatic label.
#' @param ylab Optional y-axis label overriding the automatic label.
#' @param axes Logical; whether to draw axes (default `TRUE`).
#' @param cex.axis Numeric size multiplier for axes and axis labels (default `1`).
#' @param cex Global size multiplier applied to plot elements (default `1`).
#' @param group Optional factor (or vector coercible to factor) of length `nrow(x)`
#'   grouping observations by colour (and symbol if multiple `pch` values are supplied
#'   without `group2`).
#' @param group2 Optional second grouping factor controlling point `pch` symbols.
#'   Requires `group`.
#' @param group3 Optional third grouping factor controlling point size.
#'   Requires `group`. The `group3` legend is condensed to show minimum, middle,
#'   and maximum symbol sizes. For `survey::svyprcomp()` objects, can also be
#'   set to `"scaled"` to scale point sizes by the survey weights
#'   (`sqrt(w / max(w))`, which does not require `group` and displays a
#'   condensed min/median/max weight legend when legends are enabled).
#' @param labd Optional character vector of observation labels. Use `NA` to
#'   omit individual labels.
#' @param lab.args Named list of label styling options:
#'   * `col`: Colour for observation labels (default `"black"`).
#'   * `cex`: Numeric size for observation labels (default `0.5`).
#'   * `font`: Font face for observation labels (`1` = normal, `2` = bold,
#'     `3` = italic, `4` = bold italic; default `1`).
#'   * `rotate`: Logical; rotate observation and variable labels radially (default `FALSE`).
#'   * `ran.adj`: Logical; randomly adjust label alignment when `rotate = TRUE` (default `FALSE`).
#'   * `valign`: Vertical alignment for rotated variable labels (`0`, `0.5`, or `1`; default `0.5`).
#'   * `pos`: Placement for unrotated observation labels: `"auto"` (deterministic
#'     radial quadrant, default), `"random"`, or integer position(s) in `1:4`.
#' @param variables Logical; whether to plot variable arrows (default `TRUE`).
#' @param circle.eq Logical; draw the circle of equilibrium contribution when
#'   `scale = 0` (default `FALSE`).
#' @param var.args Named list of variable-arrow styling options:
#'   * `vname`: Optional vector of variable names (defaults to `rownames(y)`).
#'   * `which`: Optional character vector of variable names, integer indices, or
#'     `"circle.eq"` (requires `circle.eq = TRUE` and `scale = 0`) to subset plotted variables.
#'   * `expand`: Numeric multiplier for variable arrow lengths (default `1`).
#'   * `lab.offset`: Numeric multiplier nudging variable labels past arrow tips (default `1.05`).
#'   * `arrow.len`: Length of arrow heads in inches; `0` suppresses heads (default `0.15`).
#'   * `lwd`: Line width for variable arrows (default `2`).
#'   * `col`: Colour for variable arrows and circle of equilibrium (default `"#d12631"`).
#'   * `col.lab`: Colour for variable labels (default `"black"`).
#'   * `cex.lab`: Numeric size for variable labels (default `1`).
#'   * `font.lab`: Font face for variable labels (default `1`).
#'   * `axes`: Logical; whether to draw top and right variable axes (default `TRUE`).
#' @param chull Logical; draw convex hulls around `group` clusters (default `FALSE`).
#' @param ellipse Logical; draw ellipses around `group` clusters (default `FALSE`).
#' @param centroids Controls group centroid display when `group` is supplied
#'   (default `FALSE`). Accepts `FALSE`, `TRUE` (diamond markers + labels),
#'   `"points"` (markers only), or `"labels"` (labels only).
#' @param spider Logical; draw segments connecting each observation to its
#'   `group` centroid (default `FALSE`).
#' @param ellipse.args Named list of hull and ellipse styling options (also
#'   accepted via the alias `chull.args`):
#'   * `method`: Ellipse method: `"cov"` (bivariate normal confidence ellipse,
#'     default) or `"t"` (small-sample F / Hotelling's T-squared ellipse).
#'   * `level`: Confidence level in `(0, 1)` for ellipses (default `0.95`).
#'   * `autolim`: Logical; automatically expand plot limits if ellipses exceed them (default `TRUE`).
#'   * `nofill`: Logical; if `TRUE`, draw only outlines for ellipses and convex hulls (default `FALSE`).
#'   * `alpha`: Fill opacity in `[0, 1]` for ellipses and convex hulls (default `0.25`).
#'   * `border.alpha`: Border opacity in `[0, 1]` for ellipses, convex hulls, and spider lines (default `0.9`).
#'   * `lwd`: Line width for ellipses and convex hulls (default `2`).
#' @param legend Controls legend display (default `FALSE`). Accepts `FALSE`
#'   (no legend), `TRUE` or `"all"` (legends for all supplied groups), or an
#'   integer/numeric vector of group numbers in `1:3` (e.g. `1`, `c(1, 2)`, `3`)
#'   to show only selected group legends in the given order. Requires `group`
#'   (except when using `group3 = "scaled"` with `survey::svyprcomp()`).
#' @param leg.args Named list of legend styling options:
#'   * `cex`: Numeric size multiplier for legends (default `0.75`).
#'   * `horiz`: Logical; arrange multiple group legends horizontally (default `FALSE`).
#'   * `x`: Optional x-coordinate for legend placement.
#'   * `y`: Optional y-coordinate for legend placement.
#'   * `title`: Optional character vector of legend title(s). Can match the
#'     number of shown legends positionally, or supply titles for all groups
#'     (indexed by group number). Defaults to `"Group"`, `"Group 2"`,
#'     `"Group 3"` (or `"Weight"` for `survey::svyprcomp()` weights).
#' @param f1 Optional zero-argument function (or name) called before plotting points.
#' @param f2 Optional zero-argument function (or name) called before plotting variable arrows.
#' @param ... Additional graphical arguments passed to [graphics::plot()], plus
#'   `chull.args` (alias for `ellipse.args`) and deprecated v1.x arguments.
#'
#' @section Notes:
#' `group`, `group2`, and `group3` preserve existing factor level order, so
#' setting `factor(x, levels = ...)` beforehand controls the order of legend
#' items as well as colour, symbol, and smallest-to-largest `group3` symbol
#' size mapping.
#'
#' @section Backwards Compatibility:
#' `betterbiplot` 2.0.0 streamlined the v1.x argument list by grouping fine-grained
#' styling into `lab.args`, `var.args`, `ellipse.args` (or `chull.args`), and
#' `leg.args`, and replacing `pc1`/`pc2` with `pc = c(1, 2)`. Legacy v1.x
#' argument names (such as `col.v`, `whichv`, `expand`, `lab_rotation`,
#' `angle`, `se`, `title.leg`, `cex.leg`, `lpx`, `f1.a`, `pc1`, `pc2`) are still
#' accepted with a deprecation warning.
#'
#' @return Invisibly returns an object of class `"bb_biplot"` containing
#'   `scores`, `loadings`, `lim`, `limr`, `pc`, `pvar`, `group_cols`, and
#'   `centroids` (for layering custom base R annotations).
#'
#' @examples
#' p <- stats::prcomp(iris[-5], scale. = TRUE)
#' gr <- factor(iris[, 5])
#' bp <- bb_biplot(p, scale = 0, group = gr, col = "Okabe-Ito",
#'                 circle.eq = TRUE, ellipse = TRUE,
#'                 ellipse.args = list(method = "cov", level = 0.95),
#'                 centroids = TRUE, spider = TRUE,
#'                 legend = TRUE, leg.args = list(title = "Species"),
#'                 lab.args = list(rotate = TRUE))
#'
#' @export
bb_biplot <- function(x, ...) UseMethod("bb_biplot")

#' @rdname bb_biplot
#' @export
bb_biplot.default <- function(
    x, y, xsd, pc = c(1, 2), limx, grid = TRUE, col, pch = 21, cex.pt = 1,
    alpha = 0.85, xlab, ylab, axes = TRUE, cex.axis = 1, cex = 1,
    group, group2, group3,
    labd, lab.args = list(),
    variables = TRUE, circle.eq = FALSE, var.args = list(),
    chull = FALSE, ellipse = FALSE, centroids = FALSE, spider = FALSE,
    ellipse.args = list(),
    legend = FALSE, leg.args = list(),
    f1, f2, ...
) {
    # Resolve pc, chull.args alias, and deprecated v1.x arguments
    pc_res <- .bb_resolve_pc(pc, list(...))
    pc <- pc_res$pc

    leg_res <- .bb_extract_legacy_args(
        pc_res$dots,
        lab.args = lab.args,
        var.args = var.args,
        ellipse.args = ellipse.args,
        leg.args = leg.args
    )
    plot_dots <- leg_res$dots
    f1_a <- leg_res$f1_a
    f2_a <- leg_res$f2_a

    lab_opt <- .bb_merge_args(
        list(col = "black", cex = 0.5, font = 1, rotate = FALSE, ran.adj = FALSE,
             valign = 0.5, pos = "auto"),
        leg_res$lab.args, "lab.args"
    )
    var_opt <- .bb_merge_args(
        list(vname = NULL, which = NULL, expand = 1, lab.offset = 1.05, arrow.len = 0.15,
             lwd = 2, col = "#d12631", col.lab = "black", cex.lab = 1, font.lab = 1, axes = TRUE),
        leg_res$var.args, "var.args"
    )
    ell_opt <- .bb_merge_args(
        list(method = "cov", level = 0.95, autolim = TRUE,
             nofill = FALSE, alpha = 0.25, border.alpha = 0.9, lwd = 2),
        leg_res$ellipse.args, "ellipse.args"
    )
    leg_opt <- .bb_merge_args(
        list(cex = 0.75, horiz = FALSE, x = NULL, y = NULL, title = NULL),
        leg_res$leg.args, "leg.args"
    )

    if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
        stop("Error! alpha must be a numeric value in [0, 1]")
    }
    if (!(ell_opt$method %in% c("cov", "t"))) {
        stop("Error! ellipse.args$method should be either \"cov\" or \"t\"")
    }
    if (!is.numeric(ell_opt$level) || length(ell_opt$level) != 1 || is.na(ell_opt$level) ||
        ell_opt$level <= 0 || ell_opt$level >= 1) {
        stop("Error! ellipse.args$level must be a numeric value strictly between 0 and 1")
    }
    if (!is.numeric(ell_opt$alpha) || length(ell_opt$alpha) != 1 || is.na(ell_opt$alpha) ||
        ell_opt$alpha < 0 || ell_opt$alpha > 1) {
        stop("Error! ellipse.args$alpha must be a numeric value in [0, 1]")
    }
    if (!is.numeric(ell_opt$border.alpha) || length(ell_opt$border.alpha) != 1 ||
        is.na(ell_opt$border.alpha) || ell_opt$border.alpha < 0 || ell_opt$border.alpha > 1) {
        stop("Error! ellipse.args$border.alpha must be a numeric value in [0, 1]")
    }
    if (!(isFALSE(centroids) || isTRUE(centroids) ||
          identical(centroids, "points") || identical(centroids, "labels"))) {
        stop("Error! centroids must be FALSE, TRUE, \"points\", or \"labels\"")
    }
    if (!identical(lab_opt$pos, "auto") && !identical(lab_opt$pos, "random")) {
        if (!is.numeric(lab_opt$pos) || anyNA(lab_opt$pos) ||
            any(lab_opt$pos != as.integer(lab_opt$pos)) || !all(lab_opt$pos %in% 1:4)) {
            stop("Error! lab.args$pos must be \"auto\", \"random\", or integer values in 1:4")
        }
    }

    has_xsd <- !missing(xsd) && !is.null(xsd)
    has_limx <- !missing(limx) && !is.null(limx)
    has_group <- !missing(group) && !is.null(group)
    has_group2 <- !missing(group2) && !is.null(group2)
    has_group3 <- !missing(group3) && !is.null(group3)
    has_col <- !missing(col) && !is.null(col)
    has_pch <- !missing(pch) && !is.null(pch)
    if (is.null(pch)) pch <- 21
    has_xlab <- !missing(xlab) && !is.null(xlab)
    has_ylab <- !missing(ylab) && !is.null(ylab)
    has_labd <- !missing(labd) && !is.null(labd)
    has_f1 <- !missing(f1) && !is.null(f1)
    has_f2 <- !missing(f2) && !is.null(f2)

    if (ncol(x) > 2) x <- x[, pc, drop = FALSE]
    if (ncol(y) > 2) y <- y[, pc, drop = FALSE]

    n <- nrow(x)
    nv <- nrow(y)

    if (has_xsd && sum(xsd^2) > 0) {
        pvar <- xsd^2 / sum(xsd^2)
    } else {
        pvar <- NULL
    }

    # Plot limits
    al <- ceiling(max(abs(range(x, na.rm = TRUE))) * 10) / 10
    lim <- c(-al, al)
    alr <- ceiling(max(abs(range(y, na.rm = TRUE))) * 10) / 10
    limr <- c(-alr, alr)

    if (has_limx) {
        lim <- lim * limx
        limr <- limr * limx
    }

    # Check groups are factors
    is_svy_g3 <- has_group3 && inherits(group3, "bb_svy_weights")
    if (has_group && !is.factor(group)) group <- as.factor(group)
    if (has_group && length(group) != n) stop("Error! group must be same length as original data used in PCA")
    if (has_group2 && !is.factor(group2)) group2 <- as.factor(group2)
    if (has_group3 && !is_svy_g3 && !is.factor(group3)) group3 <- as.factor(group3)

    # Convex hulls
    if (isTRUE(chull)) {
        if (!has_group) stop("Error! group must be supplied to draw convex hulls")
        hulls <- lapply(seq_len(nlevels(group)), function(x1) {
            idx <- which(as.integer(group) == x1)
            if (length(idx) == 0) return(NULL)
            sub_x <- as.matrix(x[idx, 1:2, drop = FALSE])
            sub_x <- sub_x[is.finite(sub_x[, 1]) & is.finite(sub_x[, 2]), , drop = FALSE]
            if (nrow(sub_x) == 0) return(NULL)
            sub_x[chull(sub_x[, 1], sub_x[, 2]), , drop = FALSE]
        })
    }

    # Ellipses
    if (isTRUE(ellipse)) {
        if (!has_group) stop("Error! group must be supplied to draw ellipses")
        e <- lapply(seq_len(nlevels(group)), function(x1) {
            idx <- which(as.integer(group) == x1)
            if (length(idx) == 0) return(NULL)
            .bb_ellipse(x = x[idx, 1:2, drop = FALSE], method = ell_opt$method, level = ell_opt$level)
        })
        if (isTRUE(ell_opt$autolim)) {
            e_valid <- Filter(Negate(is.null), e)
            if (length(e_valid) > 0) {
                ec <- do.call(rbind, e_valid)
                ecm <- max(abs(ec), na.rm = TRUE)
                if (is.finite(ecm) && max(lim) > 0 && ecm > max(lim)) {
                    pc_mult <- ceiling((ecm / max(lim)) * 10) / 10
                    lim <- lim * pc_mult
                    limr <- limr * pc_mult
                }
            }
        }
    }

    # Circle of equilibrium contribution
    if (isTRUE(circle.eq)) {
        cr <- sqrt(2 / nv) * var_opt$expand
        if ((cr / limr[2]) < 0.1) {
            message("Info: Small circle of equilibrium. Did you set scale=0?")
            circle.eq <- FALSE
        } else {
            t_seq <- seq(0, 2 * pi, 0.1)
            ceq <- data.frame(x = cr * cos(t_seq), y = cr * sin(t_seq))
        }
    }

    # Label rotation
    if (isTRUE(lab_opt$rotate)) {
        if (has_labd) {
            lsrt <- .bb_lab_rotation(x = x, valign = lab_opt$valign, ran_adj = lab_opt$ran.adj)
        }
        if (isTRUE(variables)) {
            lvsrt <- .bb_lab_rotation(x = y, valign = lab_opt$valign, ran_adj = FALSE)
        }
    }

    # Point styling by group
    p3 <- 1
    if (has_group) {
        grc <- 1

        if (has_col) {
            if (is.character(col) && length(col) == 1 && col %in% palette.pals()) {
                col <- palette.colors(nlevels(group), palette = col, recycle = TRUE)
            } else if (is.character(col) && length(col) == 1 && col %in% hcl.pals()) {
                col <- hcl.colors(nlevels(group), palette = col)
            } else if (length(col) < nlevels(group)) {
                message("Warning: Number of colours does not match number of groups")
                col <- rep(col, length.out = nlevels(group))
            }
        } else {
            col <- hcl.colors(nlevels(group), "Dark 3")
        }

        group_cols <- col
        col <- col[group]

        if (has_group2) {
            grc <- grc + 1
            if (length(group2) != n) stop("Error! group2 must be same length as original data used in PCA")
            if (!has_pch || length(pch) < nlevels(group2)) {
                pch <- seq(21, length.out = nlevels(group2))
                pch <- ifelse(pch > 25, pch - 26, pch)
                if (length(pch) > 25) pch <- ifelse(pch > 25, pch + 6, pch)
            }
            pch1 <- pch[group2]
        } else {
            if (length(pch) < nlevels(group)) pch <- rep(pch, length.out = nlevels(group))
            pch1 <- pch[group]
        }

        if (has_group3) {
            grc <- grc + 1
            if (!is_svy_g3) {
                if (length(group3) != n) stop("Error! group3 must be same length as original data used in PCA")
                k3 <- nlevels(group3)
                g3i <- 2 / k3
                pchm <- seq(from = 1, by = g3i, length.out = k3)
                p3 <- pchm[group3]
                g3_idx <- unique(c(1L, as.integer(ceiling(k3 / 2)), k3))
                g3_leg_lab <- paste(levels(group3))[g3_idx]
                g3_leg_mult <- pchm[g3_idx]
            }
        }
    } else {
        if (has_group2) message("Warning: group argument required when using group2.")
        if (has_group3 && !is_svy_g3) {
            message("Warning: group argument required when using group3.")
        }
        if (length(pch) > 1) message("Warning: Multiple pch symbols supplied, but no groups.")
        pch1 <- pch
        if (has_col) {
            if (length(col) > 1) message("Warning: Multiple colours supplied, but no groups.")
        } else {
            col <- "#4050ad"
        }
    }

    if (is_svy_g3) {
        w <- as.numeric(group3)
        if (length(w) != n) stop("Error! group3 must be same length as original data used in PCA")
        w_max <- max(w, na.rm = TRUE)
        if (!is.finite(w_max) || w_max <= 0) {
            stop("Error! survey weights must contain positive finite values")
        }
        p3 <- pmax(0.2, sqrt(w / w_max))
        w_breaks <- unique(c(min(w, na.rm = TRUE), median(w, na.rm = TRUE), w_max))
        g3_leg_lab <- format(signif(w_breaks, 3), trim = TRUE)
        g3_leg_mult <- pmax(0.2, sqrt(w_breaks / w_max))
    }

    avail_groups <- if (has_group) {
        c(1L, if (has_group2) 2L, if (has_group3) 3L)
    } else if (is_svy_g3) {
        3L
    } else {
        integer(0)
    }
    if (isFALSE(legend) || is.null(legend) || length(legend) == 0) {
        leg_groups <- integer(0)
    } else if (isTRUE(legend) || identical(legend, "all")) {
        if (length(avail_groups) == 0) {
            message("Info: legend only works when group specified.")
            leg_groups <- integer(0)
        } else {
            leg_groups <- avail_groups
        }
    } else if (is.numeric(legend)) {
        if (anyNA(legend) || any(legend != as.integer(legend)) || !all(legend %in% 1:3)) {
            stop("Error! legend must be FALSE, TRUE, \"all\", or integer values in 1:3")
        }
        if (length(avail_groups) == 0) {
            message("Info: legend only works when group specified.")
            leg_groups <- integer(0)
        } else {
            req_groups <- unique(as.integer(legend))
            missing_groups <- setdiff(req_groups, avail_groups)
            if (length(missing_groups) > 0) {
                message("Info: Requested legend group(s) not supplied: ",
                        paste(missing_groups, collapse = ", "))
            }
            leg_groups <- req_groups[req_groups %in% avail_groups]
        }
    } else {
        stop("Error! legend must be FALSE, TRUE, \"all\", or integer values in 1:3")
    }

    colf <- col
    if (length(pch1) == 1 && length(col) > 1) {
        if (pch1 %in% 21:25) col <- "black"
    } else {
        if (length(pch1) > 1 && length(col) > 1 && length(pch1) < length(col)) {
            pch1 <- rep(pch1, length.out = length(col))
        }
        pchw <- which(pch1 %in% 21:25)
        if (length(pchw) > 0) {
            if (length(col) < length(pch1)) col <- rep(col, length.out = length(pch1))
            col[pchw] <- "black"
        }
    }

    alp <- if (isFALSE(ell_opt$nofill)) c(ell_opt$alpha, ell_opt$border.alpha) else c(0, ell_opt$border.alpha)

    if (!isFALSE(centroids) && !has_group) {
        stop("Error! group must be supplied to draw centroids")
    }
    if (isTRUE(spider) && !has_group) {
        stop("Error! group must be supplied to draw spider lines")
    }

    cen_mat <- NULL
    if (has_group) {
        cen_mat <- matrix(NA_real_, nrow = nlevels(group), ncol = 2,
                          dimnames = list(levels(group), colnames(x)))
        for (i in seq_len(nlevels(group))) {
            idx <- which(as.integer(group) == i)
            if (length(idx) > 0) {
                cen_mat[i, ] <- colMeans(x[idx, , drop = FALSE], na.rm = TRUE)
            }
        }
    }

    if (!has_xlab) {
        xlab <- if (!is.null(pvar)) {
            paste0("PC", pc[1], " (", round(pvar[pc[1]] * 100, 1), "%)")
        } else {
            paste0("PC", pc[1])
        }
    }
    if (!has_ylab) {
        ylab <- if (!is.null(pvar)) {
            paste0("PC", pc[2], " (", round(pvar[pc[2]] * 100, 1), "%)")
        } else {
            paste0("PC", pc[2])
        }
    }

    # Main plot
    if (length(leg_groups) > 0 && (is.null(leg_opt$x) || leg_opt$x > max(lim))) {
        op <- par(pty = "s", mar = c(5, 4, 4, 7))
    } else {
        op <- par(pty = "s")
    }
    on.exit(par(op), add = TRUE)

    do.call(plot, c(
        list(x = x[, 1], y = x[, 2], type = "n", asp = 1, xlim = lim, ylim = lim,
             xlab = xlab, ylab = ylab, axes = FALSE, cex.lab = cex.axis * cex),
        plot_dots
    ))
    if (isTRUE(axes)) {
        axis(1, lwd = 0, lwd.ticks = 1 * (cex * 0.5), cex.axis = cex.axis * cex)
        axis(2, lwd = 0, lwd.ticks = 1 * (cex * 0.5), cex.axis = cex.axis * cex, las = 1)
        box(lwd = 1 * (cex * 0.5))
    }

    if (has_f1) {
        f1 <- match.fun(f1)
        do.call(f1, if (!is.null(f1_a)) f1_a else list())
    }

    if (isTRUE(chull)) {
        for (i in seq_len(nlevels(group))) {
            if (!is.null(hulls[[i]])) {
                polygon(x = hulls[[i]][, 1], y = hulls[[i]][, 2],
                        col = adjustcolor(group_cols[i], alp[1]),
                        border = adjustcolor(group_cols[i], alp[2]),
                        lwd = ell_opt$lwd * cex)
            }
        }
    }

    if (isTRUE(ellipse)) {
        for (i in seq_len(nlevels(group))) {
            if (!is.null(e[[i]])) {
                polygon(x = e[[i]]$x, y = e[[i]]$y,
                        col = adjustcolor(group_cols[i], alp[1]),
                        border = adjustcolor(group_cols[i], alp[2]),
                        lwd = ell_opt$lwd * cex)
            }
        }
    }

    if (isTRUE(grid)) {
        abline(v = 0, h = 0, lty = 2, lwd = 1 * cex, col = adjustcolor("black", 0.75))
    }

    if (isTRUE(spider)) {
        for (i in seq_len(nlevels(group))) {
            if (anyNA(cen_mat[i, ])) next
            idx <- which(as.integer(group) == i)
            if (length(idx) > 0) {
                segments(x0 = cen_mat[i, 1], y0 = cen_mat[i, 2],
                         x1 = x[idx, 1], y1 = x[idx, 2],
                         col = adjustcolor(group_cols[i], ell_opt$border.alpha * 0.6),
                         lwd = max(0.75, ell_opt$lwd * 0.5) * cex)
            }
        }
    }

    points(x[, 1], x[, 2], pch = pch1, col = adjustcolor(col, alpha),
           bg = adjustcolor(colf, alpha), cex = cex.pt * p3 * cex)

    if (!isFALSE(centroids)) {
        cex_pt_cen <- if (length(cex.pt) > 1) 1 else cex.pt
        for (i in seq_len(nlevels(group))) {
            if (anyNA(cen_mat[i, ])) next
            if (isTRUE(centroids) || identical(centroids, "points")) {
                points(cen_mat[i, 1], cen_mat[i, 2], pch = 23,
                       bg = group_cols[i], col = "black",
                       lwd = 1.5 * cex, cex = cex_pt_cen * 1.5 * cex)
            }
            if (isTRUE(centroids) || identical(centroids, "labels")) {
                text(cen_mat[i, 1], cen_mat[i, 2], labels = levels(group)[i],
                     font = 2, col = "black",
                     pos = if (isTRUE(centroids)) 3 else NULL,
                     cex = cex.axis * cex, xpd = TRUE)
            }
        }
    }

    if (has_labd) {
        if (length(labd) != nrow(x)) {
            message("Info: Number of labels does not match number of data observations")
        }
        if (isTRUE(lab_opt$rotate)) {
            for (i in seq_len(nrow(x))) {
                text(x[i, 1], x[i, 2], labels = labd[i], font = lab_opt$font,
                     col = lab_opt$col, srt = lsrt[i, 1], adj = c(lsrt[i, 2], lsrt[i, 3]),
                     cex = lab_opt$cex * cex, xpd = TRUE)
            }
        } else {
            adjld <- if (identical(lab_opt$pos, "auto")) {
                pos_auto <- ifelse(abs(x[, 1]) >= abs(x[, 2]),
                                   ifelse(x[, 1] >= 0, 4L, 2L),
                                   ifelse(x[, 2] >= 0, 3L, 1L))
                pos_auto[is.na(pos_auto)] <- 1L
                pos_auto
            } else if (identical(lab_opt$pos, "random")) {
                sample.int(4L, nrow(x), replace = TRUE)
            } else {
                as.integer(lab_opt$pos)
            }
            text(x[, 1], x[, 2], labels = labd, font = lab_opt$font,
                 col = lab_opt$col, pos = adjld, cex = lab_opt$cex * cex, xpd = TRUE)
        }
    }

    # Legend
    if (length(leg_groups) > 0) {
        default_titles <- c("Group", "Group 2", if (is_svy_g3) "Weight" else "Group 3")
        user_titles <- leg_opt$title
        leg_titles <- if (is.null(user_titles)) {
            default_titles[leg_groups]
        } else if (length(user_titles) == length(leg_groups)) {
            user_titles
        } else if (length(user_titles) >= max(leg_groups)) {
            user_titles[leg_groups]
        } else {
            rep(user_titles, length.out = length(leg_groups))
        }

        lpx <- if (!is.null(leg_opt$x)) leg_opt$x else max(lim) * 1.38
        lpy <- if (!is.null(leg_opt$y)) leg_opt$y else ifelse(length(leg_groups) > 1, max(lim), 0)

        cex_pt_leg <- if (length(cex.pt) > 1) 1 else cex.pt

        for (k in seq_along(leg_groups)) {
            g <- leg_groups[k]
            if (k == 1) {
                y_pos <- lpy
            } else if (isTRUE(leg_opt$horiz)) {
                lpx <- lp$rect$left + (lp$rect$w * 1.5)
                y_pos <- lpy
            } else {
                lpy <- lp$text$y[length(lp$text$y)]
                y_pos <- lpy - (lp1 * 1.25)
            }

            if (g == 1L) {
                if (grc == 1) {
                    g1p <- pch
                    g1c <- group_cols
                    g1c[which(g1p %in% 21:25)] <- "black"
                } else {
                    g1p <- 22
                    g1c <- "black"
                }
                lp <- legend(x = lpx, y = y_pos, yjust = 1, xjust = 0,
                             legend = paste(levels(group)), bty = "n", col = g1c,
                             pt.bg = group_cols, pch = g1p, pt.cex = leg_opt$cex * 2 * cex,
                             xpd = NA, cex = leg_opt$cex * cex, text.width = 1,
                             title = leg_titles[k], title.font = 2, title.adj = 0)
            } else if (g == 2L) {
                lp <- legend(x = lpx, y = y_pos, yjust = 1, xjust = 0,
                             legend = paste(levels(group2)), bty = "n", col = "black",
                             pt.bg = adjustcolor("black", 0.25), pch = pch,
                             pt.cex = cex_pt_leg * cex, xpd = NA, cex = leg_opt$cex * cex,
                             text.width = 1, title = leg_titles[k], title.font = 2,
                             title.adj = 0)
            } else if (g == 3L) {
                lp <- legend(x = lpx, y = y_pos, yjust = 1, xjust = 0,
                             legend = g3_leg_lab, bty = "n", col = "black",
                             pt.bg = adjustcolor("black", 0.25), pch = 21,
                             pt.cex = cex_pt_leg * g3_leg_mult * cex, xpd = NA,
                             cex = leg_opt$cex * cex, text.width = 1,
                             title = leg_titles[k], title.font = 2,
                             title.adj = 0, y.intersp = 1.25)
            }

            lp1 <- if (length(lp$text$y) > 1) {
                abs(lp$text$y[1] - lp$text$y[2])
            } else {
                lp$rect$h / 2
            }
        }
    }

    # Variables overlay
    if (isTRUE(variables)) {
        vname <- if (!is.null(var_opt$vname)) var_opt$vname else rownames(y)
        if (is.null(vname)) vname <- seq_len(nrow(y))

        if (!is.null(var_opt$which)) {
            whichv <- var_opt$which
            if (is.character(whichv)) {
                if ("circle.eq" %in% whichv) {
                    if (!isTRUE(circle.eq)) {
                        stop("Error! circle.eq must be TRUE (with scale = 0) when which = \"circle.eq\"")
                    }
                    vmatch <- which(sqrt((y[, 1] * var_opt$expand)^2 + (y[, 2] * var_opt$expand)^2) > cr)
                } else {
                    vmatch <- which(vname %in% whichv)
                    if (length(vmatch) > 0 && length(vmatch) != length(whichv)) {
                        message("Info: Some variables do not match variable names")
                    }
                }
            } else {
                vmatch <- whichv
            }
            if (length(vmatch) == 0) {
                stop("Error! Supplied variables do not match any variable names")
            }
            y <- y[vmatch, , drop = FALSE]
            vname <- vname[vmatch]
            if (isTRUE(lab_opt$rotate)) lvsrt <- lvsrt[vmatch, , drop = FALSE]
        }

        par(new = TRUE)
        plot(y[, 1], y[, 2], type = "n", asp = 1, xlim = limr, ylim = limr,
             ann = FALSE, axes = FALSE)
        if (isTRUE(var_opt$axes) && !isFALSE(axes)) {
            axis(3, lwd = 0, lwd.ticks = 1 * (cex * 0.5), cex.axis = cex.axis * cex)
            axis(4, lwd = 0, lwd.ticks = 1 * (cex * 0.5), cex.axis = cex.axis * cex, las = 1)
        }

        if (has_f2) {
            f2 <- match.fun(f2)
            do.call(f2, if (!is.null(f2_a)) f2_a else list())
        }

        if (isTRUE(circle.eq)) {
            polygon(x = ceq[, 1], y = ceq[, 2], border = var_opt$col,
                    lwd = ifelse(var_opt$lwd <= 2, 1.25, (var_opt$lwd / 2) * 1.25) * cex)
        }

        arrows(x0 = 0, x1 = y[, 1] * var_opt$expand, y0 = 0, y1 = y[, 2] * var_opt$expand,
               col = var_opt$col, length = var_opt$arrow.len * cex, lwd = var_opt$lwd * cex, xpd = TRUE)

        if (isTRUE(lab_opt$rotate)) {
            for (i in seq_len(nrow(y))) {
                text(y[i, 1] * var_opt$expand * var_opt$lab.offset,
                     y[i, 2] * var_opt$expand * var_opt$lab.offset,
                     labels = vname[i],
                     font = var_opt$font.lab, col = var_opt$col.lab, srt = lvsrt[i, 1],
                     adj = c(lvsrt[i, 2], lvsrt[i, 3]), cex = var_opt$cex.lab * cex, xpd = NA)
            }
        } else {
            text(y[, 1] * var_opt$expand * var_opt$lab.offset,
                 y[, 2] * var_opt$expand * var_opt$lab.offset,
                 labels = vname,
                 font = var_opt$font.lab, col = var_opt$col.lab, cex = var_opt$cex.lab * cex, xpd = TRUE)
        }
    }

    invisible(structure(
        list(
            scores = x,
            loadings = if (isTRUE(variables)) y else NULL,
            lim = lim,
            limr = limr,
            pc = pc,
            pvar = pvar,
            group_cols = if (has_group) group_cols else colf,
            centroids = cen_mat
        ),
        class = "bb_biplot"
    ))
}
