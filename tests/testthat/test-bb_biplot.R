test_that("bb_biplot works with prcomp and princomp objects", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p1 <- stats::prcomp(iris[-5], scale. = TRUE)
    expect_invisible(bb_biplot(p1))
    expect_invisible(bb_biplot(p1, scale = 0, pc.biplot = TRUE, varimax.rotate = TRUE))
    expect_invisible(bb_biplot(p1, scale = 0, varimax.rotate = TRUE,
                               varimax.args = list(normalize = FALSE, eps = 1e-6)))
    expect_error(bb_biplot(p1, varimax.rotate = TRUE,
                           varimax.args = list(bad_key = 1)), "Unknown option")

    p2 <- stats::princomp(iris[-5], cor = TRUE)
    expect_invisible(bb_biplot(p2))
    expect_invisible(bb_biplot(p2, pc = c(2, 3), scale = 0))
})

test_that("bb_biplot.default works with matrices, without xsd, and with ... args", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    expect_invisible(bb_biplot(p$x, p$rotation, main = "Custom Title"))
    expect_invisible(bb_biplot(p$x, p$rotation, xsd = p$sdev, variables = FALSE))
})

test_that("grouping, legends, hulls, and ellipses work", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    gr <- factor(iris$Species)
    set.seed(1)
    gr2 <- factor(sample(1:3, nrow(iris), replace = TRUE))
    gr3 <- factor(sample(1:3, nrow(iris), replace = TRUE))

    expect_invisible(bb_biplot(p, scale = 0, group = gr, chull = TRUE, legend = TRUE))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, chull = TRUE,
                               chull.args = list(nofill = TRUE, lwd = 1.5)))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, ellipse = TRUE,
                               ellipse.args = list(method = "cov")))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, ellipse = TRUE,
                               ellipse.args = list(method = "t", level = 0.90, nofill = TRUE)))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3,
                               legend = TRUE, title.leg = c("Sp", "G2", "G3"),
                               leg.args = list(horiz = TRUE)))
})

test_that("circle.eq, var.args$which, and lab.args$rotate work", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    expect_invisible(bb_biplot(p, scale = 0, circle.eq = TRUE,
                               var.args = list(which = "circle.eq")))
    expect_invisible(bb_biplot(p, scale = 0, circle.eq = TRUE,
                               var.args = list(which = "circle.eq", expand = 2)))
    expect_invisible(bb_biplot(p, scale = 0, var.args = list(which = "Petal.Length"),
                               lab.args = list(rotate = TRUE)))
    expect_invisible(bb_biplot(p, scale = 0, var.args = list(which = c(1, 3)),
                               lab.args = list(rotate = TRUE)))
    expect_invisible(bb_biplot(p, scale = 0, lab.args = list(rotate = TRUE),
                               labd = rownames(iris)))
})

test_that("f1 and f2 callbacks work with default arguments including '...'", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    expect_invisible(bb_biplot(p, scale = 0, f1 = graphics::grid, f2 = graphics::grid))
})

test_that("edge cases for groups, legends, ellipses, and label rotation work", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)

    # group2 with fewer levels than group using default pch
    gr5 <- factor(rep(1:5, length.out = nrow(iris)))
    gr2 <- factor(rep(1:2, length.out = nrow(iris)))
    expect_invisible(bb_biplot(p, scale = 0, group = gr5, group2 = gr2, legend = TRUE))

    # Single-level group with group2 and legend
    gr1 <- factor(rep("all", nrow(iris)))
    expect_invisible(bb_biplot(p, scale = 0, group = gr1, group2 = gr2, legend = TRUE))

    # Unused factor level in group with chull and ellipse
    gr_unused <- factor(iris$Species, levels = c(levels(iris$Species), "empty"))
    expect_invisible(bb_biplot(p, scale = 0, group = gr_unused, chull = TRUE, ellipse = TRUE))

    # Constant x coordinates in .bb_ellipse
    df_const <- data.frame(PC1 = c(1, 1, 1), PC2 = c(0, 1, 2))
    el_cov <- betterbiplot:::.bb_ellipse(df_const, method = "cov", level = 0.95)
    el_t <- betterbiplot:::.bb_ellipse(df_const, method = "t", level = 0.95)
    expect_false(anyNA(el_cov$x))
    expect_false(anyNA(el_cov$y))
    expect_false(anyNA(el_t$x))
    expect_false(anyNA(el_t$y))
    expect_gt(diff(range(el_t$y)), diff(range(el_cov$y)))

    # Origin and negative x-axis in .bb_lab_rotation
    coords <- matrix(c(0, 0, -1, 0), ncol = 2, byrow = TRUE)
    rot <- betterbiplot:::.bb_lab_rotation(coords, valign = 0.5, ran_adj = FALSE)
    expect_false(anyNA(rot))
    expect_equal(unname(rot[1, 1]), 0)
    expect_equal(unname(rot[2, 2]), 1)
})

test_that("selective group legends work", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    gr <- factor(iris$Species)
    set.seed(1)
    gr2 <- factor(sample(1:3, nrow(iris), replace = TRUE))
    gr3 <- factor(sample(1:3, nrow(iris), replace = TRUE))

    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3, legend = "all"))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3, legend = 1))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3, legend = 2))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3, legend = 3))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3, legend = c(1, 3)))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3, legend = c(2, 3)))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3,
                               legend = c(3, 1), leg.args = list(horiz = TRUE)))

    # title.leg positional and group-index matching
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3,
                               legend = 3, title.leg = "Size"))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3,
                               legend = c(1, 3), title.leg = c("Sp", "Size")))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group2 = gr2, group3 = gr3,
                               legend = 2, title.leg = c("Sp", "Shape", "Size")))

    # Messages for missing groups
    expect_message(bb_biplot(p, scale = 0, group = gr, legend = 3), "not supplied")
    expect_message(bb_biplot(p, scale = 0, legend = 1), "when group specified")
    expect_message(bb_biplot(p, scale = 0, legend = "all"), "when group specified")

    # Errors for invalid legend values
    expect_error(bb_biplot(p, scale = 0, group = gr, legend = 4), "integer values in 1:3")
    expect_error(bb_biplot(p, scale = 0, group = gr, legend = 1.5), "integer values in 1:3")
    expect_error(bb_biplot(p, scale = 0, group = gr, legend = "invalid"), "integer values in 1:3")
})

test_that("group3 legends condense to min, mid, and max", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    gr <- factor(iris$Species)
    gr3_many <- factor(rep(c(0.5, 1, 2, 2.5, 4, 6, 7), length.out = nrow(iris)))
    gr3_two <- factor(rep(c("low", "high"), length.out = nrow(iris)))
    gr3_one <- factor(rep("only", nrow(iris)))

    expect_invisible(bb_biplot(p, scale = 0, group = gr, group3 = gr3_many, legend = TRUE))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group3 = gr3_many, legend = 3))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group3 = gr3_two, legend = 3))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, group3 = gr3_one, legend = 3))
})

test_that("invalid arguments and option lists raise clear errors", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    gr <- factor(iris$Species)

    expect_error(bb_biplot(p, chull = TRUE), "group must be supplied")
    expect_error(bb_biplot(p, ellipse = TRUE), "group must be supplied")
    expect_error(bb_biplot(p, group = gr, ellipse = TRUE,
                           ellipse.args = list(method = "invalid")), "cov.*t")
    expect_error(bb_biplot(p, group = gr[1:10]), "same length")
    expect_error(bb_biplot(p, var.args = list(which = "nonexistent")), "do not match")
    expect_error(bb_biplot(p, var.args = list(which = "circle.eq")), "circle.eq must be TRUE")

    # Option list validation
    expect_error(bb_biplot(p, lab.args = list(bad_key = 1)), "Unknown option")
    expect_error(bb_biplot(p, var.args = list(bad_key = 1)), "Unknown option")
    expect_error(bb_biplot(p, ellipse.args = list(bad_key = 1)), "Unknown option")
    expect_error(bb_biplot(p, chull.args = list(bad_key = 1)), "Unknown option")
    expect_error(bb_biplot(p, leg.args = list(bad_key = 1)), "Unknown option")
    expect_error(bb_biplot(p, lab.args = "not_a_list"), "named list")
    expect_error(bb_biplot(p, pc = 1), "length 2")
})

test_that("covariance ellipses, centroids, spider lines, palettes, alpha, lab.offset, and return object work", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    gr <- factor(iris$Species)

    # Covariance ellipses & return object
    res <- bb_biplot(p, scale = 0, group = gr, col = "Okabe-Ito", alpha = 0.6,
                     ellipse = TRUE, ellipse.args = list(method = "cov", level = 0.90,
                                                         alpha = 0.15, border.alpha = 0.8),
                     centroids = TRUE, spider = TRUE,
                     var.args = list(lab.offset = 1.15),
                     labd = rownames(iris), lab.args = list(pos = "auto"))
    expect_s3_class(res, "bb_biplot")
    expect_equal(nrow(res$scores), nrow(iris))
    expect_equal(ncol(res$loadings), 2)
    expect_length(res$lim, 2)
    expect_length(res$limr, 2)
    expect_equal(res$pc, c(1L, 2L))
    expect_equal(nrow(res$centroids), nlevels(gr))
    expect_length(res$group_cols, nlevels(gr))

    # hcl.pals palette, centroids variants, and lab.args$pos variants
    expect_invisible(bb_biplot(p, scale = 0, group = gr, col = "Viridis",
                               centroids = "points", labd = rownames(iris),
                               lab.args = list(pos = 3)))
    expect_invisible(bb_biplot(p, scale = 0, group = gr, centroids = "labels",
                               labd = rownames(iris), lab.args = list(pos = "random")))

    # Unused level (< 3 points) in covariance ellipse & centroids
    gr_unused <- factor(iris$Species, levels = c(levels(iris$Species), "empty"))
    expect_invisible(bb_biplot(p, scale = 0, group = gr_unused, ellipse = TRUE,
                               ellipse.args = list(method = "cov"),
                               centroids = TRUE, spider = TRUE))

    # Validation errors
    expect_error(bb_biplot(p, alpha = 1.5), "alpha must be")
    expect_error(bb_biplot(p, group = gr, ellipse = TRUE,
                           ellipse.args = list(method = "bad")), "cov.*t")
    expect_error(bb_biplot(p, group = gr, ellipse = TRUE,
                           ellipse.args = list(method = "cov", level = 1.2)), "between 0 and 1")
    expect_error(bb_biplot(p, group = gr, centroids = "bad"), "centroids must be")
    expect_error(bb_biplot(p, centroids = TRUE), "group must be supplied")
    expect_error(bb_biplot(p, spider = TRUE), "group must be supplied")
    expect_error(bb_biplot(p, labd = rownames(iris), lab.args = list(pos = 5)), "pos must be")

    # Unlabelled matrix y with rotated labels, NA score in x, and missing PCA scores
    y_nonames <- unname(p$rotation)
    expect_invisible(bb_biplot(p$x, y_nonames, lab.args = list(rotate = TRUE)))

    x_na <- p$x
    x_na[1, 1] <- NA_real_
    expect_invisible(bb_biplot(x_na, p$rotation))

    p_noret <- stats::prcomp(iris[-5], scale. = TRUE, retx = FALSE)
    expect_error(bb_biplot(p_noret), "retx = TRUE")

    p_nosco <- stats::princomp(iris[-5], cor = TRUE, scores = FALSE)
    expect_error(bb_biplot(p_nosco), "scores = TRUE")
})
