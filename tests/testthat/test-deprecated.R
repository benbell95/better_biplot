test_that("deprecated v1.x arguments work via soft deprecation shim", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    gr <- factor(iris$Species)

    expect_warning(bb_biplot(p, pc1 = 2, pc2 = 3, scale = 0), "deprecated")
    expect_warning(bb_biplot(p, scale = 0, lab_rotation = TRUE, labd = rownames(iris),
                             col.labd = "grey30", cex.labd = 0.6, font.labd = 2,
                             ran_adj = FALSE, valign = 0.5), "deprecated")
    expect_warning(bb_biplot(p, scale = 0, circle.eq = TRUE, whichv = "circle.eq",
                             expand = 1.1, arrow.len = 0.1, lwd.v = 1.5, col.v = "navy",
                             col.labv = "darkblue", cex.labv = 0.9, font.labv = 3,
                             axes.v = FALSE), "deprecated")
    expect_warning(bb_biplot(p, scale = 0, group = gr, ellipse = TRUE,
                             angle = "atan", autolim = TRUE, nofill = TRUE,
                             lwd.e = 1.5, se = 1), "deprecated")
    expect_warning(bb_biplot(p, scale = 0, group = gr, legend = TRUE,
                             title.leg = "Species", cex.leg = 0.8,
                             horiz = FALSE, lpx = 2, lpy = 2), "deprecated")
    expect_warning(bb_biplot(p, scale = 0, f1 = graphics::grid, f1.a = list(col = "grey"),
                             f2 = graphics::grid, f2.a = list(col = "grey")), "deprecated")
})
