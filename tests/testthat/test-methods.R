test_that("bb_biplot.rda works with vegan::rda PCA objects", {
    skip_if_not_installed("vegan")
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    r <- vegan::rda(iris[-5], scale = TRUE)
    gr <- factor(iris$Species)

    expect_invisible(bb_biplot(r))
    expect_invisible(bb_biplot(r, scale = 0, group = gr, circle.eq = TRUE, legend = TRUE))
    expect_invisible(bb_biplot(r, scale = 1, pc.biplot = TRUE, varimax.rotate = TRUE))

    # Constrained RDA should error cleanly
    r_cca <- vegan::rda(iris[-5] ~ Species, data = iris, scale = TRUE)
    expect_error(bb_biplot(r_cca), "unconstrained PCA")
})

test_that("bb_biplot.svyprcomp works with survey::svyprcomp objects", {
    skip_if_not_installed("survey")
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    env <- new.env(parent = emptyenv())
    utils::data("api", package = "survey", envir = env)
    dclus2 <- survey::svydesign(
        id = ~dnum + snum,
        fpc = ~fpc1 + fpc2,
        data = env$apiclus2
    )
    sp <- survey::svyprcomp(~api99 + api00 + ell, design = dclus2, scale = TRUE, scores = TRUE)
    gr <- factor(env$apiclus2$stype)
    gr2 <- factor(rep(1:2, length.out = nrow(env$apiclus2)))

    expect_invisible(bb_biplot(sp))
    expect_invisible(bb_biplot(sp, scale = 0, group3 = "scaled"))
    expect_invisible(bb_biplot(sp, scale = 0, group3 = "scaled", legend = TRUE))
    expect_invisible(bb_biplot(sp, scale = 0, group = gr, group3 = "scaled", legend = TRUE))
    expect_invisible(bb_biplot(sp, scale = 0, group = gr, group2 = gr2, group3 = "scaled", legend = TRUE))
    expect_invisible(bb_biplot(sp, scale = 0, group = gr, group3 = "scaled", legend = 3))
    expect_invisible(bb_biplot(sp, scale = 0, group = gr, group3 = gr2, legend = TRUE))

    # Missing scores = TRUE should error cleanly
    sp_noscores <- survey::svyprcomp(~api99 + api00 + ell, design = dclus2, scale = TRUE, scores = FALSE)
    expect_error(bb_biplot(sp_noscores), "scores = TRUE")

    # Missing weights with group3 = "scaled" should error cleanly
    sp_noweights <- sp
    sp_noweights$weights <- NULL
    expect_error(bb_biplot(sp_noweights, group3 = "scaled"), "weights")
})

test_that("bb_biplot.PCA works with FactoMineR::PCA-structured objects", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    mock_pca <- structure(
        list(
            eig = cbind(eigenvalue = p$sdev^2),
            ind = list(coord = p$x),
            var = list(coord = t(t(p$rotation) * p$sdev))
        ),
        class = c("PCA", "list")
    )

    expect_invisible(bb_biplot(mock_pca))
    expect_invisible(bb_biplot(mock_pca, scale = 0, group = factor(iris$Species), circle.eq = TRUE))
    expect_error(bb_biplot(structure(list(), class = "PCA")), "FactoMineR::PCA")
})

test_that("bb_biplot.pca works with ade4::dudi.pca-structured objects", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    mock_dudi <- structure(
        list(
            eig = p$sdev^2,
            li = as.data.frame(p$x),
            c1 = as.data.frame(p$rotation),
            co = as.data.frame(t(t(p$rotation) * p$sdev))
        ),
        class = c("pca", "dudi")
    )

    expect_invisible(bb_biplot(mock_dudi))
    expect_invisible(bb_biplot(mock_dudi, scale = 0, group = factor(iris$Species), ellipse = TRUE))
    expect_error(bb_biplot(structure(list(), class = c("pca", "dudi"))), "ade4::dudi.pca")
})

test_that("bb_biplot.principal works with psych::principal-structured objects", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    loadings_mat <- structure(t(t(p$rotation) * p$sdev), class = "loadings")
    mock_principal <- structure(
        list(
            values = p$sdev^2,
            scores = scale(p$x),
            loadings = loadings_mat
        ),
        class = c("psych", "principal")
    )

    expect_invisible(bb_biplot(mock_principal))
    expect_invisible(bb_biplot(mock_principal, group = factor(iris$Species), legend = TRUE))
    expect_error(
        bb_biplot(structure(list(loadings = loadings_mat), class = c("psych", "principal"))),
        "scores = TRUE"
    )
})

test_that("methods handle zero-variance PCs and rotated principal Vaccounted cleanly", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    iris_collinear <- iris[-5]
    iris_collinear$dup <- iris_collinear$Sepal.Length
    p_zero <- stats::prcomp(iris_collinear, scale. = TRUE)
    sdev_zero <- p_zero$sdev
    sdev_zero[length(sdev_zero)] <- 0

    mock_pca_zero <- structure(
        list(
            eig = cbind(eigenvalue = sdev_zero^2),
            ind = list(coord = p_zero$x),
            var = list(coord = t(t(p_zero$rotation) * sdev_zero))
        ),
        class = c("PCA", "list")
    )
    res_pca_zero <- bb_biplot(mock_pca_zero, pc = c(4, 5), scale = 0)
    expect_false(anyNA(res_pca_zero$scores))
    expect_false(anyNA(res_pca_zero$loadings))

    mock_dudi_zero <- structure(
        list(
            eig = sdev_zero^2,
            li = as.data.frame(p_zero$x),
            co = as.data.frame(t(t(p_zero$rotation) * sdev_zero))
        ),
        class = c("pca", "dudi")
    )
    res_dudi_zero <- bb_biplot(mock_dudi_zero, pc = c(4, 5), scale = 0)
    expect_false(anyNA(res_dudi_zero$scores))
    expect_false(anyNA(res_dudi_zero$loadings))

    p <- stats::prcomp(iris[-5], scale. = TRUE)
    loadings_mat <- structure(t(t(p$rotation) * p$sdev), class = "loadings")
    v_acc <- rbind(
        "SS loadings" = c(2.0, 1.5),
        "Proportion Var" = c(0.50, 0.375)
    )
    mock_principal_vac <- structure(
        list(
            values = p$sdev^2,
            Vaccounted = v_acc,
            scores = scale(p$x),
            loadings = loadings_mat
        ),
        class = c("psych", "principal")
    )
    res_vac <- bb_biplot(mock_principal_vac)
    expect_equal(res_vac$pvar[1:2], c(0.50, 0.375), tolerance = 1e-12)
})
