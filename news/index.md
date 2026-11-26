# Changelog

## betterbiplot 2.0.1

### Changes & Improvements

- Moved `title.leg` into `leg.args = list(title = ...)` (top-level
  `title.leg` remains accepted with a deprecation warning).
- Reorganised `R/` source files (`R/bb_biplot.R`, `R/methods.R`,
  `R/utils.R`) and simplified S3 methods to use `pc` directly.
- Removed the `methods` import by replacing `hasArg()` with explicit
  `!missing(...) && !is.null(...)` checks.
- Updated
  [`bb_biplot.principal()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
  to use rotated variance proportions from
  `Vaccounted["Proportion Var", ]` when available.
- Updated GitHub Actions workflows to `actions/checkout@v6` (Node.js 24)
  and moved README SVG figures into `man/figures/` so they display on
  the `pkgdown` site.

### Bug Fixes

- **Varimax rotation consistency:** Fixed `varimax.rotate = TRUE` to
  rotate both observation scores and variable loadings consistently
  (before `scale`/`pc.biplot` scaling), preserving inner products and
  updating axis `%` variance labels.
- **`autolim` ellipse limit expansion:** Fixed the limit multiplier
  (`ceiling((ecm / max(lim)) * 10) / 10`) so large ellipses (such as
  small-sample Hotelling’s T-squared ellipses) are no longer clipped.
- **Convex hulls and ellipses on `data.frame` and `NA` inputs:** Indexed
  `x` directly by `as.integer(group)` and filtered non-finite
  coordinates, fixing silent failures when `x` is a `data.frame` and
  crashes when `x` or `group` contains `NA`.
- **`labd` with `pos = "auto"` and `NA` coordinates:** Added a fallback
  for `NA` positions so [`text()`](https://rdrr.io/r/graphics/text.html)
  no longer throws `invalid 'pos' value`.
- **Ungrouped `group_cols` and vector `col` with `pch = 21`:** Returned
  the fill colour (`colf`) rather than `"black"` in `res$group_cols`
  when ungrouped, and applied black outlines consistently across all
  observations when passing per-observation `col` with `pch = 21:25`.
- **Zero-variance PCs:** Guarded against division by zero in
  `.bb_scale_pca()`,
  [`bb_biplot.PCA()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md),
  [`bb_biplot.pca()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md),
  and `pvar` so collinear/zero-variance components no longer produce
  `NaN`/`Inf`.
- **Input validation & explicit `NULL` arguments:** Validated `pc` as
  positive integers and handled explicit `NULL` optional arguments
  cleanly as omitted.

## betterbiplot 2.0.0

### Package Conversion

- Converted the standalone `r/better_biplot.r` script into the
  `betterbiplot` R package with proper `DESCRIPTION`, `NAMESPACE`,
  roxygen2 documentation
  ([`?bb_biplot`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)),
  and `testthat` unit tests.
- Retained a stub `r/better_biplot.r` (excluded from package builds via
  `.Rbuildignore`) that prints a migration notice when
  [`source()`](https://rdrr.io/r/base/source.html)d directly from
  GitHub.
- Updated `README.md` with package installation instructions
  (`remotes::install_github("benbell95/better_biplot")`) and fixed minor
  argument typos (`group3`, `f1.a`).

### New Features

- Added S3 methods for additional PCA objects without adding hard
  dependencies:

  - [`bb_biplot.rda()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
    for unconstrained
    [`vegan::rda()`](https://vegandevs.github.io/vegan/reference/cca.html)
    PCA models (converted to base R `scale` for full compatibility with
    [`prcomp()`](https://rdrr.io/r/stats/prcomp.html) /
    [`biplot()`](https://rdrr.io/r/stats/biplot.html)).
  - [`bb_biplot.svyprcomp()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
    for
    [`survey::svyprcomp()`](https://rdrr.io/pkg/survey/man/svyprcomp.html)
    objects, supporting `group3 = "scaled"` to scale point sizes by
    survey weights (`sqrt(w / max(w))`).
  - [`bb_biplot.PCA()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
    for `FactoMineR::PCA()` objects.
  - [`bb_biplot.pca()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
    for `ade4::dudi.pca()` objects.
  - [`bb_biplot.principal()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
    for `psych::principal()` objects.

- Added selective group legends: `legend` now accepts `TRUE` / `"all"`
  (all supplied groups) or integer values in `1:3` (e.g. `1`, `c(1, 3)`,
  `3`) to display only selected group legends, with positional or
  group-indexed `title.leg` matching.

- Condensed `group3` legends to always show minimum, middle, and maximum
  symbol sizes (with `y.intersp = 1.25`), and wired
  `survey::svyprcomp(..., group3 = "scaled")` to also display a
  condensed min/median/max weight legend when legends are enabled.

- Streamlined
  [`bb_biplot()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
  arguments by grouping fine-grained styling into `lab.args`,
  `var.args`, `ellipse.args` (also accepted via alias `chull.args`), and
  `leg.args`, and combining `pc1`/`pc2` into `pc = c(1, 2)` (while
  keeping `title.leg` and `circle.eq` top-level).

- Added a soft deprecation shim for all v1.x flat argument names (`pc1`,
  `pc2`, `col.v`, `whichv`, `expand`, `lab_rotation`, `angle`, `se`,
  `cex.leg`, `lpx`, `f1.a`, etc.) so existing scripts continue to run
  with deprecation warnings.

- Replaced the v1.x parametric ellipse approximation (`angle`, `se`)
  with statistical covariance confidence ellipses
  (`ellipse.args = list(method = "cov", level = 0.95)`, default) and
  small-sample Hotelling’s T-squared ellipses (`method = "t"`).

- Added group centroids (`centroids = TRUE`, `"points"`, or `"labels"`)
  and radial spider lines (`spider = TRUE`).

- Added `var.args = list(lab.offset = 1.05)` to nudge variable labels
  past arrow tips.

- Added palette name shorthand in `col` (any
  [`grDevices::palette.pals()`](https://rdrr.io/r/grDevices/palette.html)
  or
  [`grDevices::hcl.pals()`](https://rdrr.io/r/grDevices/palettes.html)
  palette, e.g. `col = "Okabe-Ito"` or `"Viridis"`).

- Added configurable point opacity (`alpha = 0.85`) and hull/ellipse
  opacity (`ellipse.args = list(alpha = 0.25, border.alpha = 0.9)`).

- Made unrotated observation label placement deterministic by default
  (`lab.args = list(pos = "auto")`, also supporting `"random"` or
  `1:4`).

- [`bb_biplot()`](https://benbell95.github.io/better_biplot/reference/bb_biplot.md)
  now invisibly returns a `"bb_biplot"` object (`scores`, `loadings`,
  `lim`, `limr`, `pc`, `pvar`, `group_cols`, `centroids`) for downstream
  base R annotations.

- Added `varimax.args` (`normalize = TRUE`, `eps = 1e-5`) to customise
  [`stats::varimax()`](https://rdrr.io/r/stats/varimax.html) when
  `varimax.rotate = TRUE`.

- Added package vignette
  [`vignette("betterbiplot")`](https://benbell95.github.io/better_biplot/articles/betterbiplot.md)
  with rendered examples across scaling, grouping, overlays, legends,
  [`vegan::rda()`](https://vegandevs.github.io/vegan/reference/cca.html),
  [`survey::svyprcomp()`](https://rdrr.io/pkg/survey/man/svyprcomp.html),
  and post-plot extensions.

### Bug Fixes

- **Automatic `xlab`/`ylab` generation:** Stopped forwarding `...` into
  [`paste0()`](https://rdrr.io/r/base/paste.html) (which caused errors
  when graphical parameters such as `main` were passed) and handled
  missing/`NULL` `xsd` cleanly instead of producing empty labels.
- **`whichv` with `lab_rotation = TRUE`:** Subsetted variable rotation
  angles (`lvsrt`) alongside `y` so subsetted variable labels receive
  their own rotation angles instead of the first rows’ angles.
- **`whichv = "circle.eq"`:**
  - Explicitly required `circle.eq = TRUE` instead of leaving `match`
    unassigned (which previously collided with
    [`base::match`](https://rdrr.io/r/base/match.html)).
  - Checked the 2D vector radius
    `sqrt((y[, 1] * expand)^2 + (y[, 2] * expand)^2) > cr` instead of
    only `abs(y[, 2]) > cr`, keeping `whichv = "circle.eq"` consistent
    when `expand != 1`.
- **Single-variable `whichv` subsetting:** Used
  `y[vmatch, , drop = FALSE]` to preserve the 2-column matrix structure
  cleanly.
- **Rotated observation labels:** Fixed
  `adj = c(lsrt[i, 2], lsrt[i, 3])` (previously passed column 2 twice,
  ignoring `valign` and `ran_adj` vertical alignment).
- **`f1` and `f2` callbacks:** Replaced `as.list(args(f1))` with
  [`list()`](https://rdrr.io/r/base/list.html) when `f1.a`/`f2.a` are
  omitted, fixing crashes when calling functions that accept `...` (such
  as [`graphics::grid()`](https://rdrr.io/r/graphics/grid.html)).
- **`group2` symbol assignment:** Avoided expanding default `pch` for
  `group` before checking `group2`, ensuring distinct symbols are
  generated for `group2` even when `nlevels(group) >= nlevels(group2)`.
  Added length/group validation for `group2` and `group3`.
- **Single-level `group` legends:** Handled `nlevels(group) == 1` when
  computing vertical legend spacing so subsequent `group2`/`group3`
  legends do not receive `NA` y-coordinates.
- **Ellipses and convex hulls edge cases:** Handled unused factor levels
  and constant-x groups cleanly in ellipses and convex hulls.
- **Label rotation edge cases:** Avoided `NaN` angles at `(0, 0)` and
  applied left-side text justification (`ad1 = 1`) for all `x < 0`
  (including the negative x-axis where `y == 0`).
- **Modern R checks & namespace hygiene:** Replaced `class(x) == ...` /
  `class(x) != ...` with
  [`is.factor()`](https://rdrr.io/r/base/factor.html) /
  [`is.character()`](https://rdrr.io/r/base/character.html), avoided
  shadowing [`base::c()`](https://rdrr.io/r/base/c.html) and
  [`base::match()`](https://rdrr.io/r/base/match.html), and explicitly
  imported all non-`base` functions (`graphics`, `grDevices`, `methods`,
  `stats`).
