# Ben's Better Biplot

A script to create beautiful looking biplots using base R. Designed to work with `prcomp()` and `princomp()` output for principal component analysis (pca), but can also work with any other data. The function is compatible with the base `biplot()` function, but offers significant enhancements including new features and full customisation of the biplot.

## Features

* Quickly create beautiful biplots from pca analysis.
* Uses base R - no need to install other packages.
* Compatible* with base biplot() function.
* Full customisation of the plot.
* Style the plot points based on groups (up to three), or style manually.
* Add legends (automatically generated).
* Add ellipses or convex hulls to grouped data.
* Add circle of equilibrium contribution.
* Plays nicely with layout(), mfrow(), par() and other plot functions.

## Get the code

You can load the code directly into R from Github using the following, which will always be the latest version:

```
# Load bb_biplot script from Github
source("https://raw.githubusercontent.com/benbell95/better_biplot/refs/heads/main/r/better_biplot.r")
```

Or download the [latest release](https://github.com/benbell95/better_biplot/releases/latest) to use locally. E.g.

```
# Load bb_biplot script from local
source("local script location/better_biplot.r")
```

## Example

Run pca analysis and plot the results. Scale the plot, style the observations by the species, include a legend and "circle of equilbrium" on the plot.

```
# pca
p <- prcomp(iris[-5], scale=TRUE)
# Group
gr <- as.factor(iris[,5])
# Plot
bb_biplot(p, scale=0, group=gr, lab_rotation=TRUE, circle.eq=TRUE, legend=TRUE, title.leg="Species")
```

![Better biplot example plot](https://github.com/benbell95/better_biplot/blob/90368041debbf2bd0c6c9dd9d4890c29dd510e29/better_biplot_example.svg)

## Compatibility

`bb_biplot()` is designed to be compatible with the output from base R `prcomp()` and `princomp()` functions, and create a biplot. It does not do the pca itself, rather it plots the results. For the plot, data is scaled in the same way as the base R `biplot()` function, and it also shares several arguments, so existing code should work with this function with only minor changes. 

A comparison between the arguments is shown in the table below, and a comparison between the plot output can be seen here: [bens_better_biplot_compare.pdf](https://github.com/benbell95/better_biplot/blob/7e23acd0966ee6800f6569325fa828901e9103dc/bens_better_biplot_compare.pdf).

| biplot() | bb_biplot() | Details |
| --- | --- | --- |
| x, y | x, y | x should be `prcomp()` or `princomp()` object (or see section below). |
| var.axes | variables | Logical argument to plot the second set of points as arrows. |
| col | col | Colours are generated automatically, or can be specified. |
| cex | cex.* | Size of data labels, variable labels, and data symbols can all be specified individually. |
| xlabs, ylabs | xlab, ylab | Axes labels (titles) are generated automatically, or can be specified. Only applies the main plot axes (e.g. bottom x axis, and left y axis). |
| expand | expand | Adjust the size of the variable arrows. |
| arrow.len | arrow.len | Adjust the size of the arrow heads, or suppress them. |
| xlim, ylim | limx | Plot limits for the x and y axes are generated automatically, but can be increased or decreased. |
| ... | ... | Additional arguments passed to the main plot. |

See [arguments section](#Arguments) for a full list of available arguments and an explanation of what they do / how they work.

Data is scaled for plotting when using `prcomp()` and `princomp()` objects. If you usually use `vegan` for pca analysis and biplots (but want to try out base R and this script), note that the scale value is reversed, e.g. scale = 0 in vegan equates to scale = 1 in base R and `bb_biplot()`, see table below.

| biplot() | bb_biplot() | vegan biplot.rda() |
| --- | --- | --- |
| scale = 0 | scale = 0 | scale = 1 |
| scale = 1 | scale = 1 | scale = 0 |

`bb_biplot()` can also scale the second set of data (variables) using `pc.biplot=TRUE` in the same way as `biplot()`. 

Additionally, you can apply a further rotation to the main data observations using `varimax.rotate=TRUE`. Currently, this uses the base R `varimax()` function with its default settings (see `?varimax` in R for details), and only applies to ` prcomp()` or `princomp()` objects.

You can also plot other data which was not created by `prcomp()` or `princomp()` as x and y values. `x` should be a matrix of at least two columns for the data (e.g. "scores"), and `y` should be a matrix of the variables (e.g. the rotation or loadings). Standard deviation values can also be supplied using `xsd` which are used to calculate the proportion of variance in the plot labels. This data will not be scaled by `bb_biplot()` and should be scaled before plotting.

`bb_biplot()` is a generic function, so it can be extended to plot other pca (or any other analysis object requiring a biplot) objects. See the source code for details.

## To biplot (or not to biplot?)

By default, the reduced dataset principal components (e.g. the data observations) are plotted, with the variables (loadings) plotted on top for comparison, hence the term "biplot". But, you can disable the second plot of variables using `variables=FALSE`. This can be useful if you have lots of variables, and you do not want the plot covered in variable arrows and/or labels!

You can also control which variables are plotted using the `whichv` argument, where you can specify either the name (exact match) of the variable, or the index position. Additionally, if plotting the circle of equilibrium, you can specify that only variable arrows which extend beyond the circle are plotted. Examples:

```
# Disable plotting of variables
bb_biplot(p, scale=0, group=gr, variables=FALSE)
# Plot variables, but only certain ones:
bb_biplot(p, scale=0, group=gr, whichv="Petal.Length")
# Plot variables, but only those extending beyond circle of equilibrium:
bb_biplot(p, scale=0, group=gr, whichv="circle.eq", circle.eq=TRUE)
```

## Grouping data 

`bb_biplot()` allows you to use up to three groups to style the plot points, with the arguments `group`, `group2`, and `group3`. These are optional, and you can style the plot without using groups, but using groups makes styling easier and more intuitive. 

Grouping data also allows [ellipses and/or convex hulls](#Ellipses-and-Convex-hulls) to be plotted, as well as allowing plot [legends](#Legends) to be added.

Group objects should be a factor, with a length that matches the number of data observations. For example, in the `iris` dataset, the fifth column "Species" is a grouping of the data, and this can be converted to a factor. `gr <- as.factor(iris[,5])`. If an object is supplied that is not a factor, the function will attempt to convert it.

`group` controls the colour of the plot points. However, it can also control the plot symbols if these are additionally specified, and `group2` is not used. Example:

```
# Group 
# Changes colour of the plot points depending on the group
bb_biplot(p, scale=0, group=gr)
# Additionally change the plot symbols for each group
bb_biplot(p, scale=0, group=gr, pch=c(21, 22, 23))
```

Colours are generated automatically, or you can specify them using the `col` argument.

`group2` controls the plot symbol (pch) for the plot points within each group, as specified by `group` and this must also be specified for `group2` to work. Example:

```
# Fake group for use in example
gr2 <- sample.int(6L, length(iris[,5]), replace=TRUE)
# Group 2
bb_biplot(p, scale=0, group=gr, group2=gr2)
```

Plot symbols are generated automatically, or you can specify them using the `pch` argument. The number of specified symbols should be the same as the number of groupings of data used in group2. A guide to plotting symbols in R is available on my blog: [Quick guide to pch symbols in R](https://www.benjaminbell.co.uk/2018/02/quick-guide-to-pch-symbols-in-r.html)

`group3` controls the size of the plot symbol (pch) for the plot points within each group, as specified by `group` and this must also be specified for `group3` to work, but `group2` is not required. Example:

```
# Fake group for use in example
gr3 <- sample.int(4L, length(iris[,5]), replace=TRUE)
# Group 3
bb_biplot(p, scale=0, group=gr, group2=gr2, group=gr3)
```

Plot symbol sizes are generated automatically when using `group3` based on the number of groupings of data used in group3. You can also control the symbol size manually by specifying `cex.pt` argument, however, `group3` will override this setting, so it should not be used if you wish to do it manually.

## Ellipses and Convex hulls

You can add either ellipses or convex hulls (or both!) to the biplot when `group` is specified, where the ellipse or convex hull will be drawn around the data points that make up the groupings. 

To add convex hulls to the plot, use:

```
# Convex hulls
bb_biplot(p, scale=0, group=gr, chull=TRUE)
```

The colour is based on the group, but you can customise the outer line width, and/or whether it uses a fill (background) colour, see [arguments section](#Arguments) for details.

To add ellipses, use:

```
# Ellipses
bb_biplot(p, scale=0, group=gr, ellipse=TRUE)
```

The colour is based on the group, but you can customise the outer line width, and/or whether it uses a fill (background) colour, see [arguments section](#Arguments) for details. Ellipses are automatically generated using the [parametric representation formula](https://en.wikipedia.org/wiki/Ellipse#Parametric_representation), and this is not always perfect. You can change how the angles are calculated by changing the `angle` argument to `angle="atan"` from `angle="lm"` which can sometimes improve the ellipse. You can also scale the ellipse using the `se` argument with is a multiplier.

## Legends

You can add a legend, or multiple legends to the plot based on the number of groups used. `group` is required for legends to work. These are generated automatically, and will match the style of each group, and legend text is taken from the `group` arguments (i.e. the different factors). Legend titles are also generated automatically, although these are simply "Group" followed by the number. You can specify your own titles for each group using the `title.leg` argument. Examples:

```
# Legends: 1 group
bb_biplot(p, scale=0, group=gr, legend=TRUE, title.leg="Species")
# 2 groups
bb_biplot(p, scale=0, group=gr, group2=gr2, legend=TRUE, title.leg=c("Species", "second title"))
# 3 groups
bb_biplot(p, scale=0, group=gr, group2=gr2, group3=gr3, legend=TRUE, title.leg=c("Species", "second title", "third"))
```

By default, the legends are placed in the right side margin. If there is only 1 group, the legend is vertically centred, but if more than 1 group, the legends are placed at the top of the y axis. 

You can change the location of the legend(s) by using the `lpx` and/or `lpy` arguments, to specify the absolute position using the main plot coordinates (even when variables are plotted), for example, to move the legends inside the plot area. You can also place the legend in the margins by extrapolating the likely coordinates, however, if using the bottom, left, or top margin, you must also change the plot margins using `par(mar=c())` function before calling `bb_biplot()`. When using the right margin, this should happen automatically. You can also horizontally align the different legends, rather than vertically align using `horiz=TRUE`. This only works with multiple groups. Examples:

```
# Legend inside plot (top left)
bb_biplot(p, scale=0, group=gr, group2=gr2, legend=TRUE, horiz=TRUE, lpx=-3.2)
# Legend outside plot (bottom margin)
par(mar=c(12,4,4,4))
bb_biplot(p, scale=0, group=gr, group2=gr2, legend=TRUE, horiz=TRUE, lpx=-3.2, lpy=-6.5)
```

You may need to resize the plot window and/or re-run the code to show this correctly.

## Full customisation

You can customise pretty much every aspect of the plot, beyond what is mentioned in this readme file. Refer to further examples or the arguments section below for full details. 

For example, try using using the `lab_rotation=TRUE` argument to change how labels are plotted.

You can also add other plot elements to the biplot using the various plot functions available in R. 

For example, use `text()` to add text annotations. Note that the coordinates to use relate to the variable plot (if plotted) i.e. top and right sides, or the main plot (if variables are not plotted) i.e. bottom and left sides. You can easily find out what coordinates to use by typing `par("usr")` into the R console.

## Further examples and help

More information and detailed examples on how to use this function are available on my blog (coming soon!): [benjaminbell.co.uk](https://www.benjaminbell.co.uk).

You can also call help after loading the script using `bb_info()` in the R console, and/or `bb_info(arg=TRUE)` to see a list of arguments and their default values (or see list below).

## Arguments

Details for the available arguments are shown in the table below. The list of arguments is long, but the function uses sensible defaults, and does a lot of styling automatically. You can create good biplots without delving into too many of these options, but they are there if needed, and offer full customisation of the plot.


| Argument | Details |
| -------------- | ------------------ |
| x | Observations. Matrix with at least 2 columns. Automatically determined for prcomp() / princomp() objects. |
| y | Variables. Matrix with at least 2 columns. Automatically determined for prcomp() / princomp() objects. |
| xsd | Optional. Vector of standard deviation values for the variables. Automatically determined for prcomp() / princomp() objects. (Used for automatic axis labels). |
| pc1 | The first principal component to plot (default = 1). Used when x has more than 2 columns. |
| pc2 | The second principal component to plot (default = 2). Used when y has more than 2 columns. |
| scale | Scale data between 0 and 1 for plotting (default = 1). This works the same way as R base biplot() function. Applies only to prcomp() or princomp() objects. |
| pc.biplot | Logical. Additional scaling of variables (default = FALSE). This works the same way as R base biplot() function.Applies only to prcomp() or princomp() objects. |
| varimax.rotate | Logical. Additionally apply varimax rotation (uses default settings, see ?varimax for help) (default = FALSE). Applies only to prcomp() or princomp() objects. |
| limx | xlim and ylim are generated automatically, but you can increase or decrease using a multiplier value (e.g. 1.5) |
| grid | Logical. Draw a grid through center x and y axes (default = TRUE). |
| col | Optional. Colours for the plot / groups. Automatically generates a palette if not supplied. |
| pch | Optional. Plotting symbols to use for the observations (also see group 2) (default symbol = 21). |
| cex.pt | Changes the size of the plotting symbols. |
| xlab | Overwrite the default x axis label. |
| ylab | Overwrite the default y axis label. |
| axes | Logical. Plot axes (This controls all axes including variables) (default = TRUE). |
| group | Optional. Groups for the data. Should be a factor, with length that matches original data. This affects the colour of the points and also allows for ellipses or convex hulls to be drawn. This will also affect the plot symbol used if multiple pch supplied, unless group2 is also specified where symbols will then relate to the second group. |
| group2 | Optional. Second grouping of data. Should be a factor, with length that matches original data. This affects the pch symbol only (e.g. multiple symbols within a single grouping). group must also be specified. |
| group3 | Optional. Third grouping of data. Should be a factor, with length that matches original data. This affects the size of the pch symbol (size automatically determined based on number of groups). group must also be specified (but does not need group2). |
| labd | Optional. A vector of labels for the data observations, this should be the same length as the data. Use NULL values to omit some labels. |
| col.labd | Colour for the labels (default = "black"). |
| cex.labd | Size of text labels.   |
| font.labd | Font type for data labels (1 = normal, 2 = bold, 3 = italic, 4 = bold italic).   |
| variables | Logical. Whether to plot the variables or not (default = TRUE). |
| vname | Optional. Vector of names for the variables. Uses rownames(y) by default (if exists), but, vname will override if both are supplied. |
| whichv | Specify which variables to plot. Either character vector, where names should exactly match the variables, or integer values to match index position. If specified as "circle.eq", it will only plot variables extending beyond the circle of equilibrium contribution, overriding any other variables - must also specify circle.eq=TRUE and scale=0 to work. |
| expand | Scale the variable arrows, if the arrows are too large or too small, change this value. Multiplier: values above 1 increase, while below 1 decrease the size. |
| arrow.len | Length of arrow head. Use 0 to suppress.  |
| lwd.v | Line width for arrows (default = 2). |
| col.v | Colour of the arrows (default = "#d12631").  |
| col.labv | Colour of the arrow labels (default = "black"). |
| cex.labv | Size of the text labels for the variables. |
| font.labv | Font type for labels (1 = normal, 2 = bold, 3 = italic, 4 = bold italic). |
| axes.v | Logical. Plot variable axes (default = TRUE).   |
| lab_rotation | Logical. Rotate labels for the data observations and variables (default = FALSE). Useful when you have lots of data points and labels. This may be slow if there are lots of labels to plot. |
| ran_adj | Logical. Alter the adjustment position of the labels for the data observations. |
| valign | Alignment of labels for the variables (accepted values = 0, 0.5, 1). |
| circle.eq | Logical. Plot circle of equilibrium contribution (default = FALSE). Only works when scale = 0. Takes colour from col.v value, and line width scales in proportion to lwd.v value.|
| chull | Logical. Add convex hulls to the plot - group must be supplied for this to work (default = FALSE). |
| ellipse | Logical. Draw ellipses around the data observations based on the group (default = FALSE). |
| angle | Method to determine ellipse angle, either "lm" or "atan" (default = "lm"). |
| autolim | Logical. Try to increase plot limits with large ellipses. Can also use limx to specify value (default = TRUE). |
| nofill | Logical. By default, ellipses and convex hulls have a background colour (default = FALSE). Change to TRUE to only show border outlines. |
| lwd.e | Line width for the ellipses and/or convex hulls (default = 2). |
| se  | Scale the size of the ellipses. Multiplier: values above 1 increase, while below 1 decrease the size (default = 0.7). |
| legend | Logical. Add a legend to plot. Only works when group is specified (default = FALSE). |
| title.leg | Titles for legend groups, should be a string with length of 1, 2 or 3 (depending on number of groups specified). Default is "Group", "Group 2" and "Group 3". |
| cex.leg | Size of legend (multiplier). |
| horiz | Logical. Plot legends next to each other - only works when multiple groups used (default = FALSE). Need to additionally specify lpx value to show properly (i.e. move legend position). Can also specify lpy for further control. |
| lpx | Optional. Absolute position of legend on x axis (data observations). |
| lpy | Optional. Absolute position of legend on y axis (data observations). |
| ... | Additional arguments passed to plot() (main plot only). |
