# Ben's better biplot!

A script to create beautiful looking biplots using base R, without the need to use or install any additional packages. Designed to work with `prcomp()` and `princomp()` output for principal component analysis (pca), but can also work with other data. The function is compatible with the base `biplot()` function, but offers significant enhancements including new features and full customisation of the biplot.

## Get the code

You can load the code directly into R from Github using the following:

```
# Load biplot code in R from Github
source("https://raw.githubusercontent.com/benbell95/better_biplot/refs/heads/main/r/better_biplot.r")
```

## Example

Run pca analysis and plot the results. Scale the plot, style the observations by the species, include a legend and "circle of equilbrium" on the plot.

```
# pca
p <- prcomp(iris[-5], scale=TRUE)
# Group
gr <- as.factor(iris[,5])
# Plot
bb_biplot(p, group=gr, scale=0, lab_rotation=TRUE, 
    circle.eq=TRUE, legend=TRUE, title.leg="Species")
```

![bb_biplot](https://github.com/user-attachments/assets/53c0e0d3-856b-4b85-9ef7-3db8f6b549e5)

## Details

# Features

* Quickly create beautiful biplots from pca analysis.
* Uses base R - no need to install other packages.
* Compatible with base biplot() function.
* Full customisation of the plot.
* Style the plot using up to three groups, to change colour, symbol, and size.
* Add ellipses or convex hulls to grouped data.
* Add circle of equilibrium contribution.
* Plays nicely with layout(), mfrow(), par() and other plot functions.

# Compatibility

`bb_biplot()` is designed to take the output from base R `prcomp()` or `princomp()` functions and create a biplot i.e. it does not do the pca itself, rather it plots the results. For the plot, data is scaled in the same way as the base R `biplot()` function, and it also shares several arguments, so existing code should work with this function with only minor changes. A comparison between the arguments is shown in the table below.

| biplot() | bb_biplot() | Details |
| --- | --- | --- |
| x, y | x | x should be `prcomp()` or `princomp()` object (or see section below) |
| var.axes | variables | Logical argument to plot the second set of points as arrows |
| pc.biplot | pc.biplot | Logical argument to apply additional scaling to second plot |
| col | col | Colours are generated automatically, or can be specified |
| cex | cex.* | Size of data labels, variable labels, and data symbols can all be specified individually |
| xlabs | xlab | Axes labels (titles) are generated automatically, or can be specified. Only applies the main plot (e.g. bottom x axis) |
| ylabs | ylab | Axes labels (titles) are generated automatically, or can be specified. Only applies the main plot (e.g. right y axis) |
| expand | expand | Adjust the size of the variable arrows |
| arrow.len | arrow.len | Adjust the size of the arrow heads, or suppress them |
| xlim, ylim | limx | Plot limits for the x and y axes are generated automatically, but can be increased or decreased |
| ... | ... | Additional arguments passed to the main plot |

See Arguments section for a full list of available functions and an explanation of what they do / how they work.

You can also plot other data which was not created by `prcomp()` or `princomp()`, this should be supplied as a list object, with the first item containing a matrix of at least two columns for the data (e.g. "scores"). The second list item should also contain a matrix of the variables (e.g. the rotation or loadings). A third list item can also be supplied, and this should contain the sd values. This data will not be scaled by `bb_biplot()` regardless of the scale setting, and should be scaled before plotting.

If you use 'vegan' for pca analysis and biplots, note that the scale value is reversed, e.g. scale = 0 in vegan equates to scale = 1 in base R and `bb_biplot()`, see table below.

| biplot() | bb_biplot() | 'vegan' biplot.rda() |
| --- | --- | --- |
| scale = 0 | scale = 0 | scale = 1 |
| scale = 1 | scale = 1 | scale = 0 |
