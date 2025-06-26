# Ben's better biplot!

A script to create beautiful looking biplots using base R, without the need to use or install any additional packages. Designed to work with prcomp() and princomp() output for principal component analysis (pca), but can also work with any other data. The function is compatible with the base biplot() function, but offers significant enhancements (new features) and full customisation of the biplot.

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
bb_biplot(p, group=gr, scale=0, lab_rotation=TRUE, circle.eq=TRUE, legend=TRUE, title.leg="Species")
```

![bb_biplot](https://github.com/user-attachments/assets/53c0e0d3-856b-4b85-9ef7-3db8f6b549e5)

