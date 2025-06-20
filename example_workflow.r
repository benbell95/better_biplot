# Example code for readme file

# set wd
setwd("~/code/r/graphics/better_biplot")


# Local copy
source("~/code/r/graphics/better_biplot/r/better_biplot.r")

## Test with iris
p <- prcomp(iris[-5])

gr <- as.factor(iris[,5])

layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p, group=gr, scale=0, main="bens biplot scale =0", lab_rotation=TRUE, labd=1:200)
bb_biplot(p, group=gr, scale=1, main="bens biplot scale =1", lab_rotation=TRUE, labd=1:200)



p1 <- prcomp(USArrests)

layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p1, cols="blue", scale=0, main="bens biplot scale =0", lab_rotation=TRUE, labd=rownames(USArrests))
bb_biplot(p1, cols="blue", scale=1, main="bens biplot scale =1", lab_rotation=TRUE, labd=rownames(USArrests))


layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p1, cols="blue", scale=1, main="bens biplot scale =0", lab_rotation=TRUE, labd=rownames(USArrests), ran_adj=TRUE)
bb_biplot(p1, cols="blue", scale=1, main="bens biplot scale =1", lab_rotation=TRUE, labd=rownames(USArrests), ran_adj=FALSE)




bb_biplot(p, group=gr, scale=0)


layout(matrix(1:4, ncol=2, byrow=TRUE))
bb_biplot(p, group=gr, scale=0, main="bens biplot scale =0")
bb_biplot(p, group=gr, scale=1, main="bens biplot scale =1")

# base r biplot

#layout(matrix(1:2, ncol=2))
biplot(p ,scale=0, main="base R biplot scale =0")
biplot(p, scale=1, main="base R biplot scale =1")





## Use vegan for comparison
# biplot.rda function [use plot()]
library(vegan)
pv <- rda(iris[-5])
plot(pv)


# Compare vegan plot to biplot
layout(matrix(1:8, ncol=4, byrow=TRUE))
plot(pv, main="vegan biplot.rda scale=0", scaling =0)
plot(pv, main="vegan biplot.rda scale=1", scaling =1)
plot(pv, main="vegan biplot.rda scale=2", scaling =2)
plot(pv, main="vegan biplot.rda scale=3", scaling =3)

