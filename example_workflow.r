# Example code for readme file

# set wd
setwd("~/code/r/graphics/better_biplot")


# Local copy
source("~/code/r/graphics/better_biplot/r/better_biplot.r")



## Test with iris
p <- prcomp(iris[-5])
prin <- princomp(iris[-5])

gr <- as.factor(iris[,5])
gr2 <- sample.int(6L, length(iris[,5]), replace=TRUE)


layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p, group=gr, pch=c(21:23), scale=0, lab_rotation=TRUE, cex.pt=1)
bb_biplot(p, group=gr, scale=0, lab_rotation=TRUE, cex.pt=1)


.bb_info()
bb_biplot(p, group=iris[,5])




# Local copy
source("~/code/r/graphics/better_biplot/r/better_biplot.r")


bb_biplot(p, group=gr, scale=0, ellipse=TRUE, vectors=FALSE, angle="lm")


bb_biplot(p, group=gr, scale=0, ellipse=TRUE, angle="atan")

bb_biplot(p, group=gr, scale=1, ellipse=TRUE, vectors=FALSE)
bb_biplot(p, group=gr, scale=1, ellipse=TRUE, vectors=FALSE)

bb_biplot(p, group=gr, scale=1, ellipse=TRUE, vectors=FALSE, angle="atan")
bb_biplot(p, group=gr, scale=1, ellipse=TRUE, vectors=FALSE, angle="atan")


bb_biplot(p, group=gr, scale=1, lab_rotation=TRUE, ellipse=TRUE)




bb_biplot(p, group=gr, group2=gr2, group3=gr, scale=0, labd=1:200, legend=TRUE, ellipse=TRUE)






layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p, group=gr, scale=1, lab_rotation=TRUE, ellipse=TRUE, limx=1.6)
bb_biplot(p, group=gr2, scale=1, lab_rotation=TRUE, ellipse=TRUE, limx=1.6)





layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p, group=gr, scale=0, main="bens biplot scale =0", lab_rotation=TRUE, labd=1:200, font.labd=2, font.labv=3)
bb_biplot(p, group=gr, scale=1, main="bens biplot scale =1", lab_rotation=FALSE, labd=1:200, font.labd=4, font.labv=3)


layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p, group=gr, scale=0, main="bens biplot scale =0", lab_rotation=TRUE)
bb_biplot(p1, group=gr, scale=1, main="bens biplot scale =1", lab_rotation=TRUE)

p1 <- prcomp(USArrests)
prin1 <- princomp(USArrests)

s1 <- summary(p1)
sprin1 <- summary(prin1)


layout(matrix(1:4, ncol=2, byrow=TRUE))
bb_biplot(p1, scale=0, lab_rotation=TRUE, expand=0.5)
bb_biplot(prin1, scale=0, lab_rotation=TRUE, expand=1.5)

biplot(p1 ,scale=0, expand=0.5)
biplot(prin1, scale=0, expand=1.5)


layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p1, col="blue", scale=0, main="bens biplot scale =0", lab_rotation=TRUE, labd=rownames(USArrests), lwdv=5, colvt="green")
bb_biplot(p1, col="blue", scale=1, lab_rotation=TRUE, labd=rownames(USArrests), vector.axes=FALSE)

layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p1, col="blue", scale=1, main="bens biplot scale =0", lab_rotation=TRUE, labd=rownames(USArrests), ran_adj=TRUE)
bb_biplot(p1, col="blue", scale=1, main="bens biplot scale =1", lab_rotation=TRUE, labd=rownames(USArrests), ran_adj=FALSE)



bb_biplot(p1, cex.pt=2, col=rainbow(10), scale=1, main="bens biplot scale =1", lab_rotation=TRUE, pch=c(1:5, 21:25))

# Problem is when more colours than pch

layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p1, cex.pt=2, col="gold", scale=1, main="bens biplot scale =1", lab_rotation=TRUE, legend=TRUE)

bb_biplot(p1, cex.pt=2, col="gold", scale=1, main="bens biplot scale =1", lab_rotation=TRUE, pch=21:25)


layout(matrix(1:2, ncol=2, byrow=TRUE))
bb_biplot(p1, cex.pt=2, scale=1, main="bens biplot scale =1", lab_rotation=TRUE, pch=15:25)

bb_biplot(p1, cex.pt=2, col=rainbow(50), pch=15:25, scale=1, main="bens biplot scale =1", lab_rotation=TRUE)


summary(p1)
s1$n.obs

bb_biplot(p, group=gr, scale=0)


layout(matrix(1:4, ncol=2, byrow=TRUE))
bb_biplot(p, group=gr, scale=0, main="bens biplot scale =0")
bb_biplot(p, group=gr, scale=1, main="bens biplot scale =1")

# base r biplot

layout(matrix(1:2, ncol=2))
biplot(p ,scale=0, main="base R biplot scale =0")
biplot(p, scale=1, main="base R biplot scale =1")


layout(matrix(1:2, ncol=2))
biplot(p ,scale=0, main="base R biplot scale =0", pc.biplot=TRUE)
biplot(p, scale=0, main="base R biplot scale =0")






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

