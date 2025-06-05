##########################################################################
### Ben's better biplot
### Plot function for creating better looking biplots in R (using base R)
### Copyright 2025, Benjamin Bell.
### Code and updates: https://github.com/benbell95/better_biplot
### Blog: https://www.benjaminbell.co.uk

### This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

### This function has been tested with R v4.4.3. Use at your own risk.

### See README.md file for usage.


# To do: Add ability to add legend
# More customisation of plot elements


# Run PCA in R using prcomp()
# Save output to an object, and use bb_biplot to produce a better looking biplot compared to the default base R biplot() function, e.g. bb_biplot(mypca)
# All code in function is based on base R - does not require any other packages or ggplot2.
# Data is scaled based on base R biplot() function. e.g. https://github.com/wch/r-source/blob/trunk/src/library/stats/R/biplot.R

# See Readme file for details of usage, or visit my blog for examples.


### Ben's better biplot
### Arguments (options)
# x = the pca (prcomp) object
# pc1 = the first principal component to plot (default 1, but can change)
# pc2 = the second principal component to plot (default 2, but can change)
# scale = scale data between 0 and 1 (as per base biplot() function)
# limx = x and ylim are generated automatically, you can increase/decrease using a multiplier value (e.g. 1.5)
# grid (logical) = draw a grid through center x and y axes.
# cols = Optional. Colours for the groups. (can be omitted to use default colour palette)
# pch = plotting symbols to use (also see group 2)
# group = Optional, but recommended! Groups for the data, should be a factor. (e.g. could be the sample, or species), length must match length of original data. This affects the colour of the points and allows for ellipses or convex hulls to be drawn. Will also affect the symbol used (if multiple pch supplied), unless group2 is also specified.
# group2 = Optional. Second grouping of data - this only applies to the pch symbol used to plot the data.
# labd = Optional. A vector of labels for the data observations.
# lab_rotation (logical). = Rotate labels for the data observations. This can be slow if lots of labels.
# vectors (logical) = Whether to draw vectors or not
# whichv = Optional. Choose which vectors to draw (default is all of them). This can be a character vector, which should exactly match the names of the vectors you want to draw, or a numeric vector to match index position.
# sv = scale vectors - change value if they appear off the biplot.
# colv = Colour of the vector arrows
# colvt = colour of the vector arrow labels
# chull (logical) = Add convex hulls to the groups.
# ellipse (logical) = Draw ellipses around the data observations based on group.
# refine_ellipse (logical) = Checks slope angle of data and flips ellipse if meets certain criteria. Change to FALSE if weird results, or takes too long to plot.
# cex.pt = changes size of points.
# ... = additional arguments (applied to plot())

bb_biplot <- function(x, pc1=1, pc2=2, scale=1, limx, grid=TRUE, cols, pch=21, group, group2, labd, colld="black", cex.labd=1, lab_rotation=FALSE, vectors=TRUE, whichv, sv=1, colv="red", colvt="black", chull=FALSE, ellipse=FALSE, refine_ellipse=TRUE, se=1, cex.pt=1, cex.labv=1, ...) {
    ########################################
    ### Scale data as per R base biplot() function
    n <- nrow(x$x)
    lam <- (x$sdev[c(pc1, pc2)] * sqrt(n)) ^ scale
    rot <- t(t(x$rotation[, c(pc1, pc2)]) * lam)
    sco <- t(t(x$x[, c(pc1, pc2)]) / lam)
    ########################################
    ### Plot data / config
    # Get x/ylim values
    al <- ceiling(max(abs(range(sco)))*10) /10
    lim <- c(al*-1, al)
    # If plotting ellipses, and scale > 0, the limits will be too small, so always increase (unless limx is specified)
    if(ellipse==TRUE & scale > 0 & hasArg(limx) == FALSE) {
        limx <- 1.5
    }
    # Increase x/y lim values?
    if(hasArg(limx) == TRUE) {
        lim <- lim * limx
    }
    # Convex hulls
    if(chull==TRUE) {
        # Needs group to work (and must be a factor)
        if(!hasArg(group)) stop("Error! group must be supplied to draw convex hulls")
        # Add group (factor) data to scores and convert to data.frame
        sco <- cbind(sco, group) |> as.data.frame()
        # Get the outer points (for plotting)
        op <- lapply(1:nlevels(group), \(x) chull(sco[sco$group == x,]))
    }
    # Ellipsis
    if(ellipse==TRUE) {
        # Needs group to work (and must be a factor)
        if(!hasArg(group)) stop("Error! group must be supplied to draw ellipses")
        # Add group (factor) data to scores and convert to data.frame
        sco <- cbind(sco, group) |> as.data.frame()
        # Call function to calculate ellipses
        e <- lapply(1:nlevels(group), \(x) .bb_ellipse(x=sco[sco$group == x,], se=se, refine_ellipse=refine_ellipse))
    }
    # Label rotation
    if(lab_rotation==TRUE) {
        lsrt <- .bb_lab_rotation(x=sco)
    }
    ### Colour/style points by group (factor)
    # Overly complex due to outline and fill colours
    # Check if group exists
    if(hasArg(group) == TRUE) {
        # Check length of group is correct
        if(length(group)!=n) stop("Error! group must be same length as original data used in PCA")
        # Check group is a factor and/or convert
        if(class(group)!="factor") {
            group <- as.factor(group)
            message(paste0("Info: group converted to factor with ", nlevels(group), " levels"))
        }
        # Check if colors supplied
        if(hasArg(cols)==TRUE) {
            # Check length of colors 
            if(length(cols) < nlevels(group)) message("Warning: Number of colours does not match number of groups")
            c <- cols
        } else {
            # Generate default colour palette
            c <- hcl.colors(nlevels(group), "Dark 3")
        }
        # Check pch (plotting symbols)
        # Set pch if group supplied and group 2 doesn't exist
        if(length(pch) < nlevels(group) & !hasArg(group2)) {
            pch <- rep(21, times=nlevels(group))
        }       
        # If group 2 is available
        if(hasArg(group2)) {
            if(length(pch) < nlevels(group2)) {
            # Create symbols based on number of groups
            pch <- seq(21, length.out=nlevels(group2))
            # Check if any values above 25 (reset to 1)
            pch <- ifelse(pch > 25, pch - 25, pch)
            }
        }       
        # Check which pch symbol used, to correctly apply colours
        ### ******** This needs a rework/rethink to account fo multiple symbols *******
        if(as.integer(pch[1]) > 20) {
            # Symbol has a fill (bg)
            col1 <- rep("black", times=nlevels(group))
            col2 <- c
        } else {
            # Symbol is solid, or outline (col)
            col1 <- c
            col2 <- c
        }
    } else { 
        # No group supplied, but colour supplied
        if(hasArg(cols)==TRUE) {
            col1 <- cols
            col2 <- cols
        } else {
        # No group or colour supplied
        message("Info: No group argument supplied, using default colour for plot")
        col1 <- "black" 
        col2 <- "black"
        }
    }
    # Axis labels
    s <- summary(x)
    xla <- paste("PC", pc1, " (", round(s$importance[,pc1][2]*100, 1), "%)", sep = "")
    yla <- paste("PC", pc2, " (", round(s$importance[,pc2][2]*100, 1), "%)", sep = "")
    ########################################
    ### Main Plot
    # Initial blank plot (plot order matters!)
    plot(sco[,1], sco[,2], type="n", asp=1, xlim=lim, ylim=lim, xlab=xla, ylab=yla, ...)
    axis(3)
    axis(4)
    # Plot convex hulls
    if(chull==TRUE) {
        for(i in 1:nlevels(group)) {
            polygon(sco[sco$group == i,][op[[i]],], col=adjustcolor(c[i], 0.25), border=adjustcolor(c[i], 0.9), lwd=2)            
        }
    }
    # Plot ellipses
    if(ellipse==TRUE) {
        for(i in 1:nlevels(group)) {
            polygon(x=e[[i]]$x, y=e[[i]]$y, col=adjustcolor(c[i], 0.25), border=adjustcolor(c[i], 0.5))
        } 
    }
    # Plot grid
    if(grid==TRUE) {
        abline(v=0, lty=2, lwd=1, col=adjustcolor("black", 0.75))
        abline(h=0, lty=2, lwd=1, col=adjustcolor("black", 0.75))
    }
    # Add data points (pca scores)
    if(hasArg(group2)) {
        points(sco[,1], sco[,2], pch=pch[group2], col=adjustcolor(col1[group], 0.85), bg=adjustcolor(col2[group], 0.85), cex=cex.pt)
    } else {
        points(sco[,1], sco[,2], pch=pch[group], col=adjustcolor(col1[group], 0.85), bg=adjustcolor(col2[group], 0.85), cex=cex.pt)
    }
    # Add labels
    if(hasArg(labd)) {
        # Check length matches data
        if(length(labd) != length(sco[,1])) message("Info: Number of labels does not match number of data observations")
        if(lab_rotation==TRUE) {
            # rotation is not vectorised, so must run in a loop (Can be slow if lots of labels)
            for(i in 1:length(sco[,1])) {
                text(sco[,1][i], sco[,2][i], labels=labd[i], col=colld, srt=lsrt[,1][i], adj=c(lsrt[,2][i], lsrt[,2][i]), cex=cex.labd, xpd=TRUE) 
            }
        } else {
            # Random positions
            adjld <- sample.int(4L, length(sco[,1]), replace=TRUE) 
            text(sco[,1], sco[,2], labels=labd, col=colld, pos=adjld, cex=cex.labd, xpd=TRUE)
        }
    }
    ### Add vectors (arrows)
    if(vectors==TRUE) {
        if(hasArg(whichv)) {
            # Subset vectors
            # Get vector names
            vname <- rownames(x$rotation)
            # Match and subset
            if(class(whichv)=="character") {
                match <- which(vname %in% whichv)
            } else {match <- whichv}
            # Check if any matches
            if(length(match)==0) stop("Error! Supplied vectors do not match any vector labels")
            if(length(match)>0 & length(match)!=length(whichv)) message("Info: Some vectors do not match vector labels")
            # Subset
            rot <- rot[match,]
            # If only 1 match, object is converted to numeric, so convert back (otherwise breaks rest of function)
            # *** Hacky - make this more elegant solution *** 
            if(length(rot)==2) {
                rot <- data.frame(PC1=rot[1], PC2=rot[2])
                rownames(rot) <- vname[match]
            }
        } 
        # Get vector co-ordinates
        vx <- rot[,1]
        vy <- rot[,2]
        # Add
        arrows(x0=0, x1=vx*sv, y0=0, y1=vy*sv, col=colv, length=0.15, lwd=2, xpd=TRUE)
        # Labels
        adjv<- sample.int(4L, length(vx), replace=TRUE) 
        text(vx*sv, vy*sv, labels=rownames(rot), col=colvt, pos=adjv, cex=cex.labv, xpd=TRUE)
    }
}


### Label rotation function
# x = matrix of xy values for data (scores)
.bb_lab_rotation <- function(x, alter_dir=TRUE) {
    # Get orders of original data (by y)
    x_or <- order(x[,2])
    # Create rotation sequence for all data
    # Round down to ensure enough labels are made
    x_r <- seq(-90, 90, by=floor(abs(180 / length(x[,1]))*10) / 10)
    # trim extra
    x_r <- x_r[1:length(x[,1])]
    # re-order rotation sequences
    x_r1 <- x_r[order(x_or)]
    # Create adjustment values
    ad <- rep(0, times=length(x_r1))
    if(alter_dir==TRUE) {
        # sequence
        as <- seq(1, length(x[,1]), 2)
        # Alter values
        x_r1[as] <- x_r1[as] * -1
        ad[as] <- 1
    }
    # Return data
    return(cbind(x_r1, ad))
}  

### Ellipse function
# x should be 2 column df with x/y
# se = scale factor for ellipse
.bb_ellipse <- function(x, se, refine_ellipse) {
    # Parametric representation
    # (x, y) = (a cos t, b sin t), 0 < t < 2pi
    # Adapted from: https://stackoverflow.com/questions/41820683/how-to-plot-ellipse-given-a-general-equation-in-r
    # and: https://stats.stackexchange.com/questions/9898/how-to-plot-an-ellipse-from-eigenvalues-and-eigenvectors-in-r

    # Lengths of ellipse major/minor axis
    len <- apply(x[1:2], 2, range)
    len <- len[2,] - len[1,]
    # Lengths of ellipse axes
    a <- len[which(max(len)==len)] * se # Major (longest)
    b <- len[which(min(len)==len)] * se # Minor
    # Get centres
    # This is just the mean of the co-ordinates
    cen <- apply(x, 2, mean)
    xc <- cen[1] # x centre
    yc <- cen[2] # y centre
    # Angle of ellipse major axis to x axis
    phi <- atan(a / b)
    # Plot points
    t <- seq(0, 2 * pi, 0.1) 
    x1 <- xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
    y1 <- yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi)
    # Refine
    if(refine_ellipse==TRUE) {
        # Check slope of data, and flip ellipse if negative between -0 and -1
        cslm <- lm(x[,1] ~ x[,2], data=x)
        if(coef(cslm)[2] < 0 & coef(cslm)[2]> -1) {
            x1 <- xc + (a*cos(t)*cos(phi) * -1) - b*sin(t)*sin(phi)
        }
    }
    # return data for plotting
    return(data.frame(x=x1, y=y1))
}
