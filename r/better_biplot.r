##########################################################################
### Ben's Better Biplot
### Plot function for creating better looking biplots in R (using base R)
### Copyright 2025, Benjamin Bell.
### Code and updates: https://github.com/benbell95/better_biplot
### Blog: https://www.benjaminbell.co.uk
### Contact details via blog or github. Please send comments, suggestions, bug reports.

### This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

### This function has been tested with R v4.5.0. Use at your own risk.

### Plots a biplot of pca data created with base R prcomp() function. Does not require any other packages or ggplot2. 
### Should play nicely with layout(), and other plot functions, e.g. adding additional elements after plotting.

### See README.md file for usage. Arguments listed below. Full guide also available on my blog (coming soon).

# To do: 

# Look at the other princomp (base R) scaling option - compare to vegan
# Check Wills email, add the addtional functions he requested

# Add ability to add legend
# More customisation of plot elements
# Add no fill option for ellipse and chull
# Add checks for prcomp object
# Check compatability with princomp
# Change cols argument to col - and then fix function



# *** Look at the pca papers, need to add rotation to the data?

# Run PCA in R using prcomp()
# Save output to an object, and use bb_biplot to produce a better looking biplot compared to the default base R biplot() function, e.g. bb_biplot(mypca)
# All code in function is based on base R - does not require any other packages or ggplot2.
# Data is scaled based on base R biplot() function. e.g. https://github.com/wch/r-source/blob/trunk/src/library/stats/R/biplot.R


### Arguments / details

# x             = The pca (prcomp) object. You must run prcomp() on your data before using this function.
# pc1           = The first principal component to plot (default 1).
# pc2           = The second principal component to plot (default 2).

# scale         = Scale data between 0 and 1 for plotting. This works the same way as R base biplot() function.
# limx          = x and ylim are generated automatically, but you can increase or decrease using a multiplier value (e.g. 1.5)
# grid          = (logical). Draw a grid through center x and y axes.
# cols          = Colours for the plot / groups. Optional. Uses default palette if not supplied. 
# pch           = Plotting symbols to use for the observations (also see group 2)
# cex.pt        = Changes the size of the plotting symbols.

# group         = Groups for the data. Optional. Should be a factor, with length that matches original data. This affects the colour of the points  and also allows for ellipses or convex hulls to be drawn. This will also affect the plot symbol used if multiple pch supplied, unless group2 is also specified.
# group2        = Second grouping of data. Optional. This affects the pch symbol only (e.g. multiple symbols within a single grouping). 

# labd          = A vector of labels for the data observations, this should be the same length as the data. Use NULL values to omit some labels. Optional. 
# colld         = Colours for the labels. Default is black.
# cex.labd      = Size of text labels.

# vectors       = (logical). Whether to plot the vectors or not.
# whichv        = By default, all vectors are plotted, but you can instead specify which to plot. This can be a character vector, where names should exactly match the names of the vectors you want to plot, or a numeric vector to match index position.
# sv            = Scale the vectors, if the vectors are too large or too small, change this value. Multiplier - values above 1 increase, while below 1 decrease the size.
# colv          = Colour of the vector arrows.
# colvt         = Colour of the vector arrow labels.
# cex.labv      = Size of the text labels for the vectors.

# lab_rotation  = (logical). Rotate labels for the data observations and vectors. Useful when you have lots of data points and labels. This may be slow if there are lots of labels to plot.
# ran_adj     = (logical). Alter (adjustment) position of the labels for the data observations.
# valign        = Alignment of labels for the vectors (accepted values = 0, 0.5 (default), 1)

# chull         = (logical). Add convex hulls to the plot - group must be supplied for this to work.
# ellipse       = (logical). Draw ellipses around the data observations based on the group. 
# refine_ellipse = (logical). Checks the slope angle of the grouped data and flips ellipse if meets certain criteria. Change to TRUE if some of the ellipses seem like they should be the other way around.
# se            = Scale the size of the ellipses. Multiplier - values above 1 increase, while below 1 decrease the size.

# ...           = Additional arguments passed to plot() [main plot only].


### Ben's Better biplot function
bb_biplot <- function(x, pc1=1, pc2=2, scale=1, pc.biplot=FALSE, limx, grid=TRUE, cols, pch=21, cex.pt=1, group, group2, labd, colld="black", cex.labd=1, vectors=TRUE, whichv, sv=1, colv="red", colvt="black", cex.labv=1, lab_rotation=FALSE, ran_adj=FALSE, valign=0.5, chull=FALSE, ellipse=FALSE, refine_ellipse=FALSE, se=1, ...) {
    # Add check that x is a prcomp object [or princomp]
    # not tested with princomp()
    xclass <- class(x)

    #if(scale < 0 || scale > 1) warning("scale should be a value between 0 and 1]")
    
    ########################################
    ### Scale data
    # Code modified from R base biplot() function code [Copyright (C) 1995-2012 The R Core Team]
    # prcomp() scaling
    if(xclass=="prcomp") {
        n <- nrow(x$x)
        if(scale != 0) {
            lam <- (x$sdev[c(pc1, pc2)] * sqrt(n)) ^ scale 
        } else {
            lam <- 1
        }
        #if(pc.biplot==TRUE) lam <- lam / sqrt(n)      
        sco <- t(t(x$x[, c(pc1, pc2)]) / lam)
        rot <- t(t(x$rotation[, c(pc1, pc2)]) * lam)
    }

    ########################################
    ### Plot data / config
    # Get x/ylim values
    # Scores
    al <- ceiling(max(abs(range(sco)))*10) /10
    lim <- c(al*-1, al)
    # Rotation
    alr <- ceiling(max(abs(range(rot)))*10) /10
    limr <- c(alr*-1, alr)
    # If plotting ellipses, and scale > 0, the limits will be too small, so always increase (unless limx is specified)
    if(ellipse==TRUE & scale > 0 & hasArg(limx) == FALSE) {
        limx <- 1.5
    }
    # Increase x/y lim values?
    if(hasArg(limx) == TRUE) {
        lim <- lim * limx
        limr <- limr * limx
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
        lsrt <- .bb_lab_rotation(x=sco, valign=valign, ran_adj=ran_adj)
        lvsrt <- .bb_lab_rotation(x=rot, valign=valign, ran_adj=FALSE)
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
    op <- par(pty="s")
    on.exit(par(op))
    plot(sco[,1], sco[,2], type="n", asp=1, xlim=lim, ylim=lim, xlab=xla, ylab=yla, ...)
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
            # srt rotation is not vectorised, so must run in a loop (Can be slow if lots of labels)
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
        # New plot window
        par(new=TRUE)
        plot(rot[,1], rot[,2], type="n", asp=1, xlim=limr, ylim=limr, ann=FALSE, axes=FALSE)
        axis(3)
        axis(4)
        # Get vector co-ordinates
        vx <- rot[,1]
        vy <- rot[,2]
        # Add
        arrows(x0=0, x1=vx*sv, y0=0, y1=vy*sv, col=colv, length=0.15, lwd=2, xpd=TRUE)
        # Labels
        if(lab_rotation==TRUE) {
            # srt rotation is not vectorised, so must run in a loop (Can be slow if lots of labels)
            for(i in 1:length(rot[,1])) {
                text(rot[,1][i], rot[,2][i], labels=rownames(rot)[i], col=colvt, srt=lvsrt[,1][i], adj=c(lvsrt[,2][i], lvsrt[,3][i]), cex=cex.labv, xpd=TRUE) 
            }
        } else {
            # Random positions
            adjv<- sample.int(4L, length(vx), replace=TRUE) 
            text(vx*sv, vy*sv, labels=rownames(rot), col=colvt, pos=adjv, cex=cex.labv, xpd=TRUE)
        }
    }
    ### Legend
    # Everybody needs a legend
}


### Label rotation function
# x = matrix of xy values for data (scores)
### Remove this!
.bb_lab_rotation_OLD <- function(x, alter_dir) {
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

### Vector label rotation function
# x = matrix of xy values for vector arrows (x$rotation)
.bb_lab_rotation <- function(x, valign, ran_adj) {
    # Angle in degrees
    x_an <- atan(x[,2] / x[,1]) * 180/pi
    # text adj value
    ad1 <- rep(0, times=nrow(x))
    ad2 <- rep(valign, times=nrow(x))
    # Adjust based on quadrant
    bl <- which(x[,1] < 0 & x[,2] < 0)
    tl <- which(x[,1] < 0 & x[,2] > 0)
    # bottom left
    if(length(bl)>0) {
       ad1[bl] <- 1
    }
    # top left
    if(length(tl)>0) {
        ad1[tl] <- 1
    }
    # Random adjustment
    if(ran_adj==TRUE) {
        ad1 <- sample(c(0, 0.5, 1), length(ad1), replace=TRUE) 
        ad2 <- sample(c(0, 0.5, 1), length(ad2), replace=TRUE) 
    }
    # Combine and return data
    df <- cbind(x_an, ad1, ad2)
    return(df)
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
        # ******** Check this with more data sets - needs refinement ********
        cslm <- lm(x[,1] ~ x[,2], data=x)
        if(coef(cslm)[2] < 0 & coef(cslm)[2]> -1) {
            x1 <- xc + (a*cos(t)*cos(phi) * -1) - b*sin(t)*sin(phi)
        }
    }
    # return data for plotting
    return(data.frame(x=x1, y=y1))
}
