##########################################################################
### Ben's Better Biplot
### Plot function for creating better looking biplots in R (using base R)
### Copyright 2025, Benjamin Bell.
### Code and updates: https://github.com/benbell95/better_biplot
### Blog: https://www.benjaminbell.co.uk
### Contact details via blog or github. Please send comments, suggestions, bug reports.

### This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

### This function has been tested with R v4.5.0. Use at your own risk.

### Plots a biplot of pca data created with base R prcomp() or princomp() function. 
### Replicates (+ enhances) functionality of biplot() - applying the same data scaling methods.
### All code is base R - does not require any other packages or ggplot2 to run. 
### Should play nicely with layout(), and other plot functions, e.g. adding additional elements after plotting, changing par() etc.

### See README.md file for usage. Arguments listed below. Full guide also available on my blog (coming soon).

# To do: 
# Check Wills email, add the additional functions he requested

### Arguments / details

# x             = The pca (prcomp/princomp) object. You must run prcomp() or princomp() on your data before using this function.
# pc1           = The first principal component to plot [default = 1].
# pc2           = The second principal component to plot [default = 2].

# scale         = Scale data between 0 and 1 for plotting. This works the same way as R base biplot() function.
# pc.biplot     = Logical. Additional scaling of vectors. This works the same way as R base biplot() function.

# limx          = xlim and ylim are generated automatically, but you can increase or decrease using a multiplier value (e.g. 1.5)
# grid          = Logical. Draw a grid through center x and y axes.
# col           = Optional. Colours for the plot / groups. Automatically generates a palette if not supplied. 
# pch           = Optional. Plotting symbols to use for the observations (also see group 2)
# cex.pt        = Changes the size of the plotting symbols.
# xlab          = Overwrite the default x axis label.
# ylab          = Overwrite the default y axis label.  
# axes          = Logical. Plot axes (his controls all axes including vectors) [default = TRUE].

# group         = Optional. Groups for the data. Should be a factor, with length that matches original data. This affects the colour of the points and also allows for ellipses or convex hulls to be drawn. This will also affect the plot symbol used if multiple pch supplied, unless group2 is also specified.
# group2        = Optional. Second grouping of data. Should be a factor, with length that matches original data. This affects the pch symbol only (e.g. multiple symbols within a single grouping). group must also be specified.
# group3        = Optional. Third grouping of data. Should be a factor, with length that matches original data. This affects the size of the pch symbol. group must also be specified (but does not need group2).

# labd          = Optional. A vector of labels for the data observations, this should be the same length as the data. Use NULL values to omit some labels.
# colld         = Colour for the labels [default = black].
# cex.labd      = Size of text labels.
# font.labd     = Font type for data labels [1 = normal, 2 = bold, 3 = italic, 4 = bold italic].

# vectors       = Logical. Whether to plot the vectors or not.
# whichv        = By default, all vectors are plotted, but you can instead specify which to plot. This can be a character vector, where names should exactly match the names of the vectors you want to plot, or a numeric vector to match index position.
# expand        = Scale the vectors, if the vectors are too large or too small, change this value. Multiplier: values above 1 increase, while below 1 decrease the size.
# arrow.len     = Length of arrow head. Use 0 to suppress.
# lwdv          = Line width for vector arrows [default = 2].
# colv          = Colour of the vector arrows.
# colvt         = Colour of the vector arrow labels.
# cex.labv      = Size of the text labels for the vectors.
# font.labv     = Font type for vector labels [1 = normal, 2 = bold, 3 = italic, 4 = bold italic].
# axes.v       = Logical. Plot vector axes [default = TRUE].

# lab_rotation  = Logical. Rotate labels for the data observations and vectors. Useful when you have lots of data points and labels. This may be slow if there are lots of labels to plot.
# ran_adj     = Logical. Alter (adjustment) position of the labels for the data observations.
# valign        = Alignment of labels for the vectors [accepted values = 0, 0.5, 1]

# chull         = Logical. Add convex hulls to the plot - group must be supplied for this to work.
# ellipse       = Logical. Draw ellipses around the data observations based on the group. 
# angle         = Method to determine ellipse angle, either "lm" or "atan".
# autolim       = Logical. Try to increase plot limits with large ellipses. Can also use limx to specify value.
# nofill        = Logical. By default, ellipses and convex hulls have a background colour. Change to false to only show outline.
# lwde          = Line width for the ellipses and/or convex hulls [default = 2].
# se            = Scale the size of the ellipses. Multiplier: values above 1 increase, while below 1 decrease the size.

# legend        = Logical. Add a legend to plot. Only works when group is specified.
# title.leg     = Titles for legend groups, should be a string with length of 1, 2 or 3 (depending on number of groups specified). Default is "Group", "Group 2" and "Group 3".
# cex.leg       = Size of legend (multiplier).
# horiz         = Logical. Plot legends next to each other. Applies when multiple groups uses. Need to specify lpx value to show properly (i.e. move legend position). Can also specify lpy for further control.
# lpx           = Optional. Absolute position of legend on x axis.
# lpy           = Optional. Absolute position of legend on y axis.

# ...           = Additional arguments passed to plot() [main plot only].

########################################
### Ben's Better biplot function
bb_biplot <- function(x, pc1=1, pc2=2, scale=1, pc.biplot=FALSE, limx, grid=TRUE, col, pch=21, cex.pt=1, xlab, ylab, axes=TRUE, group, group2, group3, labd, colld="black", cex.labd=0.5, font.labd=1, vectors=TRUE, whichv, expand=1, arrow.len=0.15, lwdv=2, colv="red", colvt="black", cex.labv=1, font.labv=1, axes.v=TRUE, lab_rotation=FALSE, ran_adj=FALSE, valign=0.5, chull=FALSE, ellipse=FALSE, angle="lm", autolim=TRUE, nofill=FALSE, lwde=2, se=0.7, legend=FALSE, title.leg, cex.leg=0.75, horiz=FALSE, lpx, lpy, ...) {
    ########################################
    ### Scale data
    # Code modified from R base biplot() function code [Copyright (C) 1995-2012 The R Core Team]
    xclass <- class(x)
    # prcomp() scaling
    if(xclass=="prcomp") {
        n <- nrow(x$x)
        if(scale != 0) {lam <- (x$sdev[c(pc1, pc2)] * sqrt(n)) ^ scale} else {lam <- 1 }
        if(pc.biplot==TRUE) lam <- lam / sqrt(n)      
        sco <- t(t(x$x[, c(pc1, pc2)]) / lam)
        rot <- t(t(x$rotation[, c(pc1, pc2)]) * lam)
    }
    # princomp() scaling
    if(xclass=="princomp") {
        n <- x$n.obs
        if(scale != 0) {lam <- (x$sdev[c(pc1, pc2)] * sqrt(n)) ^ scale} else {lam <- 1 }
        if(pc.biplot==TRUE) lam <- lam / sqrt(n)      
        sco <- t(t(x$scores[, c(pc1, pc2)]) / lam)
        rot <- t(t(x$loadings[, c(pc1, pc2)]) * lam)
    }
    # Proportion of variance
    pvar <- x$sdev ^ 2 / sum(x$sdev ^ 2)
    ########################################
    ### Plot data / config
    # Get x/ylim values
    # Scores
    al <- ceiling(max(abs(range(sco)))*10) / 10
    lim <- c(al*-1, al)
    # Rotation
    alr <- ceiling(max(abs(range(rot)))*10) / 10
    limr <- c(alr*-1, alr)
    # Increase x/y lim values?
    if(hasArg(limx) == TRUE) {
        lim <- lim * limx
        limr <- limr * limx
    }
    # Check group(s) are factors
    if(hasArg(group)==TRUE && class(group)!="factor") {group <- as.factor(group)}
    if(hasArg(group2)==TRUE && class(group2)!="factor") {group2 <- as.factor(group2)}
    if(hasArg(group3)==TRUE && class(group3)!="factor") {group3 <- as.factor(group3)}
    # Convex hulls
    if(chull==TRUE) {
        # Needs group to work (and must be a factor)
        if(!hasArg(group)) stop("Error! group must be supplied to draw convex hulls")
        # Add group (factor) data to scores and convert to data.frame
        sco <- cbind(sco, group) |> as.data.frame()
        # Get the outer points (for plotting)
        opin <- lapply(1:nlevels(group), \(x) chull(sco[sco$group == x,]))
    }
    # Ellipsis
    if(ellipse==TRUE) {
        # Needs group to work (and must be a factor)
        if(!hasArg(group)) stop("Error! group must be supplied to draw ellipses")
        # Check angle/se values
        if(angle %in% c("lm","atan")==FALSE) {stop("Error! angle should be specified as either \"lm\" or \"atan\"")}
        if(angle=="atan" && scale==1 && se < 1) {message("Info: You may need to increase se value (e.g. se=1) to increase ellipse size")}
        # Add group (factor) data to scores and convert to data.frame
        sco <- cbind(sco, group) |> as.data.frame()
        # Call function to calculate ellipses
        e <- lapply(1:nlevels(group), \(x) .bb_ellipse(x=sco[sco$group == x,], se=se, angle=angle))
        # Check if ellipses exceed plot limits + auto increase
        if(autolim== TRUE) {
            ec <- do.call(rbind, e)
            ecm <- max(abs(ec))
            if(ecm > max(lim)) {
                # Increase limits
                pc <- (ecm - max(lim)) / ((ecm + max(lim)) / 2) + 1
                pc <- ceiling(pc*10) / 10
                lim <- lim * pc
                limr <- limr * pc        
            }
        }
    }
    # Label rotation
    if(lab_rotation==TRUE) {
        lsrt <- .bb_lab_rotation(x=sco, valign=valign, ran_adj=ran_adj)
        lvsrt <- .bb_lab_rotation(x=rot, valign=valign, ran_adj=FALSE)
    }
    # pch multiplier (default value if no group 3)
    p3 <- 1
    ### Colour/style points by group (factor)
    # !! Overly complex due to outline and fill colours and pch combinations and size of group/group2
    # Check if group exists
    if(hasArg(group) == TRUE) {
        # Group counter (for legend)
        grc <- 1
        # Check length of group is correct
        if(length(group)!=n) stop("Error! group must be same length as original data used in PCA")
        # Check number of pch symbols
        if(length(pch) < nlevels(group)) {pch <- rep(pch, times=nlevels(group))} 
        # Check if colors supplied
        if(hasArg(col)==TRUE) {
            # Check length of colors 
            if(length(col) < nlevels(group)) {
                message("Warning: Number of colours does not match number of groups")
                col <- rep(col, times=nlevels(group)/length(col))
            }
        } else {
            # Generate default colour palette
            col <- hcl.colors(nlevels(group), "Dark 3")
        }
        # colour/pch based on group
        c <- col
        col <- col[group]
        pch1 <- pch[group]  
        # If group 2 is available
        if(hasArg(group2)) {
            # Group counter (for legend)
            grc <- grc + 1
            if(length(pch) < nlevels(group2)) {
                # Create symbols based on number of groups (starting at 21)
                pch <- seq(21, length.out=nlevels(group2))
                # Check if any values above 25 and minus 26 (pch starts at 0)
                pch <- ifelse(pch > 25, pch - 26, pch)
                # In the unlikely event that there are more than 25 symbols...
                if(length(pch)>25) {pch <- ifelse(pch > 25, pch + 6, pch)} 
            }
            pch1 <- pch[group2]
        }
        # If group 3 is available
        if(hasArg(group3)) {
            # Group counter (for legend)
            grc <- grc + 1
            # Create a multiplier based on number of groups.
            g3i <- 2 / nlevels(group3)
            pchm <- seq(from=1, by=g3i, length.out=nlevels(group3))
            p3 <- pchm[group3] 
        }
    # No group supplied                          
    } else {
        # No group, but group 2 supplied:
        if(hasArg(group2)) {message("Warning: group argument required when using group2.")}
        # Give warning message if multiple pch supplied with no group
        if(length(pch)>1) {message("Warning: Multiple pch symbols supplied, but no groups.")}  
        pch1 <- pch
        # No group but legend = TRUE
        if(legend==TRUE) {
            {message("Info: legend only works when group specified.")}
            legend <- FALSE
        }
        # No group supplied, but colour supplied
        if(hasArg(col)==TRUE) {
            # Give warning message if multiple colours supplied with no group
            if(length(col)>1) {message("Warning: Multiple colours supplied, but no groups.")}  
        } else {
        # No group or colour supplied
        col <- "black"
        }   
    }
    # Copy colours for fill symbols (21:25)
    colf <- col
    # Check if any symbols between 21-25, and change outline colour
    pchw <- which(pch1 %in% 21:25)
    if(length(pchw)>0) {
        if(length(col) < length(pch1)) {col <- rep(col, times=length(pch1)/length(col))}
        col[pchw] <- "black"
    } 
    # Set alpha level
    if(nofill==FALSE) {alp <- c(0.25, 0.9)} else {alp <- c(0, 1)}
    # Axis labels
    if(hasArg(xlab) == FALSE) {xlab <- paste("PC", pc1, " (", round(pvar[pc1]*100, 1), "%)", sep = "")}
    if(hasArg(ylab) == FALSE) {ylab <- paste("PC", pc2, " (", round(pvar[pc2]*100, 1), "%)", sep = "")}
    ########################################
    ### Main Plot
    # par
    if(legend==TRUE && !hasArg(lpx) || legend==TRUE && hasArg(lpx) && lpx > max(lim)) {
            op <- par(pty="s", mar=c(5,4,4,7))
        } else {
        op <- par(pty="s")
    }
    on.exit(par(op))
    # Initial blank plot (plot order matters!)
    plot(sco[,1], sco[,2], type="n", asp=1, xlim=lim, ylim=lim, xlab=xlab, ylab=ylab, axes=FALSE, ...)
    if(axes==TRUE) {
        axis(1, lwd=0, lwd.ticks=1)
        axis(2, lwd=0, lwd.ticks=1)
        box()
    }
    # Plot convex hulls
    if(chull==TRUE) {
        for(i in 1:nlevels(group)) {
            polygon(sco[sco$group == i,][opin[[i]],], col=adjustcolor(c[i], alp[1]), border=adjustcolor(c[i], alp[2]), lwd=lwde)       
        }
    }
    # Plot ellipses
    if(ellipse==TRUE) {
        for(i in 1:nlevels(group)) {
            polygon(x=e[[i]]$x, y=e[[i]]$y, col=adjustcolor(c[i], alp[1]), border=adjustcolor(c[i], alp[2]), lwd=lwde)
        }
    }
    # Plot grid
    if(grid==TRUE) {abline(v=0, h=0, lty=2, lwd=1, col=adjustcolor("black", 0.75))}
    # Add data points (pca scores)
    points(sco[,1], sco[,2], pch=pch1, col=adjustcolor(col, 0.85), bg=adjustcolor(colf, 0.85), cex=cex.pt*p3)
    # Add labels
    if(hasArg(labd)) {
        # Check length matches data
        if(length(labd) != length(sco[,1])) message("Info: Number of labels does not match number of data observations")
        if(lab_rotation==TRUE) {
            # srt rotation is not vectorised, so must run in a loop (Can be slow if lots of labels)
            for(i in 1:length(sco[,1])) {
                text(sco[,1][i], sco[,2][i], labels=labd[i], font=font.labd, col=colld, srt=lsrt[,1][i], adj=c(lsrt[,2][i], lsrt[,2][i]), cex=cex.labd, xpd=TRUE) 
            }
        } else {
            # Random positions
            adjld <- sample.int(4L, length(sco[,1]), replace=TRUE) 
            text(sco[,1], sco[,2], labels=labd, font=font.labd, col=colld, pos=adjld, cex=cex.labd, xpd=TRUE)
        }
    }
    ### Legend
    if(legend==TRUE) {
        if(!hasArg(title.leg)) {title.leg <- c("Group", "Group 2", "Group 3") }
        if(!hasArg(lpx)) {lpx <- max(lim)*1.37}
        if(!hasArg(lpy)) {lpy <- ifelse(grc > 1, max(lim), 0)}              
        # Group 1
        lp <- legend(x=lpx, y=lpy, yjust=1, xjust=0, legend=paste(levels(group)), bty="n", col="black", pt.bg=c, pch=22, pt.cex=cex.leg*2, xpd=NA, cex=cex.leg, text.width=1, title=title.leg[1], title.font=2, title.adj=0)
        # Position
        lp1 <- abs(lp$text$y[1] - lp$text$y[2])  
        # Group 2
        if(hasArg(group2)) {
            # Adjust positions
            if(horiz==TRUE) {
                lpx <- lp$rect$left + (lp$rect$w * 1.5)
                lp1 <- 0
            } else {
                lpy <- lp$text$y[length(lp$text$y)]
            }
            # Add legend
            lp <- legend(x=lpx, y=lpy - (lp1 * 1.25), yjust=1, xjust=0, legend=paste(levels(group2)), bty="n", col="black", pt.bg=adjustcolor("black", 0.25), pch=pch, pt.cex=cex.pt, xpd=NA, cex=cex.leg, text.width=1, title=title.leg[2], title.font=2, title.adj=0)
        }
        # Group 3
        if(hasArg(group3)) {
            # Adjust positions
            if(horiz==TRUE) {
                lpx <- lp$rect$left + (lp$rect$w * 1.5)
                lp1 <- 0
            } else {
                lpy <- lp$text$y[length(lp$text$y)]
            }
            # Add legend
            lp <- legend(x=lpx, y=lpy - (lp1 * 1.25), yjust=1, xjust=0, legend=paste(levels(group3)), bty="n", col="black", pt.bg=adjustcolor("black", 0.25), pch=21, pt.cex=cex.pt*pchm, xpd=NA, cex=cex.leg, text.width=1, title=title.leg[length(title.leg)], title.font=2, title.adj=0)
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
            if(length(rot)==2) {
                rot <- data.frame(PC1=rot[1], PC2=rot[2])
                rownames(rot) <- vname[match]
            }
        } 
        # New plot window
        par(new=TRUE)
        plot(rot[,1], rot[,2], type="n", asp=1, xlim=limr, ylim=limr, ann=FALSE, axes=FALSE)
        if(axes.v==TRUE && axes!=FALSE) {
            axis(3, lwd=0, lwd.ticks=1)
            axis(4, lwd=0, lwd.ticks=1)
        }
        # Draw arrows
        arrows(x0=0, x1=rot[,1]*expand, y0=0, y1=rot[,2]*expand, col=colv, length=arrow.len, lwd=lwdv, xpd=TRUE)
        # Labels
        if(lab_rotation==TRUE) {
            # srt rotation is not vectorised, so must run in a loop (Can be slow if lots of labels)
            for(i in 1:length(rot[,1])) {
                text(rot[,1][i]*expand, rot[,2][i]*expand, labels=rownames(rot)[i], font=font.labv, col=colvt, srt=lvsrt[,1][i], adj=c(lvsrt[,2][i], lvsrt[,3][i]), cex=cex.labv, xpd=NA) 
            }
        } else {
             text(rot[,1]*expand, rot[,2]*expand, labels=rownames(rot), font=font.labv, col=colvt, cex=cex.labv, xpd=TRUE)
        }
    }
}

### Label rotation function
# x = matrix of xy values for vector arrows (x$rotation)
.bb_lab_rotation <- function(x, valign, ran_adj) {
    # Angle in degrees
    x_an <- atan(x[,2] / x[,1]) * 180 / pi
    # text adj value
    ad1 <- rep(0, times=nrow(x))
    ad2 <- rep(valign, times=nrow(x))
    # Adjust based on quadrant
    bl <- which(x[,1] < 0 & x[,2] < 0)
    tl <- which(x[,1] < 0 & x[,2] > 0)
    # bottom left
    if(length(bl)>0) {ad1[bl] <- 1}
    # top left
    if(length(tl)>0) {ad1[tl] <- 1}
    # Random adjustment
    if(ran_adj==TRUE) {
        ad1 <- sample(c(0, 0.5, 1), length(ad1), replace=TRUE) 
        ad2 <- sample(c(0, 0.5, 1), length(ad2), replace=TRUE) 
    }
    # Combine and return data
    return(cbind(x_an, ad1, ad2))
}

### Ellipse function
# x should be 2 column df with x/y
# se = scale factor for ellipse
.bb_ellipse <- function(x, angle, se) {
    # Parametric representation
    # Coordinates of min/max on each axis 
    x1 <- x[which.max(x[,1]),]
    x2 <- x[which.min(x[,1]),]
    y1 <- x[which.max(x[,2]),]
    y2 <- x[which.min(x[,2]),]
    # x axis length
    # (x1[,1] = xa, x1[,2] = ya, x2[,1] = xb, x2[,2] = yb)
    # sqrt((xa - xb)^2 + (ya - yb)^2)
    xl <- sqrt((x1[,1] - x2[,1])^2 + (x1[,2] - x2[,2])^2)
    # y axis length
    yl <- sqrt((y1[,1] - y2[,1])^2 + (y1[,2] - y2[,2])^2)
    # combine
    len <- c(xl, yl)
    # Major/minor axis
    a <- len[which.max(len)] # Major axis (biggest)
    b <- len[which.min(len)] # Minor axis (smallest)
    # Get centres (this is just the mean of the co-ordinates)
    cen <- apply(x[,1:2], 2, mean)
    # Angle 
    if(angle=="lm") {
        # Angle using lm (works better?)
        lr <- lm(x[,2] ~ x[,1])
        an <- atan(coef(lr)[2])        
    } else if(angle=="atan") {
        # atan
        an <- atan2(b, a)    
    }
    # Scale ellipse
    a <- a * se
    b <- b * se
    # Plot points
    t <- seq(0, 2 * pi, 0.1) 
    xp <- cen[1] + a * cos(t) * cos(an) - b * sin(t) * sin(an)
    yp <- cen[2] + a * cos(t) * sin(an) + b * sin(t) * cos(an)
    # return data for plotting
    return(data.frame(x=xp, y=yp))
}