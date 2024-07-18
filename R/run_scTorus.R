#' run_scTorus
#'
#' This function generate plotting objects for making cell cycle trajectory plots.
#' It will output a scatterplot showing single cells arranged in a circle and a heatmap
#' and trajectory profiles showing the up- and downs of various cell cycle protine markers.
#' @param x A data frame with cell IDs as rows and markers as columns
#' @param genes.to.highlight Gene names to highlight using measured mean intensity values
#' @param method Dimension reduction method to use. Options are "lmds", "isoMDS", "CMD"
#' @param control.ROIs ROIs for cells in the control group
#' @param treated.ROIs ROIs for cells in the treated group
#' @param smoothFactor A integer values used to control the smoothness of the trajectory
#' @return A list of ggplot2 plotting objects, including a scatterplot showing cells arranged as circles,
#' a heatmap showing the marker intensity along a calculated pseudo time, a line plot showing the dynamics
#' of intensity along the pseudo time
#' @export
#'
run_scTorus <- function(x, genes.to.highlight = NULL, method = c("lmds", "isoMDS", "CMD", "phate"),
                        control.ROIs = NULL, treated.ROIs = NULL, smoothFactor = 5, reference.marker = NULL){
    
    # Transform data to make it closer to a normal distribution
    x.transformed <- log2(x + 1)
    x.scaled <- scale(x.transformed)
    
    #x.scaled[is.nan(x.scaled)] <- 0

    plotList <- list()
    
    if (method == 'lmds'){
        results <- lmds(
            x.scaled,
            distance_method = "pearson",
            num_landmarks = 500
        )
    } else if (method %in% c("isoMDS", "CMD")){
        # down sample the data if exceeds 10000 cells to speed up the code running
        if (dim(x.scaled)[1] > 3000){
            to.keep <- sample(1:nrow(x.scaled), 3000, replace = FALSE)
            x.scaled <- x.scaled[to.keep, ]
        }
        
        corMatrix <- cor(t(x.scaled), use="everything", method = "pearson")
        data.for.ccd <- 0.5 - abs(corMatrix)/2
        
        if (method == 'isoMDS'){
            results <- isoMDS(data.for.ccd, k = 2)
            results <- as.data.frame(results)
        } else if (method == 'CMD'){
            scaledCmd <- cmdscale(data.for.ccd, k = 2)
        }
    } else if (method == 'phate'){
        results <- phate(x.scaled, knn=3, decay=100, t=10)
    }
    
    results <- as.data.frame(results)
    colnames(results)[1:2] <- c("CMD1", "CMD2")
    
    
    # Fit a circle for the data and center the data on the origin of the fitted circle
    f = fitSS(results[, 1:2])
    results$x1 <- results$CMD1 - f$par[1]
    results$y1 <- results$CMD2 - f$par[2]
    
    # Remove the cells near the center of the origin. These cells are considered as 'noise'
    # the radius of the fitted circle is in f$par[3]. We will remove cells within the 1/4 of the radius
    #results$keep <- sqrt(results$x1**2 + results$y1**2) > f$par[3] * 1/2
    # results <- results[results$keep == TRUE, ]
    
    # Calculate angle of each data point in radians
    results$angle <- atan2(results$y1, results$x1)
    
    # convert to continous radian angle from 0 to 2pi, counterclockwise
    results$angle2 <-  ifelse(results$angle >= 0, results$angle, 3.14*2 + results$angle)
    
    # Join with scaled intensity data for making trajectory heatmap and pseudo-time profiles
    results <- merge(results, x.scaled, by = 0)
    rownames(results) <- results$Row.names
    results$Row.names <- NULL
    
    # Order by angle2
    results <- results[order(results$angle2), ]
    
    # add a coreID column for deconvoluting cells from control and treated
    results$ROI <- sub('.im3:\\d+', '', rownames(results))
    
    # Use reference marker to calculate a rotation angle so that the reference marker are plotted
    # at the 12 o'clock position. We will do this by generating trajectory profiles for all the 
    # markers and identify the smoothed profile peak location for marker from user input. Then
    # calculate a rotation angle and use this angle to calculate the new x and y coordinates
    # of the cells.
    
    data.for.heatmap <- results[, colnames(x.scaled)]
    
    exp.cell.cycle.smoothed <- apply(data.for.heatmap, 2, smoothByWindow)
    exp.cell.cycle.smoothed <- as.data.frame(exp.cell.cycle.smoothed)
    
    for (i in 1:smoothFactor){
      exp.cell.cycle.smoothed <- apply(exp.cell.cycle.smoothed, 2, smoothByWindow)
      exp.cell.cycle.smoothed <- as.data.frame(exp.cell.cycle.smoothed)
    }
    
    # Get the peak of the smoothed profile for the reference marker
    if (!is.null(reference.marker)){
      
      peak <- which(exp.cell.cycle.smoothed[[reference.marker]] == max(exp.cell.cycle.smoothed[[reference.marker]], na.rm = TRUE))
      peak.angle <- results[peak, ]$angle2
      
      if (peak.angle >= 0 & peak.angle < pi/2){
        amount.to.rotate <- 1.5*pi + peak.angle 
      } else if (peak.angle >= pi/2 & peak.angle < pi){
        amount.to.rotate <- peak.angle - pi/2
      } else if (peak.angle >= pi & peak.angle < 3/2 * pi){
        amount.to.rotate <- peak.angle - pi 
      } else {
        amount.to.rotate <- peak.angle - 1.5*pi 
      }
      
      results$x2 <- results$x1 * cos(amount.to.rotate) + results$y1 * sin(amount.to.rotate)
      results$y2 <- -results$x1 * sin(amount.to.rotate) + results$y1 * cos(amount.to.rotate)
      results$x1 <- results$x2
      results$y1 <- results$y2
      results$x2 <- NULL
      results$y2 <- NULL
      
      # shift the data points to start with the peak position of reference marker
      data.for.heatmap <- rbind(data.for.heatmap[peak:nrow(data.for.heatmap), ], data.for.heatmap[1:peak,])
    }
    
    if (!is.null(control.ROIs) & !is.null(treated.ROIs)){
      
      # Add ROI column for splitting cells by treatment
      data.for.heatmap$ROI <- sub('.im3:\\d+', '', rownames(data.for.heatmap))
      
      # Split into treated and control and plot separately
      data.for.heatmap.treated <- data.for.heatmap
      
      # Replace missing values with zero
      data.for.heatmap.treated[data.for.heatmap.treated$ROI %in% control.ROIs, 1:(ncol(data.for.heatmap.treated) -1)] <- 0
      data.for.heatmap.treated$ROI <- NULL
      
      data.for.heatmap.ctrl <- data.for.heatmap
      
      # Replace missing values with zero
      data.for.heatmap.ctrl[data.for.heatmap.ctrl$ROI %in% treated.ROIs, 1:(ncol(data.for.heatmap.ctrl) - 1) ] <- 0
      data.for.heatmap.ctrl$ROI <- NULL
      data.for.heatmap$ROI <- NULL
      
      # make smoothed trajectory line plots
      exp.cell.cycle.smoothed.treated <- apply(data.for.heatmap.treated, 2, smoothByWindow)
      exp.cell.cycle.smoothed.treated <- as.data.frame(exp.cell.cycle.smoothed.treated)
      for (i in 1:smoothFactor){
        exp.cell.cycle.smoothed.treated <- apply(exp.cell.cycle.smoothed.treated, 2, smoothByWindow)
        exp.cell.cycle.smoothed.treated <- as.data.frame(exp.cell.cycle.smoothed.treated)
      }
      
      exp.cell.cycle.smoothed.ctrl <- apply(data.for.heatmap.ctrl, 2, smoothByWindow)
      exp.cell.cycle.smoothed.ctrl <- as.data.frame(exp.cell.cycle.smoothed.ctrl)
      for (i in 1:smoothFactor){
        exp.cell.cycle.smoothed.ctrl <- apply(exp.cell.cycle.smoothed.ctrl, 2, smoothByWindow)
        exp.cell.cycle.smoothed.ctrl <- as.data.frame(exp.cell.cycle.smoothed.ctrl)
      }
      
      suppressMessages(
        ph.t <- Heatmap(t(data.for.heatmap.treated),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_names = F,
                        name = 'Z-Score of\nIntensity'
        )
      )
      suppressMessages(
        ph.c <- Heatmap(t(data.for.heatmap.ctrl),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_names = F,
                        name = 'Z-Score of\nIntensity'
        )
      )
      plotList[['heatmap.treated']] <- ph.t
      plotList[['heatmap.ctrl']] <- ph.c
      
    } else {
      suppressMessages(
        ph <- Heatmap(t(data.for.heatmap),
                      cluster_rows = F,
                      cluster_columns = F,
                      show_column_names = F,
                      name = 'Z-Score of\nIntensity'
        )
      )
      plotList[['heatmap']] <- ph
    }
    
    # add a pseudo time column for ordering the cells
    if (!is.null(control.ROIs) & !is.null(treated.ROIs)){
      exp.cell.cycle.smoothed.treated$time <- 1:nrow(exp.cell.cycle.smoothed.treated)
      dd = melt(exp.cell.cycle.smoothed.treated, id=c("time"))
      colnames(dd) <- c('Pseudo Time', 'Gene', 'Z-Score of Intensity')
      
      # make lines to have similar colors for each gene
      pt.t <- ggplot(dd) + geom_line(aes(x=`Pseudo Time`, y=`Z-Score of Intensity`, colour=Gene, linetype=Gene)) + labs(title=NULL) +
        theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 6)) + theme_classic()
      plotList[['traject.treated']] <- pt.t
      
      exp.cell.cycle.smoothed.ctrl$time <- 1:nrow(exp.cell.cycle.smoothed.ctrl)
      dd = melt(exp.cell.cycle.smoothed.ctrl, id=c("time"))
      colnames(dd) <- c('Pseudo Time', 'Gene', 'Z-Score of Intensity')
      
      # make lines to have similar colors for each gene
      pt.c <- ggplot(dd) + geom_line(aes(x=`Pseudo Time`, y=`Z-Score of Intensity`, colour=Gene, linetype=Gene)) + labs(title=NULL) +
        theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 6)) + theme_classic()
      plotList[['traject.control']] <- pt.c
      
    } else {
      exp.cell.cycle.smoothed$time <- 1:nrow(exp.cell.cycle.smoothed)
      dd = melt(exp.cell.cycle.smoothed, id=c("time"))
      colnames(dd) <- c('Pseudo Time', 'Gene', 'Z-Score of Intensity')
      
      # make lines to have similar colors for each gene
      pt <- ggplot(dd) + geom_line(aes(x=`Pseudo Time`, y=`Z-Score of Intensity`, colour=Gene, linetype=Gene)) + labs(title=NULL) +
        theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 6)) + theme_classic()
      plotList[['traject']] <- pt
    }
    
    # Join with original intensity data for marker expression level highlighting
    results <- merge(results, x, by = 0)
    rownames(results) <- results$Row.names
    results$Row.names <- NULL
    
    plotList2 <- list()
    
    if (!is.null(control.ROIs) & !is.null(treated.ROIs)){
        
        cells.ctrl <- results %>% filter(ROI %in% control.ROIs)
        cells.treated <- results %>% filter(ROI %in% treated.ROIs)

        # Use the same number of cells in the two conditions
        if (nrow(cells.ctrl) > nrow(cells.treated)){
          cells.ctrl <- cells.ctrl[sample(1:nrow(cells.ctrl), nrow(cells.treated)), ]
        } else {
          cells.treated <- cells.treated[sample(1:nrow(cells.treated), nrow(cells.ctrl)), ]
        }
        
        dsList <- list()
        dsList[['Ctrl']] <- cells.ctrl
        dsList[['Treated']] <- cells.treated
        
        for (ds in names(dsList)){
            data.for.circle <- dsList[[ds]]
            
            for (marker in colnames(x)){
              marker <- paste(marker, 'y', sep = '.')
              data.for.circle$Intensity <- log10(data.for.circle[[marker]] + 1)
              
              p <- ggplot(data.for.circle, aes(x = x1, y = y1, colour = Intensity)) + geom_point(alpha = 1, size = 0.5) +
                labs(x = "Comp_1", y = "Comp_2")
              
              p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_rect(fill = "white"))
              p <- p + scale_colour_viridis(option = "H", name = "Intensity")
                p <- p +  ggtitle(paste(ds, sub('.y', '', marker), sep = ', ')) + theme(plot.title = element_text(hjust = 0.5))
                
                plotList2[[paste(ds, sub('.y', '', marker))]] <- p
            }
        }
    } else {
        for (marker in colnames(x)){
            data.for.circle <- results
            
            marker <- paste(marker, 'y', sep = '.')
            
            data.for.circle$Intensity <- log10(data.for.circle[[marker]] + 1)

            p <- ggplot(data.for.circle, aes(x = x1, y = y1, colour = Intensity)) + geom_point(alpha = 1, size = 0.5) +
              labs(x = "Comp_1", y = "Comp_2")
            
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_rect(fill = "white"))
            p <- p + scale_colour_viridis(option = "H", name = "Intensity")
            p <- p +  ggtitle(sub('.y', '', marker)) + theme(plot.title = element_text(hjust = 0.5))
            
            plotList2[[marker]] <- p
        }
    }
    
    plotList[['Circle']] <- plotList2
    
    # calculate inter octile variance
    results$octile <- 8
    results$octile[results$angle2 < pi * 7/4] <- 7
    results$octile[results$angle2 < pi * 3/2] <- 6
    results$octile[results$angle2 < pi * 5/4] <- 5
    results$octile[results$angle2 < pi] <- 4
    results$octile[results$angle2 < pi *3/4] <- 3
    results$octile[results$angle2 < pi/2] <- 2
    results$octile[results$angle2 < pi/4] <- 1
    
    paramList <- list()
    
    iov <- sd(table(results$octil))/mean(table(results$octil)) * 100
    paramList[['iov']] <- iov
    
    # calculate CFD (circle fit distance)
    results$distance <- abs(sqrt(results$x1^2 + results$y1^2) - f$par[3])
    paramList[['cfd']] <- sum(results$distance)/nrow(results) * 100 # scale up to per 100 cells
    paramList[['count']] <- nrow(results)
    
    res <- list(plot = plotList, param = paramList)
    
    res
}

