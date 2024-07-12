# run_scTorus
#
# This function generate plotting objects for making cell cycle trajectory plots.
# It will output a scatterplot showing single cells arranged in a circle and a heatmap
# and trajectory profiles showing the up- and downs of various cell cycle protine markers.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

run_scTorus <- function(x, genes.to.highlight = NULL,
                        correlation.plot = FALSE, method = c("lmds", "isoMDS", "CMD"),
                        control.ROIs = NULL, treated.ROIs = NULL, smoothFactor = 5){

    # Transform data to make it closer to a normal distribution
    x.transformed <- log2(x + 1)

    # Remove cells exceed 4 standard deviations. These are considered as outliers
    tmp <- x.transformed
    for (marker in colnames(x.transformed)){
        intensity <- x.transformed[[marker]]
        sd_ <- sd(intensity)
        tmp <- tmp %>% filter(tmp[[marker]] < 4*sd_)
    }

    x.transformed <- tmp
    x.scaled <- scale(x.transformed)
    x.scaled[is.nan(x.scaled)] <- 0
    x.scaled <- as.data.frame(x.scaled)

    plotList <- list()

    if (method == 'lmds'){
        results <- lmds(
            as.matrix(x.scaled),
            distance_method = "pearson",
            num_landmarks = 500
        )

        results <- as.data.frame(results)
        colnames(results) <- c("CMD1", "CMD2", "CMD3")
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
            colnames(results) <- c("CMD1", "CMD2", "Stress")
        } else if (method == 'CMD'){
            scaledCmd <- cmdscale(data.for.ccd, k = 2)
            results <- as.data.frame(results)
            colnames(results) <- c("CMD1", "CMD2")
        }
    }

    plotList2 <- list()
    # add a coreID column for deconvoluting cells from control and treated
    results$ROI <- sub('.im3:\\d+', '', rownames(results))

    results <- merge(results, x, by = 0)
    rownames(results) <- results$Row.names
    results$Row.names <- NULL

    if (!is.null(control.ROIs) & !is.null(treated.ROIs)){

        cells.ctrl <- results %>% filter(ROI %in% control.ROIs)
        cells.treated <- results %>% filter(ROI %in% treated.ROIs)

        dsList <- list()
        dsList[['Ctrl']] <- cells.ctrl
        dsList[['Treated']] <- cells.treated



        for (ds in names(dsList)){
            data.for.circle <- dsList[[ds]]

            for (marker in colnames(x)){
                data.for.circle$Intensity <- data.for.circle[[marker]]

                p <- ggplot(data.for.circle, aes(x = CMD1, y = CMD2, colour = Intensity)) + geom_point(alpha = 1, size = 0.5) +
                    labs(x = "Comp_1", y = "Comp_2")
                p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "white"))
                p <- p + scale_colour_viridis(option = "H", name = "Intensity")
                p <- p +  ggtitle(paste(ds, marker, sep = ', ')) + theme(plot.title = element_text(hjust = 0.5))

                plotList2[[paste(ds, marker)]] <- p
            }

        }
    } else {
        for (marker in colnames(x)){
            data.for.circle <- results
            data.for.circle$Intensity <- data.for.circle[[marker]]

            p <- ggplot(data.for.circle, aes(x = CMD1, y = CMD2, colour = Intensity)) + geom_point(alpha = 1, size = 0.5) +
                labs(x = "Comp_1", y = "Comp_2")
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_rect(fill = "white"))
            p <- p + scale_colour_viridis(option = "H", name = "Intensity")
            p <- p +  ggtitle(paste(ds, marker, sep = ', ')) + theme(plot.title = element_text(hjust = 0.5))

            plotList2[[marker]] <- p
        }
    }

    plotList[['Circle']] <- plotList2

    #
    # arrange the cells based on cell cycle stages
    #


    # center the data on the origin of the fitted circle
    f = fitSS(results[, 1:2])
    results$x1 <- results$CMD1 - f$par[1]
    results$y1 <- results$CMD2 - f$par[2]

    # Remove the cells near the center of the origin. These cells are considered as 'noise'
    # the radius of the fitted circle is in f$par[3]. We will remove cells within the 1/4 of the radius
    results$keep <- sqrt(results$x1**2 + results$y1**2) > f$par[3] * 1/2

    results <- results[results$keep == TRUE, ]


    results$angle <- atan2(results$y1, results$x1)

    results$angle2 <-  ifelse(results$angle >= 0, results$angle, 3.14*2 + results$angle)

    results$angle <- NULL # no langer needed

    results <- results[order(results$angle2), ]

    # Heatmap
    data.for.heatmap <- x.scaled[match(rownames(results), rownames(x.scaled)), ]


    suppressMessages(
        ph <- Heatmap(t(data.for.heatmap),
                      cluster_rows = F,
                      cluster_columns = F,
                      show_column_names = F,
                      name = 'Zscore of\nIntensity'
        )
    )

    plotList[['heatmap']] <- ph

    # make smoothed trajectory line plots
    exp.cell.cycle.smoothed <- apply(data.for.heatmap, 2, smoothByWindow)
    exp.cell.cycle.smoothed <- as.data.frame(exp.cell.cycle.smoothed)
    for (i in 1:smoothFactor){
        exp.cell.cycle.smoothed <- apply(exp.cell.cycle.smoothed, 2, smoothByWindow)
        exp.cell.cycle.smoothed <- as.data.frame(exp.cell.cycle.smoothed)
    }

    # add a pseudo time column for ordering the cells
    exp.cell.cycle.smoothed$time <- 1:nrow(exp.cell.cycle.smoothed)

    dd = melt(exp.cell.cycle.smoothed, id=c("time"))

    colnames(dd) <- c('Pseudo Time', 'Gene', 'Z-Score of Expression')

    # make lines to have similar colors for each gene
    pt <- ggplot(dd) + geom_line(aes(x=`Pseudo Time`, y=`Z-Score of Expression`, colour=Gene, linetype=Gene)) + labs(title=NULL) +
        theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 6)) + theme_classic()


    plotList[['traject']] <- pt

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

