# utils
#
# Utility functions for scTorus

smoothByWindow <- function(x){
	    res <- NULL
    
    # pad 50 data points on each end of x
    y <- c(x[(length(x) - 50):length(x)], x, x[1:50])
	    
	    res <- lapply(1:length(x), function(x) mean(y[x:(x+100)]))
	    res <- unlist(res)
		    
		    res
}

smoothByLowess <- function(x){
    res <- lowess(1:length(x), x, f = 0.05)
    res[[2]]
}

fitSS <- function(xy,
                  a0=mean(xy[,1]),
                  b0=mean(xy[,2]),
                  r0 = mean(sqrt((xy[,1]-a0)^2 + (xy[,2]-b0)^2)),
                  ...){
    SS <- function(abr){
        sum((abr[3] - sqrt((xy[,1]-abr[1])^2 + (xy[,2]-abr[2])^2))^2)
    }
    optim(c(a0,b0,r0), SS, ...)
}

circlexy <- function(xyr, n=180){
    theta = seq(0,2*pi,len=n)
    cbind(xyr[1] + xyr[3]*cos(theta),
          xyr[2] + xyr[3]*sin(theta)
    )
}

make.hist <- function(x, y, log2.scale = FALSE){

    histList <- list()

    # log2 transform data
    if (log2.scale){
        data.for.hist <- as.data.frame(lapply((x + 1), log2))
    } else {
        data.for.hist <- x
    }
    for (marker in names(data.for.hist)){
        p <- ggplot(data.for.hist, aes_string(x = marker)) + geom_histogram(color="black", fill="white", binwidth = 0.02)

        histList[[marker]] <- p
    }

    p <- wrap_plots(plots = histList, ncol = y)

    p
}

