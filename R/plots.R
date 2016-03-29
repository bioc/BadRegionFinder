#plot summary variant
plotSummary <- function(threshold1, threshold2, percentage1, percentage2,
                        badCoverageSummary, output){
    coverage_temp <- as.data.frame(badCoverageSummary)[,c(1:3,6:length(as.data.frame(badCoverageSummary)[1,]))]
    xdim <- sum(coverage_temp[,3]-coverage_temp[,2])
    dif <- coverage_temp[,3]-coverage_temp[,2]+1
    ticks <- c(0)
    ticks1 <- c()
    labels <- c(coverage_temp[1,5])
    pointer <- 2
    if(length(coverage_temp[,1]) > 1){
        for(i in 2:length(coverage_temp[,1])){
            if((is.na(coverage_temp[i,5]) && !is.na(labels[pointer-1]))
               || (!is.na(coverage_temp[i,5]) && is.na(labels[pointer-1]))
               || !is.na(coverage_temp[i,5]) && !is.na(labels[pointer-1])
               && coverage_temp[i,5] != labels[pointer-1]){
                labels[pointer] <- coverage_temp[i,5]
                ticks[pointer] <- sum(dif[1:(i-1)])
                ticks1[pointer-1] <- ticks[pointer-1]+(ticks[pointer]-ticks[pointer-1])/2
                pointer <- pointer+1
            }
        } 
    }
    ticks[length(ticks)+1] <- sum(dif)
    ticks1[length(ticks1)+1] <- ticks[length(ticks)-1]+(ticks[length(ticks)]-ticks[length(ticks)-1])/2
    ticks1[length(ticks1)+1] <- ticks1[length(ticks1)]
    labels[length(labels)+1] <- ""
    if(output != ""){
        png(filename = paste(output, "/CoverageQuality_Summary.png", sep=""),
            width = 2000, height = 1200)
    }
    par(mar = c(5,5,3,1))
    plot(1, type="n", xlim=c(0,xdim), ylim=c(0,6), xaxt='n', xlab="Genes",
         ylab="Coverage Quality", main="Summarized Coverage Quality")
    for(i in 1:length(coverage_temp[,1])){
        if(i == 1){
            if(coverage_temp[i,4] == 1){
                rect(xleft = 1, xright = dif[1], 
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1),border=NA, col = "red")
            }
            if(coverage_temp[i,4] == 3){
                rect(xleft = 1, xright = dif[1], 
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "darkgoldenrod1")
            }
            if(coverage_temp[i,4] == 5){
                rect(xleft = 1, xright = dif[1],
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "green")
            } 
            if(coverage_temp[i,4] == 0){
                rect(xleft = 1, xright = dif[1],
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "black")
            }
            if(coverage_temp[i,4] == 2){
                rect(xleft = 1, xright = dif[1],
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "gray40")
            }
            if(coverage_temp[i,4] == 4){
                rect(xleft = 1, xright = dif[1],
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "gray65")
            } 
        }
        if(i > 1){
            if(coverage_temp[i,4] == 1){
                rect(xleft = sum(dif[1:(i-1)]), xright = sum(dif[1:i]),
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "red")
            }
            if(coverage_temp[i,4] == 3){
                rect(xleft = sum(dif[1:(i-1)]), xright = sum(dif[1:i]),
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "darkgoldenrod1")
            }
            if(coverage_temp[i,4] == 5){
                rect(xleft = sum(dif[1:(i-1)]), xright = sum(dif[1:i]),
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border=NA,
                     col = "green")
            }
            if(coverage_temp[i,4] == 0){
                rect(xleft = sum(dif[1:(i-1)]), xright = sum(dif[1:i]),
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border=NA,
                     col = "black")
            }
            if(coverage_temp[i,4] == 2){
                rect(xleft = sum(dif[1:(i-1)]), xright = sum(dif[1:i]),
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border=NA,
                     col = "gray40")
            }
            if(coverage_temp[i,4] == 4){
                rect(xleft = sum(dif[1:(i-1)]), xright = sum(dif[1:i]),
                     ytop = (coverage_temp[i,4]+0.1),
                     ybottom = (coverage_temp[i,4]-0.1), border = NA,
                     col = "gray65")
            } 
        }   
    }
    for(i in 1:length(ticks)){
        lines(rep(ticks[i],2), c(0,5.15), lty=2)
    }
    axis(1, at = ticks1, labels = labels, tick = FALSE)
    axis(1, at = ticks, labels = FALSE, tick = TRUE)
    legend(0, 6, legend = c(paste("Bad Region on Target Region (<",
                                  percentage1*100, "% with coverage >=",
                                  threshold1, "x)", sep=""),
                            paste("Acceptable Region on Target Region (<",
                                  percentage2*100, "% with coverage >",
                                  threshold2, "x, but >=", percentage1*100,
                                  "% with coverage >=", threshold1, "x)", sep=""),
                            paste("Good Region on Target Region (>=",
                                  percentage2*100, "% with coverage >=",
                                  threshold2,"x)", sep="")), bty="n",
           fill = c("red", "goldenrod1", "green"))
    legend(0.5*xdim,6,legend=c(paste("Bad Region off Target Region (<",
                                     percentage1*100, "% with coverage >=",
                                     threshold1, "x)", sep=""),
                               paste("Acceptable Region off Target Region (<",
                                     percentage2*100, "% with coverage >",
                                     threshold2, "x, but >=", percentage1*100,
                                     "% with coverage >=", threshold1, "x)", 
                                     sep=""),
                               paste("Good Region off Target Region (>=",
                                     percentage2*100, "% with coverage >=",
                                     threshold2, "x)", sep="")), bty="n",
           fill=c("black", "gray45", "gray65"))
    if(output != ""){
        dev.off()   
    }
    return()
}

#plot detailed variant
plotDetailed <- function(threshold1, threshold2, percentage1, percentage2,
                         coverage_indicators_temp, output){
    xdim <- 0
    ydim <- 0
    for(i in 1:length(coverage_indicators_temp)){
        if(!is.na(as.data.frame(coverage_indicators_temp[[i]])[1,1])){
            coverage_temp <- as.data.frame(coverage_indicators_temp[[i]])[,c(1:3,6:length(as.data.frame(coverage_indicators_temp[[i]])[1,]))]
            if(length(coverage_temp) > 0 && !is.na(coverage_temp[1,1])){
                xdim <- xdim+length(coverage_temp[,1])
                if(length(coverage_temp[,1]) > 0
                   && max(apply(coverage_temp[,4:(length(coverage_temp[1,])-4)],1,median)) > ydim){
                    ydim <- max(apply(coverage_temp[,4:(length(coverage_temp[1,])-4)],
                                      1, median))
                }
            }
        }   
    }
    first <- 0
    found <- FALSE
    for(i in 1:length(coverage_indicators_temp)){
        if(found == FALSE 
           && !is.na(as.data.frame(coverage_indicators_temp[[i]])[1,1])){
            first <- i
            found <- TRUE
        }
    }
    ticks <- c(0,0)
    ticks1 <- c()
    coverage_temp1 <- as.data.frame(coverage_indicators_temp[[first]])[,c(1:3,6:length(as.data.frame(coverage_indicators_temp[[first]])[1,]))]
    labels <- c(coverage_temp1[1,length(coverage_temp1[1,])-1])
    pointer <- 2
    for(i in first:length(coverage_indicators_temp)){
        message("Analyzing Chromosome ", i)
        if(!is.na(as.data.frame(coverage_indicators_temp[[i]])[1,1])){
            coverage_temp <- as.data.frame(coverage_indicators_temp[[i]])[,c(1:3,6:length(as.data.frame(coverage_indicators_temp[[i]])[1,]))]
            for(j in 1:length(coverage_temp[,1])){
                if(!is.na(coverage_temp[j,length(coverage_temp[1,])-1])
                   && !is.na(labels[pointer-1])
                   && coverage_temp[j,length(coverage_temp[1,])-1] == labels[pointer-1]){
                    ticks[pointer] <- ticks[pointer]+1
                    ticks1[pointer-1] <- ticks[pointer-1]+(ticks[pointer]-ticks[pointer-1])/2
                }
                if(!is.na(coverage_temp[j,length(coverage_temp[1,])-1])
                   && !is.na(labels[pointer-1])
                   && coverage_temp[j,length(coverage_temp[1,])-1] != labels[pointer-1]){
                    pointer <- pointer+1
                    labels[pointer-1] <- coverage_temp[j,length(coverage_temp[1,])-1]
                    ticks[pointer] <- ticks[pointer-1]+1
                    ticks1[pointer-1] <- ticks[pointer-1]+(ticks[pointer]-ticks[pointer-1])/2
                    
                }
                if(is.na(coverage_temp[j,length(coverage_temp[1,])-1])
                   && (is.na(labels[pointer-1])||labels[pointer-1] == "")){
                    ticks[pointer] <- ticks[pointer]+1
                    ticks1[pointer-1] <- ticks[pointer-1]+(ticks[pointer]-ticks[pointer-1])/2
                }
                if(is.na(coverage_temp[j,length(coverage_temp[1,])-1])
                   && !is.na(labels[pointer-1]) && labels[pointer-1]!=""){
                    pointer <- pointer+1
                    labels[pointer-1] <- ""
                    ticks[pointer] <- ticks[pointer-1]+1
                    ticks1[pointer-1] <- ticks[pointer-1]+(ticks[pointer]-ticks[pointer-1])/2 
                }
            } 
        }
    }
    ticks[length(ticks)+1] <- xdim
    ticks1[length(ticks1)+1] <- ticks[length(ticks)-1]+(ticks[length(ticks)]-ticks[length(ticks)-1])/2
    ticks1[length(ticks1)+1] <- ticks1[length(ticks1)]
    labels[length(labels)+1] <- ""
    labels[length(labels)+1] <- ""
    if(output != ""){
        png(filename = paste(output, "/CoverageQuality_Details.png", sep=""),
            width = 2000, height = 1200)   
    }
    par(mar = c(5,5,3,1))
    plot(1, type="n", xlim=c(0,xdim), ylim=c(0,ydim), xaxt='n', xlab="Genes",
         ylab="Median Coverage", main="Detailed Coverage Quality")
    alreadyVisited <- 0
    last <- 0
    for(i in first:length(coverage_indicators_temp)){
        if(!is.na(as.data.frame(coverage_indicators_temp[[i]])[1,1])){
            coverage_temp <- as.data.frame(coverage_indicators_temp[[i]])[,c(1:3,6:length(as.data.frame(coverage_indicators_temp[[i]])[1,]))]
            for(j in 1:length(coverage_temp[,1])){
                if(j == 1 && i == first){
                    last <- apply(coverage_temp[j,4:(length(coverage_temp[1,])-4)],
                                  1, median)
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 1){
                        rect(xleft = alreadyVisited, 
                             xright = (alreadyVisited+1),
                             ytop = last, ybottom = 0, col = "red", border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 3){
                        rect(xleft = alreadyVisited, 
                             xright = (alreadyVisited+1),
                             ytop = last, ybottom = 0, col = "darkgoldenrod1",
                             border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 5){
                        rect(xleft = alreadyVisited, 
                             xright = (alreadyVisited+1),
                             ytop = last, ybottom = 0, col = "green",
                             border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 0){
                        rect(xleft = alreadyVisited, 
                             xright = (alreadyVisited+1),
                             ytop = last, ybottom = 0, col = "black",
                             border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 2){
                        rect(xleft = alreadyVisited, 
                             xright = (alreadyVisited+1),
                             ytop = last, ybottom = 0, col = "gray45",
                             border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 4){
                        rect(xleft = alreadyVisited,
                             xright = (alreadyVisited+1),
                             ytop = last, ybottom = 0, col = "gray65", 
                             border = NA)
                    }
                    alreadyVisited <- alreadyVisited+1
                }
                if(j > 1 || i > 1){
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 1){
                        rect(xleft = alreadyVisited, 
                             xright = (alreadyVisited+1),
                             ytop = apply(coverage_temp[j,4:(length(coverage_temp[1,])-4)],
                                          1, median),
                             ybottom = 0, col = "red", border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 3){
                        rect(xleft = alreadyVisited, xright = (alreadyVisited+1),
                             ytop = apply(coverage_temp[j,4:(length(coverage_temp[1,])-4)],
                                          1, median),
                             ybottom = 0, col = "darkgoldenrod1", border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 5){
                        rect(xleft = alreadyVisited, 
                             xright = (alreadyVisited+1),
                             ytop = apply(coverage_temp[j,4:(length(coverage_temp[1,])-4)],
                                          1, median),
                             ybottom = 0, col = "green", border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 0){
                        rect(xleft = alreadyVisited,
                             xright = (alreadyVisited+1),
                             ytop = apply(coverage_temp[j,4:(length(coverage_temp[1,])-4)],
                                          1, median),
                             ybottom = 0, col = "black", border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 2){
                        rect(xleft = alreadyVisited,
                             xright = (alreadyVisited+1),
                             ytop = apply(coverage_temp[j,4:(length(coverage_temp[1,])-4)],
                                          1,
                                          median),
                             ybottom = 0, col = "gray45", border = NA)
                    }
                    if(coverage_temp[j,length(coverage_temp[1,])-2] == 4){
                        rect(xleft = alreadyVisited,xright = (alreadyVisited+1),
                             ytop = apply(coverage_temp[j,4:(length(coverage_temp[1,])-4)],
                                          1, median),
                             ybottom = 0, col = "gray65", border = NA)
                    }
                    alreadyVisited <- alreadyVisited+1   
                }
            }
        }   
    }
    for(i in 1:length(ticks)){
        lines(rep(ticks[i],2), c(0,ydim*0.9), lty=2)
    }
    axis(1, at = ticks1, labels = labels, tick = FALSE)
    axis(1, at = ticks, labels = FALSE, tick = TRUE)
    legend(0, ydim, legend = c(paste("Bad Region on Target Region (<",
                                     percentage1*100, "% with coverage >",
                                     threshold1, "x)", sep=""),
                               paste("Acceptable Region on Target Region (<",
                                     percentage2*100, "% with coverage >",
                                     threshold2, "x, but >=", percentage1*100,
                                     "% with coverage >=", threshold1, "x)", 
                                     sep=""),
                               paste("Good Region on Target Region (>=",
                                     percentage2*100, "% with coverage >=",
                                     threshold2, "x)", sep="")),
           bty = "n", fill = c("red", "goldenrod1", "green"))
    legend(0.5*xdim, ydim, legend = c(paste("Bad Region off Target Region (<",
                                            percentage1*100, "% with coverage >",
                                            threshold1, "x)", sep=""),
                                      paste("Acceptable Region off Target Region (<",
                                            percentage2*100, "% with coverage >",
                                            threshold2, "x, but >=", percentage1*100,
                                            "% with coverage >=", threshold1, "x)",
                                            sep=""),
                                      paste("Good Region off Target Region (>=",
                                            percentage2*100, "% with coverage >=",
                                            threshold2, "x)", sep="")),
           bty = "n", fill = c("black", "gray45", "gray65"))
    if(output != ""){
        dev.off()   
    }
}
#plot badCoverageOverview as barplot:
plotSummaryGenes <- function(threshold1, threshold2, percentage1, percentage2,
                             badCoverageGenes, output){
    onOffTarget <- matrix(rep(0,3*2*length(badCoverageGenes[,1])), ncol=3)
    for(i in 1:length(badCoverageGenes[,1])){
        onOffTarget[i*2-1,] <- c(badCoverageGenes[i,3], 
                                 badCoverageGenes[i,5], badCoverageGenes[i,7])
        onOffTarget[i*2,] <- c(badCoverageGenes[i,4],
                               badCoverageGenes[i,6], badCoverageGenes[i,8])
    }
    if(output != ""){
        png(filename = paste(output, "/CoverageQuality_Summary_Genes.png", 
                             sep=""), width = 2000, height = 1200)
    }
    barplot(height = t(onOffTarget*100), width = 2,
            xlim = c(0,length(onOffTarget[,1])*3),
            ylim = c(0,110), col=c("red","goldenrod1","green"),
            space = c(1,0), main = "Summarized Coverage Quality",
            xlab = "Genes", ylab = "Coverage Quality (%)")
    barplot(height = t(onOffTarget*100)[,seq(1,length(onOffTarget[,1]),2)],
            width = 2, xlim = c(0,length(onOffTarget[,1])*3), ylim = c(0,110),
            col = c("black", "gray45", "gray65"),
            space = c(1,rep(2,length(onOffTarget[,1])-1)), add = TRUE)
    axis(1, at = sort(c(seq(2, length(onOffTarget[,1])*3, 6),
                        seq(6, length(onOffTarget[,1])*3, 6))), labels = FALSE)
    axis(1, at = seq(4, length(onOffTarget[,1])*3, 6),
         labels = badCoverageGenes[,1], tick = FALSE)
    legend(0, 110, legend = c(paste("Bad Region on Target Region (<",
                                    percentage1*100, "% with coverage >",
                                    threshold1, "x)", sep = ""),
                              paste("Acceptable Region on Target Region (<",
                                    percentage2*100, "% with coverage >", threshold2,
                                    "x, but >=", percentage1*100, "% with coverage >=",
                                    threshold1, "x)", sep = ""),
                              paste("Good Region on Target Region (>=",
                                    percentage2*100, "% with coverage >=",
                                    threshold2, "x)", sep="")), bty = "n",
           fill = c("red", "goldenrod1", "green"))
    legend(length(onOffTarget[,1])*1.5,110,
           legend = c(paste("Bad Region off Target Region (<", percentage1*100,
                            "% with coverage >", threshold1, "x)", sep = ""),
                      paste("Acceptable Region off Target Region (<",
                            percentage2*100, "% with coverage >", threshold2,
                            "x, but >=", percentage1*100, "% with coverage >=",
                            threshold1, "x)", sep=""),
                      paste("Good Region off Target Region (>=",
                            percentage2*100, "% with coverage >=", threshold2,
                            "x)", sep = "")),
           bty = "n", fill = c("black", "gray45", "gray65"))
    if(output != ""){
        dev.off()   
    }
}