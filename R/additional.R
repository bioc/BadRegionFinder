###weitere funktionen:
#user-defined quaniles
determineQuantiles <- function(coverage_summary, quantiles, output){
    coverage_summary2 <- list()
    for(i in 1:length(coverage_summary)){
        if(!is.na(as.data.frame(coverage_summary[[i]])[1,1])){
            coverage_temp <- as.data.frame(coverage_summary[[i]])[,c(1:3,6:length(as.data.frame(coverage_summary[[i]])[1,]))]
            coverage_summary2[[i]] <- coverage_temp[,1:3]
            for(j in 1:length(quantiles)){
                coverage_summary2[[i]][,3+j] <- apply(coverage_temp[,4:(length(coverage_temp[1,])-1)],
                                                      1,quantile,probs=quantiles[j],na.rm=TRUE)
            }
            coverage_summary2[[i]][,length(coverage_summary2[[i]][1,])+1] <- coverage_temp[,length(coverage_temp[1,])]
            names(coverage_summary2[[i]]) <- c("Chr", "Start", "End",
                                               quantiles, "OnTarget")
            if(output != ""){
                write.table(coverage_summary2[[i]],
                            paste(output, "/Quantiles_chr", i, ".txt", sep=""),
                            row.names=FALSE, sep="\t", quote=FALSE)   
            }
            coverage_summary2[[i]] <- makeGRangesFromDataFrame(coverage_summary2[[i]],
                                                               start.field="start",
                                                               end.field="end",
                                                               keep.extra.columns=TRUE)
        }
        if(is.na(as.data.frame(coverage_summary[[i]])[1,1])){
            coverage_summary2[[i]] <- NA
        }
    }
    return(coverage_summary2)
}