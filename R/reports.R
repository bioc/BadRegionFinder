#summary variant
reportBadRegionsSummary <- function(threshold1, threshold2, percentage1,
                                    percentage2, coverage_indicators,
                                    mart, output){
    badCoverage <- list()
    for(i in 1:length(coverage_indicators)){
        badCoverage[[i]] <- data.frame(Chr=NA, start=NA, end=NA,
                                       QualityMarker=NA, Gene=NA, GeneID=NA)
    }
    for(i in 1:length(coverage_indicators)){
        message("Analyzing Chromosome ", i)
        if(!is.na(as.data.frame(coverage_indicators[[i]])[1,1])){
            coverage_temp <- as.data.frame(coverage_indicators[[i]])[,c(1:3,6:length(as.data.frame(coverage_indicators[[i]])[1,]))]
            line_counter <- 1
            line_counter2 <- 1        
            while(line_counter2 <= length(coverage_temp[,1])){
                badCoverage[[i]][line_counter,1] <- as.character(coverage_temp[line_counter2,1])
                badCoverage[[i]][line_counter,2] <- coverage_temp[line_counter2,2]
                badCoverage[[i]][line_counter,3] <- coverage_temp[line_counter2,2]
                badCoverage[[i]][line_counter,4] <- coverage_temp[line_counter2,length(coverage_temp[1,])]
                line_counter2 <- line_counter2+1
                while(line_counter2 <= length(coverage_temp[,1])
                      && badCoverage[[i]][line_counter,4] == coverage_temp[line_counter2,length(coverage_temp[1,])]
                      && (coverage_temp[line_counter2-1,2]+1) == coverage_temp[line_counter2,2]){
                    badCoverage[[i]][line_counter,3] <- badCoverage[[i]][line_counter,3]+1
                    line_counter2 <- line_counter2+1
                }
                if(line_counter2 <= length(coverage_temp[,1])
                   && (badCoverage[[i]][line_counter,4] != coverage_temp[line_counter2,length(coverage_temp[1,])]
                       || (coverage_temp[line_counter2-1,2]+1) != coverage_temp[line_counter2,2])){
                    if(coverage_temp[line_counter2-1,2] != coverage_temp[line_counter2-1,3]){
                        badCoverage[[i]][line_counter,3] <- coverage_temp[line_counter2-1,3]
                    }
                    line_counter <- line_counter+1
                }
            }           
            if(!is.na(badCoverage[[i]][1,1])){
                if(is.character(mart)){
                    mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                   host="grch37.ensembl.org",
                                   path="/biomart/martservice",
                                   dataset="hsapiens_gene_ensembl")
                }
                g2 <- data.frame(chr=NA, start=NA, end=NA,gene=NA,
                                 ensembl_gene=NA)
                chromosome <- badCoverage[[i]][1,1]
                g3 <- getBM(attributes=c("chromosome_name", "start_position",
                                         "end_position", "external_gene_name",
                                         "ensembl_gene_id"),
                            filters="chromosome_name", values=chromosome, mart)
                for(k in 1:length(badCoverage[[i]][,1])){
                    distance_best <- 999999999
                    for(j in 1:length(g3[,1])){
                        if(g3[j,2] <= badCoverage[[i]][k,2]
                           && g3[j,3] >= badCoverage[[i]][k,3]){
                            distance_comp <- (badCoverage[[i]][k,2]-g3[j,2])+(g3[j,3]-badCoverage[[i]][k,3])
                            if(distance_comp < distance_best){
                                badCoverage[[i]][k,5] <- g3[j,4]
                                badCoverage[[i]][k,6] <- g3[j,5]
                                distance_best <- distance_comp
                            }
                        }
                    } 
                    if(is.na(badCoverage[[i]][k,5])){
                        for(j in 1:length(g3[,1])){
                            if(g3[j,2]<=badCoverage[[i]][k,2]
                               && g3[j,3] >= badCoverage[[i]][k,2]
                               || g3[j,2] <= badCoverage[[i]][k,3]
                               && g3[j,3] >= badCoverage[[i]][k,3]
                               || g3[j,2] > badCoverage[[i]][k,2]
                               && g3[j,3] < badCoverage[[i]][k,3]){
                                badCoverage[[i]][k,5] <- g3[j,4]
                                badCoverage[[i]][k,6] <- g3[j,5]
                            }
                        } 
                    }
                }
            } 
        }
    }    
    first <- 0
    found <- FALSE
    for(i in 1:length(badCoverage)){
        if(found == FALSE && !is.na(badCoverage[[i]][1,1])){
            first <- i
            found <- TRUE
        }
    }
    badCoverageSummary <- badCoverage[[first]]
    if(first < length(badCoverage)){
        for(i in (first+1):length(badCoverage)){
            if(!is.na(badCoverage[[i]][1,1])){
                badCoverageSummary <- rbind(badCoverageSummary,badCoverage[[i]])
            }
        }
    }
    if(output != ""){
        write.table(badCoverageSummary,
                    paste(output, "/BadCoverageSummary", threshold1, ";",
                          percentage1, ";", threshold2, ";", percentage2,
                          ".txt", sep=""),
                    row.names=FALSE, sep="\t", quote=FALSE)   
    }
    badCoverageSummary <- makeGRangesFromDataFrame(badCoverageSummary,
                                                   start.field="start",
                                                   end.field="end",
                                                   keep.extra.columns=TRUE)
    return(badCoverageSummary)
}

#detailed variant
reportBadRegionsDetailed <- function(threshold1, threshold2, percentage1, 
                                     percentage2, coverage_indicators, mart,
                                     samples, output){
    coverage_indicators_temp <- list()
    for(i in 1:length(coverage_indicators)){
        coverage_indicators_temp[[i]] <- data.frame()
    }
    for(i in 1:length(coverage_indicators_temp)){
        message("Analyzing Chromosome ", i)
        if(!is.na(as.data.frame(coverage_indicators[[i]])[1,1])){
            coverage_temp <- as.data.frame(coverage_indicators[[i]])[,c(1:3,6:length(as.data.frame(coverage_indicators[[i]])[1,]))]
            coverage_indicators_temp[[i]] <- coverage_temp[coverage_temp[,3]-coverage_temp[2]==0,]
            coverage_indicators_temp[[i]] <- cbind(coverage_indicators_temp[[i]],
                                                   c(rep(NA,length(coverage_indicators_temp[[i]][,1]))))
            coverage_indicators_temp[[i]] <- cbind(coverage_indicators_temp[[i]],
                                                   c(rep(NA,length(coverage_indicators_temp[[i]][,1]))))
            if(!is.na(coverage_indicators_temp[[i]][1,1])){
                if(is.character(mart)){
                    mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                   host="grch37.ensembl.org",
                                   path="/biomart/martservice",
                                   dataset="hsapiens_gene_ensembl")
                }
                g2 <- data.frame(chr=NA,start=NA,end=NA,gene=NA,ensembl_gene=NA)
                chromosome <- coverage_indicators_temp[[i]][1,1]
                g3 <- getBM(attributes=c("chromosome_name", "start_position",
                                         "end_position", "external_gene_name",
                                         "ensembl_gene_id"),
                            filters="chromosome_name", values=chromosome, mart)
                k <- 1
                while(k <= length(coverage_indicators_temp[[i]][,1])){
                    end <- 0
                    distance_best <- 999999999
                    for(j in 1:length(g3[,1])){
                        if(g3[j,2] <= coverage_indicators_temp[[i]][k,2]
                           && g3[j,3] >= coverage_indicators_temp[[i]][k,2]){
                            distance_comp <- (coverage_indicators_temp[[i]][k,2]-g3[j,2])+(g3[j,3]-coverage_indicators_temp[[i]][k,2])
                            if(distance_comp < distance_best){
                                coverage_indicators_temp[[i]][k,length(coverage_indicators_temp[[i]][1,])-1] <- g3[j,4]
                                coverage_indicators_temp[[i]][k,length(coverage_indicators_temp[[i]][1,])] <- g3[j,5]
                                end <- g3[j,3]
                                distance_best <- distance_comp
                            }
                        }
                    } 
                    if(is.na(coverage_indicators_temp[[i]][k,length(coverage_indicators_temp[[i]][1,])])){
                        k <- k+1
                    }
                    if(!is.na(coverage_indicators_temp[[i]][k,length(coverage_indicators_temp[[i]][1,])])){
                        k <- k+1
                        while(k < length(coverage_indicators_temp[[i]][,1])
                              && end >= coverage_indicators_temp[[i]][k,2]){
                            coverage_indicators_temp[[i]][k,length(coverage_indicators_temp[[i]][1,])-1] <- coverage_indicators_temp[[i]][k-1,length(coverage_indicators_temp[[i]][1,])-1]
                            coverage_indicators_temp[[i]][k,length(coverage_indicators_temp[[i]][1,])] <- coverage_indicators_temp[[i]][k-1,length(coverage_indicators_temp[[i]][1,])]
                            k <- k+1
                        }     
                    }
                }
                names(coverage_indicators_temp[[i]]) <- c("Chr", "Start", "End",
                                                          as.character(samples[,1]),
                                                          "targetBases",
                                                          "QualityMarker",
                                                          "Gene", "GeneID")
                if(output != ""){
                    write.table(coverage_indicators_temp[[i]],
                                paste(output, "/BadCoverageChr", i, ";",
                                      threshold1, ";", percentage2, ";", 
                                      threshold2, ";", percentage2, ".txt", 
                                      sep=""), row.names=FALSE, sep="\t",
                                quote=FALSE)
                }
                coverage_indicators_temp[[i]] <- makeGRangesFromDataFrame(coverage_indicators_temp[[i]],
                                                                          start.field="start",
                                                                          end.field="end",
                                                                          keep.extra.columns=TRUE)
            }
        }
        if(is.na(as.data.frame(coverage_indicators[[i]])[1,1])){
            coverage_indicators_temp[[i]] <- NA
        }
    }
    return(coverage_indicators_temp)
}

#summary for genes
reportBadRegionsGenes <- function(threshold1, threshold2, percentage1, 
                                  percentage2, badCoverageSummary, output){
    badCoverageOverview<-data.frame(Gene=NA,GeneID=NA,
                                    BadRegion_offTarget=0,
                                    BadRegion_onTarget=0,
                                    AcceptableRegion_offTarget=0,
                                    AcceptableRegion_onTarget=0,
                                    GoodRegions_offTarget=0,
                                    GoodRegions_onTarget=0)
    coverage_temp <- as.data.frame(badCoverageSummary)[,c(1:3,6:length(as.data.frame(badCoverageSummary)[1,]))]
    pointer <- 1
    badCoverageOverview[1,1] <- as.character(coverage_temp[1,5])
    badCoverageOverview[1,2] <- as.character(coverage_temp[1,6])
    badCoverageOverview[1,3+coverage_temp[1,4]] <- coverage_temp[1,3]-coverage_temp[1,2]+1
    if(length(coverage_temp[,1]) > 1){
        for(i in 2:length(coverage_temp[,1])){
            if((is.na(coverage_temp[i,5]) && is.na(coverage_temp[i-1,5]))
               || !is.na(coverage_temp[i,5]) && !is.na(coverage_temp[i-1,5])
               && coverage_temp[i,5] == coverage_temp[i-1,5]){
                if(coverage_temp[i,4] == 0){
                    badCoverageOverview[pointer,3] <- badCoverageOverview[pointer,3]+coverage_temp[i,3]-coverage_temp[i,2]+1
                }
                if(coverage_temp[i,4] == 1){
                    badCoverageOverview[pointer,4] <- badCoverageOverview[pointer,4]+coverage_temp[i,3]-coverage_temp[i,2]+1
                }
                if(coverage_temp[i,4] == 2){
                    badCoverageOverview[pointer,5] <- badCoverageOverview[pointer,5]+coverage_temp[i,3]-coverage_temp[i,2]+1
                }
                if(coverage_temp[i,4] == 3){
                    badCoverageOverview[pointer,6] <- badCoverageOverview[pointer,6]+coverage_temp[i,3]-coverage_temp[i,2]+1
                }
                if(coverage_temp[i,4] == 4){
                    badCoverageOverview[pointer,7] <- badCoverageOverview[pointer,7]+coverage_temp[i,3]-coverage_temp[i,2]+1
                }
                if(coverage_temp[i,4] == 5){
                    badCoverageOverview[pointer,8] <- badCoverageOverview[pointer,8]+coverage_temp[i,3]-coverage_temp[i,2]+1
                }        
            }
            if((is.na(coverage_temp[i,5]) && !is.na(coverage_temp[i-1,5]))
               || (!is.na(coverage_temp[i,5]) && is.na(coverage_temp[i-1,5]))
               || !is.na(coverage_temp[i,5]) && !is.na(coverage_temp[i-1,5])
               && coverage_temp[i,5] != coverage_temp[i-1,5]){
                pointer <- pointer+1
                badCoverageOverview[pointer,1] <- as.character(coverage_temp[i,5])
                badCoverageOverview[pointer,2] <- as.character(coverage_temp[i,6])
                if(coverage_temp[i,4] == 0){
                    badCoverageOverview[pointer,3] <- coverage_temp[i,3]-coverage_temp[i,2]+1
                    badCoverageOverview[pointer,4] <- 0
                    badCoverageOverview[pointer,5] <- 0
                    badCoverageOverview[pointer,6] <- 0
                    badCoverageOverview[pointer,7] <- 0
                    badCoverageOverview[pointer,8] <- 0            
                }
                if(coverage_temp[i,4] == 1){
                    badCoverageOverview[pointer,4] <- coverage_temp[i,3]-coverage_temp[i,2]+1
                    badCoverageOverview[pointer,3] <- 0
                    badCoverageOverview[pointer,5] <- 0
                    badCoverageOverview[pointer,6] <- 0
                    badCoverageOverview[pointer,7] <- 0
                    badCoverageOverview[pointer,8] <- 0  
                }
                if(coverage_temp[i,4] == 2){
                    badCoverageOverview[pointer,5] <- coverage_temp[i,3]-coverage_temp[i,2]+1
                    badCoverageOverview[pointer,3] <- 0
                    badCoverageOverview[pointer,4] <- 0
                    badCoverageOverview[pointer,6] <- 0
                    badCoverageOverview[pointer,7] <- 0
                    badCoverageOverview[pointer,8] <- 0  
                }
                if(coverage_temp[i,4] == 3){
                    badCoverageOverview[pointer,6] <- coverage_temp[i,3]-coverage_temp[i,2]+1
                    badCoverageOverview[pointer,3] <- 0
                    badCoverageOverview[pointer,4] <- 0
                    badCoverageOverview[pointer,5] <- 0
                    badCoverageOverview[pointer,7] <- 0
                    badCoverageOverview[pointer,8] <- 0  
                }
                if(coverage_temp[i,4] == 4){
                    badCoverageOverview[pointer,7] <- coverage_temp[i,3]-coverage_temp[i,2]+1
                    badCoverageOverview[pointer,3] <- 0
                    badCoverageOverview[pointer,4] <- 0
                    badCoverageOverview[pointer,5] <- 0
                    badCoverageOverview[pointer,6] <- 0
                    badCoverageOverview[pointer,8] <- 0  
                }
                if(coverage_temp[i,4] == 5){
                    badCoverageOverview[pointer,8] <- coverage_temp[i,3]-coverage_temp[i,2]+1
                    badCoverageOverview[pointer,3] <- 0
                    badCoverageOverview[pointer,4] <- 0
                    badCoverageOverview[pointer,5] <- 0
                    badCoverageOverview[pointer,6] <- 0
                    badCoverageOverview[pointer,7] <- 0  
                } 
            }
        }
    }    
    summe <- rowSums(badCoverageOverview[,3:8])
    for(i in 1:length(badCoverageOverview[,1])){
        for(j in 3:8){
            badCoverageOverview[i,j] <- badCoverageOverview[i,j]/summe[i]
        }
    }
    if(output != ""){
        write.table(badCoverageOverview,
                    paste(output, "/BadCoverageGenes", threshold1, ";",
                          percentage1, ";", threshold2, ";", percentage2,
                          ".txt", sep=""), row.names=FALSE, sep="\t",
                    quote=FALSE)   
    }
    return(badCoverageOverview)
}
