###sum up coverage on one chromosome
sumUpCoverage <- function(coverage_samples, targetRegions_coverage, samples, 
                          chr, flag){
    summary <- data.frame()
    #first line
    summary[1,1] <- chr
    summary[1,2] <- 1
    if(runValue(targetRegions_coverage)[1] == 0){
        start <- runLength(targetRegions_coverage)[1]
    }
    for(i in 1:length(coverage_samples)){
        if(runValue(coverage_samples[[i]])[1] == 0 
           && runLength(coverage_samples[[i]])[1] < start){
            start <- runLength(coverage_samples[[i]])[1]
        }
    }
    summary[1,3] <- start
    for(j in 1:length(samples[,1])){
        summary[1,3+j] <- 0
    }
    summary[1,4+length(samples[,1])] <- 0
    names(summary) <- c("Chr", "start", "end", as.character(samples[,1]), 
                        "targetBases")
    
    #remove first uncovered bases
    for(i in 1:length(coverage_samples)){
        runLength(coverage_samples[[i]])[1] <- runLength(coverage_samples[[i]])[1]-start
    }
    runLength(targetRegions_coverage)[1] <- runLength(targetRegions_coverage)[1]-start
    
    coverage_line <- 2
    while(length(coverage_samples[[1]]) > 0){
        summary[coverage_line,1] <- chr
        summary[coverage_line,2] <- summary[coverage_line-1,3]+1
        summary[coverage_line,3] <- summary[coverage_line-1,3]+1
        for(j in 1:length(samples[,1])){
            summary[coverage_line,3+j] <- runValue(coverage_samples[[j]])[1]
            runLength(coverage_samples[[j]])[1] <- runLength(coverage_samples[[j]])[1]-1
        }
        if(is.na(runValue(targetRegions_coverage)[1]) == TRUE){
            summary[coverage_line,4+length(samples[,1])] <- 0
        }
        if(is.na(runValue(targetRegions_coverage)[1]) == FALSE){
            summary[coverage_line,4+length(samples[,1])] <- runValue(targetRegions_coverage)[1]
            runLength(targetRegions_coverage)[1] <- runLength(targetRegions_coverage)[1]-1
        }
        coverage_line <- coverage_line+1
        
        if(sum(summary[coverage_line-1,4:(4+length(samples[,1]))]) == 0){
            rest_coverage <- 0
            rest_regions <- 0
            for(i in 1:length(coverage_samples)){
                rest_regions <- length(runLength(coverage_samples[[i]]))+rest_regions
                rest_coverage <- sum(runValue(coverage_samples[[i]]))+rest_coverage
            }
            if(is.na(runValue(targetRegions_coverage)[1]) == TRUE
               && rest_regions == length(coverage_samples)
               && rest_coverage == 0){
                summary[(coverage_line-1),3] <- summary[(coverage_line-1),3]+length(coverage_samples[[1]])
                runLength(coverage_samples[[1]])[1] <- runLength(coverage_samples[[1]])[1]-runLength(coverage_samples[[1]])[1]
            }
            if(is.na(runValue(targetRegions_coverage)[1]) == FALSE
               || rest_regions != length(coverage_samples)
               || rest_coverage != 0){
                no_coverage <- runLength(coverage_samples[[1]])[1]
                if(length(coverage_samples) > 1){
                    for(i in 2:length(coverage_samples)){
                        if(runValue(coverage_samples[[i]])[1] == 0
                           && runLength(coverage_samples[[i]])[1] < no_coverage){
                            no_coverage <- runLength(coverage_samples[[i]])[1]
                        }
                    }  
                }
                if(is.na(runValue(targetRegions_coverage)[1]) == FALSE
                   && runLength(targetRegions_coverage)[1] < no_coverage){
                    no_coverage <- runLength(targetRegions_coverage)[1]
                }
                
                summary[(coverage_line-1),3] <- summary[coverage_line-1,3]+no_coverage
                for(j in 1:length(samples[,1])){
                    runLength(coverage_samples[[j]])[1] <- runLength(coverage_samples[[j]])[1]-no_coverage
                }
                if(is.na(runValue(targetRegions_coverage)[1]) == FALSE){
                    runLength(targetRegions_coverage)[1] <- runLength(targetRegions_coverage)[1]-no_coverage
                }
                summary<-summary[1:(coverage_line-1),]
            }
            
        }
    }
    if(flag == FALSE){
        return(summary)
    }
    if(flag == TRUE){
        return(summary[summary[,4+length(samples[,1])]==1,])
    }
}

###sum up coverage on one chromosome (not targeted)
sumUpCoverageNoTarget <- function(coverage_samples, samples, chr){
    summary <- data.frame()
    #first line
    summary[1,1] <- chr
    summary[1,2] <- 1
    if(runValue(coverage_samples[[1]])[1] == 0){
        start <- runLength(coverage_samples[[1]])[1]
    }
    if(length(coverage_samples) > 1){
        for(i in 2:length(coverage_samples)){
            if(runValue(coverage_samples[[i]])[1] == 0
               && runLength(coverage_samples[[i]])[1] < start){
                start <- runLength(coverage_samples[[i]])[1]
            }
        }  
    }
    summary[1,3] <- start
    for(j in 1:length(samples[,1])){
        summary[1,3+j] <- 0
    }
    summary[1,4+length(samples[,1])] <- 0
    names(summary) <- c("Chr", "start", "end", as.character(samples[,1]),
                        "targetBases")
    
    #remove first uncovered bases
    for(i in 1:length(coverage_samples)){
        runLength(coverage_samples[[i]])[1] <- runLength(coverage_samples[[i]])[1]-start
    }
    
    coverage_line <- 2
    while(length(coverage_samples[[1]]) > 0){
        summary[coverage_line,1] <- chr
        summary[coverage_line,2] <- summary[coverage_line-1,3]+1
        summary[coverage_line,3] <- summary[coverage_line-1,3]+1
        for(j in 1:length(samples[,1])){
            summary[coverage_line,3+j] <- runValue(coverage_samples[[j]])[1]
            runLength(coverage_samples[[j]])[1] <- runLength(coverage_samples[[j]])[1]-1
        }
        summary[coverage_line,4+length(samples[,1])] <- 0
        coverage_line <- coverage_line+1
        
        if(sum(summary[coverage_line-1,4:(4+length(samples[,1]))]) == 0){
            rest_coverage <- 0
            rest_regions <- 0
            for(i in 1:length(coverage_samples)){
                rest_regions <- length(runLength(coverage_samples[[i]]))+rest_regions
                rest_coverage <- sum(runValue(coverage_samples[[i]]))+rest_coverage
            }
            if(rest_regions == length(coverage_samples)
               && rest_coverage == 0){
                summary[(coverage_line-1),3] <- summary[(coverage_line-1),3]+length(coverage_samples[[1]])
                runLength(coverage_samples[[1]])[1] <- runLength(coverage_samples[[1]])[1]-runLength(coverage_samples[[1]])[1]
            }
            if(rest_regions != length(coverage_samples)
               || rest_coverage != 0){
                no_coverage <- runLength(coverage_samples[[1]])[1]
                if(length(coverage_samples) > 1){
                    for(i in 2:length(coverage_samples)){
                        if(runValue(coverage_samples[[i]])[1] == 0
                           && runLength(coverage_samples[[i]])[1] < no_coverage){
                            no_coverage <- runLength(coverage_samples[[i]])[1]
                        }
                    }  
                }
                
                summary[(coverage_line-1),3] <- summary[coverage_line-1,3]+no_coverage
                for(j in 1:length(samples[,1])){
                    runLength(coverage_samples[[j]])[1] <- runLength(coverage_samples[[j]])[1]-no_coverage
                }
                summary <- summary[1:(coverage_line-1),]
            }
            
        }
    }
    return(summary)
}


#scans the whole genome, but if output on the target regions only is desired, the
#flag TRonly may be set to TRUE:
determineCoverage <- function(samples, bam_input, targetRegions, output, TRonly){
    coverage_summary <- list()
    coverage_temp <- list()
    
    message("Determine Coverage")
    message("Sample ", samples[1,1])
    if(is.character(bam_input) == TRUE){
        coverage_temp[[1]] <- coverage(paste(bam_input, "/", samples[1,1], 
                                             ".bam", sep=""))
    }
    if(is.character(bam_input) == FALSE){
        coverage_temp[[1]] <- coverage(path(bam_input)[1])
    } 
    if(length(samples[,1]) > 1){
        for(k in 2:length(samples[,1])){
            message("Sample ", samples[k,1])
            if(is.character(bam_input) == TRUE){
                coverage_temp[[k]] <- coverage(paste(bam_input, "/", 
                                                     samples[k,1], ".bam",
                                                     sep=""))
            }
            if(is.character(bam_input) == FALSE){
                coverage_temp[[k]] <- coverage(path(bam_input)[1])
            } 
        }
    }
    
    message("Determine target bases")
    if(is.data.frame(targetRegions) == TRUE){
        names(targetRegions) <- c("Chr", "start", "end")
        targetRegions_temp <- makeGRangesFromDataFrame(targetRegions, 
                                                       start.field="start",
                                                       end.field="end")
        targetRegions <- targetRegions_temp
        for(i in 1:length(targetRegions)){
            start(targetRegions[i]) <- start(targetRegions[i])+1
        }
    }
    targetRegions_cov <- as.list(coverage(targetRegions))
    
    message("Combine Information")
    target <- 1
    for(i in 1:length(coverage_temp[[1]])){
        message("Chromosome: ", names(coverage_temp[[1]])[i])
        covered <- FALSE
        targeted <- FALSE
        for(j in 1:length(coverage_temp)){
            if(nrun(coverage_temp[[j]][[i]]) > 1
               || (nrun(coverage_temp[[j]][[i]]) == 1
                   && runValue(coverage_temp[[j]][[i]])[1] > 0)){
                covered <- TRUE
            }
            if(!is.na(names(targetRegions_cov)[target])
               && names(targetRegions_cov)[target] == names(coverage_temp[[j]])[i]){
                targeted <- TRUE
            }
        }
        if(covered == TRUE && targeted == TRUE
           || covered == FALSE && targeted == TRUE){
            coverage_samples <- list()
            for(k in 1:length(coverage_temp)){
                coverage_samples[[k]] <- coverage_temp[[k]][[i]]
            }
            coverage_summary[[i]] <- sumUpCoverage(coverage_samples,
                                                   targetRegions_cov[[target]],
                                                   samples,
                                                   names(coverage_temp[[1]][i]),
                                                   TRonly)
            target <- target+1
            if(output != ""){
                write.table(coverage_summary[[i]],
                            paste(output, "/Summary_chr",
                                  names(coverage_temp[[1]])[i], ".txt", sep=""),
                            row.names=FALSE, sep="\t", quote=FALSE)   
            }
            coverage_summary[[i]] <- makeGRangesFromDataFrame(coverage_summary[[i]],
                                                              start.field="start",
                                                              end.field="end",
                                                              keep.extra.columns=TRUE)
        }
        if(covered == TRUE && targeted == FALSE && TRonly == FALSE){
            coverage_samples <- list()
            for(k in 1:length(coverage_temp)){
                coverage_samples[[k]] <- coverage_temp[[k]][[i]]
            }
            coverage_summary[[i]] <- sumUpCoverageNoTarget(coverage_samples,
                                                           samples,
                                                           names(coverage_temp[[1]][i]))
            if(output != ""){
                write.table(coverage_summary[[i]],
                            paste(output, "/Summary_chr",
                                  names(coverage_temp[[1]])[i], ".txt", sep=""),
                            row.names=FALSE, sep="\t", quote=FALSE)   
            }
            coverage_summary[[i]] <- makeGRangesFromDataFrame(coverage_summary[[i]],
                                                              start.field="start",
                                                              end.field="end",
                                                              keep.extra.columns=TRUE)
        }
        if(covered == FALSE && targeted == FALSE && TRonly == FALSE){
            coverage_summary[[i]] <- data.frame()
            coverage_summary[[i]][1,1] <- names(coverage_temp[[1]][i])
            coverage_summary[[i]][1,2] <- 1
            coverage_summary[[i]][1,3] <- runLength(coverage_temp[[1]][i])[[1]]
            for(j in 1:length(samples[,1])){
                coverage_summary[[i]][1,3+j] <- 0
            }
            coverage_summary[[i]][1,4+length(samples[,1])] <- 0
            names(coverage_summary[[i]]) <- c("Chr", "start", "end",
                                              as.character(samples[,1]),
                                              "targetBases")
            if(output != ""){
                write.table(coverage_summary[[i]],
                            paste(output, "/Summary_chr",
                                  names(coverage_temp[[1]])[i], ".txt", sep=""),
                            row.names=FALSE, sep="\t", quote=FALSE)   
            }
            coverage_summary[[i]] <- makeGRangesFromDataFrame(coverage_summary[[i]],
                                                              start.field="start",
                                                              end.field="end",
                                                              keep.extra.columns=TRUE)
        }
        if(covered == TRUE && targeted == FALSE && TRonly == TRUE
           || covered == FALSE && targeted == FALSE && TRonly == TRUE){
            coverage_summary[[i]] <- data.frame()
            coverage_summary[[i]][1,1] <- names(coverage_temp[[1]][i])
            coverage_summary[[i]][1,2] <- 0
            coverage_summary[[i]][1,3] <- 0
            for(j in 1:length(samples[,1])){
                coverage_summary[[i]][1,3+j] <- 0
            }
            coverage_summary[[i]][1,4+length(samples[,1])] <- 0
            names(coverage_summary[[i]]) <- c("Chr", "start", "end",
                                              as.character(samples[,1]),
                                              "targetBases")
            coverage_summary[[i]] <- makeGRangesFromDataFrame(coverage_summary[[i]],
                                                              start.field="start",
                                                              end.field="end",
                                                              keep.extra.columns=TRUE)
        }
    }
    coverage_summary_list <- GRangesList("Chromosome 1" = coverage_summary[[1]],
                                         "Chromosome 2" = coverage_summary[[2]],
                                         "Chromosome 3" = coverage_summary[[3]],
                                         "Chromosome 4" = coverage_summary[[4]],
                                         "Chromosome 5" = coverage_summary[[5]],
                                         "Chromosome 6" = coverage_summary[[6]],
                                         "Chromosome 7" = coverage_summary[[7]],
                                         "Chromosome 8" = coverage_summary[[8]],
                                         "Chromosome 9" = coverage_summary[[9]],
                                         "Chromosome 10" = coverage_summary[[10]],
                                         "Chromosome 11" = coverage_summary[[11]],
                                         "Chromosome 12" = coverage_summary[[12]],
                                         "Chromosome 13" = coverage_summary[[13]],
                                         "Chromosome 14" = coverage_summary[[14]],
                                         "Chromosome 15" = coverage_summary[[15]],
                                         "Chromosome 16" = coverage_summary[[16]],
                                         "Chromosome 17" = coverage_summary[[17]],
                                         "Chromosome 18" = coverage_summary[[18]],
                                         "Chromosome 19" = coverage_summary[[19]],
                                         "Chromosome 20" = coverage_summary[[20]],
                                         "Chromosome 21" = coverage_summary[[21]],
                                         "Chromosome 22" = coverage_summary[[22]],
                                         "Chromosome X" = coverage_summary[[23]],
                                         "Chromosome Y" = coverage_summary[[24]],
                                         "Chromosome MT" = coverage_summary[[25]])
    return(coverage_summary_list)
}


#separate into no coverage (0 and 1), bad region (2 and 3) 
#and good region (4 and 5) (works for whole genome and target only)
determineCoverageQuality <- function(threshold1, threshold2, percentage1, 
                                     percentage2, coverage_summary){
    coverage_indicators <- list()
    for(i in 1:length(coverage_summary)){
        message("Analyzing Chromosome ", i)
        if(as.data.frame(coverage_summary[[i]])[1,2] != 0){
            coverage_temp <- as.data.frame(coverage_summary[[i]])[,c(1:3,6:length(as.data.frame(coverage_summary[[i]])[1,]))]
            coverage_indicators[[i]] <- data.frame(coverage_temp, indicator=NA)
            for(j in 1:length(coverage_temp[,1])){
                if(threshold1 > threshold2){
                    message("Threshold1 is greater than threshold2. This won't lead to useful results!")
                }
                bad <- sum(as.numeric(coverage_temp[j,c(4:(length(coverage_temp)-1))] < threshold1))
                if(sum(as.numeric(coverage_temp[j,c(4:(length(coverage_temp)-1))] >= threshold1)) > 0){
                    acceptable <- sum(as.numeric(coverage_temp[j,c(FALSE, FALSE, FALSE, coverage_temp[j,c(4:(length(coverage_temp)-1))] >= threshold1, FALSE)] < threshold2))
                }
                if(sum(as.numeric(coverage_temp[j,c(4:(length(coverage_temp)-1))] >= threshold1)) == 0){
                    acceptable <- 0
                }
                good <- sum(as.numeric(coverage_temp[j,c(4:(length(coverage_temp)-1))] >= threshold2))
                flag <- FALSE
                
                if(percentage2 == 0){
                    message("It is not wise to set percentage2 to zero. Please increase the value to obtain usefull results!")
                }
                #good region: at least percentage2 percent of samples with good coverage
                if(good >= percentage2*(length(coverage_temp[j,])-4)){
                    coverage_indicators[[i]][j,length(coverage_temp[j,])+1] <- 4+as.numeric(coverage_temp[j,length(coverage_temp)]>0)
                    flag <- TRUE
                }
                #acceptable region: at least percentage1 percent of samples with good or acceptable coverage
                if(flag == FALSE && (acceptable+good) >= percentage1*(length(coverage_temp[j,])-4)){
                    coverage_indicators[[i]][j,length(coverage_temp[j,])+1] <- 2+as.numeric(coverage_temp[j,length(coverage_temp)]>0)
                    flag <- TRUE
                }
                #bad region: not even percentage1 percent of samples with good or acceptable coverage
                if(flag == FALSE && (acceptable+good) < percentage1*(length(coverage_temp[j,])-4)){
                    coverage_indicators[[i]][j,length(coverage_temp[j,])+1] <- 0+as.numeric(coverage_temp[j,length(coverage_temp)]>0)
                    flag <- TRUE
                }
                
            }  
            coverage_indicators[[i]] <- makeGRangesFromDataFrame(coverage_indicators[[i]],
                                                                 start.field="start",
                                                                 end.field="end",
                                                                 keep.extra.columns=TRUE)
        }
        if(as.data.frame(coverage_summary[[i]])[1,2] == 0){
            coverage_indicators[[i]] <- NA
        }
    }
    return(coverage_indicators)
}

#select other regions of interest:
#available as bed-file 
determineRegionsOfInterest <- function(regionsOfInterest, coverage_indicators){
    coverage_indicators_2 <- list()
    for(i in 1:length(coverage_indicators)){
        coverage_indicators_2[[i]] <- data.frame()
    }
    pointer_target <- 1
    if(is.data.frame(regionsOfInterest) == FALSE){
        regionsOfInterest <- as.data.frame(regionsOfInterest)
    }
    for(i in 1:length(coverage_indicators)){  
        message("Analyze Chromosome ", i)
        pointer_all <- 1
        pointer_new <- 1 
        if(!is.na(as.data.frame(coverage_indicators[[i]])[1,1])){
            coverage_temp <- as.data.frame(coverage_indicators[[i]])[,c(1:3,6:length(as.data.frame(coverage_indicators[[i]])[1,]))]
            while(length(coverage_temp[,1]) >= pointer_all
                  && pointer_target <= length(regionsOfInterest[,1])
                  && as.character(regionsOfInterest[pointer_target,1]) == coverage_temp[pointer_all,1]
                  && regionsOfInterest[pointer_target,2] <= regionsOfInterest[pointer_target,3]){
                if(pointer_target <= length(regionsOfInterest[,1])
                   && regionsOfInterest[pointer_target,2] <= coverage_temp[pointer_all,3]
                   && regionsOfInterest[pointer_target,2] <= regionsOfInterest[pointer_target,3]){
                    coverage_indicators_2[[i]][pointer_new,1] <- as.character(regionsOfInterest[pointer_target,1])
                    coverage_indicators_2[[i]][pointer_new,2] <- regionsOfInterest[pointer_target,2]
                    coverage_indicators_2[[i]][pointer_new,3] <- regionsOfInterest[pointer_target,2]
                    for(j in 4:length(coverage_temp[pointer_all,])){
                        coverage_indicators_2[[i]][pointer_new,j] <- coverage_temp[pointer_all,j]
                    }
                    pointer_new <- pointer_new+1
                    if(pointer_target<=length(regionsOfInterest[,1])
                       &&regionsOfInterest[pointer_target,2] == coverage_temp[pointer_all,3]){
                        pointer_all <- pointer_all+1
                    }
                    if(pointer_target <= length(regionsOfInterest[,1])
                       && length(coverage_temp[,1]) >= pointer_all
                       && regionsOfInterest[pointer_target,2] < coverage_temp[pointer_all,3]
                       && regionsOfInterest[pointer_target,2] <= regionsOfInterest[pointer_target,3]){
                        regionsOfInterest[pointer_target,2] <- regionsOfInterest[pointer_target,2]+1
                    }
                    if(pointer_target <= length(regionsOfInterest[,1])
                       && length(coverage_temp[,1]) >= pointer_all
                       && regionsOfInterest[pointer_target,2] > regionsOfInterest[pointer_target,3]){
                        pointer_target <- pointer_target+1
                    }
                }
                if(length(coverage_temp[,1]) >= pointer_all
                   && pointer_target <= length(regionsOfInterest[,1])
                   && regionsOfInterest[pointer_target,2] > coverage_temp[pointer_all,3]){
                    pointer_all <- pointer_all+1
                }
            }
            if(length(coverage_indicators_2[[i]]) > 0){
                names(coverage_indicators_2[[i]]) <- names(coverage_temp)
                coverage_indicators_2[[i]] <- makeGRangesFromDataFrame(coverage_indicators_2[[i]],
                                                                       start.field="start",
                                                                       end.field="end",
                                                                       keep.extra.columns=TRUE)
            }
            if(length(coverage_indicators_2[[i]]) == 0){
                coverage_indicators_2[[i]] <- NA
            } 
        }
        if(is.na(as.data.frame(coverage_indicators[[i]])[1,1])){
            coverage_indicators_2[[i]] <- NA
        } 
    }
    return(coverage_indicators_2)
}