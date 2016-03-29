\name{determineRegionsOfInterest}
\alias{determineRegionsOfInterest}
\alias{determine Regions Of Interest}
\title{
Determines the regions of interest
}
\description{
The function \code{determineRegionsOfInterest} serves to select the coverage information (including the coverage of all samples and - depending on the input object - their assigned quality value) of one or more subsets of regions. 
}
\usage{
determineRegionsOfInterest(regionsOfInterest, coverage_indicators)
}
\arguments{
  \item{regionsOfInterest}{
Data frame- or GRanges object containing the regions of interest (if data frame: chromosome: first column, start position: second column and end position: third column).
}
  \item{coverage_indicators}{
List object, return value of function \code{determineCoverageQuality}.
}
}
\details{
Every chromosome is analyzed individually. For every base defined in \code{regionsOfInterest} the previously determined coverage information is written out. The function thereby serves to select special subsets of regions, e.g. targeted and untargeted regions when using \code{determineCoverage} with \code{TRonly=FALSE} in the first place or particular targeted regions when using \code{TRonly=TRUE}.

It is not recommended to use \code{TRonly=TRUE}, but to select regions off target using \code{determineRegionsOfInterest}. In this case, no coverage will be registered for all bases off target, as this information was not saved during the step of \code{determineCoverage}.

If \code{regionsOfInterest} is identical compared to \code{targetRegions} and the function \code{determineCoverage} with \code{TRonly=TRUE} has been used in the first place, the selection of regions that is returned is not changed at all.  
}
\value{
A list is returned. Every component contains the coverage information of one chromosomeas as a GRanges object. The metadata columns contain information on the concrete coverage of each sample at a specific position. Furthermore, the column 'TargetBases' contains information on whether the considered region or position contains target bases (value 1) or not (value 0). The column 'indicator' contains information on the coverage quality of the corresponding region/position (0: bad region off target; 1: bad region on target; 2: acceptable region off target; 3: acceptable region on target; 4: good region off target; 5: good region on target).

If a chromosome has not been targeted and/or not covered by any sample, but defined in \code{regionsOfInterest}, the component is "NA".
}

\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
}
\examples{
threshold1 <- 20
threshold2 <- 100
percentage1 <- 0.80
percentage2 <- 0.90
sample_file <- system.file("extdata", "SampleNames.txt", 
                           package = "BadRegionFinder")
samples <- read.table(sample_file)
bam_input <- system.file("extdata", package = "BadRegionFinder")
output <- system.file("extdata", package = "BadRegionFinder")
target_regions <- system.file("extdata", "targetRegions.bed",
                              package = "BadRegionFinder")
targetRegions <- read.table(target_regions, header = FALSE,
                            stringsAsFactors = FALSE)

\donttest{
coverage_summary <- determineCoverage(samples, bam_input, targetRegions, output,
                                      TRonly = FALSE)
coverage_indicators <- determineCoverageQuality(threshold1, threshold2,
                                                percentage1, percentage2,
                                                coverage_summary)
coverage_indicators_2 <- determineRegionsOfInterest(targetRegions,
                                                    coverage_indicators)
}
}

