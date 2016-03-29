\name{determineCoverageQuality}
\alias{determineCoverageQuality}
\alias{determine Coverage Quality}
\title{
Classifies the determined coverage
}
\description{
The previously determined coverage (using \code{determineCoverage} with \code{TRonly = TRUE} or \code{TRonly = FALSE}) for all samples gets combined to be classified into six categories: bad coverage off target, bad coverage on target, acceptable coverage off target, acceptable coverage on target, good coverage off target, good coverage on target. These categories are user-defined.
}
\usage{
determineCoverageQuality(threshold1, threshold2, percentage1, percentage2, 
                         coverage_summary)
}
\arguments{
  \item{threshold1}{
Integer, threshold defining the number of reads that have to be registered for a sample that its coverage is classified as acceptable.
}
  \item{threshold2}{
Integer, threshold defining the number of reads that have to be registered for a sample that its coverage is classified as good. To obtain useful results, \code{threshold2} has to be greater than \code{threshold1}. 
}
  \item{percentage1}{
Float, defining the percentage of samples that have to feature a coverage of at least \code{threshold1} so that the position is classified as acceptably covered.
}
  \item{percentage2}{
Float, defining the percentage of samples that have to feature a coverage of at least \code{threshold2} so that the position is classified as well covered. To obtain useful results, \code{percentage2} should be greater than zero.
}
  \item{coverage_summary}{
GRangesList object, return value of function \code{determineCoverage}.
}
}
\details{
Every chromosome is analyzed individually. First, the coverage of each sample is categorized according to \code{threshold1} and \code{threshold2} into three different categories:

bad coverage: less than \code{threshold1} reads

acceptable coverage: at least \code{threshold1}, but less than \code{threshold2} reads

good coverage: at least \code{threshold2} reads

Subsequently this information gets combined with the defined precentages to obtain a numerically coded quality value for each region saved in the previously created list object \code{coverage_summary}:

0: off target; not even \code{percentage1} percent of all samples have a good or acceptable coverage (bad region)

1: on target; not even \code{percentage1} percent of all samples have a good or acceptable coverage (bad region)

2: off target; at least \code{percentage1} percent of all samples have a good or acceptable coverage, but less than \code{percentage2} percent of all samples have a good coverage (acceptable region)

3: on target; at least \code{percentage1} percent of all samples have a good or acceptable coverage, but less than \code{percentage2} percent of all samples have a good coverage (acceptable region)

4: off target; at least \code{percentage2} percent of all samples have a good coverage (good region)

5: on target; at least \code{percentage2} percent of all samples have a good coverage (good region)
}
\value{
A list is returned. Every component contains the coverage information of one chromosome as a GRanges object. The metadata columns contain information on the concrete coverage of each sample at a specific position. Furthermore, the column 'TargetBases' contains information on whether the considered region or position contains target bases (value 1) or not (value 0). The column 'indicator' contains information on the coverage quality of the corresponding region/position.

If a chromosome is not covered and was not targeted as well, the GRanges object solely contains a single line, considering a whole chromosome if \code{TRonly=FALSE}. If \code{TRonly=TRUE} the corresponding component is "NA".
}

\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
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
}                                              
}

