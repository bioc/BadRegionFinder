\name{plotSummary}
\alias{plotSummary}
\alias{plot Summary}
\title{
Plots a summary of the coverage quality
}
\description{
The function \code{plotSummary} provides a possibility to visualize the output of \code{reportBadRegionsSummary}. A line graph is returned, visualizing the number of bases that fall into each category of coverage quality. Furthermore, information on the genes located in these regions is included.
}
\usage{
plotSummary(threshold1, threshold2, percentage1, percentage2, 
            badCoverageSummary, output)
}
\arguments{
  \item{threshold1}{
Integer, threshold defining the number of reads that have to be registered for a sample that its coverage is classified as acceptable.
}
  \item{threshold2}{
Integer, threshold defining the number of reads that have to be registered for a sample that its coverage is classified as good.
}
  \item{percentage1}{
Float, defining the percentage of samples that have to feature a coverage of at least \code{threshold1} so that the position is classified as acceptably covered.
}
  \item{percentage2}{
Float, defining the percentage of samples that have to feature a coverage of at least \code{threshold2} so that the position is classified as well covered.
}
  \item{badCoverageSummary}{
GRanges object, return value of function \code{reportBadRegionsSummary}.
}
  \item{output}{
The folder to write the output file into. If this argument is an empty string, the plot is printed on the screen.
}
}
\details{
The function \code{plotSummary} serves to summarize the previously determined coverage quality in a visual way. 

On the y axis the coverage quality is coded. The different categories are color coded as well as height coded. As numbers from 0 to 5 were previously assigned to the different categories, thick lines are now drawn at the height of the category. Furthermore, the categories are color coded in the following way: red - bad region on target; yellow - acceptable region on target; green - good region on target; black - bad region off target; dark gray - acceptable region off target; light gray - good region off target.

On the x axis the detected genes are printed. Wherever a new region covering a new gene is registered, a dashed line is drawn. 
}
\value{
No value is returned.
}
\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
}
\examples{
library("BSgenome.Hsapiens.UCSC.hg19")
\donttest{
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

coverage_summary <- determineCoverage(samples, bam_input, targetRegions, output,
                                      TRonly = TRUE)
coverage_indicators <- determineCoverageQuality(threshold1, threshold2,
                                                percentage1, percentage2,
                                                coverage_summary)
badCoverageSummary <- reportBadRegionsSummary(threshold1, threshold2,
                                              percentage1, percentage2,
                                              coverage_indicators, "", output)
plotSummary(threshold1, threshold2, percentage1, percentage2, 
            badCoverageSummary, output)
}
}

