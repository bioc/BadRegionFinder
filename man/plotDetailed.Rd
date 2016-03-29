\name{plotDetailed}
\alias{plotDetailed}
\alias{plot Detailed}
\title{
Plots a more detailed overview of the coverage quality
}
\description{
The function \code{plotDetailed} provides a possibility to visualize the output of \code{reportBadRegionsDetailed}. A line graph is returned, visualizing the median coverage at each base that was chosen to be in the region of interest. For each base, the category of coverage quality is color coded. Furthermore, information on the genes that are located at the positions analyzed is included.
}
\usage{
plotDetailed(threshold1, threshold2, percentage1, percentage2, 
             coverage_indicators_temp, output)
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
  \item{coverage_indicators_temp}{
List object, return value of function \code{reportBadRegionsDetailed}.
}
  \item{output}{
The folder to write the output file into. If this argument is an empty string, the plot is printed on the screen.
}
}
\details{
The function \code{plotDetailed} serves to summarize the previously determined coverage quality in a visual way, including additionally information on the mean coverage over all samples at every position.

On the y axis the median coverage over all samples is coded. Every position is considered individually.

On the x axis the detected genes are printed. Wherever a new region covering a new gene is registered, a dashed line is drawn. 

Yet, additionally to the mere median coverage, the corresponding coverage quality at each position is also included in the plot. The different categories of coverage quality are color coded in the following way: red - bad region on target; yellow - acceptable region on target; green - good region on target; black - bad region off target; dark gray - acceptable region off target; light gray - good region off target.
}
\value{
No value is returned.
}

\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
}
\examples{
library("BSgenome.Hsapiens.UCSC.hg19")
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
                                      TRonly = FALSE)
coverage_indicators <- determineCoverageQuality(threshold1, threshold2,
                                                percentage1, percentage2,
                                                coverage_summary)
coverage_indicators_temp <- reportBadRegionsDetailed(threshold1, threshold2,
                                                     percentage1, percentage2,
                                                     coverage_indicators, "",
                                                     samples, output)
plotDetailed(threshold1, threshold2, percentage1, percentage2,
             coverage_indicators_temp, output)
}

