\name{plotSummaryGenes}
\alias{plotSummaryGenes}
\alias{plot Summary Genes}
\title{
Plots a summary of the coverage quality concerning the genes only
}
\description{
The function \code{plotSummaryGenes} provides a possibility to visualize the output of \code{reportBadRegionsGenes}. A barplot is returned, visualizing the percent of each gene that falls into each category of coverage quality. The plot thereby serves to quickly distinguish well from bad covered genes.
}
\usage{
plotSummaryGenes(threshold1, threshold2, percentage1, percentage2, 
                 badCoverageGenes, output)
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
  \item{badCoverageGenes}{
Data frame object, return value of function \code{reportBadRegionsGenes}.
}
  \item{output}{
The folder to write the output file into. If this argument is an empty string, the plot is printed on the screen.
}
}
\details{
The function \code{plotSummaryGenes} serves to summarize the previously determined coverage quality in a visual way concerning the genes only. 

For every gene either one or two stacked bars are plotted. If a gene is covered, but it was not originally targeted, a bar is plotted containing the following color code: black - bad region off target; dark gray - acceptable region off target; light gray - good region off target. If a gene was originally targeted, a bar is plotted containing the following color code: red - bad region on target; yellow - acceptable region on target; green - good region on target.

}
\value{
No value is returned.
}

\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{determineQuantiles}}
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
                                      TRonly = TRUE)
coverage_indicators <- determineCoverageQuality(threshold1, threshold2,
                                                percentage1, percentage2,
                                                coverage_summary)
badCoverageSummary <- reportBadRegionsSummary(threshold1, threshold2,
                                              percentage1, percentage2,
                                              coverage_indicators, "", output)
badCoverageGenes <- reportBadRegionsGenes(threshold1, threshold2, percentage1,
                                          percentage2, badCoverageSummary, 
                                          output)
plotSummaryGenes(threshold1, threshold2, percentage1, percentage2,
                 badCoverageGenes, output)             
}

