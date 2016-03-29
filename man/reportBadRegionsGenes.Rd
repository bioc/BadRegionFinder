\name{reportBadRegionsGenes}
\alias{reportBadRegionsGenes}
\alias{report Bad Regions Genes}
\title{
Sums up the coverage quality on a gene basis
}
\description{
The function \code{reportBadRegionsGenes} creates a summary report considering the coverage quality on a genewise level. Taking the output of \code{reportBadRegionsSummary} as an input, the coverage quality for every previously identified gene is reported.
}
\usage{
reportBadRegionsGenes(threshold1, threshold2, percentage1, percentage2, 
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
Data frame object, return value of function \code{reportBadRegionsSummary}.
}
  \item{output}{
The folder to write the output file into. If \code{output} is just an empty string, no output file is written out.
}
}
\details{
To gain an overview of the coverage quality of each targeted/covered gene, a summary file may be created by the function \code{reportBadRegionsGenes}. The function takes the output of \code{reportBadRegionsSummary} as an input.

All regions covering the same gene are summed up in the following way: 

The number of bases falling into each quality category is summed up. Thereby, regions which were orignially targeted may easily be separated from those which were not, as targeted regions always feature an uneven number characterizing their coverage quality. If a region is broader than the detected gene, but the quality category is the same for the whole region, the whole region is assigned to the gene.

If no gene is reported in the input file, the coverage quality is summed up for a gene named "NA".

The output file is saved as: "BadCoverageGenes\code{threshold1};\code{percentage1};\code{threshold2};\code{percentage2}.txt". The output file may be visualized using \code{plotSummaryGenes}.
}
\value{
A data frame object is returned. The first column contains the name and the geneID of the gene. The following columns contain the percentage of bases falling into the following categories: bad region off traget, bad region on target, acceptable region off target, acceptable region on target, good region off target, good region on target.
}
\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
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
                                              coverage_indicators, "",output)
badCoverageGenes <- reportBadRegionsGenes(threshold1, threshold2, percentage1,
                                          percentage2, badCoverageSummary, output)
}
}
