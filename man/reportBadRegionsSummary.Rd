\name{reportBadRegionsSummary}
\alias{reportBadRegionsSummary}
\alias{report Bad Regions Summary}
\title{
Sums up the coverage quality
}
\description{
The function \code{reportBadRegionsSummary} creates a summary report containing all regions of interest, their coverage quality and the corresponding gene (name and geneID).
}
\usage{
reportBadRegionsSummary(threshold1, threshold2, percentage1, percentage2, 
                        coverage_indicators, mart, output)
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
  \item{coverage_indicators}{
List object, return value of function \code{determineCoverageQuality} or \code{determineRegionsOfInterest}.
}
  \item{mart}{
mart as defined in the manual for package 'biomaRt'. If the human genome (hg19) shall be used, an empty string may be provided and the mart is automatically generated.
}
  \item{output}{
The folder to write the output file into. If \code{output} is just an empty string, no output file is written out.
}
}
\details{
To gain an overview of the coverage quality, a summary file may be created by the function \code{reportBadRegionsSummary}. The function may either take information on the whole genome (output from \code{determineCoverage} with \code{TRonly=FALSE}, processed using \code{determineCoverageQuality}) as an input, or information on the target regions (output from \code{determineCoverage} with \code{TRonly=TRUE}, processed using \code{determineCoverageQuality}), or information on a selection of regions of interest (output from \code{determineRegionsOfInterest}). 

Wherever subsequent bases feature the same coverage quality, the region gets summed up. Although it is not directly reported whether a region contains on or off target bases, this information can be gained from the coverage quality: all bases off target feature an even number characterizing the coverage quality; all bases on target feature an uneven number characterizing the coverage quality.

For each summed up region the gene that is most likely to be targeted by the original experiment gets reported using biomaRt. If no gene can be found, "NA" is saved for the corresponding region. If not all bases in the summed up region cover a gene, the gene gets reported for the whole region nonetheless.

The output file is saved as: "BadCoverageSummary\code{threshold1};\code{percentage1};\code{threshold2};\code{percentage2}.txt". The output file may be visualized using \code{plotSummary}.
}
\value{
A GRanges object is returned. It represents a summary of the those adjacent regions that feature the same base quality. In the metadata columns the coverage quality of the region, the name and the geneID of the gene that is located in the corresponding region is saved.
}
\references{
More information on the R/Bioconductor package 'biomaRt' may be found at:

http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html
}
\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
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
badCoverageSummary <- reportBadRegionsSummary(threshold1, threshold2, percentage1,
                                              percentage2, coverage_indicators,
                                              "", output)
}
}

