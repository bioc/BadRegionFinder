\name{reportBadRegionsDetailed}
\alias{reportBadRegionsDetailed}
\alias{report Bad Regions Detailed}
\title{
Gives a detailed report on the coverage quality
}
\description{
The function \code{reportBadRegionsDetailed} creates a detailed report containing all regions of interest (basewise), the coverage of each sample at the corresponding positions, the indicator whether the bases were originally targeted, their coverage quality and the corresponding gene (name and geneID).
}
\usage{
reportBadRegionsDetailed(threshold1, threshold2, percentage1, percentage2, 
                         coverage_indicators, mart, samples, output)
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
  \item{samples}{
Data frame object containing the names of the samples to be analyzed (in one column).
}
  \item{output}{
The folder to write the output files into. If \code{output} is just an empty string, no output file is written out.
}
}
\details{
To gain more detailed information of the coverage quality, a file for every chromosome to be analyzed may be created by the function \code{reportBadRegionsDetailed}. The function may either take information on the whole genome (output from \code{determineCoverage} with \code{TRonly=FALSE}, processed using \code{determineCoverageQuality}) as an input, or information on the target regions (output from \code{determineCoverage} with \code{TRonly=TRUE}, processed using \code{determineCoverageQuality}), or information on a selection of regions of interest (output from \code{determineRegionsOfInterest}). 

Different from the summed-up variant \code{reportBadRegionsSummary}, information on every single base of interest gets reported (except for completely uncovered and untargeted regions, which are summed up). For every base its position, the coverage of each sample, information on whether this base was originally targeted (value 1) or not (value 0), the coverage quality and the most likely gene (name and geneID) that was targeted by the original experiment get reported. Information on the gene names and the geneIDs results from biomaRt. If no gene can be found at a position, "NA" is reported for the corresponding base.

The output files are saved as: "BadCoverageChromosome<chromosomename>;\code{threshold1};\code{percentage1};\code{threshold2};\code{percentage2}.txt". The output file may be visualized using \code{plotDetailed}.
}
\value{
A list is returned. Every component contains the coverage information of one chromosome as a GRanges object. The metadata columns contain information on the concrete coverage of each sample at a specific position. Furthermore, the column 'TargetBases' contains information on whether the considered region or position contains target bases (value 1) or not (value 0). The column 'indicator' contains information on the coverage quality of the corresponding region/position (0: bad region off target; 1: bad region on target; 2: acceptable region off target; 3: acceptable region on target; 4: good region off target; 5: good region on target). Furthermore, the name and the geneID of the gene that is located at the corresponding position is saved.

If a chromosome is not covered and was not targeted as well, the component is "NA".
}
\references{
More information on the R/Bioconductor package 'biomaRt' may be found at:

http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html
}
\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}
\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
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
coverage_indicators_temp <- reportBadRegionsDetailed(threshold1, threshold2,
                                                     percentage1, percentage2,
                                                     coverage_indicators, "",
                                                     samples, output)
}
}

