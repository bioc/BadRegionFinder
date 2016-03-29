\name{determineQuantiles}
\alias{determineQuantiles}
\alias{determine Quantiles}
\title{
Determines basewise user-defined quantiles
}
\description{
The function \code{determineQuantiles} provides a possibility to determine user-definable quantiles for every base previously analyzed. Thereby, the quantiles are determined over all samples.
}
\usage{
determineQuantiles(coverage_summary, quantiles, output)
}
\arguments{
  \item{coverage_summary}{
List object, return value of function \code{determineCoverage}.
}
  \item{quantiles}{
Vector determining the quantiles to be calculated.
}
  \item{output}{
The folder to write the output files into.
}
}
\details{
The function \code{determineQuantiles} serves to determine a set of user-defined quantiles at each position over all samples. Every single base is analyzed, except for the case when the bases were not originally targeted and if no coverage is detected by any of the samples. In this case the corresponding region is summed up.

Files get written out in the form: "Quantiles_chr<chromosomename>.txt".
}
\value{
A list is returned. Every component contains the coverage information of one chromosome as a GRanges object. The metadata columns contain information on the coverage according to the previously defined quantiles. Furthermore, the column 'TargetBases' contains information on whether the considered region or position contains target bases (value 1) or not (value 0).

If a chromosome has not been targeted and/or not covered by any sample, but defined in \code{regionsOfInterest}, the component is "NA".
}

\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}

\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}
}
\examples{
library("BSgenome.Hsapiens.UCSC.hg19")
\donttest{
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
quantiles <- c(0.25, 0.5, 0.75)
coverage_summary2 <- determineQuantiles(coverage_summary, quantiles, output)
}
}