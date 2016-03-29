\name{determineCoverage}
\alias{determineCoverage}
\title{
Determines the coverage (recommended for whole-genome analyses)
}
\description{
\bold{BadRegionFinder} performs a coverage analysis of various samples at a time. The first, essential step of the analysis pipeline -- the coverage determination -- is performed by the function \code{determineCoverage}. Thereby, the whole genome is scanned and wherever a covered base is registered or an originally targeted base is detected, detailed information concerning this position is written out.
}
\usage{
determineCoverage(samples, bam_input, targetRegions, output, TRonly)
}
\arguments{
  \item{samples}{
Data frame object containing the names of the samples to be analyzed (in one column).
}
  \item{bam_input}{
Folder containing the alignment data in bam- and bai format or BamFileList.
}
  \item{targetRegions}{
Data frame- or GRanges object containing the target regions to be analyzed (chromosome: first column, start position: second column and end position: third column).
}
  \item{output}{
The folder to write the output files into. If an empty string is provided, no files are written out.
}
  \item{TRonly}{
Boolean, indicating whether the coverage of the whole genome should be analyzed and reported (\code{FALSE}) or the coverage of the target regions only (\code{TRUE}).
}
}
\details{
The coverage which is determined by the function \code{determineCoverage} contains different steps:

For every sample that is defined in \code{samples}, the coverage is determined using the function \code{coverage} ("Determine Coverage"). To combine information on the coverage with information on whether a set of bases were originally targeted by some sequencing experiment, the \code{targetRegions} get processed ("Determine target bases"). Finally, the information gets combined ("Combine information"): Those positions where no sample shows any coverage and no target base is registered, are summed up. All other positions are reported basewise.

Files get written out in the form: "Summary_chr<chromosomename>.txt".

As sequencing does often not mean whole-genome- or whole-exome sequencing, but targeted sequencing, the function \code{determienCoverage} contains a switch: \code{TRonly}. In case misaligned reads in a targeted sequencing experiment shall be analyzed, it is advisable to set \code{TRonly} to \code{FALSE}. Yet, if only the coverage of the targeted regions are of interest, it is advisable to set \code{TRonly} to \code{TRUE}.
}
\value{
A GRangesList is returned. Every GRanges object contains the coverage information of one chromosome. The metadata columns contain information on the concrete coverage of each sample at a specific position. Furthermore, the column 'TargetBases' contains information on whether the considered region or position contains target bases (value 1) or not (value 0). A region cannot contain both as two regions would be defined in this case. 

If a chromosome is not covered and was not targeted as well, the GRanges object solely contains a single line, considering a whole chromosome if \code{TRonly=FALSE}. If \code{TRonly=TRUE} the starting and end position of the corresponding chromosome is set to zero.
}
\author{
Sarah Sandmann <sarah.sandmann@uni-muenster.de>
}
\seealso{
\code{\link{BadRegionFinder}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
}
\examples{
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
}
}