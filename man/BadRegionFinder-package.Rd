\name{BadRegionFinder-package}
\alias{BadRegionFinder-package}
\alias{BadRegionFinder}
\alias{Coverage Classifier}
\alias{Bad Coverage}
\docType{package}
\title{
\packageTitle{BadRegionFinder}
}
\description{
\packageDescription{BadRegionFinder}
}
\details{

\packageDESCRIPTION{BadRegionFinder}
In the use case of targeted sequencing it is most important to design the set of used primers in a way that the targeted regions are sequenced with a sufficient coverage. Yet, due to e.g. high GC-content the aimed at coverage may not always be obtained. Thus, a tool performing a detailed coverage analysis comparing many samples at a time -- and not considering all available samples individually -- appears to be most useful. Furthermore, with regards to reads mapping off target, it seems helpful to have a tool for investigating those regions, which show a relatively high coverage, but which were not originally targeted.

BadRegionFinder is a package for classifying a selection of regions or the whole genome into the user-definable categories of bad, acceptable and good coverage in any sequence alignment data available as bam files. Various visual and textual types of output are available including detailed output files considering every base that is or should be covered and an overview file considering the coverage of the different genes that were targeted.

\packageIndices{BadRegionFinder}
The package contains a function performing the coverage determination - \code{determineCoverage} (switch for whole-genome- and target region analyses). The actual classification of the coverage is performed by the function \code{determineCoverageQuality}. If any subsets of regions are of interest, these may be selected by the function \code{determineRegionsOfInterest}.

There are three different forms of textual reports available: a summary variant (\code{reportBadRegionsSummary}), a detailed variant (\code{reportBadRegionsDetailed}) and a summary variant focussing on the coverage of the genes (\code{reportBadRegionsGenes}).

Furthermore, there exist three different forms of visual reports: a summary variant (\code{plotSummary}), a detailed variant (\code{plotDetailed}) and a summary variant visualizing the coverage of the genes as a barplot (\code{plotSummaryGenes}).

Additionally, BadRegionFinder may be used to determine user-definable, basewise quantiles over all samples at any position (\code{determineQuantiles}).
}
\author{
\packageAuthor{BadRegionFinder}

Maintainer: \packageMaintainer{BadRegionFinder}
}
\references{
More information on the bam format can be found at: \url{http://samtools.github.io/hts-specs/SAMv1.pdf}
}
\keyword{ package }
\seealso{
\code{\link{determineCoverage}}, \code{\link{determineCoverageQuality}}, \code{\link{determineRegionsOfInterest}}, \code{\link{reportBadRegionsSummary}}, \code{\link{reportBadRegionsDetailed}}, \code{\link{reportBadRegionsGenes}}, \code{\link{plotSummary}}, \code{\link{plotDetailed}}, \code{\link{plotSummaryGenes}}, \code{\link{determineQuantiles}}
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

coverage_summary <- determineCoverage(samples, bam_input, targetRegions, 
                                      output, TRonly = FALSE)
coverage_indicators <- determineCoverageQuality(threshold1, threshold2,
                                                percentage1, percentage2,
                                                coverage_summary)
badCoverageSummary <- reportBadRegionsSummary(threshold1, threshold2, 
                                              percentage1, percentage2,
                                              coverage_indicators, "", output)
coverage_indicators_temp <- reportBadRegionsDetailed(threshold1, threshold2,
                                                     percentage1, percentage2,
                                                     coverage_indicators, "",
                                                     samples, output)
badCoverageOverview <- reportBadRegionsGenes(threshold1, threshold2, percentage1,
                                            percentage2, badCoverageSummary,
                                            output)

plotSummary(threshold1, threshold2, percentage1, percentage2,
            badCoverageSummary, output)
plotDetailed(threshold1, threshold2, percentage1, percentage2,
             coverage_indicators_temp, output)
plotSummaryGenes(threshold1, threshold2, percentage1, percentage2,
                 badCoverageOverview, output)

quantiles <- c(0.5)
coverage_summary2 <- determineQuantiles(coverage_summary, quantiles, output)
}
}
