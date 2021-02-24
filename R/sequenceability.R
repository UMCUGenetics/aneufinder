

#' Procedure to determine the sequenceability factors
#'
#' @param binpath.uncorrected  The uncorrected binned data from an Aneufinder run.
#' @param bins                 The output of either 'fixedWidthBins' or 'variableWidthBins'.
#'
#' @return A vector containing the sequenceability factors.
#'
#' @export
determineSequenceabilityFactors <- function (binpath.uncorrected, bins) {
    binned.files  <- list.files(binpath.uncorrected, pattern=".*RData$", full.names=TRUE)
    numberOfCells <- length(binned.files)
    numberOfBins  <- length(bins[[1]][[1]])

    totals        <- array(NA, dim = c(numberOfCells, numberOfBins))
    for (cellIndex in 1:numberOfCells)
    {
        file <- binned.files[cellIndex]
        cell <- get(load(file))
        for (binIndex in 1:numberOfBins)
        {
            totals[cellIndex,binIndex] <- cell[[1]]$counts[binIndex]
        }
        rm(cell)
    }

    medianPerBin <- numeric(numberOfBins)
    for (binIndex in 1:numberOfBins)
    {
        medianPerBin[binIndex] <- median(totals[,binIndex])
    }

    averageBinCount <- mean(medianPerBin)
    sequenceability.factors <- averageBinCount / medianPerBin

    ## When there are no reads in a bin, we'd get an infinite factor.
    ## But we want to ignore these cases, and therefore we set a factor
    ## of 1.0 for these bins.
    sequenceability.factors[which(!is.finite(sequenceability.factors))] <- 1.0

    return (sequenceability.factors)
}
