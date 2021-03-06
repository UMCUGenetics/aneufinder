% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hotspotter.R
\name{findHotspots}
\alias{findHotspots}
\title{Find breakpoint hotspots}
\usage{
findHotspots(models, bw, pval = 0.05, spacing.bp = 5000, filename = NULL)
}
\arguments{
\item{models}{A list of \code{\link{GRanges-class}} or \code{\link{aneuHMM}} objects or a character vector with files that contain such objects.}

\item{bw}{Bandwidth used for kernel density estimation (see \code{\link[stats]{density}}).}

\item{pval}{P-value cutoff for hotspots.}

\item{spacing.bp}{Spacing of datapoints for KDE in basepairs.}

\item{filename}{Will write hotspot coordinates and densities to the specified file. Endings "_breakpoint-hotspots.bed.gz" and "_breakpoint-densities.wig.gz" will be appended to \code{filename}.}
}
\value{
A list of \code{\link{GRanges-class}} objects containing 1) coordinates of hotspots and 2) p-values within the hotspot.
}
\description{
Find breakpoint hotspots with kernel density estimation (KDE).
}
\details{
\code{findHotspots} uses \code{\link[stats]{density}} to perform a KDE. A p-value is calculated by comparing the density profile of the genomic events with the density profile of a randomly subsampled set of genomic events. Due to this random sampling, the result can vary for each function call, most likely for hotspots whose p-value is close to the specified \code{pval}.
}
