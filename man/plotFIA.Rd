\name{plotFIA}
\alias{plotFIA}
\title{Static and animated plots of spatial FIA summaries}
\description{
Produce static and temporally animated plots of spatial FIA summaries. Requires that summaries be returned as spatial objects (e.g. specify \code{returnSpatial = TRUE} when computing summaries using \code{\link{tpa}}). For animated plots, also requires that multiple reporting years be present in the summary data (animations iterate through time). Specify \code{savePath} and \code{fileName} to save plots (animations saved as .gif files).
}

\usage{
plotFIA(data, fillVar, animate = FALSE, title = NULL, colOption = 'viridis',
        lineCol = "gray30", lineWidth =1, minYear = 2005, direction = 1,
        alpha = .9, transform = "identity", text.size = 1, text.font = '',
        lab.width = 1, legend.height = 1, legend.width = 1, device = NULL,
        savePath = NULL, fileName = NULL, ...)
           }

\arguments{
  \item{data}{sf object; spatial FIA summary produced from other \code{rFIA} functions (e.g. \code{\link{tpa}}, \code{\link{biomass}}, etc.).}

  \item{fillVar}{character; variable contained in \code{data} which will be used to fill polygons. Variable name must be quoted.}

  \item{animate}{logical; if TRUE, produces temporally animated plots (iterates over years).}

  \item{title}{character; plot title.}

  \item{colOption}{character; one of: "viridis" (default), "magma", "inferno", "plasma", or "cividis".}

  \item{lineCol}{character; polygon outline color.}

  \item{lineWidth}{numeric; scalar for polygon outline width. Specify \code{lineWidth = 0} for no outline.}

  \item{minYear}{numeric; earliest year to be included in animation. FIA data is sparse in years prior to 2005 and estimates are unlikely to be available.}

  \item{direction}{numeric; sets the order of colors in the scale. If 1, the default, colors are ordered from darkest to lightest. If -1, the order of colors is reversed.}

  \item{alpha}{numeric; alpha transparency, a number in [0,1], see argument alpha in \code{hsv}.}

  \item{transform}{character; transformations to apply to plotted variable \code{fillVar}. Options include: "asn", "atanh", "boxcox", "exp", "identity", "log", "log10", "log1p", "log2", "logit", "reciprocal", "reverse", "sqrt".}

  \item{text.size}{numeric; scalar for text size (e.g. text.size = 2 would be twice the default size).}

  \item{text.font}{character; font family. Choose from: 'Short', 'Canonical', 'mono', 'Courier', 'sans', 'Helvetica', 'serif', 'Times', 'AvantGarde', 'Bookman', 'Helvetica-Narrow', 'NewCenturySchoolbook', 'Palatino', 'URWGothic',  'URWBookman', 'NimbusMon', 'URWHelvetica', 'NimbusSan', 'NimbusSanCond', 'CenturySch', 'URWPalladio', 'URWTimes', or 'NimbusRom'.}

  \item{lab.width}{numeric; scalar for legend title width. This value controls text wrapping in title.}

  \item{legend.height}{numeric; scalar for legend height.}

  \item{legend.width}{numeric; scalar for legend width.}

  \item{savePath}{character; path to save plot to (combined with fileName).}

  \item{fileName}{character; file name to create on disk.}

  \item{...}{other arguments passed on to \code{\link{scale_fill_viridis_c}}.}
}


\author{
Hunter Stanke and Andrew Finley
}

\examples{
\dontrun{
## Compute abundance estimates for live stems in Rhode Island
## for all available inventory years, summarized by counties
tpaRI <- tpa(fiaRI, polys = countiesRI, returnSpatial = TRUE)

## Produce animated plot
plotFIA(tpaRI, fillVar = 'TPA', animate = TRUE, title = 'Abundance (TPA)')

## With a square root transform
plotFIA(tpaRI, fillVar = 'TPA', animate = TRUE, title = 'Abundance (TPA)', transform = 'sqrt')


## Same as above, but for static plots (most recent subset from RI)
tpaMR <- tpa(clipFIA(fiaRI), polys = countiesRI, returnSpatial = TRUE)

## Produce animated plot
plotFIA(tpaMR, fillVar = 'TPA', animate = FALSE, title = 'Abundance (TPA)')

## With a square root transform
plotFIA(tpaMR, fillVar = 'TPA', animate = FALSE, title = 'Abundance (TPA)', transform = 'sqrt')
}
}

