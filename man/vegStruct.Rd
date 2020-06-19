\name{vegStruct}
\alias{vegStruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate vegation cover by canopy layer with the FIADB}
\description{Produces estimates of vegetation cover by canopy layer and species growth form from the Forest Inventory and Analysis Database. Estimates can be produced for regions defined within the FIA Database (e.g. counties), at the plot level, or within user-defined areal units. All estimates are returned by species although can be grouped by other variables defined in the FIADB. If multiple reporting years (EVALIDs) are included in the data, estimates will be output as a time series. If multiple states are represented by the data, estimates will be output for the full region (all area combined), unless specified otherwise (e.g. \code{grpBy = STATECD}). Easy options to implement parallel processing.}

\usage{
vegStruct(db, grpBy = NULL, polys = NULL,
          returnSpatial = FALSE, landType = "forest",
          method = "TI", lambda = 0.5,
          areaDomain = NULL, byPlot = FALSE,
          totals = FALSE, nCores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{\code{FIA.Database} or \code{Remote.FIA.Database} object produced from \code{\link{readFIA}}. If a \code{Remote.FIA.Database}, data will be read in and processed state-by-state to conserve RAM.}

 \item{grpBy}{variables from PLOT, COND, or TREE tables to group estimates by (NOT quoted). Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(ECOSUBCD, OWNGRPCD)}.}

  \item{polys}{\code{sp} or \code{sf} Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit.}

  \item{returnSpatial}{logical; if TRUE, return sf spatial object (polys must also be specified).}

  \item{landType}{character ('forest' or 'timber'); Type of land which estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{method}{character; Method used for annual panel combination (see details). One of: 'TI' (temporally indifferent), 'annual' (annual panels), 'SMA' (simple moving average), 'LMA' (linear moving average), or 'EMA' (exponential moving average)}

  \item{lambda}{numeric (0,1); if \code{method == 'EMA'}, the decay parameter used to define weighting scheme for annual panels. Low values place higher weight on more recent panels, and vice versa. Specify a vector of values to compute estimates using mulitple wieghting schemes, and use \code{plotFIA} with \code{grp} set to \code{lambda} to produce moving average ribbon plots.}

  \item{areaDomain}{logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (e.g. within 1 mile of improved road: \code{RDDISTCD \%in\% c(1:6)}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations (population totals not returned).}

  \item{totals}{logical; if TRUE, returns population estimates in addition to \% coverage.}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}
}



\details{
Estimation of attributes follows the procedures documented in Bechtold and Patterson (2005). Specifically, percent areal coverage is computed using a sample-based ratio-of-means estimator of total coverage area / total land area within the domain of interest.

Users may specify alternatives to the 'Temporally Indifferent' method for panel combination using the \code{method} argument. Estimators are currently implemented for the temporally indifferent (TI) method (default, used by EVALIDator), simple moving average (SMA), exponential moving average (EMA), and annual estimation (no panel combination). Users may modify the temporal weighting scheme employed by the exponential moving average using the \code{yrs} argument. See Chapter 5 of Bechtold & Patterson (2005) for more details.

Stratified random sampling techniques are most often employed to compute estimates in recent inventories, although double sampling and simple random sampling may be employed for early inventories. Estimates are adjusted for non-response bias by assuming attributes of non-response plot locations to be equal to the mean of other plots included within thier respective stratum or population.

Forest land must be at least 10-percent stocked by trees of any size, including land that formerly had such tree cover and that will be naturally or artificially regenerated. Forest land includes transition zones, such as areas between heavily forested and nonforested lands that are at least 10-percent stocked with trees and forest areas adjacent to urban and builtup lands. The minimum area for classification of forest land is 1 acre and 120 feet wide measured stem-to-stem from the outer-most edge. Unimproved roads and trails, streams, and clearings in forest areas are classified as forest if less than 120 feet wide. Timber land is a subset of forest land that is producing or is capable of producing crops of industrial wood and not withdrawn from timber utilization by statute or administrative regulation. (Note: Areas qualifying as timberland are capable of producing at least 20 cubic feet per acre per year of industrial wood in natural stands. Currently inaccessible and inoperable areas are NOT included).

Easy, efficient parallelization is implemented with the \code{\link{parallel}} package. Users must only specify the \code{nCores} argument with a value greater than 1 in order to implement parallel processing on their machines. Parallel implementation is achieved using a snow type cluster on any Windows OS, and with multicore forking on any Unix OS (Linux, Mac). Implementing parallel processing may substantially decrease free memory during processing, particularly on Windows OS. Thus, users should be cautious when running in parallel, and consider implementing serial processing for this task if computational resources are limited (\code{nCores = 1}).
}


\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, values of areal coverage are returned for each plot (sq. ft.). All variables with names ending in \code{SE}, represent the estimate of sampling error (\%) of the variable.

\itemize{
      \item{\strong{YEAR}: reporting year associated with estimates}
      \item{\strong{LAYER}: canopy layer}
      \item{\strong{GROWTH_HABIT}: species growth habit}
      \item{\strong{COVER_PCT}: estimate of percent areal coverage of the growth habit within the canopy layer}
      \item{\strong{COVER_AREA}: estimate of areal coverage of the growth habit within the canopy layer (acres) }
      \item{\strong{AREA}: estimate of total land area (acres) }
      \item{\strong{nPlots_VEG}: number of non-zero plots used to compute areal coverage estimates}
      \item{\strong{nPlots_AREA}: number of non-zero plots used to compute land area estimates}
      }
}

\author{
Hunter Stanke and Andrew Finley
}

\references{
FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}
}

\note{All sampling error estimates are returned as percentages, and represent ~68\% confidence (1 standard deviation). To compute sampling error percent with 95\% confidence, multiply by 1.96.}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{invasive}}, \code{\link{dwm}}
}


\examples{
## Load data from the rFIA package
data(fiaRI)
data(countiesRI)

## Most recents subset
fiaRI_mr <- clipFIA(fiaRI)
}

