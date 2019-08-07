\name{invasive}
\alias{invasive}

\title{Estimate invasive species coverage from FIADB
}

\description{
Produces estimates of areal coverage of invasive species from the Forest Inventory and Analysis Database. Estimates can be produced for regions defined within the FIA Database (ex. counties), at the plot level, or within user-defined areal units. All estimates are returned by species although can be grouped by other variables defined in the FIADB. If multiple reporting years (EVALIDs) are included in the data, estimates will be output as a time series. If multiple states are represented by the data, estimates will be output for the full region (all area combined), unless specified otherwise (ex. \code{grpBy = 'STATECD'}). Easy options to implement parallel processing.}

\usage{
invasive(db, grpBy = NULL, polys = NULL, returnSpatial = FALSE, landType = "forest",
         areaDomain = NULL, byPlot = FALSE, totals = FALSE, SE = TRUE,
         progress = TRUE, nCores = 1)
}

\arguments{
  \item{db}{\code{FIA.Database} object produced from \code{\link{readFIA}}. Function requires that PLOT, INVASIVE_SUBPLOT_SPP, COND, POP_PLOT_STRATUM_ASSGN, POP_ESTN_UNIT, POP_EVAL, POP_STRATUM, POP_EVAL_TYP, POP_EVAL_GRP tables exist in FIA.Database object.}

 \item{grpBy}{character; variables from PLOT, COND, or TREE tables to group estimates by. Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c('ECOSUBCD', 'OWNGRPCD')}.}

  \item{polys}{sp or sf Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit.}

  \item{returnSpatial}{logical; if TRUE, return sf spatial object (polys must also be specified).}

  \item{landType}{character ('forest' or 'timber'); Type of land which estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{areaDomain}{Logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (ex. within 1 mile of improved road: \code{RDDISTCD \%in\% c(1:6)}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations (population totals not returned).}

  \item{totals}{logical; if TRUE, returns population estimates in addition to \% coverage.}

  \item{SE}{logical; if TRUE, returns estimates with samping error (approx. 5x faster without returning samping errors)}

  \item{progress}{logical; if TRUE, a dynamic progress bar is displayed. Approximately 30\% faster when FALSE.}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}
}

\details{
Estimation of attributes follows the procedures documented in Bechtold & Patterson (2005). Specifically, percent areal coverage is computed using a sample-based ratio-of-means estimator of total invasive coverage area / total land area within the domain of interest. Estimates of areal coverage of individual invasive species should NOT be summed to produce estimates of areal coverage by ALL invasive species, as areal coverage by species is not mutually exclusive (multiple species my occur in the same area). Current FIA data collection protocols do not allow for the unbiased estimation of areal coverage by all invasive species.

Stratified random sampling techniques are most often employed to compute estimates in recent inventories, although double sampling and simple random sampling may be employed for early inventories. Estimates are adjusted for non-response bias by assuming attributes of non-response plot locations to be equal to the mean of other plots included within thier respective stratum or population.

Forest land must be at least 10-percent stocked by trees of any size, including land that formerly had such tree cover and that will be naturally or artificially regenerated. Forest land includes transition zones, such as areas between heavily forested and nonforested lands that are at least 10-percent stocked with trees and forest areas adjacent to urban and builtup lands. The minimum area for classification of forest land is 1 acre and 120 feet wide measured stem-to-stem from the outer-most edge. Unimproved roads and trails, streams, and clearings in forest areas are classified as forest if less than 120 feet wide. Timber land is a subset of forest land that is producing or is capable of producing crops of industrial wood and not withdrawn from timber utilization by statute or administrative regulation. (Note: Areas qualifying as timberland are capable of producing at least 20 cubic feet per acre per year of industrial wood in natural stands. Currently inaccessible and inoperable areas are NOT included).

Easy, efficient parallelization is implemented with the \code{\link{parallel}} package. Users must only specify the \code{nCores} argument with a value greater than 1 in order to implement parallel processing on their machines. Parallel implementation is achieved using a snow type cluster on any Windows OS, and with multicore forking on any Unix OS (Linux, Mac). Implementing parallel processing may substantially decrease free memory during processing, particularly on Windows OS. Thus, users should be cautious when running in parallel, and consider implementing serial processing for this task (\code{nCores = 1}).
}

\note{All sampling error estimates are returned as percentages, and represent ~68\% confidence (1 standard deviation). To compute sampling error percent with 95\% confidence, multiply by 1.96.}

\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, values of areal coverage are returned for each plot. All variables with names ending in \code{SE}, represent the estimate of sampling error (\%) of the variable.

\itemize{
      \item{\strong{YEAR}: reporting year associated with estimates}
      \item{\strong{SYMBOL}: unique species ID from NRCS Plant Reference Guide}
      \item{\strong{SCIENTIFIC_NAME}: scientific name of the species}
      \item{\strong{COMMON_NAME}: common name of the species}
      \item{\strong{COVER_PCT}: estimate of percent areal coverage of the species }
      \item{\strong{COVER_AREA}: estimate of areal coverage of the species (acres) }
      \item{\strong{AREA}: estimate of total land area (acres) }
      \item{\strong{nPlots_INV}: number of non-zero plots used to compute invasive coverage estimates}
      \item{\strong{nPlots_AREA}: number of non-zero plots used to compute land area estimates}
      }
}

\author{
Hunter Stanke & Andrew Finley, Michigan State University
}

\references{
FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}
}

\seealso{
\code{\link{dwm}}, \code{\link{tpa}}, \code{\link{clipFIA}}
}

\examples{
\dontrun{

## Load data from rFIADB package
library(rFIA)
data(fiaRI)
data(countiesRI)


## Most recent estimates on forest land
invasive(db = clipFIA(fiaRI),
    landType = 'forest')

## Same as above, but implemented in parallel (much quicker)
detectCores() # 4 cores available, we will take 3
invasive(db = clipFIA(fiaRI),
         landType = 'forest',
         nCores = 3)

## Same as above, but return population totals as well (total areal coverage, total forest area)
invasive(db = clipFIA(fiaRI, mostRecent = TRUE),
    landType = 'forest',
    totals = TRUE)

## Estimates for all available inventories (time-series) on forestland,
##  with no sampling error estimates (~5x faster)
invasive(db = fiaRI, SE = FALSE)

## Most recent estimates grouped by stand age on forest land
# Make a categorical variable which represents stand age (grouped by 10 yr intervals)
fiaRI$COND$STAND_AGE <- makeClasses(fiaRI$COND$STDAGE, interval = 10)
invasive(db = clipFIA(fiaRI, mostRecent = TRUE),
         grpBy = 'STAND_AGE')

## Estimates on forested mesic sites (all available inventories)
invasive(fiaRI,
         areaDomain = PHYSCLCD \%in\% c(21:29)) # Mesic Physiographic classes

## Most recent estimates on forest land grouped by user-defined areal units
ctSF <- invasive(clipFIA(fiaRI, mostRecent = TRUE),
                 polys = countiesRI,
                 returnSpatial = TRUE)
plot(ctSF) # Plot multiple variables simultaneously
plotFIA(ctSF[ctSF$SYMBOL == 'ROMU',], 'COVER_PCT') # Plot Multiflora rose coverage
}
}

