\name{dwm}
\alias{dwm}

\title{Estimate volume, biomass, and carbon stocks of down woody material (fuels) from FIADB
}

\description{
Produces estimates of down woody material stocks on a per acre basis from the Forest Inventory and Analysis Database (FIADB), along with population totals for each variable. Estimates are returned by fuel class (1HR, 10HR, 100HR, 1000HR, Piles) for application in fuels management. Estimates can be produced for regions defined within the FIA Database (e.g. counties), at the plot level, or within user-defined areal units. If multiple reporting years (EVALIDs) are included in the data, estimates will be output as a time series. If multiple states are represented by the data, estimates will be output for the full region (all area combined), unless specified otherwise (e.g. \code{grpBy = STATECD}). Easy options to implement parallel processing.
}

\usage{
dwm(db, grpBy = NULL, polys = NULL, returnSpatial = FALSE, landType = 'forest',
    areaDomain = NULL, byPlot = FALSE, totals = FALSE, tidy = TRUE,
    SE = TRUE, progress = TRUE, nCores = 1)
}

\arguments{
  \item{db}{\code{FIA.Database} object produced from \code{\link{readFIA}}. Function requires that PLOT, COND_DWM_CALC, COND, POP_PLOT_STRATUM_ASSGN, POP_ESTN_UNIT, POP_EVAL, POP_STRATUM, POP_EVAL_TYP, POP_EVAL_GRP tables exist in FIA.Database object.}

  \item{grpBy}{variables from PLOT, COND, or TREE tables to group estimates by (NOT quoted). Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(ECOSUBCD, OWNGRPCD)}.}

  \item{polys}{\code{sp} or \code{sf} Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit.}

  \item{returnSpatial}{logical; if TRUE, return sf spatial object (polys must also be specified).}

  \item{landType}{character ('forest' or 'timber'); Type of land which estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{areaDomain}{Logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (e.g. within 1 mile of improved road: \code{RDDISTCD \%in\% c(1:6)}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations (population totals not returned).}

  \item{totals}{logical; if TRUE, returns population estimates in addition to ratios.}

  \item{tidy}{logical; if TRUE, returns estimates grouped by fuel type, rather than including individual columns for each fuel type (For use in tidyverse packages, e.g. \code{ggplot2}, \code{\link{dplyr}}). Not recommended when returning spatial objects (\code{returnSpatial = TRUE}), for consistency with shapefile data structures. }

  \item{SE}{logical; if TRUE, returns estimates with samping error (approx. 5x faster without returning samping errors)}

  \item{progress}{logical; if TRUE, a dynamic progress bar is displayed (approx. 30\% faster when FALSE.}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}
}
\details{
Estimation of attributes follows the procedures documented in Bechtold and Patterson (2005). Specifically, per acre estimates are computed using a sample-based ratio-of-means estimator of total volume (biomass or carbon) / total land area within the domain of interest.

As defined by FIA, down woody material includes dead organic materials (resulting from plant mortality and leaf turnover) and fuel complexes of live shrubs and herbs. To maintain relevance for forest fuels management, we report estimates grouped by fuel lag-time classes. Specifically, we report estimates for 1HR fuels (small, fine woody debris), 10HR fuels (medium, fine woody debris), 100HR fuels (large, fine woody debris), 1000HR fuels (coarse woody debris), and slash piles. See Woodall and Monleon (2007) for definitions of fuel lag-time classes and for details on sampling and estimation procedures.

Stratified random sampling techniques are most often employed to compute estimates in recent inventories, although double sampling and simple random sampling may be employed for early inventories. Estimates are adjusted for non-response bias by assuming attributes of non-response plot locations to be equal to the mean of other plots included within thier respective stratum or population.

Forest land must be at least 10-percent stocked by trees of any size, including land that formerly had such tree cover and that will be naturally or artificially regenerated. Forest land includes transition zones, such as areas between heavily forested and nonforested lands that are at least 10-percent stocked with trees and forest areas adjacent to urban and builtup lands. The minimum area for classification of forest land is 1 acre and 120 feet wide measured stem-to-stem from the outer-most edge. Unimproved roads and trails, streams, and clearings in forest areas are classified as forest if less than 120 feet wide. Timber land is a subset of forest land that is producing or is capable of producing crops of industrial wood and not withdrawn from timber utilization by statute or administrative regulation. (Note: Areas qualifying as timberland are capable of producing at least 20 cubic feet per acre per year of industrial wood in natural stands. Currently inaccessible and inoperable areas are NOT included).

Easy, efficient parallelization is implemented with the \code{\link{parallel}} package. Users must only specify the \code{nCores} argument with a value greater than 1 in order to implement parallel processing on their machines. Parallel implementation is achieved using a snow type cluster on any Windows OS, and with multicore forking on any Unix OS (Linux, Mac). Implementing parallel processing may substantially decrease free memory during processing, particularly on Windows OS. Thus, users should be cautious when running in parallel, and consider implementing serial processing for this task if computational resources are limited (\code{nCores = 1}).
}

\note{All sampling error estimates are returned as percentages, and represent ~68\% confidence (1 standard deviation). To compute sampling error percent with 95\% confidence, multiply by 1.96.}

\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, values are returned for each plot. All variables with names ending in \code{SE}, represent the estimate of sampling error (\%) of the variable. All variables with names ending in \code{TOTAL} represent the population total of the variable.

\itemize{
      \item{\strong{YEAR}: reporting year associated with estimates}
      \item{\strong{VOL_ACRE}: estimate of mean volume per acre of dwm (cu.ft/acre)}
      \item{\strong{BIO_ACRE}: estimate of mean biomass per acre of dwm (tons/acre)}
      \item{\strong{CARB_ACRE}: estimate of mean carbon mass per acre of dwm (tons/acre)}
      \item{\strong{nPlots}: number of non-zero plots used to compute estimates}
      }
}

\author{
Hunter Stanke and Andrew Finley
}

\references{
FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}

Woodall, C.; Monleon, V.J., eds. 2007. Sampling Protocol, Estimation, and Analysis Procedures for the Down Woody Materials Indicator of the FIA Program. Gen. Tech. Rep. NRS - 22. ewtown Square, PA: U.S. Department of Agriculture, Forest Service, Northern Research Station. \url{https://www.nrs.fs.fed.us/pubs/gtr/gtr_nrs22.pdf}
}

\seealso{
\code{\link{tpa}}, \code{\link{biomass}}
}


\examples{
\dontrun{

## Load data from rFIA package
library(rFIA)
data(fiaRI)
data(countiesRI)

## Estimates of all forestland, over time
dwm(fiaRI)

## Same as above, but output contains seperate column for each structural stage,
##   rathern than grouping variable
dwm(fiaRI,
    tidy = FALSE)

## Same as above, but implemented in parallel (much quicker)
detectCores() # 4 cores available, we will take 3
dwm(fiaRI,
    tidy = FALSE,
    nCores = 3)

## Estimates of all forestland by owner group (most recent subset)
dwm(clipFIA(fiaRI, mostRecent = TRUE),
    grpBy = OWNGRPCD)

## Estimates of all forestland on mesic sites (most recent)
dwm(clipFIA(fiaRI, mostRecent = TRUE),
    areaDomain = PHYSCLCD \%in\% 21:29)

## Estimates of all forestland by county and return
#   return spatial object
dwmSF <- dwm(clipFIA(fiaRI, mostRecent = TRUE),
            polys = countiesRI,
            returnSpatial = TRUE,
            tidy = FALSE)
plot(dwmSF)
plotFIA(dwmSF, BIO_ACRE) # TOTAL BIOMASS / ACRE (tons)
}
}

