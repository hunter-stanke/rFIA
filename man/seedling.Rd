\name{seedling}
\alias{seedling}

\title{Estimate seedling abundance per acre from FIADB}

\description{
Produces seedling (< 1 inch DBH) tree per acre (TPA) estimates from FIA data, along with population totals. Estimates can be produced for regions defined within the FIA Database (e.g. counties), at the plot level, or within user-defined areal units. Options to group estimates by species and other variables defined in the FIADB. If multiple reporting years (EVALIDs) are included in the data, estimates will be output as a time series. If multiple states are represented by the data, estimates will be output for the full region (all area combined), unless specified otherwise (e.g. \code{grpBy = STATECD}). Easy options to implement parallel processing.
}


\usage{
seedling(db, grpBy = NULL, polys = NULL, returnSpatial = FALSE,
         bySpecies = FALSE, landType = "forest", treeDomain = NULL,
         areaDomain = NULL, totals = FALSE, byPlot = FALSE,
         SE = TRUE, nCores = 1)
}

\arguments{
  \item{db}{\code{FIA.Database} object produced from \code{\link{readFIA}}. Function requires that PLOT, TREE, COND, POP_PLOT_STRATUM_ASSGN, POP_ESTN_UNIT, POP_EVAL, POP_STRATUM, POP_EVAL_TYP, POP_EVAL_GRP tables exist in FIA.Database object.}

\item{grpBy}{variables from PLOT, COND, or SEEDLING tables to group estimates by (NOT quoted). Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(ECOSUBCD, OWNGRPCD)}.}

  \item{polys}{\code{sp} or \code{sf} Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit.}

  \item{returnSpatial}{logical; if TRUE, return sf spatial object (polys must also be specified).}

  \item{bySpecies}{logical; if TRUE, returns estimates grouped by species.}


  \item{landType}{character ('forest' or 'timber'); Type of land which estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{treeDomain}{logical predicates defined in terms of the variables in PLOT, SEEDLING, and/or COND tables. Used to define the type of trees for which estimates will be produced (e.g. white pine: \code{SPCD == 129}). Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only trees where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{areaDomain}{logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (e.g. within 1 mile of improved road: \code{RDDISTCD \%in\% c(1:6)}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{totals}{logical; if TRUE, return population estimates (e.g. total area, total trees) along with ratio estimates (e.g. mean trees per acre).}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations (population totals not returned).}

  \item{SE}{logical; if TRUE, returns estimates with samping error (approx. 5x faster without returning samping errors)}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}
}


\details{
Estimation of attributes follows the procedures documented in Bechtold and Patterson (2005). Specifically, TPA is computed using a sample-based ratio-of-means estimator of total seedlings / total land area within the domain of interest. Percentages of live TPA in the domain of interest are represented as the total number of trees of a particular type (e.g., white pine) / total number of trees (live, all species) within the region. The total populations used to compute these percentages will not change by changing treeType, but will vary if the user specifies an areaDomain or treeDomain.

Stratified random sampling techniques are most often employed to compute estimates in recent inventories, although double sampling and simple random sampling may be employed for early inventories. Estimates are adjusted for non-response bias by assuming attributes of non-response plot locations to be equal to the mean of other plots included within thier respective stratum or population.

Forest land must be at least 10-percent stocked by trees of any size, including land that formerly had such tree cover and that will be naturally or artificially regenerated. Forest land includes transition zones, such as areas between heavily forested and nonforested lands that are at least 10-percent stocked with trees and forest areas adjacent to urban and builtup lands. The minimum area for classification of forest land is 1 acre and 120 feet wide measured stem-to-stem from the outer-most edge. Unimproved roads and trails, streams, and clearings in forest areas are classified as forest if less than 120 feet wide. Timber land is a subset of forest land that is producing or is capable of producing crops of industrial wood and not withdrawn from timber utilization by statute or administrative regulation. (Note: Areas qualifying as timberland are capable of producing at least 20 cubic feet per acre per year of industrial wood in natural stands. Currently inaccessible and inoperable areas are NOT included).

Easy, efficient parallelization is implemented with the \code{\link{parallel}} package. Users must only specify the \code{nCores} argument with a value greater than 1 in order to implement parallel processing on their machines. Parallel implementation is achieved using a snow type cluster on any Windows OS, and with multicore forking on any Unix OS (Linux, Mac). Implementing parallel processing may substantially decrease free memory during processing, particularly on Windows OS. Thus, users should be cautious when running in parallel, and consider implementing serial processing for this task if computational resources are limited (\code{nCores = 1}).
}

\note{All sampling error estimates are returned as percentages, and represent ~68\% confidence (1 standard deviation). To compute sampling error percent with 95\% confidence, multiply by 1.96.}

\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, values are returned for each plot. All variables with names ending in \code{SE}, represent the estimate of sampling error (\%) of the variable. All variables with names ending in \code{TOTAL} represent the population total of the variable.

\itemize{
      \item{\strong{YEAR}: reporting year associated with estimates}
      \item{\strong{TPA}: estimate of mean trees per acre}
      \item{\strong{TPA_PERC}: estimate of mean proportion of live trees falling within the domain of interest, with respect to trees per acre}
      \item{\strong{nPlots_SEEDLING}: number of non-zero plots used to compute tpa estimates}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tpa}}, \code{\link{growMort}}, \code{\link{biomass}}
}

\examples{
## Load data from the rFIA package
data(fiaRI)
data(countiesRI)

## Most recents subset
fiaRI_mr <- clipFIA(fiaRI)


## Most recent estimates on timber land by species
seedling(db = fiaRI_mr,
    landType = 'timber')
\donttest{

## Estimates for white pine on forested mesic sites (all available inventories)
seedling(fiaRI_mr,
    treeDomain = SPCD == 129, # Species code for white pine
    areaDomain = PHYSCLCD \%in\% 21:29) # Mesic Physiographic classes

## Most recent estimates grouped by stand age on forest land
# Make a categorical variable which represents stand age (grouped by 10 yr intervals)
fiaRI_mr$COND$STAND_AGE <- makeClasses(fiaRI_mr$COND$STDAGE, interval = 10)
seedling(db = fiaRI_mr,
    grpBy = STAND_AGE)

## Most recent estimates for live stems on forest land by species
seedling(db = fiaRI_mr,
    landType = 'forest',
    bySpecies = TRUE)

## Same as above, but implemented in parallel (much quicker)
parallel::detectCores(logical = FALSE) # 4 cores available, we will take 2
seedling(db = fiaRI_mr,
    landType = 'forest',
    bySpecies = TRUE,
    nCores = 2)


## Most recent estimates for all stems on forest land grouped by user-defined areal units
ctSF <- seedling(fiaRI_mr,
            polys = countiesRI,
            returnSpatial = TRUE)
plot(ctSF) # Plot multiple variables simultaneously
plotFIA(ctSF, TPA) # Plot of TPA with color scale
}
}