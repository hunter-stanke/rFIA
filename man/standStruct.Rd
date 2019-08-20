\name{standStruct}
\alias{standStruct}

\title{Estimate forest structural stage distribution from FIADB
}

\description{
Estimates the stand structural stage distribution of an area of forest/ timberland from FIA data. Estimates can be produced for regions defined within the FIA Database (e.g. counties), at the plot level, or within user-defined areal units. If multiple reporting years (EVALIDs) are included in the data, estimates will be output as a time series. Easy options to implement parallel processing. Stand structural stage is classified for each stand (condition) using a method similar to that of Frelich and Lorimer (1991) but substitute basal area for exposed crown area (see Details, References). If multiple states are represented by the data, estimates will be output for the full region (all area combined), unless specified otherwise (e.g. \code{grpBy = STATECD}).
}

\usage{
standStruct(db, grpBy = NULL, polys = NULL, returnSpatial = FALSE,
            landType = 'forest', areaDomain = NULL,
            byPlot = FALSE, totals = FALSE, tidy = TRUE, SE = TRUE,
            progress = TRUE, nCores = 1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{db}{\code{FIA.Database} object produced from \code{\link{readFIA}}. Function requires that PLOT, TREE, COND, POP_PLOT_STRATUM_ASSGN, POP_ESTN_UNIT, POP_EVAL, POP_STRATUM, POP_EVAL_TYP, POP_EVAL_GRP tables exist in FIA.Database object.}

   \item{grpBy}{variables from PLOT, COND, or TREE tables to group estimates by (NOT quoted). Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(ECOSUBCD, OWNGRPCD)}.}

  \item{polys}{\code{sp} or \code{sf} Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit.}

  \item{returnSpatial}{logical; if TRUE, return sf spatial object (polys must also be specified).}

\item{landType}{character ('forest' or 'timber'); Type of land which estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{areaDomain}{logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (e.g. within 1 mile of improved road: \code{RDDISTCD \%in\% c(1:6)}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations (population totals not returned).}

  \item{totals}{logical; if TRUE, returns population estimates in addition to ratios.}

  \item{tidy}{logical; if TRUE, returns estimates grouped by structural stage, rather than including individual columns for each stand structural stage (For use in tidyverse packages, e.g. \code{ggplot2}, \code{\link{dplyr}}). Not recommended when returning spatial objects (\code{returnSpatial = TRUE}), for consistency with shapefile data structures. }

  \item{SE}{logical; if TRUE, returns estimates with samping error (approx. 5x faster without returning samping errors)}

  \item{progress}{logical; if TRUE, a dynamic progress bar is displayed (approx. 30\% faster when FALSE.}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}
}

\details{
Estimation of attributes follows the procedures documented in Bechtold and Patterson (2005). Specifically, the percent land area occupied by forest in each stand structural stage are computed using a sample-based ratio-of-means estimator of total area in structural stage / total land area within the domain of interest. Stand structural stage is classified based on the relative basal area of canopy stems in various size classes (defined below). Only stems which are identified on-site dominant, subdominant, or intermdediate crown-classes are used to classify stand structural stage.

\strong{Diameter Classes}
\itemize{
\item \emph{Pole}: 11 - 25.9 cm
\item \emph{Mature}: 26 - 45.9 cm
\item \emph{Large}: 46+ cm
}

\strong{Structural Stage Classification}
\itemize{
\item \emph{Pole Stage}: > 67\% BA in pole and mature classes, with more BA in pole than mature.
\item \emph{Mature Stage}: > 67\% BA in pole and mature classes, with more BA in mature than pole OR > 67\% BA in mature and large classes, with more BA in mature.
\item \emph{Late-Successional Stage}: > 67\% BA in mature and large classes, with more in large
\item \emph{Mosiac}: Any plot not meeting above criteria.
}

Stratified random sampling techniques are most often employed to compute estimates in recent inventories, although double sampling and simple random sampling may be employed for early inventories. Estimates are adjusted for non-response bias by assuming attributes of non-response plot locations to be equal to the mean of other plots included within thier respective stratum or population.

Forest land must be at least 10-percent stocked by trees of any size, including land that formerly had such tree cover and that will be naturally or artificially regenerated. Forest land includes transition zones, such as areas between heavily forested and nonforested lands that are at least 10-percent stocked with trees and forest areas adjacent to urban and builtup lands. The minimum area for classification of forest land is 1 acre and 120 feet wide measured stem-to-stem from the outer-most edge. Unimproved roads and trails, streams, and clearings in forest areas are classified as forest if less than 120 feet wide. Timber land is a subset of forest land that is producing or is capable of producing crops of industrial wood and not withdrawn from timber utilization by statute or administrative regulation. (Note: Areas qualifying as timberland are capable of producing at least 20 cubic feet per acre per year of industrial wood in natural stands. Currently inaccessible and inoperable areas are NOT included).

Easy, efficient parallelization is implemented with the \code{\link{parallel}} package. Users must only specify the \code{nCores} argument with a value greater than 1 in order to implement parallel processing on their machines. Parallel implementation is achieved using a snow type cluster on any Windows OS, and with multicore forking on any Unix OS (Linux, Mac). Implementing parallel processing may substantially decrease free memory during processing, particularly on Windows OS. Thus, users should be cautious when running in parallel, and consider implementing serial processing for this task if computational resources are limited (\code{nCores = 1}).
}

\note{All sampling error estimates are returned as percentages, and represent ~68\% confidence (1 standard deviation). To compute sampling error percent with 95\% confidence, multiply by 1.96.}

\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, values of areal coverage are returned for each plot. All variables with names ending in \code{SE}, represent the estimate of sampling error (\%) of the variable.

\itemize{
      \item{\strong{STAGE}: Stand structural stage.}
      \item{\strong{PERC_AREA}: \% land area in each structural stage.}
      \item{\strong{TOTAL_AREA}: Total land area used to produce estimates (forest or timberland) (acres).}
      }
}


\references{
FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}

Frelich, L. E., and Lorimer, C. G. (1991). Natural Disturbance Regimes in Hemlock-Hardwood Forests of the Upper Great Lakes Region. Ecological Monographs, 61(2), 145-164. doi:10.2307/1943005

Goodell, L., and Faber-Langendoen, D. (2007). Development of stand structural stage indices to characterize forest condition in Upstate New York. Forest Ecology and Management, 249(3), 158-170. doi:10.1016/j.foreco.2007.04.052
}

\author{
Hunter Stanke and Andrew Finley
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tpa}}, \code{\link{diversity}}
}
\examples{
\dontrun{

## Load data from rFIA package
library(rFIA)
data(fiaRI)
data(countiesRI)

## Calculate structural stage distribution of all forestland, over time
standStruct(fiaRI)

## Same as above, but implemented in parallel
detectCores() # 4 cores available, we will take 3
standStruct(fiaRI, nCores = 3)

## Same as above, but output contains seperate column for each structural stage,
##   rathern than grouping variable
standStruct(fiaRI, tidy = FALSE)

## Calculate structural stage distribution of all forestland by owner group, over time
standStruct(clipFIA(fiaRI, mostRecent = TRUE),
            grpBy = c(OWNGRPCD))

## Calculate structural stage distribution of all forestland on xeric sites, over time
standStruct(clipFIA(fiaRI, mostRecent = TRUE),
            areaDomain = PHYSCLCD \%in\% c(11:19))

## Calculate structural stage distribution of all forestland by county and return
#   return spatial object
sdSF <- standStruct(clipFIA(fiaRI, mostRecent = TRUE),
            polys = countiesRI,
            returnSpatial = TRUE,
            tidy = FALSE)
plot(sdSF)
plotFIA(sdSF, POLE_PERC)
}
}
