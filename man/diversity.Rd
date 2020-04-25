\name{diversity}
\alias{diversity}
\title{Estimate diversity from FIADB}

\description{
Produces estimates of diversity from FIA data. Returns shannon's index, shannon's equitability, and richness for alpha (mean/SE of stands), beta, and gamma diversity. Default behavior estimates species diversity, using \code{TPA} as a state variable and \code{SPCD} to groups of individuals. Estimates can be produced for regions defined within the FIA Database (e.g. counties), at the plot level, or within user-defined areal units. Options to group estimates by size class and other variables defined in the FIADB. If multiple reporting years (EVALIDs) are included in the data, estimates will be output as a time series. If multiple states are represented by the data, estimates will be output for the full region (all area combined), unless specified otherwise (e.g. \code{grpBy = STATECD}). Easy options to implement parallel processing.
}

\usage{
diversity(db, grpBy = NULL, polys = NULL, returnSpatial = FALSE, bySizeClass = FALSE,
          landType = 'forest', treeType = 'live', method = 'TI', lambda = .5,
          stateVar = TPA_UNADJ, grpVar = SPCD, treeDomain = NULL,
          areaDomain = NULL, byPlot = FALSE, totals = FALSE, nCores = 1)
}

\arguments{
  \item{db}{\code{FIA.Database} object produced from \code{\link{readFIA}}. Function requires that PLOT, TREE, COND, POP_PLOT_STRATUM_ASSGN, POP_ESTN_UNIT, POP_EVAL, POP_STRATUM, POP_EVAL_TYP, POP_EVAL_GRP tables exist in FIA.Database object.}

  \item{grpBy}{variables from PLOT, COND, or TREE tables to group estimates by (NOT quoted). Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(ECOSUBCD, OWNGRPCD)}.}

  \item{polys}{\code{sp} or \code{sf} Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit.}

  \item{returnSpatial}{logical; if TRUE, return sf spatial object (polys must also be specified).}

   \item{bySizeClass}{logical; if TRUE, returns estimates grouped by size class (default 2-inch intervals, see \code{\link{makeClasses}} to compute other size class intervals).}

  \item{landType}{character ('forest' or 'timber'); Type of land which estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{treeType}{character ('all', 'live', 'dead', or 'gs'); Type of tree which estimates will be produced for. All (default) includes all stems, live and dead, greater than 1 in. DBH. Live/Dead includes all stems greater than 1 in. DBH which are live or dead (leaning less than 45 degrees), respectively. GS (growing-stock) includes live stems greater than 5 in. DBH which contain at least one 8 ft merchantable log.}

  \item{method}{character; Method used for annual panel combination (see details). One of: 'TI' (temporally indifferent), 'annual' (annual panels), 'SMA' (simple moving average), 'LMA' (linear moving average), or 'EMA' (exponential moving average)}

  \item{lambda}{numeric (0,1); if \code{method == 'EMA'}, the decay parameter used to define weighting scheme for annual panels. Low values place higher weight on more recent panels, and vice versa. Specify a vector of values to compute estimates using mulitple wieghting schemes, and use \code{plotFIA} with \code{grp} set to \code{lambda} to produce moving average ribbon plots.}

  \item{stateVar}{variable from TREE table to use as state variable (NOT quoted). Default, \code{TPA_UNADJ}. Try, \code{DRYBIO_AG} for aboveground biomass, \code{pi*(DIA/2)^2} for basal area, or others.}

  \item{grpVar}{factor, variable from TREE table to define individual groups (NOT quoted). Default, \code{SPCD}. Try, \code{SPGRPCD} for species group, \code{makeClasses(db$TREE$DIA, interval = 2)} for diameter class, or others.}

  \item{treeDomain}{logical predicates defined in terms of the variables in PLOT, TREE, and/or COND tables. Used to define the type of trees for which estimates will be produced (e.g. DBH greater than 20 inches: \code{DIA > 20}, Dominant/Co-dominant crowns only: \code{CCLCD \%in\% c(2,3))}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only trees where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{areaDomain}{logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (e.g. within 1 mile of improved road: \code{RDDISTCD \%in\% c(1:6)}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations (population totals not returned).}

  \item{totals}{logical; if TRUE, returns population estimates in addition to ratios.}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}
}


\details{
Procedures for computing diversity indices are outlined in Hill (1973) and Shannon (1948), and estimation of mean/ sampling error follows the procedures documented in Bechtold and Patterson (2005).

Alpha-level indices are computed as the mean (SE) diversity of a stand or condition as defined by FIA. Specifically, alpha diversity is estimated using a sample-based ratio-of-means estimator of stand diversity (e.g. Richness) * land area of stand / total land area within the domain of interest. Thus estimates of alpha diversity within a stand are weighted by the area which that stand represents. Gamma-level diversity is computed as regional indicies, pooling all plot data. Beta diversity is computed as gamma diversity - alpha diversity, and thus represents the excess of regional diversity with respect to local diversity.

Users may specify alternatives to the 'Temporally Indifferent' method for panel combination using the \code{method} argument. Estimators are currently implemented for the temporally indifferent (TI) method (default, used by EVALIDator), simple moving average (SMA), exponential moving average (EMA), and annual estimation (no panel combination). Users may modify the temporal weighting scheme employed by the exponential moving average using the \code{yrs} argument. See Chapter 5 of Bechtold & Patterson (2005) for more details.

Stratified random sampling techniques are most often employed to compute estimates in recent inventories, although double sampling and simple random sampling may be employed for early inventories. Estimates are adjusted for non-response bias by assuming attributes of non-response plot locations to be equal to the mean of other plots included within thier respective stratum or population.

Forest land must be at least 10-percent stocked by trees of any size, including land that formerly had such tree cover and that will be naturally or artificially regenerated. Forest land includes transition zones, such as areas between heavily forested and nonforested lands that are at least 10-percent stocked with trees and forest areas adjacent to urban and builtup lands. The minimum area for classification of forest land is 1 acre and 120 feet wide measured stem-to-stem from the outer-most edge. Unimproved roads and trails, streams, and clearings in forest areas are classified as forest if less than 120 feet wide. Timber land is a subset of forest land that is producing or is capable of producing crops of industrial wood and not withdrawn from timber utilization by statute or administrative regulation. (Note: Areas qualifying as timberland are capable of producing at least 20 cubic feet per acre per year of industrial wood in natural stands. Currently inaccessible and inoperable areas are NOT included).

Easy, efficient parallelization is implemented with the \code{\link{parallel}} package. Users must only specify the \code{nCores} argument with a value greater than 1 in order to implement parallel processing on their machines. Parallel implementation is achieved using a snow type cluster on any Windows OS, and with multicore forking on any Unix OS (Linux, Mac). Implementing parallel processing may substantially decrease free memory during processing, particularly on Windows OS. Thus, users should be cautious when running in parallel, and consider implementing serial processing for this task if computational resources are limited (\code{nCores = 1}).
}

\note{All sampling error estimates are returned as percentages, and represent ~68\% confidence (1 standard deviation). To compute sampling error percent with 95\% confidence, multiply by 1.96.}

\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, indices are returned for each plot. All variables with names ending in \code{SE}, represent the estimate of sampling error (\%) of the variable.

\itemize{
      \item{\strong{H_a}: mean Shannon's Diversity Index, alpha (stand) level}
      \item{\strong{H_b}: Shannon's Diversity Index, beta (landscape) level}
      \item{\strong{H_g}: Shannon's Diversity Index, gamma (regional) level}
      \item{\strong{Eh_a}: mean Shannon's Equitability Index, alpha (stand) level}
      \item{\strong{Eh_b}: Shannon's Equitability Index, beta (landscape) level}
      \item{\strong{Eh_g}: Shannon's Equitability Index, alpha (stand) level}
      \item{\strong{S_a}: mean Species Richness, alpha (stand) level}
      \item{\strong{S_b}: Species Richness, beta (landscape) level}
      \item{\strong{S_g}: Species Richness, gamma (regional) level}
      \item{\strong{nStands}: number of stands with non-zero plots used to compute alpha diversity estimates}
      }
}

\references{
FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}

Analysis of ecological communities. (2002). United States: M G M SOFTWARE DESIGN (OR).

Hill, M. O. (1973). Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology, 54(2), 427-432. doi:10.2307/1934352.

Shannon, C. E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal, 27(3), 379-423. doi:10.1002/j.1538-7305.1948.tb01338.x.
}

\author{
Hunter Stanke and Andrew Finley
}

\seealso{
\code{\link{tpa}}, \code{\link{standStruct}}, \code{\link{invasive}}
}


\examples{
## Load data from rFIA package
data(fiaRI)
data(countiesRI)

## Make a most recent subset
fiaRI_mr <- clipFIA(fiaRI)

## Most recent estimates for live stems on forest land
diversity(db = fiaRI_mr,
          landType = 'forest',
          treeType = 'live')

\donttest{
## Most recent estimates grouped by stand age on forest land
# Make a categorical variable which represents stand age (grouped by 10 yr intervals)
fiaRI_mr$COND$STAND_AGE <- makeClasses(fiaRI_mr$COND$STDAGE, interval = 10)
diversity(db = fiaRI_mr,
          grpBy = STAND_AGE)

## Estimates for live white pine ( > 12" DBH) on forested mesic sites (all available inventories)
diversity(fiaRI,
          treeType = 'live',
          treeDomain = DIA > 12,
          areaDomain = PHYSCLCD \%in\% 21:29) # Mesic Physiographic classes

## Most recent estimates for growing-stock on timber land by species
diversity(db = fiaRI_mr,
          landType = 'timber',
          treeType = 'gs',
          bySizeClass = TRUE)

## Same as above, implemented in parallel
parallel::detectCores(logical = FALSE) # 4 cores available, we will take 2
diversity(db = fiaRI_mr,
          landType = 'timber',
          treeType = 'gs',
          bySizeClass = TRUE,
          nCores = 2)

## Most recent estimates for all stems on forest land grouped by user-defined areal units
ctSF <- diversity(clipFIA(fiaRI, mostRecent = TRUE),
                  polys = countiesRI,
                  returnSpatial = TRUE)
plot(ctSF) # Plot multiple variables simultaneously
plotFIA(ctSF, H_a) # Plot of mean Shannons Index of stands
}
}
