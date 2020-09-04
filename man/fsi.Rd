\name{fsi}
\alias{fsi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the Forest Stability Index from the FIADB}
\description{Estimate forest population performance from the FIADB using the Forest Stability Index. This function is experimental, please check back soon for updated documentation. }
\usage{
fsi(db, grpBy = NULL, polys = NULL, returnSpatial = FALSE,
    bySpecies = FALSE, bySizeClass = FALSE,
    landType = "forest", treeType = "live", method = "sma",
    lambda = 0.5, treeDomain = NULL, areaDomain = NULL,
    totals = TRUE, variance = TRUE, byPlot = FALSE,
    useSeries = FALSE, scaleBy = NULL, betas = NULL,
    returnBetas = FALSE, nCores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{\code{FIA.Database} object produced from \code{\link{readFIA}}; Function requires that PLOT, TREE, COND, POP_PLOT_STRATUM_ASSGN, POP_ESTN_UNIT, POP_EVAL, POP_STRATUM, POP_EVAL_TYP, POP_EVAL_GRP tables exist in FIA.Database object.}

\item{grpBy}{variables from PLOT, COND, or TREE tables to group estimates by (NOT quoted). Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(ECOSUBCD, OWNGRPCD)}.}

  \item{polys}{\code{sp} or \code{sf} Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit.}

  \item{returnSpatial}{logical; if TRUE, return sf spatial object (polys must also be specified).}

\item{bySpecies}{logical; if TRUE, returns estimates grouped by species.}

  \item{bySizeClass}{logical; if TRUE, returns estimates grouped by size class (2-inch intervals, see \code{\link{makeClasses}} to compute different size class intervals).}

  \item{landType}{character ('forest' or 'timber'); Type of land which estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{treeType}{character ('all', 'live', 'dead', or 'gs'); Type of tree which estimates will be produced for. All (default) includes all stems, live and dead, greater than 1 in. DBH. Live/Dead includes all stems greater than 1 in. DBH which are live or dead (leaning less than 45 degrees), respectively. GS (growing-stock) includes live stems greater than 5 in. DBH which contain at least one 8 ft merchantable log. }

  \item{method}{character; Method used for annual panel combination (see details). One of: 'TI' (temporally indifferent), 'annual' (annual panels), 'SMA' (simple moving average), 'LMA' (linear moving average), or 'EMA' (exponential moving average)}

  \item{lambda}{numeric (0,1); if \code{method == 'EMA'}, the decay parameter used to define weighting scheme for annual panels. Low values place higher weight on more recent panels, and vice versa. Specify a vector of values to compute estimates using mulitple wieghting schemes, and use \code{plotFIA} with \code{grp} set to \code{lambda} to produce moving average ribbon plots.}

  \item{treeDomain}{logical predicates defined in terms of the variables in PLOT, TREE, and/or COND tables. Used to define the type of trees for which estimates will be produced (e.g. DBH greater than 20 inches: \code{DIA > 20}, Dominant/Co-dominant crowns only: \code{CCLCD \%in\% c(2,3))}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only trees where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{areaDomain}{logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (e.g. within 1 mile of improved road: \code{RDDISTCD \%in\% c(1:6)}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{totals}{logical; if TRUE, return population estimates (e.g. total area, total biomass) along with ratio estimates (e.g. mean biomass per acre).}

  \item{variance}{logical; if TRUE, return estimates of uncertainty as variance instead of sampling error (default). Sampling error cannot be used to estimate confidence intervals.}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations (population totals not returned).}

  \item{useSeries}{logical; If TRUE, use multiple remeasurements to estimate change in TPA & BA on each plot, when available.}

  \item{scaleBy}{variables from PLOT or COND tables to use as 'random effects' in model of size-density relationships.Multiple variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(FORTYPCD, SITECLCD)}.}
  \item{betas}{more soon.}

  \item{returnBetas}{more soon.}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}

}
\details{
Please check back soon for more details
}
\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, totals are returned for each plot. All variables with names ending in \code{SE} represent the estimate of sampling error (\%) of the variable. All variables with names ending in \code{TOTAL} represent the population total of the variable.

\itemize{
      \item{\strong{YEAR}: reporting year associated with estimates}
      \item{\strong{SI}: estimate of forest stability index }
      \item{\strong{SI_STATUS}: indication of the forest stability index (i.e., decline, stable, or expand)}
      \item{\strong{SI_INT}: width of 95\% confidence interval of mean FSI}
      \item{\strong{TPA_RATE}: standardized estimate of annual change in TPA}
      \item{\strong{BAA_RATE}: standardized estimate of annual change in BAA}
}
}

\author{
Hunter Stanke and Andrew Finley
}

\references{
FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}
}

\seealso{
\code{\link{tpa}}, \code{\link{vitalRates}}, \code{\link{growMort}}
}

\note{All sampling error estimates are returned as percentages, and represent ~68\% confidence (1 standard deviation). To compute sampling error percent with 95\% confidence, multiply by 1.96.}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
## Load data from the rFIA package
data(fiaRI)
data(countiesRI)

## Most recents subset
fiaRI_mr <- clipFIA(fiaRI)

## FSI for all forestland 2018 in RI
fsi(fiaRI_mr, method = 'sma')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
