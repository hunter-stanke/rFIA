\name{customPSE}
\alias{customPSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post-stratified estimator for custom variables in FIA Data}
\description{
Produces estimates of population totals, ratios, and associated variances for custom variables using FIA's post-stratified inventories. Accepts tree- and condition-level summaries from rFIA estimator functions as input, with potential modifications to associated variables (e.g., custom allometrics can be applied to estimate tree biomass/carbon). See our website for example use cases.
}
\usage{
customPSE(db, x, xVars, xGrpBy = NULL, y = NULL,
          yVars = NULL, yGrpBy = NULL, method = "TI",
          lambda = 0.5, totals = TRUE, variance = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{\code{FIA.Database} or \code{Remote.FIA.Database} object produced from \code{\link{readFIA}} or \code{\link{getFIA}}. If a \code{Remote.FIA.Database}, data will be read in and processed state-by-state to conserve RAM (see details for an example).}

  \item{x}{data.frame; tree- or condition-list containing numerator variable(s). See details more info on producing acceptable tree- and condition-lists using \code{rFIA} estimator functions.}

  \item{xVars}{name of variable(s) in \code{x} to be treated as numerator variables, unquoted. Multiple variables should be combined with \code{c()}, e.g., \code{xVars = c(TPA, BAA)}.}

  \item{xGrpBy}{names of variables in \code{x} to be treated as grouping variables for the numerator, unquoted. Multiple variables should be combined with \code{c()}, e.g., \code{xVars = c(SPCD, OWNGRPCD)}.}

  \item{y}{data.frame; tree- or condition-list containing denominator variable. See details more info on producing acceptable tree- and condition-lists using \code{rFIA} estimator functions.}

  \item{yVars}{name of variable in \code{y} to be treated as denominator variables, unquoted. Maximum one denominator variable allowed at this time.}

  \item{yGrpBy}{names of variables in \code{y} to be treated as grouping variables for the numerator, unquoted. Multiple variables should be combined with \code{c()}, e.g., \code{xVars = c(SPCD, OWNGRPCD)}.}

  \item{method}{character; design-based estimator to use. One of:  "TI" (temporally indifferent, default), "annual" (annual), "SMA"" (simple moving average), "LMA" (linear moving average), or "EMA" (exponential moving average). See \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020} for a complete description of these estimators.}

  \item{lambda}{numeric (0,1); if \code{method = 'EMA'}, the decay parameter used to define weighting scheme for annual panels. Low values place higher weight on more recent panels, and vice versa. Specify a vector of values to compute estimates using mulitple wieghting schemes, and use \code{plotFIA} with \code{grp} set to \code{lambda} to produce moving average ribbon plots. See \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020} for examples.}

  \item{totals}{logical; if TRUE, return total population estimates (e.g. total area) along with ratio estimates (e.g. mean trees per acre).}

    \item{variance}{logical; if TRUE, return estimated variance (\code{VAR}) and sample size (\code{N}). If FALSE, return 'sampling error' (\code{SE}) as returned by EVALIDator. Note: sampling error cannot be used to construct confidence intervals.}

}


\details{

\strong{Workflow and intended use cases}


\code{customPSE} is intended to be used in combination with standard rFIA estimator functions, like \code{\link{tpa}}, \code{\link{area}}, and \code{\link{volume}}, among others. Standard rFIA estimator functions generate tree- and or condition-lists for standard variables of interest (see \code{treeList} and \code{condList} arguments in estimator functions). Users may make modifications to these standard variables, for example a variable representing tree crown area may be added to a tree-list produced by \code{\link{tpa}} (via some suite of allometrics). Users may then hand their modified tree-list to \code{customPSE} to estimate the total and proportion of forested land area in their domain of interest that is covered by tree crowns.

\code{customPSE} may be used to estimate population totals for multiple variables simultaneously (total number of trees in a region), and given a denominator variable, the associated population ratios (trees per forested acre, where forested land area is the denominator). Estimation follows the procedures documented in Bechtold and Patterson (2005) and \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020}.

Three general forms of ratio estimates may be produced: tree-tree, tree-area, and area-area ratios. For example, if tree height is specified as the numerator variable (adjusted for sampling area by multiplying by TPA), and TPA is specified as the denominator variable, a tree-tree ratio will be produced that represents the height of the average tree within a region of interest. Similarly, if stand age is specified as the numerator (adjusted for sampling area by proportionate area of the forested condition on the plot, i.e., \code{PROP_FOREST}), and the proportion area of the plot that is forested is specified as the denominator, an area-area ratio will be produced that represents average stand age within the region of interest. Tree-area ratios are more familiar, such as trees per acre, tree biomass per acre, etc, where a tree variable is specified as the numerator, and proportion of plot area occupied by forestland is the denominator. See our website for detailed examples of each of these ratio estimates.




\strong{Input requirements}


Estimation of tree variables require the following columns be present in \code{x} and/or \code{y}: PLT_CN, EVAL_TYP, SUBP, TREE, and TREE_BASIS. Similarly, estimation of area variables require the following columns be present in \code{x} and/or \code{y}: PLT_CN, EVAL_TYP, CONDID, and AREA_BASIS. Each of these required variables will be returned in tree- and condition-lists generated by standard rFIA estimator functions.

IMPORTANT: Only one of TREE_BASIS or AREA_BASIS may be present \code{x} or \code{y}, as the presence of these columns are used to determine if variables to be estimated are tree variables or area variables. Some standard rFIA estimator functions will produce tree-lists with both TREE_BASIS and AREA_BASIS listed in output, as the tree-list will contain tree variables (e.g., TPA, BAA) as well as area variables (e.g., PROP_FOREST, proportion of plot represented by the forested condition where each tree is growing)). To produce a tree-area ratio with such an output, AREA_BASIS must be removed from the data.frame specified in \code{x}, and TREE_BASIS must be removed from that specified in \code{y}.


}


\value{

\itemize{
      \item{\strong{YEAR}: reporting year associated with estimates}
      \item{\strong{*_RATIO}: population ratio estimate, where * will be replaced with the name of each numerator variable}
      \item{\strong{*_TOTAL}: population total estimate, where * will be replaced with the name of each numerator/ denominator variable}
      \item{\strong{*_RATIO_VAR}: estimated variance of the population ratio}
      \item{\strong{*_TOTAL_VAR}: estimated variance of the population total}
      \item{\strong{nPlots_x}: number of non-zero plots used to compute numerator estimates}
      \item{\strong{nPlots_y}: number of non-zero plots used to compute denominator estimates}
      \item{\strong{N}: total number of plots (including zeros) associated with each inventory}
      }
}



\references{
rFIA website: \url{https://rfia.netlify.app/}

FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}

Stanke, H., Finley, A. O., Weed, A. S., Walters, B. F., & Domke, G. M. (2020). rFIA: An R package for estimation of forest attributes with the US Forest Inventory and Analysis database. Environmental Modelling & Software, 127, 104664.
}


\author{
Hunter Stanke and Andrew Finley
}

\note{All sampling error estimates (SE) are returned as the "percent coefficient of variation" (standard deviation / mean * 100) for consistency with EVALIDator. IMPORTANT: sampling error cannot be used to construct confidence intervals. Please use \code{variance = TRUE} for that (i.e., return variance and sample size instead of sampling error).}


\seealso{\code{\link{getDesignInfo}}}


\examples{
## See our website for a more thorough suite of examples
data(fiaRI)

# Get tree-list from tpa
tree.list <- tpa(fiaRI,
                 treeList = TRUE)

# Estimate trees per acre and basal area per acre
customPSE(db = fiaRI,
          # Numerator variables
          x = dplyr::select(tree.list, -c(AREA_BASIS)),
          xVars = c(TPA, BAA),
          # Denominator variables
          y = dplyr::select(tree.list, -c(TREE_BASIS)),
          yVars = PROP_FOREST)

# Same as above, but rename variables for a clean output
customPSE(db = fiaRI,
          x = dplyr::select(tree.list, -c(AREA_BASIS)),
          # Variables can be renamed using c()
          xVars = c(NUM = TPA,
                    BA = BAA),
          y = dplyr::select(tree.list, -c(TREE_BASIS)),
          # Variables can be renamed using c()
          yVars = c(FOREST_AREA = PROP_FOREST))

# Ensure the above matches expected output
tpa(fiaRI,
    totals = TRUE,
    variance = TRUE)
}
