\name{customPSE}
\alias{customPSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generic post-stratified estimator for FIA Data}

\description{More soon.}
\usage{
customPSE(db, evalType, x, xVars, xGrpBy = NULL, y = NULL,
          yVars = NULL, yGrpBy = NULL, method = "TI",
          lambda = 0.5, totals = TRUE, variance = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{\code{FIA.Database} or \code{Remote.FIA.Database} object produced from \code{\link{readFIA}} or \code{\link{getFIA}}. If a \code{Remote.FIA.Database}, data will be read in and processed state-by-state to conserve RAM (see details for an example).}

  \item{evalType}{character ('ALL', 'CURR', 'VOL', 'GROW', 'MORT', 'REMV', 'CHNG', 'DWM', 'REGEN'). See Reference Population Evaluation Table Description Type Table in FIADB P2 User Guide (References) for descriptions of evaluation types.}

  \item{x}{data.frame containing numerator variables of interest. More soon on requirements.}
  \item{xVars}{names of variables in \code{x} to be treated as numerator variables, unquoted.}
  \item{xGrpBy}{names of variables \code{x} to be treated as grouping variables, unquoted}
  \item{y}{data.frame containing numerator variables of interest. More soon on requirements.}
  \item{yVars}{names of variables in \code{y} to be treated as numerator variables, unquoted.}
  \item{yGrpBy}{names of variables \code{y} to be treated as grouping variables, unquoted.}

  \item{method}{character; design-based estimator to use. One of:  "TI" (temporally indifferent, default), "annual" (annual), "SMA"" (simple moving average), "LMA" (linear moving average), or "EMA" (exponential moving average). See \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020} for a complete description of these estimators.}

  \item{lambda}{numeric (0,1); if \code{method = 'EMA'}, the decay parameter used to define weighting scheme for annual panels. Low values place higher weight on more recent panels, and vice versa. Specify a vector of values to compute estimates using mulitple wieghting schemes, and use \code{plotFIA} with \code{grp} set to \code{lambda} to produce moving average ribbon plots. See \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020} for examples.}

  \item{totals}{logical; if TRUE, return total population estimates (e.g. total area) along with ratio estimates (e.g. mean trees per acre).}

    \item{variance}{logical; if TRUE, return estimated variance (\code{VAR}) and sample size (\code{N}). If FALSE, return 'sampling error' (\code{SE}) as returned by EVALIDator. Note: sampling error cannot be used to construct confidence intervals.}

}


\details{
More soon.
}


\value{
More soon.
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

data(fiaRI)

}
