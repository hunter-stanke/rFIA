\name{fiaRI}
\alias{fiaRI}
\docType{data}
\title{FIADB for Rhode Island 1985 - 2017}

\description{Subset of the Forest Inventory and Analysis Database for the state of Rhode Island. Reporting years range from 1985 - 2017. Specify \code{fiaRI} as the \code{db} argument in any \code{rFIA} function to produce estimates for the state fo Rhode Island.

Download other subsets of the FIA Database from the FIA Datamart: \url{https://apps.fs.usda.gov/fia/datamart/datamart.html}. Once downloaded, unzip the directory, and read into R using \code{\link{readFIA}}.
}

\usage{data("fiaRI")}
\format{
---- FIA Database Object -----

Reporting Years:  1985 1998 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017

States:           RHODE ISLAND

Total Plots:      222

Memory Used:      35.5 Mb

Tables:           COND_DWM_CALC COND INVASIVE_SUBPLOT_SPP PLOT POP_ESTN_UNIT POP_EVAL_GRP POP_EVAL_TYP POP_EVAL POP_PLOT_STRATUM_ASSGN POP_STRATUM SUBP_COND_CHNG_MTRX SUBPLOT TREE_GRM_COMPONENT TREE_GRM_ESTN TREE
}

\examples{
data(fiaRI)
summary(fiaRI)
print(fiaRI)
}
\keyword{datasets}
