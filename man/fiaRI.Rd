\name{fiaRI}
\alias{fiaRI}
\docType{data}
\title{FIADB for Rhode Island 2013 - 2018}

\description{Subset of the Forest Inventory and Analysis Database for the state of Rhode Island. Reporting years range from 2013 - 2018. Specify \code{fiaRI} as the \code{db} argument in any \code{rFIA} function to produce estimates for the state fo Rhode Island.

Download other subsets of the FIA Database from the FIA Datamart: \url{https://apps.fs.usda.gov/fia/datamart/datamart.html}. Once downloaded, unzip the directory, and read into R using \code{\link{readFIA}}.
}

\usage{data("fiaRI")}
\format{
---- FIA Database Object -----
Reporting Years:  2013 2014 2015 2016 2017 2018
States:           RHODE ISLAND
Total Plots:      769
Memory Used:      19.5 Mb
Tables:           COND COND_DWM_CALC INVASIVE_SUBPLOT_SPP PLOT POP_ESTN_UNIT POP_EVAL POP_EVAL_GRP POP_EVAL_TYP POP_PLOT_STRATUM_ASSGN POP_STRATUM SUBPLOT TREE TREE_GRM_COMPONENT TREE_GRM_ESTN SUBP_COND_CHNG_MTRX
}

\examples{
data(fiaRI)
summary(fiaRI)
print(fiaRI)
}
\keyword{datasets}
