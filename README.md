
<!-- README.md is generated from README.Rmd. Please edit that file -->
rFIA: Unlocking the FIA Database in R
=====================================

![US Biomass](man/figures/usBiomass.png) <!-- badges: start --> <!-- badges: end -->

The goal of `rFIA` is to increase the accessibility and use of the USFS Forest Inventory and Analysis (FIA) Database by providing a user-friendly, open source platform to easily query and analyze FIA Data. Designed to accommodate a wide range of potential user objectives, `rFIA` simplifies the estimation of forest variables from the FIA Database and allows all R users (experts and newcomers alike) to unlock the flexibility and potential inherent to the Enhanced FIA design.

Specifically, `rFIA` improves accessibility to the spatio-temporal estimation capacity of the FIA Database by producing space-time indexed summaries of forest variables within user-defined population boundaries. Direct integration with other popular R packages (e.g., dplyr, sp, and sf) facilitates efficient space-time query and data summary, and supports common data representations and API design. The package implements design-based estimation procedures outlined by Bechtold & Patterson (2005), and has been validated against estimates and sampling errors produced by EVALIDator. Current development is focused on the implementation of spatially-enabled model-assisted estimators to improve population, change, and ratio estimates.

You can download subsets of the FIA Database at the FIA Datamart: <https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html>

<br>

Installation
------------

You can install the released version of `rFIA` from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rFIA")
```

Alternatively, you can install the development version from GitHub:

``` r
devtools::install_github('hunter-stanke/rFIA')
```

<br>

Functionality
-------------

| `rFIA` Function | Description                                                        |
|-----------------|--------------------------------------------------------------------|
| `area`          | Estimate land area                                                 |
| `biomass`       | Estimate volume, biomass, & carbon stocks of standing trees        |
| `clipFIA`       | Spatial & temporal queries                                         |
| `diversity`     | Estimate species diversity                                         |
| `dwm`           | Estimate volume, biomass, and carbon stocks of down woody material |
| `getFIA`        | Download FIA data, load into R, and optionally save to disk        |
| `growMort`      | Estimate recruitment, mortality, and harvest rates                 |
| `invasive`      | Estimate areal coverage of invasive species                        |
| `plotFIA`       | Produce static & animated plots of spatial FIA summaries           |
| `readFIA`       | Load FIA database into R environment                               |
| `standStruct`   | Estimate forest structural stage distributions                     |
| `tpa`           | Estimate abundance of standing trees (TPA & BAA)                   |
| `vitalRates`    | Estimate live tree growth rates                                    |

<br>

Example Usage
-------------

### *Download FIA Data and Load into R*

The first step to using `rFIA` is to download subsets of the FIA Database. The easiest way to accomplish this is using `getFIA`. Using one line of code, you can download state subsets of the FIA Database, load data into your R environment, and optionally save those data to a local directory for future use!

``` r
## Download the state subset or Connecticut (requires an internet connection)
# All data acquired from FIA Datamart: https://apps.fs.usda.gov/fia/datamart/datamart.html
ct <- getFIA(states = 'CT', dir = '/path/to/save/data')
```

By default, `getFIA` only downloads the portions of the database required to produce summaries with other `rFIA` functions (`common = TRUE`). This conserves memory on your machine and speeds download time. If you would like to download all available tables for a state, simple specify `common = FALSE` in the call to `getFIA`.

**But what if I want to load multiple states worth of FIA data into R?** No problem! Simply specify mutiple state abbreviations in the `states` argument of `getFIA` (e.g. `states = c('MI', 'IN', 'WI', 'IL'`)), and all state subsets will be downloaded and merged into a single `FIA.Database` object. This will allow you to use other `rFIA` functions to produce estimates within polygons which straddle state boundaries!

Note: given the massive size of the full FIA Database, users are cautioned to only download the subsets containing their region of interest.

**If you have previously downloaded FIA data would simply like to load into R from a local directory, use `readFIA`:**

``` r
## Load FIA Data from a local directory
db <- readFIA('/path/to/your/directory/')
```

------------------------------------------------------------------------

### *Compute Estimates of Forest Variables*

Now that you have loaded your FIA data into R, it's time to put it to work. Let's explore the basic functionality of `rFIA` with `tpa`, a function to compute tree abundance estimates (TPA, BAA, & relative abundance) from FIA data, and `fiaRI`, a subset of the FIA Database for Rhode Island including inventories from 2013-2018.

**Estimate the abundance of live trees in Rhode Island:**

``` r
library(rFIA)
## Load the Rhode Island subset of the FIADB (included w/ rFIA)
## NOTE: This object can be produced using getFIA and/or readFIA
data("fiaRI")

## Only estimates for the most recent inventory year
fiaRI_MR <- clipFIA(fiaRI, mostRecent = TRUE) ## subset the most recent data
tpaRI_MR <- tpa(fiaRI_MR)
head(tpaRI_MR)
#>   YEAR      TPA      BAA TPA_PERC BAA_PERC   TPA_SE   BAA_SE TPA_PERC_SE
#> 1 2018 426.7119 122.1143 93.22791   93.677 6.631549 3.057083     7.61641
#>   BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    4.477168         126         127

## All Inventory Years Available (i.e., returns a time series)
tpaRI <- tpa(fiaRI)

## Time Series plot
plotFIA(tpaRI, BAA, plot.title = 'Basal area per acre in Rhode Island over time')
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

**What if I want to group estimates by species? How about by size class?**

``` r
## Group estimates by species
tpaRI_species <- tpa(fiaRI_MR, bySpecies = TRUE)
head(tpaRI_species, n = 3)
#>   YEAR SPCD          COMMON_NAME        SCIENTIFIC_NAME        TPA
#> 1 2018   12           balsam fir         Abies balsamea 0.08732464
#> 2 2018   43 Atlantic white-cedar Chamaecyparis thyoides 0.24652918
#> 3 2018   68     eastern redcedar   Juniperus virginiana 1.13812339
#>          BAA   TPA_PERC   BAA_PERC    TPA_SE    BAA_SE TPA_PERC_SE
#> 1 0.02945721 0.01907867 0.02259738 114.02897 114.02897    7.700936
#> 2 0.17960008 0.05386163 0.13777578  59.08447  55.97810    7.638942
#> 3 0.13844250 0.24865690 0.10620276  64.77184  67.45948    7.643501
#>   BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    4.619978           1         127
#> 2    4.511914           3         127
#> 3    4.527460           5         127

## Group estimates by size class
## NOTE: Default 2-inch size classes, but you can make your own using makeClasses()
tpaRI_sizeClass <- tpa(fiaRI_MR, bySizeClass = TRUE)
head(tpaRI_sizeClass, n = 3)
#>   YEAR sizeClass       TPA      BAA TPA_PERC BAA_PERC    TPA_SE    BAA_SE
#> 1 2018     [1,3) 200.04626 3.703843 43.68492 3.004808  9.350683  9.599735
#> 2 2018     [3,5)  67.04905 5.636636 14.64177 4.572821 12.110857 12.699784
#> 3 2018     [5,7)  44.10363 8.579190  9.63109 6.960020  5.338930  5.392966
#>   TPA_PERC_SE BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    4.464367    1.978547          76         126
#> 2    4.465090    1.980327          46         126
#> 3    4.463748    1.976972         115         126

## Group by species and size class, and plot the distribution 
##  for the most recent inventory year
tpaRI_spsc <- tpa(fiaRI_MR, bySpecies = TRUE, bySizeClass = TRUE)
plotFIA(tpaRI_spsc, BAA, grp = COMMON_NAME, x = sizeClass,
        plot.title = 'Size-class distributions of BAA by species', 
        x.lab = '', text.size = .64, n.max = 5) # Only want the top 5 species, try n.max = -5 for bottom 5
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

**What if I want estimates for a specific type of tree (ex. greater than 12-inches DBH and in a canopy dominant or subdominant position) in specific area (ex. growing on mesic sites), and I want to group by estimates by some variable other than species or size class (ex. ownsership group)?** Easy! Each of these specifications are described in the FIA Database, and all `rFIA` functions can leverage these data to easily implement complex queries!

``` r
## grpBy specifies what to group estimates by (just like species and size class above)
## treeDomain describes the trees of interest, in terms of FIA variables 
## areaDomain, just like above,describes the land area of interest
tpaRI_own <- tpa(fiaRI_MR, 
                     grpBy = OWNGRPCD, 
                     treeDomain = DIA > 12 & CCLCD %in% c(1,2),
                     areaDomain = PHYSCLCD %in% c(20:29))
head(tpaRI_own)
#>   YEAR OWNGRPCD       TPA      BAA TPA_PERC BAA_PERC   TPA_SE   BAA_SE
#> 1 2018       30 0.8482522 3.567751      100      100 58.95933 59.06494
#> 2 2018       40 1.4920940 3.992055      100      100 25.69191 27.70682
#> 3 2018       NA 1.2884293 3.857836      100      100 23.72942 25.98661
#>   TPA_PERC_SE BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    60.69603    60.74652           3          38
#> 2    26.85686    28.61007          12          82
#> 3    24.15310    26.25434          15         118
```

**What if I want to produce estimates within my own population boundaries (within user-defined spatial zones/polygons)?** This is where things get really exciting.

``` r
## Load the county boundaries for Rhode Island
data('countiesRI') ## Load your own spatial data from shapefiles using readOGR() (rgdal)

## polys specifies the polygons (zones) where you are interested in producing estimates
## returnSpatial = TRUE indicates that the resulting estimates will be joined with the 
##    polygons we specified, thus allowing us to visualize the estimates across space
tpaRI_counties <- tpa(fiaRI_MR, polys = countiesRI, returnSpatial = TRUE)

## NOTE: Any grey polygons below simply means no FIA data was available for that region
plotFIA(tpaRI_counties, BAA) # Plotting method for spatial FIA summaries, also try 'TPA' or 'TPA_PERC'
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

**We produced a really cool time series earlier, how would I marry the spatial and temporal capacity of `rFIA` to produce estimates across user-defined polygons and through time?** Easy! Just hand `tpa` the full FIA.Database object you produced with `readFIA` (not the most recent subset produced with `clipFIA`). For stunning space-time visualizations, hand the output of `tpa` to `plotFIA`. To save the animation as a .gif file, simpy specify `fileName` (name of output file) and `savePath` (directory to save file, combined with `fileName`).

``` r
## Using the full FIA dataset, all available inventories
tpaRI_st <- tpa(fiaRI, polys = countiesRI, returnSpatial = TRUE)

## Animate the output
plotFIA(tpaRI_st, TPA, animate = TRUE, legend.title = 'Abundance (TPA)', legend.height = .8)
```

<img src="man/figures/README-unnamed-chunk-8-1.gif" width="100%" />
