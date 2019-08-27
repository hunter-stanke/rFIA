
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rFIA: Unlocking the FIA Database in R

![US Biomass](static/usBiomass.png) <!-- badges: start -->
<!-- badges: end -->

The goal of `rFIA` is to increase the accessibility and use of the USFS
Forest Inventory and Analysis (FIA) Database by providing a
user-friendly, open source platform to easily query and analyze FIA
Data. Designed to accommodate a wide range of potential user objectives,
`rFIA` simplifies the estimation of forest variables from the FIA
Database and allows all R users (experts and newcomers alike) to unlock
the flexibility and potential inherent to the Enhanced FIA design.

Specifically, `rFIA` improves accessibility to the spatio-temporal
estimation capacity of the FIA Database by producing space-time indexed
summaries of forest variables within user-defined population boundaries.
Direct integration with other popular R packages (e.g., dplyr, sp, and
sf) facilitates efficient space-time query and data summary, and
supports common data representations and API design. The package
implements design-based estimation procedures outlined by Bechtold &
Patterson (2005), and has been validated against estimates and sampling
errors produced by EVALIDator. Current development is focused on the
implementation of spatially-enabled model-assisted estimators to improve
population, change, and ratio estimates.

You can download subsets of the FIA Database at the FIA Datamart:
<https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html>

<br>

## Installation

You can install the released version of `rFIA` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rFIA")
```

Alternatively, you can install the development version from
GitHub:

``` r
devtools::install_github('hunter-stanke/rFIA')
```

<br>

## Functionality

| `rFIA` Function | Description                                                        |
| --------------- | ------------------------------------------------------------------ |
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

## Example Usage

### *Download FIA Data and Load into R*

The first step to using `rFIA` is to download subsets of the FIA
Database. The easiest way to accomplish this is using `getFIA`. Using
one line of code, you can download state subsets of the FIA Database,
load data into your R environment, and optionally save those data to a
local directory for future
use\!

``` r
## Download the state subset or Connecticut (requires an internet connection)
# All data acquired from FIA Datamart: https://apps.fs.usda.gov/fia/datamart/datamart.html
ct <- getFIA(states = 'CT', dir = '/path/to/save/data')
```

By default, `getFIA` only downloads the portions of the database
required to produce summaries with other `rFIA` functions (`common =
TRUE`). This conserves memory on your machine and speeds download time.
If you would like to download all available tables for a state, simple
specify `common = FALSE` in the call to `getFIA`.

**But what if I want to load multiple states worth of FIA data into R?**
No problem\! Simply specify mutiple state abbreviations in the `states`
argument of `getFIA` (e.g. `states = c('MI', 'IN', 'WI', 'IL'`)), and
all state subsets will be downloaded and merged into a single
`FIA.Database` object. This will allow you to use other `rFIA` functions
to produce estimates within polygons which straddle state boundaries\!

Note: given the massive size of the full FIA Database, users are
cautioned to only download the subsets containing their region of
interest.

**If you have previously downloaded FIA data would simply like to load
into R from a local directory, use `readFIA`:**

``` r
## Load FIA Data from a local directory
db <- readFIA('/path/to/your/directory/')
```

-----

### *Compute Estimates of Forest Variables*

Now that you have loaded your FIA data into R, it’s time to put it to
work. Let’s explore the basic functionality of `rFIA` with `tpa`, a
function to compute tree abundance estimates (TPA, BAA, & relative
abundance) from FIA data, and `fiaRI`, a subset of the FIA Database for
Rhode Island including all inventories up to 2017.

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
#> 1 2017 441.1381 123.2122 95.38343  94.6289 6.658465 3.012546     7.65795
#>   BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    4.486253         124         125

## All Inventory Years Available (i.e., returns a time series)
tpaRI <- tpa(fiaRI)
## 2005 marks the first reporting year of the annual inventory, prior data is sparse
tail(tpaRI[tpaRI$YEAR >= 2005,], n = 4) 
#>    YEAR      TPA      BAA TPA_PERC BAA_PERC   TPA_SE   BAA_SE TPA_PERC_SE
#> 12 2014 466.0828 119.8941 96.81215 95.01087 6.726866 3.093355    7.624014
#> 13 2015 444.4095 120.5557 96.56753 94.95990 6.397598 3.056033    7.370677
#> 14 2016 449.7928 122.8870 95.81793 95.01799 6.458278 2.940997    7.454994
#> 15 2017 441.1381 123.2122 95.38343 94.62890 6.658465 3.012546    7.657950
#>    BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 12    4.597449         121         123
#> 13    4.479756         122         124
#> 14    4.486982         124         125
#> 15    4.486253         124         125

## Time Series plot
plotFIA(tpaRI, BAA, plot.title = 'Basal area per acre in Rhode Island over time')
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

**What if I want to group estimates by species? How about by size
class?**

``` r
## Group estimates by species
tpaRI_species <- tpa(fiaRI_MR, bySpecies = TRUE)
head(tpaRI_species, n = 3)
#>   YEAR SPCD          COMMON_NAME        SCIENTIFIC_NAME        TPA
#> 1 2017   12           balsam fir         Abies balsamea 0.08543602
#> 2 2017   43 Atlantic white-cedar Chamaecyparis thyoides 0.25204548
#> 3 2017   68     eastern redcedar   Juniperus virginiana 1.14522403
#>          BAA   TPA_PERC   BAA_PERC    TPA_SE    BAA_SE TPA_PERC_SE
#> 1 0.02882013 0.01847308 0.02213431 116.12825 116.12825    7.745135
#> 2 0.18331354 0.05449759 0.14078767  58.76023  55.63951    7.680104
#> 3 0.13385013 0.24762176 0.10279899  64.41035  65.79782    7.684579
#>   BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    4.634001           1         125
#> 2    4.520513           3         125
#> 3    4.534002           5         125

## Group estimates by size class
## NOTE: Default 2-inch size classes, but you can make your own using makeClasses()
tpaRI_sizeClass <- tpa(fiaRI_MR, bySizeClass = TRUE)
head(tpaRI_sizeClass, n = 3)
#>   YEAR sizeClass       TPA      BAA  TPA_PERC BAA_PERC    TPA_SE    BAA_SE
#> 1 2017     [1,3) 205.23681 3.669717 44.579534 2.986311  9.374920  8.810176
#> 2 2017     [3,5)  69.80959 5.821516 15.163356 4.737384 12.440803 13.026966
#> 3 2017     [5,7)  45.22111 8.753295  9.822488 7.123183  5.241419  5.313978
#>   TPA_PERC_SE BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    4.557014    2.051695          75         124
#> 2    4.557813    2.053974          44         124
#> 3    4.556390    2.050508         113         124

## Group by species and size class, and plot the distribution 
##  for the most recent inventory year
tpaRI_spsc <- tpa(fiaRI_MR, bySpecies = TRUE, bySizeClass = TRUE)
plotFIA(tpaRI_spsc, BAA, grp = COMMON_NAME, x = sizeClass,
        plot.title = 'Size-class distributions of BAA by species', 
        x.lab = '', text.size = .64, n.max = 5) # Only want the top 5 species, try n.max = -5 for bottom 5
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

**What if I want estimates for a specific type of tree (ex. greater than
12-inches DBH and in a canopy dominant or subdominant position) in
specific area (ex. growing on mesic sites), and I want to group by
estimates by some variable other than species or size class (ex.
ownsership group)?** Easy\! Each of these specifications are described
in the FIA Database, and all `rFIA` functions can leverage these data to
easily implement complex
queries\!

``` r
## grpBy specifies what to group estimates by (just like species and size class above)
## treeDomain describes the trees of interest, in terms of FIA variables 
## areaDomain, just like above,describes the land area of interest
tpaRI_own <- tpa(fiaRI_MR, 
                     grpBy = OWNGRPCD, 
                     treeDomain = DIA > 12 & CCLCD %in% c(1,2),
                     areaDomain = PHYSCLCD %in% c(20:29))
head(tpaRI_own)
#>   YEAR OWNGRPCD      TPA      BAA TPA_PERC BAA_PERC   TPA_SE   BAA_SE
#> 1 2017       30 1.238332 4.783645      100      100 49.84574 49.34242
#> 2 2017       40 1.435430 3.737055      100      100 24.80016 26.94152
#> 3 2017       NA 1.372623 4.070559      100      100 22.68298 24.92343
#>   TPA_PERC_SE BAA_PERC_SE nPlots_TREE nPlots_AREA
#> 1    51.46922    50.95504           4          38
#> 2    26.04668    27.89684          12          80
#> 3    23.14604    25.23145          16         116
```

**What if I want to produce estimates within my own population
boundaries (within user-defined spatial zones/polygons)?** This is where
things get really exciting.

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

**We produced a really cool time series earlier, how would I marry the
spatial and temporal capacity of `rFIA` to produce estimates across
user-defined polygons and through time?** Easy\! Just hand `tpa` the
full FIA.Database object you produced with `readFIA` (not the most
recent subset produced with `clipFIA`). For stunning space-time
visualizations, hand the output of `tpa` to `plotFIA`. To save the
animation as a .gif file, simpy specify `fileName` (name of output file)
and `savePath` (directory to save file, combined with `fileName`).

``` r
## Using the full FIA dataset, all available inventories
tpaRI_st <- tpa(fiaRI, polys = countiesRI, returnSpatial = TRUE)

## Animate the output
plotFIA(tpaRI_st, TPA, animate = TRUE, legend.title = 'Abundance (TPA)', legend.height = .8)
```

<img src="man/figures/README-unnamed-chunk-8-1.gif" width="100%" />
