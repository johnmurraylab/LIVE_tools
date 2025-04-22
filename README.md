
# embryo_lineage_analysis_R

R based programs to analyze and visualize embryo lineage after the cell
data are extracted from 3D movies in Acetree Intended to process data
tables where:

- each table represent one embryo
- each row represents one cell at a specific time point, with columns:
  - `time`: imaging time this cell read is created
  - `cell`: cell name
  - `x`: x position
  - `y`
  - `z`
  - `blot`: reproter expression

`LineageProcessing.R` contains functions to process one single embryo
entity, and slices from embryo entities, which are the dependencies of
other codes. INcluding functions `CD_Processing.R` are functions to work
with a directory of several embryo data tables. Including fucntions to
load data tables, retrieve specific cells and lineages data from
multiple embryos, and operations that are only possible with multiple
embryos (e.g. depth correction) `drawEmb.R` contains functions that plot
all nucleus in an embryo in 3D

## Getting started

``` r
source("LineageProcessing.R")
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
source("CD_Processing.R")
source("drawEmb.R")
```

    ## Loading required package: plotly

    ## Loading required package: ggplot2

    ## 
    ## Attaching package: 'plotly'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     last_plot

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

1.  read a directory of embryo tables

``` r
embryos <- CD_In(directory = "JIM767_25/", prefix = "CD", TIME = FALSE, AuxInfo = FALSE)
```

- Optional: implement depth correction

``` r
# get the correction parameters for the model of fluoresence loss by depth 
model <- depthCorrectionParm(embryos, lineage = c("ABara", "ABalp", "E"), alignAt = "E", startT = 30, endT = 150, exc_zMin = 0.2, exc_zMax = 1, zMax = 67)
```

    ## `summarise()` has grouped output by 'time'. You can override using the
    ## `.groups` argument.

``` r
# correct the embryos' data with the model
embryos <- embryos |> 
  dataCorrection(lineage=c("ABara", "ABalp", "P1"), model = model, zMax = 67, exc_zMin = 0.2, exc_zMax = 1)
```

2.  
