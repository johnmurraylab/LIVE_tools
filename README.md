<div class="container-fluid main-container">

<div id="header">

</div>

<div id="embryo_lineage_analysis_r" class="section level1">

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
other codes. Including functions `CD_Processing.R` are functions to work
with a directory of several embryo data tables. Including fucntions to
load data tables, retrieve specific cells and lineages data from
multiple embryos, and operations that are only possible with multiple
embryos (e.g. depth correction) `drawEmb.R` contains functions that plot
all nucleus in an embryo in 3D `tree_plots.R` plots lineage trees 
with customizable styles

## Installation
Please ensure the package and its dependencies are installed **BEFORE** installing LIVE_tools

### 1. Install Dependencies

This package depends on several CRAN and Bioconductor packages.  
You can install all dependencies automatically using the following commands:

#### Install CRAN Dependencies
```r
cran_deps <- c(
  "data.table", "dplyr", "plotly", "reticulate", "ggplot2", "tidyr", "viridis" 
) # <-- replace with actual CRAN dependencies from DESCRIPTION
cran_missing <- setdiff(cran_deps, rownames(installed.packages()))
if (length(cran_missing)) {
  install.packages(cran_missing)
}
```

#### Install Bioconductor Dependencies
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") #first, make sure Bioconductor package manager is present

bioc_deps <- c(
  "ggtree" #Bioconductor packages to install
)
bioc_missing <- setdiff(bioc_deps, rownames(installed.packages())) #check missing and install
if (length(bioc_missing)) {
  BiocManager::install(bioc_missing)
}
```

#### Setup Python Dependencies
The `saveEmbImg` function depends on `plotly::save_image`, which further utilizes reticulate to call Python kaleido module and convert plotly plots (html widgets) to static images. 
The best practice is to setup a python environment specifically for this functionality because R plotly can only utilize a legacy version of kaleido (v0.1.0). 


Install python dependencies for the environment
```sh
# assume that we have an active python venv environment (you can also use conda for this purpose)
pip install numpy plotly kaleido==0.1.0 # install dependencies
```
The, you shall inform reticualte which python installation to use (see section [Configuring Python](#configuring-python), otherwise reticualte might install a cache python that do not have the proper modules.

### Install LIVEtools

Simplest way is to use install_github (assumes you have the devtools library installed)
```r
devtools::install_github("johnmurraylab/LIVE_tools")
library(LIVEtools)
```

### Alternative local download install:
Download LIVEtools release package (with file name `LIVEtools_?.?.?.????.tar.gz` where the `?`s are version numbers) from <a href="https://github.com/johnmurraylab/LIVE_tools/releases">**Releases**</a> on github repository

Then run the following command
```r
#replace package name with the actual release file name you downloaded
install.packages("full/path/to/tar.gz/file", repos = NULL, type = "source")
```

## Modules
### LineageProcessing
functions to process one single embryo entity, and slices from embryo entities, which are the dependencies of other codes. Including functions

### CD_Processing
functions to work with a directory of several embryo data tables. Including fucntions to load data tables, retrieve specific cells and lineages data from multiple embryos, and operations that are only possible with multiple embryos (e.g. depth correction)

### drawEmb
functions that plot all nucleus in an embryo in 3D

#### Configuring Python
The `saveEmbImg` function uses the [reticulate](https://rstudio.github.io/reticulate/) package
to interface with a Python instance with "kaleido" and "plotly" installed. 

Before loading the package, you should tell reticulate which Python to use. 

You can either set a permanent reticulate python installation in .Rprofile: 
```r
Sys.setenv(RETICULATE_PYTHON = "/path/to/python")
```
Or *EACH TIME BEFORE* you load LIVEtools and use saveEmbImg, you have to configure reticulate (in cause you might work with different python installations through reticualte):

1. load reticulate 
```r
library("reticulate")
```
2. specify the python environment
```r
# Use a conda environment
reticulate::use_condaenv("my_env", required = TRUE)
# Or use a python vnev
reticulate::use_virtualenv("my_venv")
```
3. Activate environment and load LIVEtools module in R
```r
reticulate::py_run_string("import sys") # have to so this to ensure consistent behavior
library(LIVEtools) # only load LIVEtools after everything above is executed
```

### tree_plots
This package requires a non-CRAN package [ggtree](https://doi.org/doi:10.18129/B9.bioc.ggtree) that won't be automatically installed
To install ggtree (offered by bioconductor), see section [Install Bioconductor Dependencies](#install-bioconductor-dependencies).
