# LIVEtools

R based package to analyze and visualize embryo lineage after the cell
data are extracted from 3D movies in StarryNite/Acetree. Allowable input data 
formats are StarryNite ".zip" files or "CD" .csv files where:

- each table represent one embryo
- the first row contains column names
- each row represents one cell at a specific time point, with required columns:
  - `time`: imaging time this cell read is created
  - `cell`: cell name
  - `x`, `y`, `z`: to define nucleus positions in 3-dimensional space
  - `blot`: reporter expression

Example input files of each format are available in the sample_data folder

`LineageProcessing.R` contains functions to process one single embryo
entity, and slices from embryo entities, which are the dependencies of
other codes. Including functions `CD_Processing.R` are functions to work
with a directory of several embryo data tables. Including fucntions to
load data tables, retrieve specific cells and lineages data from
multiple embryos, and operations that are only possible with multiple
embryos (e.g. depth correction) `drawEmb.R` contains functions that plot
all nucleus in an embryo in 3D `tree_plots.R` plots lineage trees 
with customizable styles

## Before You Begin

All installation commands below are R commands that should be run in the **Console** within [RStudio](https://posit.co/downloads/) (not the "Terminal"), 
except for the optional Python setup. 

### For Tutorial
After the installation, you can follow the examples in **GettingStarted.Rmd**, you'll need the sample data included in this repository:

1. Download or clone the entire GitHub repository from https://github.com/johnmurraylab/LIVE_tools
2. Open RStudio
3. Set the downloaded repository as working directory in **Console** with R command 

     ```r
     setwd("directory/of/local/copy")
     ```

4. Open `GettingStarted.Rmd`

### For Your Own Projects (After Completing the Tutorial)
For analyzing your own embryo data, we recommend creating a dedicated project directory to keep your work organized:

```
# Example structure for a new project:
# ~/my_embryo_project/
#   ├── data/              # Your embryo tracking data files
#       ├── group1/        # Sould organize data with sub-directories
#       ├── group2/
#       └── .../
#   ├── scripts/           # Your analysis scripts
#   ├── output/            # Generated data and plots
#   └── .Rprofile          # (Optional) Project-specific R settings
```

You can set up project-specific settings (such as python paths for image export) in a `.Rprofile` file to avoid configuring them each session. This is optional and for advanced users.

## Step 1: Install R Package Dependencies

This package depends on several CRAN and Bioconductor packages.  
You can install all dependencies automatically using the following commands:

### 1.1: Install CRAN Dependencies
```r
# List of required CRAN packages
cran_deps <- c("data.table", "dplyr", "plotly", "reticulate", "ggplot2", "tidyr", "viridis")
# Check which packages are missing and install them
cran_missing <- setdiff(cran_deps, rownames(installed.packages()))
if (length(cran_missing)) {
  install.packages(cran_missing)
}
```

### 1.2 Install Bioconductor Dependencies
```r
#first, make sure Bioconductor package manager is present
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Required Bioconductor packages
bioc_deps <- c("ggtree")
#check whats missing packages and install them
bioc_missing <- setdiff(bioc_deps, rownames(installed.packages()))
if (length(bioc_missing)) {
  BiocManager::install(bioc_missing)
}
```

## Step 2: Install LIVEtools

Choose one of the two methods below to install the LIVEtools package. 

### Method A: Install from GitHub (Recommended)

This method installs the latest version directly from GitHub.

```r
# First, install the devtools package if you don't have it
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
# Install LIVEtools from GitHub
devtools::install_github("johnmurraylab/LIVE_tools")
# Load the package to verify installation
library(LIVEtools)
```

### Method B: Install from Downloaded Release
1. Go to the [Releases page](https://github.com/johnmurraylab/LIVE_tools/releases) on the GitHub repository
2. Download LIVEtools release package (with file name `LIVEtools_?.?.?.????.tar.gz` where the `?`s are version numbers) from <a href="https://github.com/johnmurraylab/LIVE_tools/releases">**Releases**</a> on github repository
3. Then run the following commands in R Console
```r
#replace package name with the actual release file name you downloaded
install.packages("full/path/to/tar.gz", repos = NULL, type = "source")
```

## Step 3 (Optional): Setup Python for Image Export

**This step is only required if you want to use the `saveEmbImg()` function** to 
export 3D plots as static images (PNG files) in a **fully scripted workflow** 
(R markdown provides buttons to manually capture static snapshots of 3D plots). 
All other package functionality works without Python.

The `saveEmbImg()` function depends on `plotly::save_image`, which uses the *reticulate* package to call Python's *kaleido* module and convert plotly plots (HTML widgets) to static images.

### Requirements
- Python 3.8, 3.9, or 3.10 (recommended for compatibility with kaleido 0.2.1)
- Python packages: plotly, numpy, kaleido==0.2.1

### 3.1 Create a Python Virtual Environment
We recommend using a Python virtual environment to keep these dependencies isolated. This is because LIVEtools requires an older version of *kaleido* (0.2.1) that may conflict with other Python projects. 
The guide will provide the most lightweight option to manage python: [python venv](https://docs.python.org/3/library/venv.html). But *reticulate* do allow other methods (like conda) for users more proficient with python. 

**In your system Terminal** (not R Console):

1. Navigate to a directory where you would like to put your python venv directory.
2. Create a new virtual environment. 
(For this example we are creating one called `livetools_env` under home directory `~`. 
You can choose other directories and names, just remember to refer to the name you choose in the following steps)

    ```sh
    python -m venv livetools_env
    ```
   This would create a directory under the current working directory with structure:
   
    ```
    # ./livetools_env/
    #    ├── bin/
    #        ├── python        #python executable, where python codes can be executed
    #        ├── activate      #script that can activate the virtual envrionment
    #        └── .......
    #    └── ......            #Other folders and files
    ```
3. Activate the virtual environment to make changes to it. 
  On macOS/Linux:
  
    ```sh
    source livetools_env/bin/activate
    ```
  On Windows:
  
    ```powershell
    livetools_env\Scripts\activate
    ```
4. Install required Python packages

    ```bash
    pip install numpy plotly kaleido==0.2.1
    # You can now deactivate the environment
    deactivate
    ```

### 3.2 Configure Reticulate to Use Your Python Environment

After creating the Python environment, you need to tell R's *reticulate* package which Python installation to use.

**Option A: Set permanently in .Rprofile (Recommended for regular use)**

Add this line to your `.Rprofile` file (put it in your [project direcotry](#for-your-own-projects-after-completing-the-tutorial)) or modify your home directory .Rprofile):

```r
Sys.setenv(RETICULATE_PYTHON = "~/livetools_env/bin/python")
```

Replace the path with the actual path to your virtual environment's python executable (see [3.1](#31-create-a-python-virtual-environment) for venv).

**Option B: Set each time before using saveEmbImg (For occasional use)**

Each time you want to use `saveEmbImg()`, run these commands in R Console **before** loading LIVEtools:
```r
# Load reticulate
library("reticulate")

# Specify the python virtual environment
# Replace the path with your actual venv path
reticulate::use_virtualenv("~/livetools_env", required = TRUE)

# Alternatively, if using conda:
# reticulate::use_condaenv("my_env", required = TRUE)

# Activate python kernel to ensure consistent behavior
reticulate::py_run_string("import sys")

# Now load LIVEtools
library(LIVEtools)
```

### 3.3 Verify Python Setup

After completing the python configuration, test that the required modules are available:

```r
library(LIVEtools)
library(reticulate)

# Check that python modules are available
reticulate::py_module_available("kaleido")  # Should return TRUE
reticulate::py_module_available("plotly")   # Should return TRUE
reticulate::py_module_available("numpy")    # Should return TRUE

# If all return TRUE, the setup is successful!
```

### Troubleshooting Reticulate

**Problem: Wrong Python being used**

```r
reticulate::py_config()
# Shows unexpected Python path
```

**Solution:**

Clear reticulate's Python cache: 
```r
Sys.setenv(RETICULATE_PYTHON = "")
```
Restart R session. Then explicitly set the correct path before loading any packages:
```r
library(reticulate)
reticulate::use_virtualenv("~/livetools_env", required = TRUE)
reticulate::py_config()  # Verify correct Python is now being used
```

**Problem: Architecture mismatch**

Some mac users might have an x86(Intel) architecture Python despite the computer having ARM(Apple silicon) processors, which might cause issue when using venv. 
```
# Error when running reticulate::py_run_string
Error:
......
...incompatible architecture...
```

**Solution:** Install Compatible Python Specifically

*Perform these steps in Terminal*: 
Use ARM version of `brew` (not Intel) to install a stable mac architecture python
```sh
/opt/homebrew/bin/brew install python@3.10
```
Then remove your old virtual environment,
```sh
rm -rf ~/livetools_env
```
and create a new venv.
```sh
/opt/homebrew/bin/python3.10 -m venv ~/livetools_env #specifically call ARM Python
```
Finally check the architecture of the new venv
```sh
~/livetools_env/bin/python -c "import platform; print(platform.machine())"
```

## Verify Installation

Test that LIVEtools is properly installed and working:

```r
library(LIVEtools)

# Check package version
packageVersion("LIVEtools")

# Open help documentation for a main function
?readEmbryoTable

# If these commands work without errors, installation was successful!
```

## Modules
### LineageProcessing
functions to process one single embryo entity, and slices from embryo entities, which are the dependencies of other codes. 

### CD_Processing
functions to work with a directory of several embryo data tables. Including functions to load data tables, retrieve specific cells and lineages data from multiple embryos, and operations that are only possible with multiple embryos (e.g. depth correction)

### drawEmb
functions that plot all nucleus in an embryo in 3D

### tree_plots
This package requires a non-CRAN package [ggtree](https://doi.org/doi:10.18129/B9.bioc.ggtree) that won't be automatically installed
To install ggtree (offered by bioconductor), see section [Install Bioconductor Dependencies](#install-bioconductor-dependencies).
