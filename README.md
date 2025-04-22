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
other codes. INcluding functions `CD_Processing.R` are functions to work
with a directory of several embryo data tables. Including fucntions to
load data tables, retrieve specific cells and lineages data from
multiple embryos, and operations that are only possible with multiple
embryos (e.g. depth correction) `drawEmb.R` contains functions that plot
all nucleus in an embryo in 3D