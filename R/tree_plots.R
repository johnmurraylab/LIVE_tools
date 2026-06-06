#' plot_trees
#' @description
#' tree plots for a list of embryo cell tracking data
#' generate one file for each CD file in the given directory, with consistent parameters
#'
#' @details
#' This function uses [ggtree] (Bioconductor) to plot trees.
#' #' The `ggtree` package must be installed separately:
#' ```r
#' if (!requireNamespace("BiocManager", quietly = TRUE))
#'     install.packages("BiocManager")
#' BiocManager::install("ggtree")
#' ```
#' @param inDir directory of embryo cell tracking files to plot into trees
#' @param root root cell of the tree
#' @param plotExpression whether to color each branch by the average expression of the cell
#' @param lowerBound lower bound of each cell's expression to count into average expression value
#' @param upperBound upper bound of each cell's expression to count into average expression value
#' @param min_gain for color scale, the value corresponding to "minimum" color
#' @param max_gain for color scale, the value corresponding to "maximum" color
#' @param transform color scale transformation to apply
#' @param start_time only consider data before `start_time`
#' @param end_time only consider data before `end_time`
#' @param tip_lab Boolean, whether to label the tree "leaves" with cell names
#' @param node_lab Boolean, whether branch poins with cell names
#' @param outDir output directory
#' @param exp_col name of the column that correspond to expression attributes
#' @param cell_col name of the column that correspond to cell names attributes
#' @param time_col name of the column that correspond to cell time attributes
#' @param width in inch, figure width
#' @param height in inch, figure height
#' @param dpi dpi of output figure
#'
#' @return `NA`
#' @export
#' @import ggplot2
plot_trees <- function(inDir, root = "P0",
                       plotExpression = T, lowerBound = 0.4, upperBound = 0.97,
                       min_gain=100, max_gain = 10000, transform = "identity",
                       start_time = FALSE, end_time=FALSE,
                       tip_lab=TRUE,node_lab=FALSE, outDir = NaN,
                       exp_col="blot", cell_col="cell", time_col="time",
                       width = 7, height = 7, dpi = 300){
  if(is.na(outDir)){outDir = inDir}
  files <- list.files(inDir, pattern=".+\\.csv$")
  embNames <- gsub("(.*)\\.csv$", "\\1", files)
  for (i in seq_along(files)) {
    directory <- paste0(inDir,files[i])
    CD <- read.csv(directory, header=T)
    tree <- CD_tree_plot(CD, root = root, transform = transform,
                         start_time = start_time, end_time = end_time,
                         exp_col=exp_col, cell_col=cell_col, time_col=time_col,
                         plotExpression = plotExpression, min_gain=min_gain,max_gain=max_gain,
                         tip_lab=tip_lab,node_lab=node_lab)
    ggplot2::ggsave(filename = paste0(outDir, embNames[i], root,".png"), plot = tree,
           width = width, height = height, dpi = dpi, units = "in", device='png')
  }
}

GetSister <- function(x){
    if(is.na(x)) return(NA)
    if(substr(x,nchar(x),nchar(x))=="a"){
        return(paste(substr(x,1,nchar(x)-1),"p",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="d"){
        return(paste(substr(x,1,nchar(x)-1),"v",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="l"){
        return(paste(substr(x,1,nchar(x)-1),"r",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="p"){
        return(paste(substr(x,1,nchar(x)-1),"a",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="v"){
        return(paste(substr(x,1,nchar(x)-1),"d",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="r"){
        return(paste(substr(x,1,nchar(x)-1),"l",sep=""))
        }



    if(x=="AB") return("P1")
    if(x=="EMS") return("P2")
    if(x=="C") return("P3")
    if(x=="D") return("P4")
    if(x=="Z2") return("Z3")
    if(x=="E") return("MS")

    if(x=="P1") return("AB")
    if(x=="P2") return("EMS")
    if(x=="P3") return("C")
    if(x=="P4") return("D")
    if(x=="Z3") return("Z2")
    if(x=="MS") return("E")

    return(NA)
}

GetParent <- function(x){

    if(is.na(x) || x=="P") return(NA)
    if(is.element(substr(x,nchar(x),nchar(x)),c("a","p","d","v","l","r"))) return (substr(x,1,nchar(x)-1))
    if(x=="AB" || x=="P1") return("P0")
    if(x=="EMS" || x=="P2") return("P1")
    if(x=="E" || x=="MS") return("EMS")
    if(x=="C" || x=="P3") return("P2")
    if(x=="D" || x=="P4") return("P3")
    if(x=="Z2" || x=="Z3") return("P4")
    return("P")
}

##recursively create Newick format string from complete data frame (columns:
MakeNewick <- function(data, root="P0", outer=TRUE,verbose=FALSE, tiplab = TRUE, nodelab=FALSE){
    if(verbose){
        message(root)
    }
    if(data[root,"hasDaughter"]){
        daughters <- sort(rownames(data[data[,"parents"]==root,]))
        if(outer){
            return(paste("((",MakeNewick(data,daughters[1],outer=FALSE),",",MakeNewick(data,daughters[2],outer=FALSE),")",root,":",data[root,"length"],");", sep=""))
        }else{
            return(paste("(",MakeNewick(data,daughters[1],outer=FALSE),",",MakeNewick(data,daughters[2],outer=FALSE),")",root,":",data[root,"length"], sep=""))
        }

        ##        message(paste("(",MakeNewick(data,daughters[1]),",",MakeNewick(data,daughters[2]),")",root,data[root,"length"], sep=""))
    }else{
        return(paste(root,data[root,"length"],sep=":"))
        ##        message(paste(root,data[root,"length"],sep=":"))
    }
}

#' CD_tree_plot
#'
#' @description
#' tree plots for a list of embryo cell tracking data
#'
#' @param CD dataframe for the CD table
#' @param root root cell of the tree
#' @param min_gain for color scale, the value corresponding to "minimum" color
#' @param max_gain for color scale, the value corresponding to "maximum" color
#' @param transform color scale transformation to apply
#' @param start_time only consider data before `start_time`
#' @param end_time only consider data before `end_time`
#' @param plotExpression whether to color each branch by the average expression of the cell
#' @param lowerBound lower bound of each cell's expression to count into average expression value
#' @param upperBound upper bound of each cell's expression to count into average expression value
#' @param tip_lab Boolean, whether to label the tree "leaves" with cell names
#' @param node_lab Boolean, whether branch poins with cell names
#' @param exp_col name of the column that correspond to expression attributes
#' @param cell_col name of the column that correspond to cell names attributes
#' @param time_col name of the column that correspond to cell time attributes
#' @param extra_branch_l length of dummy branches to connect trees to root P0 or P if such cell were not tracked
#' @param exp_legend name of the color scale
#'
#' @return ggplot figure of the tree plot
#' @export
#' @import ggtree
#' @import ggplot2
#' @import tidyr
CD_tree_plot <- function(CD, root="P0",
                         min_gain=-1000, max_gain=FALSE, transform = "identity",
                         start_time = FALSE, end_time=FALSE,
                         plotExpression=TRUE, lowerBound = 0.4, upperBound = 0.97,
                         tip_lab=TRUE,node_lab=FALSE,
                         exp_col="blot", cell_col="cell", time_col="time",
                         extra_branch_l=15, exp_legend="Expression"){

  if (!requireNamespace("ggtree", quietly = TRUE)) {
    stop(
      "Package 'ggtree' is required for this function. It is a non-CRAN package.
      Please install it with bioconductor:
      if (!require('BiocManager', quietly = TRUE))
        install.packages('BiocManager')
      BiocManager::install('ggtree')",
      call. = FALSE
    )
  }

  if(is.numeric(start_time)){
    CD <- CD[CD$time >= start_time,]
  }
  if(is.numeric(end_time)){
    CD <- CD[CD$time <= end_time,]
    X_Max <- (end_time-start_time)+15
  }
  else{X_Max <- NA}

  if(cell_col != "cell"){
      colnames(CD)[colnames(CD) == cell_col] <- "cell"
  }
  if(time_col != "time"){
      colnames(CD)[colnames(CD) == time_col] <- "time"
  }
  CD <- CD|>drop_na(time)

  CD <- CD[order(CD$cell),]

  cells<-unique(CD$cell)
  birthTime <- aggregate(CD$time,by=list(CD$cell), FUN=min)[,2]
  names(birthTime) <- cells
  endTime <- aggregate(CD$time,by=list(CD$cell), FUN=max)[,2]
  names(endTime) <- cells
  length <- (endTime-birthTime) + 1
  parents <- sapply(cells,GetParent)
  hasParent <- cells %in% cells[parents %in% cells]
  names(hasParent) <- cells

  ##debugging code
  ##message(paste("Eal", birthTime["Eal"], endTime["Eal"], length["Eal"], "Ear", birthTime["Ear"], endTime["Ear"],length["Ear"],sep=" "))

  ##add dummy branches to connect trees to root P0 or P. Each gets the same minute/tp length
  while(sum(!hasParent & !(cells %in% c("P0", "P"))) > 0){
      message("ADDING PARENTS FOR:", sum(!hasParent), "CELLS", sep=" ")
      message(paste(cells[!hasParent], collapse=" "))
      newCells <- unique(parents[cells[!hasParent]])
      cells <- c(cells, newCells)
      length <- c(length, rep(extra_branch_l,length(newCells)))
      parents <- c(parents, sapply(newCells,GetParent))
      names(length) <- cells
      names(parents) <- cells
      hasParent <- cells %in% cells[parents %in% cells]
      names(hasParent) <- cells
  }

  hasDaughter <- cells %in% parents
  names(hasDaughter) <- cells

  ##QC
  sisters <- sapply(cells,GetSister)
  hasSister <- cells %in% sisters
  message(paste(sum(!hasSister),"cells don't have sisters - should equal one for P0"))

  tree_table <- data.frame(parents,sisters,length,hasDaughter)
  newick_tree <- MakeNewick(tree_table, root=root)
##   message(newick_tree)
  ## load into ggtree
  tree <- read.tree(text=newick_tree)

##    message(exp_col)
##    message(colnames(CD))
  if(plotExpression){
      funAvg <- function(X, lower, upper){
        lowerV <- quantile(X, lower)
        upperV <- quantile(X, upper)
        Xfiltered <- X[X>=lowerV[[1]]]
        Xfiltered <- Xfiltered[Xfiltered<=upperV[[1]]]
        return(mean(Xfiltered))
      }
      blot <- aggregate(CD[,exp_col], by=list(CD$cell), FUN=funAvg, lower=lowerBound, upper = upperBound)
      rownames(blot) <- blot[,1]
      blot[is.na(blot[,2]), 2] <- min_gain
      ##blot <- blot[-1]
      ##fixme consider NA
      blot[setdiff(union(tree$tip.label, tree$node.label),rownames(blot)),2] <- 0
      blot[,1] <- rownames(blot)
      colnames(blot) <- c("label", "blot")

      if(is.numeric(min_gain)){blot[blot[,2] < min_gain,2] <- min_gain}
      if(is.numeric(max_gain)){blot[blot[,2] > max_gain,2] <- max_gain
        message("adjusting max gain")}
      else{max_gain <- max(blot[,2])}

      g <- ggtree::ggtree(tree, ladderize=F) %<+% blot
      g <- g + xlim(0, X_Max) + ggtree::geom_tree(aes(color=blot))

      #define color scaling and labels according to transform argument
      if(transform == "identity"){
        newMin <- min_gain
        newMax <- max_gain
        midpoint <- (newMin+newMax)/2
      }
      else if(transform == "log2"){
        newMin <- log2(min_gain)
        newMax <- log2(max_gain)
        midpoint <- 2**((newMin+newMax)/2)
      }
      else if(transform == "log10"){
        newMin <- log10(min_gain)
        newMax <- log10(max_gain)
        midpoint <- 10**((newMin+newMax)/2)
      }
      else if(transform == "log"){
        newMin <- log(min_gain)
        newMax <- log(max_gain)
        midpoint <- exp((newMin+newMax)/2)
      }

      g <- g + ggplot2::scale_color_gradient2(name=exp_legend,
                                     trans = transform,
                                     low=rgb(0,0,0.9),mid=rgb(0.9,0.85,0.8),high="red",
                                     midpoint = midpoint, limits = c(min_gain,max_gain))
  }else{
      g <- ggtree::ggtree(tree, ladderize=F)+ xlim(0, X_Max)
  }

  if(tip_lab){
      g <- g + ggtree::xlim(0, X_Max*1.1) + ggtree::geom_tiplab()
  }
  if(node_lab){
      g <- g + ggtree::geom_nodelab()
  }
  g <- g + ggtree::theme_tree2()
  return(g)
}


# =============================================================================
# Continuous-resolution lineage tree plotting
# =============================================================================
#
# `plot_lineage_tree()` is a more general-purpose tree plotter than
# `CD_tree_plot()`. It supports two resolutions of value-coloring:
#
#   resolution = "cell"        -> one color per cell branch (legacy-style)
#   resolution = "timepoint"   -> per-timepoint segment coloring (continuous
#                                 along the branch; ported from the
#                                 Antigravity-era MS expression-tree code)
#
# Both modes accept the same long-format input (cell, time, value), with
# `cell_col`, `time_col`, and `value_col` configurable. Cells/timepoints
# whose value is `NA` render in a configurable neutral color (`na_color`)
# so missing data is visible without skewing the color scale.
#
# Other features ported from the MS-tree improvements:
#   * cross-division recursive smoothing (`smoothing_w`)
#   * dynamic end-time milestone (e.g. "stop when MS reaches 32 cells")
#   * percentile-based color limits, optionally scoped to a sublineage
#   * value floor/ceiling caps on the color scale
#
# Designed as a standalone plotter (no defect-pipeline coupling) — works for
# expression, deviation scores, or any per-cell-per-time numeric metric.

# Internal: cross-division recursive smoothing helper.
# Returns the value at `target_time` for `current_cell`, recursing to the
# parent if before birth or to the daughters if after the cell divides.
# All of `cell_value_map`, `parents_map`, `cell_birth`, `cell_end` are
# precomputed look-ups indexed by cell name; `daughters_map` is indexed
# by parent name.
.value_at_time <- function(current_cell, target_time, cell_value_map,
                           parents_map, daughters_map, cell_birth, cell_end) {
    val <- cell_value_map[[paste(current_cell, target_time, sep = "@")]]
    if (!is.null(val) && !is.na(val)) return(val)

    birth <- cell_birth[[current_cell]]
    end   <- cell_end[[current_cell]]
    if (is.null(birth) || is.null(end)) return(NA_real_)

    if (target_time < birth) {
        p <- parents_map[[current_cell]]
        if (is.null(p) || is.na(p)) return(NA_real_)
        return(.value_at_time(p, target_time, cell_value_map, parents_map,
                              daughters_map, cell_birth, cell_end))
    }
    if (target_time > end) {
        ds <- daughters_map[[current_cell]]
        if (is.null(ds) || length(ds) == 0) return(NA_real_)
        vals <- vapply(ds, function(d) {
            .value_at_time(d, target_time, cell_value_map, parents_map,
                           daughters_map, cell_birth, cell_end)
        }, numeric(1))
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0) return(mean(vals))
        return(NA_real_)
    }
    NA_real_
}

# Default biological-convention daughter-ordering exceptions. Sulston-style
# convention puts these daughter pairs in *non*-alphabetical order; without
# overrides, MakeNewick's `sort()` would render them alphabetically. Each
# entry maps `parent_cell -> c(bottom_daughter, top_daughter)` — the order
# they should appear in Newick (since ggtree assigns y=1 to the first leaf
# encountered in Newick = the bottom of the rendered tree).
.SULSTON_DAUGHTER_ORDER_OVERRIDES <- list(
    # EMS: E conventionally drawn above MS in Sulston lineage diagrams,
    # so MS comes first in Newick (bottom), E second (top). Within each
    # subtree, daughters fall back to alphabetical (Sulston-default).
    EMS = c("MS", "E")
)

# Internal: Newick builder that supports per-parent daughter-order overrides.
# Replaces `LIVEtools::MakeNewick` for our layout. At every internal node
# whose name appears in `overrides`, the daughters are emitted in the
# specified order; all other nodes fall back to alphabetical sort (matching
# `MakeNewick`'s behavior).
.MakeNewickWithOverrides <- function(data, root, overrides = list(), outer = TRUE) {
    if (!root %in% rownames(data)) {
        # Cell missing from tree_table — emit a leaf with no length.
        return(paste(root, "0", sep = ":"))
    }
    if (isTRUE(data[root, "hasDaughter"])) {
        raw <- rownames(data[data[, "parents"] == root & !is.na(data[, "parents"]), ])
        if (root %in% names(overrides)) {
            want <- overrides[[root]]
            if (length(want) == 2 && all(want %in% raw)) {
                daughters <- want
            } else {
                daughters <- sort(raw)
            }
        } else {
            daughters <- sort(raw)
        }
        d1 <- .MakeNewickWithOverrides(data, daughters[1], overrides, outer = FALSE)
        d2 <- .MakeNewickWithOverrides(data, daughters[2], overrides, outer = FALSE)
        if (outer) {
            paste("((", d1, ",", d2, ")", root, ":", data[root, "length"], ");", sep = "")
        } else {
            paste("(", d1, ",", d2, ")", root, ":", data[root, "length"], sep = "")
        }
    } else {
        paste(root, data[root, "length"], sep = ":")
    }
}

# Internal: build the ggtree-derived layout (tree_y per cell) plus
# division-connector and tip-label tables. Returns a list of frames.
#
# `daughter_order_overrides` is a named list mapping parent cell ->
# `c(bottom_daughter, top_daughter)`. Default is the Sulston-convention
# E/MS exception. Use `list()` to disable all overrides.
.build_tree_layout <- function(CD, root = "P0", extra_branch_l = 15,
                               daughter_order_overrides = .SULSTON_DAUGHTER_ORDER_OVERRIDES) {
    cells <- unique(CD$cell)
    birthTime <- aggregate(CD$time, by = list(CD$cell), FUN = min)[, 2]
    names(birthTime) <- cells
    endTime <- aggregate(CD$time, by = list(CD$cell), FUN = max)[, 2]
    names(endTime) <- cells
    branch_len <- (endTime - birthTime) + 1
    parents <- sapply(cells, GetParent)
    names(parents) <- cells
    hasParent <- cells %in% cells[parents %in% cells]
    names(hasParent) <- cells

    # add dummy ancestors up to root
    while (sum(!hasParent & !(cells %in% c("P0", "P"))) > 0) {
        newCells <- unique(parents[cells[!hasParent]])
        cells <- c(cells, newCells)
        branch_len <- c(branch_len, rep(extra_branch_l, length(newCells)))
        parents <- c(parents, sapply(newCells, GetParent))
        names(branch_len) <- cells
        names(parents) <- cells
        hasParent <- cells %in% cells[parents %in% cells]
        names(hasParent) <- cells
    }
    hasDaughter <- cells %in% parents
    names(hasDaughter) <- cells
    sisters <- sapply(cells, GetSister)

    tree_table <- data.frame(parents, sisters,
                             length = branch_len, hasDaughter)
    # Use our custom Newick builder so daughter ordering can be overridden
    # at specific parents (e.g. EMS) while preserving alphabetical order
    # at all other nodes.
    newick_tree <- .MakeNewickWithOverrides(
        tree_table, root = root,
        overrides = daughter_order_overrides %||% list())
    tree <- ape::read.tree(text = newick_tree)
    g <- ggtree::ggtree(tree, ladderize = FALSE)

    coords <- g$data
    coords <- coords[!is.na(coords$label), c("label", "y")]
    colnames(coords) <- c("cell", "tree_y")

    list(coords = coords, parents = parents,
         birthTime = birthTime, endTime = endTime)
}

#' Lineage tree colored by an arbitrary numeric metric
#'
#' `plot_lineage_tree()` draws a *C. elegans* Sulston lineage tree from a
#' long-format `(cell, time, value)` data frame and colors each branch by
#' the value column. Unlike [CD_tree_plot()], it supports two resolutions:
#'
#' - **`resolution = "cell"`** — one color per cell branch, where the
#'   per-cell value is collapsed from the input via `cell_aggregate`
#'   (`"mean"`, `"max"`, `"first"`, or `"last"`). Best for one-value-per-cell
#'   metrics like cell-cycle deviation or division-orientation scores.
#' - **`resolution = "timepoint"`** — per-timepoint `geom_segment` coloring,
#'   producing a continuous color gradient along each branch. Best for
#'   per-timepoint metrics like position deviation or expression
#'   trajectories. Optional cross-division smoothing via `smoothing_w`.
#'
#' Cells/timepoints whose value is `NA` are drawn in `na_color` (default
#' light grey) so missing data is visible without skewing the color scale.
#'
#' This function is the generalization of the per-timepoint expression-tree
#' code that was developed for the MS lineage trees during the LineagePhenotyping
#' refactor. See `MS_Tree_Improvements_Summary.md` for the rationale behind the
#' continuous-tracking, smoothing, and percentile-color-scale features.
#'
#' @param CD long-format data frame; must contain columns named by
#'   `cell_col`, `time_col`, and `value_col`. Other columns are ignored.
#' @param root root cell of the tree (default `"P0"`).
#' @param value_col,cell_col,time_col column names to use from `CD`.
#' @param resolution `"cell"` or `"timepoint"`.
#' @param cell_aggregate function used to collapse per-timepoint values to
#'   one value per cell when `resolution = "cell"`. One of `"mean"`, `"max"`,
#'   `"first"`, `"last"`.
#' @param smoothing_w window half-width for cross-division smoothing
#'   (only used when `resolution = "timepoint"`). 0 = no smoothing
#'   (default — appropriate for already-aggregated metrics like deviations).
#'   2 reproduces the Antigravity expression-trajectory smoothing.
#' @param lineage_subset optional sublineage to restrict the plot to (e.g.
#'   `"P1"`); passed through [grepCells()].
#' @param end_time numeric hard cap on the time axis. `NULL` = use the data's
#'   max time.
#' @param end_time_milestone alternative to `end_time`: a list with elements
#'   `lineage` (e.g. `"MS"`) and `n_cells` (e.g. `32`); the tree stops at the
#'   first time when the named lineage has at least `n_cells` distinct cells
#'   alive. If both are given, `end_time_milestone` wins.
#' @param color_low,color_high endpoints of the 2-color gradient. Used only
#'   when `colors` is `NULL`.
#' @param colors optional character vector of colors for a multi-stop
#'   gradient via [ggplot2::scale_color_gradientn()]. Useful for
#'   highlighting extremes (e.g. `c("blue","white","yellow","red")`). When
#'   set, `color_low`/`color_high` are ignored. `NULL` (default) uses the
#'   2-color gradient.
#' @param percentile_bounds two numerics in `[0, 1]` giving the lower and
#'   upper quantiles used to set the color scale limits when `value_min`
#'   / `value_max` are not supplied. `c(0.15, 0.85)` matches the MS-tree
#'   default.
#' @param value_min,value_max explicit numeric limits for the color scale.
#'   When non-`NULL`, override `percentile_bounds`. `value_max` accepts the
#'   sentinel `"data_max"` meaning "use the observed maximum but no less
#'   than `value_min + 1`"; analogously `value_min = "data_min"`. Use for
#'   defect plots where you want a fixed "no-deviation floor" (e.g. 3 µm
#'   for position deviation) and the upper end to follow the data.
#' @param value_floor,value_ceiling optional hard caps applied to the
#'   percentile-derived limits (legacy; ignored when `value_min` /
#'   `value_max` are set). Useful to prevent amplifying noise in
#'   low-signal channels (set `value_floor` to keep the upper end of the
#'   gradient above some threshold).
#' @param color_values numeric vector of breakpoints in the same units as
#'   `value_col`, corresponding to the entries of `colors`. When supplied
#'   together with `colors` it produces a piecewise-stop gradient via
#'   [ggplot2::scale_color_gradientn()] (`values = scales::rescale(...)`).
#'   Use for diverging or asymmetric scales — e.g. for CC deviation with
#'   colors `c("green","darkblue","black","magenta4","orange")` and
#'   `color_values = c(min, -5, 0, 5, max)` to flatten the ±5 neutral
#'   zone and emphasize larger excursions.
#' @param percentile_scope optional sublineage name; if given, the
#'   percentile bounds are computed from this subset only (mirrors the
#'   MS-tree P1-scoped scaling). `NULL` = use all values.
#' @param na_color color used for `NA` values (default `"grey80"`).
#' @param title plot title.
#' @param tip_lab,node_lab whether to label tree leaves and internal nodes.
#' @param truncate_to_last_data when `TRUE` (default), per cell, drops CD
#'   rows whose `time` exceeds the latest time at which the cell has a
#'   non-NA value. Trims trailing-NA "grey tails" produced by cells whose
#'   data ends before the global `end_time`.
#' @param drop_empty_cells when `TRUE`, iteratively prunes cells whose
#'   value is `NA` across every observed timepoint (no data at all).
#'   Useful for division-time / cell-cycle plots where terminal leaves
#'   without divisions have no measurement and would otherwise render
#'   as empty `na_color` branches. Default `FALSE` to preserve the
#'   classic "every observed cell drawn" behavior.
#' @param branch_width linewidth for cell-branch segments (default 3).
#' @param extra_branch_l branch length for dummy ancestor cells inserted
#'   to connect orphan tips to `root`.
#' @param value_legend label shown above the color legend.
#'
#' @return a `ggplot` object.
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
plot_lineage_tree <- function(CD, root = "P0",
                              value_col = "blot",
                              cell_col = "cell", time_col = "time",
                              resolution = c("cell", "timepoint"),
                              cell_aggregate = c("mean", "max", "first", "last"),
                              smoothing_w = 0,
                              lineage_subset = NULL,
                              end_time = NULL,
                              end_time_milestone = NULL,
                              color_low = "blue", color_high = "red",
                              colors = NULL,
                              color_values = NULL,
                              percentile_bounds = c(0.15, 0.85),
                              value_min = NULL, value_max = NULL,
                              value_floor = NULL, value_ceiling = NULL,
                              percentile_scope = NULL,
                              na_color = "grey80",
                              title = NULL,
                              tip_lab = TRUE, node_lab = FALSE,
                              tip_size = NULL,
                              plot_height_in = 10,
                              tip_size_min = 0.4,
                              tip_size_max = 4,
                              truncate_to_last_data = TRUE,
                              drop_empty_cells = FALSE,
                              branch_width = 3,
                              extra_branch_l = 15,
                              value_legend = NULL) {
    # ---- validation & required-package checks ----
    if (!requireNamespace("ggtree", quietly = TRUE)) {
        stop("Package 'ggtree' is required. Install with BiocManager::install('ggtree').")
    }
    if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required. Install with install.packages('ape').")
    }
    if (!requireNamespace("scales", quietly = TRUE)) {
        stop("Package 'scales' is required. Install with install.packages('scales').")
    }

    resolution <- match.arg(resolution)
    cell_aggregate <- match.arg(cell_aggregate)
    if (length(percentile_bounds) != 2 ||
        any(percentile_bounds < 0) || any(percentile_bounds > 1) ||
        percentile_bounds[1] >= percentile_bounds[2]) {
        stop("percentile_bounds must be two numerics in [0,1] with lo < hi")
    }
    if (is.null(value_legend)) value_legend <- value_col

    # ---- normalize column names ----
    if (!(cell_col  %in% colnames(CD))) stop("cell_col '",  cell_col,  "' not in CD")
    if (!(time_col  %in% colnames(CD))) stop("time_col '",  time_col,  "' not in CD")
    if (!(value_col %in% colnames(CD))) stop("value_col '", value_col, "' not in CD")

    if (cell_col  != "cell")  colnames(CD)[colnames(CD) == cell_col]  <- "cell"
    if (time_col  != "time")  colnames(CD)[colnames(CD) == time_col]  <- "time"
    if (value_col != "value") colnames(CD)[colnames(CD) == value_col] <- "value"

    CD <- CD |> drop_na(time)
    CD <- CD[order(CD$cell, CD$time), ]

    # ---- optional lineage subset ----
    if (!is.null(lineage_subset)) {
        CD <- grepCells(CD, lineages = lineage_subset, times = "ALL")
    }

    # ---- resolve end_time ----
    if (!is.null(end_time_milestone)) {
        lin    <- end_time_milestone$lineage
        n_min  <- end_time_milestone$n_cells
        if (is.null(lin) || is.null(n_min)) {
            stop("end_time_milestone must be a list with `lineage` and `n_cells`")
        }
        re <- paste0("^", lin)
        cells_per_t <- CD %>%
            dplyr::filter(grepl(re, .data$cell)) %>%
            dplyr::group_by(.data$time) %>%
            dplyr::summarize(n = dplyr::n_distinct(.data$cell), .groups = "drop") %>%
            dplyr::filter(.data$n >= n_min) %>%
            dplyr::arrange(.data$time)
        if (nrow(cells_per_t) > 0) {
            end_time <- cells_per_t$time[1]
            message(sprintf("[plot_lineage_tree] end_time milestone (%s>=%d) -> %s",
                            lin, n_min, end_time))
        } else {
            message(sprintf("[plot_lineage_tree] milestone (%s>=%d) not reached; using max(time)",
                            lin, n_min))
        }
    }
    if (is.null(end_time) || isFALSE(end_time)) end_time <- max(CD$time, na.rm = TRUE)
    CD <- CD[CD$time <= end_time, ]
    if (nrow(CD) == 0) stop("CD is empty after filtering; nothing to plot")

    # Drop "phantom" cells: cells whose only data lies at exactly `end_time`.
    # These are newborn stubs — drawing them produces zero-length segments
    # and tip labels at the cap that don't correspond to displayed branches.
    # (Their parents are still present and will be drawn as the real leaves.)
    cell_min_t <- tapply(CD$time, CD$cell, min)
    phantom <- names(cell_min_t)[cell_min_t >= end_time]
    if (length(phantom) > 0) {
        CD <- CD[!(CD$cell %in% phantom), ]
        if (nrow(CD) == 0) stop("CD is empty after dropping phantom cells")
    }

    # ---- truncate trailing-NA per cell (if requested) ----
    # Each cell's effective end is the last `time` at which `value` is non-NA.
    # CD rows past that point would otherwise render as `na_color` segments
    # extending well beyond the cell's actual data; trim them out so the
    # branch visually ends where measurements end.
    if (truncate_to_last_data) {
        last_data_t <- tapply(CD$time[!is.na(CD$value)],
                              CD$cell[!is.na(CD$value)], max)
        # Cells with no non-NA value at all are not in last_data_t — leave
        # them for the drop_empty_cells pass below, or keep them as fully NA.
        cells_with_data <- names(last_data_t)
        keep <- rep(TRUE, nrow(CD))
        in_data <- CD$cell %in% cells_with_data
        keep[in_data] <- CD$time[in_data] <= last_data_t[CD$cell[in_data]]
        CD <- CD[keep, , drop = FALSE]
        if (nrow(CD) == 0) stop("CD is empty after truncate_to_last_data")
    }

    # ---- prune cells whose value is NA at every timepoint (if requested) ----
    # Iterative because removing a leaf may make its (now-only) sister's parent
    # the new effective leaf — repeated passes settle quickly because each
    # iteration strictly shrinks the unique-cell set.
    if (drop_empty_cells) {
        repeat {
            non_na_per_cell <- tapply(!is.na(CD$value), CD$cell, any)
            empty <- names(non_na_per_cell)[!non_na_per_cell]
            if (length(empty) == 0) break
            CD <- CD[!(CD$cell %in% empty), , drop = FALSE]
            if (nrow(CD) == 0) {
                stop("CD is empty after drop_empty_cells (no cells with data)")
            }
        }
    }

    # ---- optional cross-division smoothing (timepoint mode only) ----
    if (resolution == "timepoint" && smoothing_w > 0) {
        cells <- unique(CD$cell)
        parents_map <- as.list(sapply(cells, GetParent))
        names(parents_map) <- cells
        daughters_map <- split(cells, unlist(parents_map))
        cell_birth <- as.list(tapply(CD$time, CD$cell, min))
        cell_end   <- as.list(tapply(CD$time, CD$cell, max))
        cell_value_map <- new.env(hash = TRUE, parent = emptyenv())
        for (i in seq_len(nrow(CD))) {
            k <- paste(CD$cell[i], CD$time[i], sep = "@")
            cell_value_map[[k]] <- CD$value[i]
        }
        cv_list <- as.list(cell_value_map)

        smoothed <- numeric(nrow(CD))
        for (i in seq_len(nrow(CD))) {
            c_cell <- CD$cell[i]; c_time <- CD$time[i]
            vals <- numeric()
            for (off in seq(-smoothing_w, smoothing_w)) {
                v <- .value_at_time(c_cell, c_time + off, cv_list,
                                    parents_map, daughters_map,
                                    cell_birth, cell_end)
                if (!is.na(v)) vals <- c(vals, v)
            }
            smoothed[i] <- if (length(vals) > 0) mean(vals) else CD$value[i]
        }
        CD$value <- smoothed
    }

    # ---- build the tree layout (tree_y per cell) ----
    layout <- .build_tree_layout(CD, root = root, extra_branch_l = extra_branch_l)
    coords <- layout$coords
    parents <- layout$parents

    # ---- per-resolution segment frame ----
    daughters_of <- function(p) names(parents)[parents == p & !is.na(parents)]
    parent_end_obs <- CD %>%
        dplyr::group_by(.data$cell) %>%
        dplyr::summarize(end_t = max(.data$time, na.rm = TRUE), .groups = "drop")
    pet_vec <- parent_end_obs$end_t; names(pet_vec) <- parent_end_obs$cell

    if (resolution == "cell") {
        agg_fn <- switch(cell_aggregate,
                         mean  = function(x) mean(x, na.rm = TRUE),
                         max   = function(x) max(x, na.rm = TRUE),
                         first = function(x) x[1],
                         last  = function(x) x[length(x)])
        per_cell <- CD %>%
            dplyr::group_by(.data$cell) %>%
            dplyr::summarize(value = agg_fn(.data$value),
                             birth = min(.data$time, na.rm = TRUE),
                             end   = max(.data$time, na.rm = TRUE),
                             .groups = "drop")
        CD_plot <- per_cell %>%
            dplyr::inner_join(coords, by = "cell") %>%
            dplyr::mutate(time_start = .data$birth, time = .data$end)
    } else {
        # timepoint mode: per-timepoint segments via lag(time)
        CD_plot <- CD %>%
            dplyr::inner_join(coords, by = "cell") %>%
            dplyr::group_by(.data$cell) %>%
            dplyr::arrange(.data$time, .by_group = TRUE) %>%
            dplyr::mutate(
                p             = parents[.data$cell],
                default_start = ifelse(is.na(.data$p) | !(.data$p %in% names(pet_vec)),
                                        dplyr::first(.data$time) - 1,
                                        pet_vec[.data$p]),
                time_start    = dplyr::lag(.data$time, default = .data$default_start[1])
            ) %>%
            dplyr::ungroup()
    }

    # ---- division connectors ----
    div_data <- data.frame(parent = character(), time = numeric(),
                           ybottom = numeric(), ytop = numeric())
    for (p in unique(CD$cell)) {
        ds <- daughters_of(p)
        ds <- ds[ds %in% coords$cell]
        if (length(ds) >= 2) {
            p_end <- max(CD$time[CD$cell == p], na.rm = TRUE)
            if (p_end < end_time) {
                y_vals <- coords$tree_y[coords$cell %in% ds]
                div_data <- rbind(div_data,
                                  data.frame(parent = p, time = p_end,
                                             ybottom = min(y_vals),
                                             ytop    = max(y_vals)))
            }
        }
    }

    # ---- tip labels (cells without drawn daughters) ----
    drawn_parents <- div_data$parent
    tip_data <- CD_plot %>%
        dplyr::group_by(.data$cell) %>%
        dplyr::summarize(time = max(.data$time, na.rm = TRUE),
                         tree_y = dplyr::first(.data$tree_y),
                         .groups = "drop") %>%
        dplyr::filter(!(.data$cell %in% drawn_parents))

    # ---- color limits ----
    if (!is.null(percentile_scope)) {
        scope_cells <- grepCells(CD, lineages = percentile_scope, times = "ALL")$cell
        scope_vals <- CD_plot$value[CD_plot$cell %in% scope_cells]
    } else {
        scope_vals <- CD_plot$value
    }
    scope_vals <- scope_vals[!is.na(scope_vals)]

    # Resolve the "data_min" / "data_max" sentinels.
    .resolve_sentinel <- function(x, fallback) {
        if (is.null(x)) return(NULL)
        if (is.character(x) && length(x) == 1) {
            if (x == "data_min") return(if (length(scope_vals) > 0) min(scope_vals) else fallback)
            if (x == "data_max") return(if (length(scope_vals) > 0) max(scope_vals) else fallback)
            stop("value_min/value_max sentinel must be 'data_min' or 'data_max'")
        }
        as.numeric(x)
    }
    vmin <- .resolve_sentinel(value_min, 0)
    vmax <- .resolve_sentinel(value_max, 1)

    if (length(scope_vals) == 0 && is.null(vmin) && is.null(vmax)) {
        clim_lo <- 0; clim_hi <- 1
    } else if (!is.null(vmin) || !is.null(vmax)) {
        # Explicit limits override percentile.
        clim_lo <- if (!is.null(vmin)) vmin else as.numeric(quantile(scope_vals, percentile_bounds[1]))
        clim_hi <- if (!is.null(vmax)) vmax else as.numeric(quantile(scope_vals, percentile_bounds[2]))
        if (!is.null(vmin) && !is.null(vmax) && clim_hi <= clim_lo) {
            clim_hi <- clim_lo + 1
        }
    } else {
        clim_lo <- as.numeric(quantile(scope_vals, percentile_bounds[1]))
        clim_hi <- as.numeric(quantile(scope_vals, percentile_bounds[2]))
        if (!is.null(value_floor))   clim_hi <- max(clim_hi, value_floor)
        if (!is.null(value_ceiling)) clim_hi <- min(clim_hi, value_ceiling)
    }
    if (clim_lo == clim_hi) clim_hi <- clim_lo + .Machine$double.eps

    # ---- build the plot ----
    g <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = div_data,
            ggplot2::aes(x = .data$time, xend = .data$time,
                         y = .data$ybottom, yend = .data$ytop),
            color = "grey50", linewidth = 1.0
        ) +
        ggplot2::geom_segment(
            data = CD_plot,
            ggplot2::aes(x = .data$time_start, xend = .data$time,
                         y = .data$tree_y, yend = .data$tree_y,
                         color = .data$value),
            lineend = "butt", linewidth = branch_width
        ) +
        (if (!is.null(colors) && length(colors) >= 2) {
            # Piecewise stops if `color_values` provided; otherwise even spacing.
            scale_values <- NULL
            if (!is.null(color_values)) {
                if (length(color_values) != length(colors)) {
                    stop("`color_values` must have the same length as `colors`")
                }
                # Clip breakpoints to the active limits, then rescale to [0,1].
                cv <- pmin(pmax(as.numeric(color_values), clim_lo), clim_hi)
                scale_values <- scales::rescale(cv, to = c(0, 1),
                                                from = c(clim_lo, clim_hi))
            }
            ggplot2::scale_color_gradientn(
                colours = colors,
                values = scale_values,
                limits = c(clim_lo, clim_hi),
                oob = scales::squish,
                na.value = na_color,
                name = value_legend
            )
        } else {
            ggplot2::scale_color_gradient(
                low = color_low, high = color_high,
                limits = c(clim_lo, clim_hi),
                oob = scales::squish,
                na.value = na_color,
                name = value_legend
            )
        }) +
        ggplot2::xlim(min(CD_plot$time_start, na.rm = TRUE),
                      max(CD_plot$time, na.rm = TRUE) + 20) +
        ggtree::theme_tree2()

    if (tip_lab) {
        # Auto-scale tip-label size so labels don't overlap. Tips share the
        # vertical panel space evenly (after tree y-layout). geom_text size
        # is in mm; ~70% of mm-per-tip avoids overlap with descender room.
        # Approx 17.78 mm per inch * panel_height_in / n_tips * 0.7.
        if (is.null(tip_size)) {
            n_tips <- nrow(tip_data)
            mm_per_tip <- (plot_height_in * 25.4) / max(n_tips, 1)
            tip_size <- pmin(tip_size_max, pmax(tip_size_min, 0.7 * mm_per_tip))
        }
        g <- g + ggplot2::geom_text(
            data = tip_data,
            ggplot2::aes(x = .data$time + 2, y = .data$tree_y, label = .data$cell),
            hjust = 0, size = tip_size
        )
    }
    if (!is.null(title)) g <- g + ggplot2::ggtitle(title)
    g
}


#' Save a list of plots to a paginated PDF
#'
#' Wraps `patchwork::wrap_plots()` + `patchwork::plot_annotation()` to
#' chunk a list of ggplots into pages of `plots_per_page` and write a
#' single multi-page PDF. Convenient for comparing many lineage trees
#' (one per gene, condition, embryo, etc.) in a single document.
#'
#' @param plots a list of ggplot objects.
#' @param file output PDF path.
#' @param ncol number of columns per page (default 2).
#' @param plots_per_page maximum plots per page (default 4).
#' @param title optional shared title applied to every page via
#'   `patchwork::plot_annotation`.
#' @param width,height PDF page dimensions in inches (default 24 x 16).
#'
#' @return invisibly returns `file` (the output path).
#' @export
paginate_tree_plots <- function(plots, file,
                                ncol = 2, plots_per_page = 4,
                                title = NULL,
                                width = 24, height = 16) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        stop("Package 'patchwork' is required. Install with install.packages('patchwork').")
    }
    if (length(plots) == 0) {
        warning("paginate_tree_plots: empty plot list; nothing written")
        return(invisible(file))
    }

    pages <- split(plots, ceiling(seq_along(plots) / plots_per_page))

    grDevices::pdf(file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (pg in pages) {
        grid <- patchwork::wrap_plots(pg, ncol = ncol)
        if (!is.null(title)) {
            grid <- grid + patchwork::plot_annotation(
                title = title,
                theme = ggplot2::theme(
                    plot.title = ggplot2::element_text(size = 20, face = "bold")))
        }
        print(grid)
    }
    invisible(file)
}

