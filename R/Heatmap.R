#' coerce one input into a NAMED list of embryo tracking data.frames
#' @param x a dataframe or named list of dataframe
#' @param label a dummy variable to identify input from main function for debug
#' @return `x` as a list of dataframe with names added as "embryoX"
#'   where X is the idenx of unnamed list element
.asEmbryoList <- function(x, label) {
  if (is.data.frame(x)) x <- list(x)
  if (!is.list(x) || !all(vapply(x, is.data.frame, logical(1)))) {
    stop("`", label, "` must be a data.frame or a list of data.frames. ",
         "data.frames cannot be stored in atomic vectors, use a list().",
         call. = FALSE)
  }
  if (length(x) == 0) stop("`", label, "` contains no data.frames.", call. = FALSE)
  nms <- names(x)
  if (is.null(nms)) nms <- rep(NA_character_, length(x))
  unnamed <- is.na(nms) | nms == ""
  nms[unnamed] <- paste0("embryo", which(unnamed))
  setNames(x, make.unique(nms))
}

#' parse data0/data1 into a named list of groups (data0 -> "control");
#' each group is a named list of embryo tracking data.frames
#' @return data1 is dataframe/list of dataframe: list(control=list(), experiment = list())
#'   data1 is list-of-list of dataframe: list(control=list(), data1_name1 = list(), data1_name2......)
.parseGroups <- function(data0, data1) {
  groups <- list(control = .asEmbryoList(data0, "data0"))
  if (is.null(data1)) return(groups)
  if (is.data.frame(data1) || # data1 is dataframe or one group only
      (is.list(data1) && all(vapply(data1, is.data.frame, logical(1))))) {
    groups$experiment <- .asEmbryoList(data1, "data1")
    return(groups)
  }
  if (!is.list(data1)) {
    stop("`data1` must be a data.frame, a list of data.frames, or a list of ",
         "lists of data.frames.", call. = FALSE)
  }
  if (length(data1) == 0) stop("`data1` is empty.", call. = FALSE)
  gnms <- names(data1)
  if (is.null(gnms)) gnms <- rep(NA_character_, length(data1))
  unnamed <- is.na(gnms) | gnms == ""
  gnms[unnamed] <- paste0("group", which(unnamed))
  groups <- c(groups,
              setNames(lapply(data1, .asEmbryoList, label = "data1"), gnms))
  setNames(groups, make.unique(names(groups)))
}

# decide which cells touching the first/last timepoint of their data.frame
# are dropped, according to a head/tail argument
# ("drop" | "keep" | character vector of exceptions to keep)
.boundaryDrop <- function(atBoundary, boundaryArg, argName) {
  if (identical(boundaryArg, "drop")) {
    atBoundary
  } else if (identical(boundaryArg, "keep")) {
    atBoundary[] <- FALSE
    atBoundary
  } else if (is.character(boundaryArg)) {
    atBoundary & !(names(atBoundary) %in% boundaryArg)
  } else {
    stop("`", argName, "` must be 'drop', 'keep', or a character vector of ",
         "cell names kept at the boundary.", call. = FALSE)
  }
}

#' summarise ONE embryo data.frame into a named numeric vector with one
#' aggregated attribute value per cell; head/tail-dropped cells carry NA
.summariseEmbryo <- function(df, attribute, agg, aggArgs, tail, head) {
  missingCols <- setdiff(c("cell", "time", attribute), names(df))
  if (length(missingCols) > 0) {
    stop("missing column(s): ", paste(missingCols, collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(df$time)) stop("`time` column must be numeric.", call. = FALSE)
  if (!is.numeric(df[[attribute]])) {
    stop("`", attribute, "` column must be numeric.", call. = FALSE)
  }
  if (nrow(df) == 0) return(setNames(numeric(0), character(0)))
  cells <- as.character(df$cell)
  time  <- as.numeric(df$time)
  vals  <- as.numeric(df[[attribute]])

  bad <- is.na(cells) | is.na(time)          # rows unusable for aggregation
  if (any(bad)) {
    warning("dropped ", sum(bad), " row(s) with NA in `cell` or `time`.",
            call. = FALSE)
    cells <- cells[!bad]; time <- time[!bad]; vals <- vals[!bad]
  }
  if (length(cells) == 0L){
    warning("no useable rows found in an embryo dataframe")
    return(setNames(numeric(0), character(0)))
  }

  timeByCell <- split(time, cells)
  cellMax <- vapply(timeByCell, max, numeric(1))
  cellMin <- vapply(timeByCell, min, numeric(1))

  dropped <- .boundaryDrop(cellMax == max(time), tail, "tail") |
    .boundaryDrop(cellMin == min(time), head, "head")

  valsByCell <- split(vals, cells)
  values <- vapply(valsByCell, function(x) {
    v <- do.call(agg, c(list(x), aggArgs))
    if (length(v) != 1L || !is.numeric(v)) {
      stop("`agg` must return a single numeric value (NA allowed) when given ",
           "the attribute values of one cell.", call. = FALSE)
    }
    as.numeric(v)
  }, numeric(1))

  values[dropped] <- NA_real_
  values
}

#' cellSet / cellOccurrence selection on the raw cells-x-embryos matrix.
#' A cell "occurs" in an embryo when it is observed AND not head/tail-dropped
#' (i.e. has a non-NA value). Returns the selected cell names.
.selectCellSet <- function(rawMat, groupOfEmbryo, cellSet, cellOccurrence) {
  observed <- !is.na(rawMat)
  countInGroup <- function(group) {
    rowSums(observed[, groupOfEmbryo == group, drop = FALSE])
  }
  if (cellSet == "union") {
    rownames(rawMat)[rowSums(observed) >= cellOccurrence]
  } else if (cellSet == "control") {
    rownames(rawMat)[countInGroup("control") >= cellOccurrence]
  } else {# intersection: must qualify in every single group (control included)
    qualifies <- lapply(unique(groupOfEmbryo), function(group) {
      countInGroup(group) >= cellOccurrence
    })
    rownames(rawMat)[Reduce(`&`, qualifies)]
  }
}

# test each data1 group against data0 ("control"), per cell, on the RAW
# per-embryo aggregated values; NA where the test is not possible
.pValues <- function(rawMat, groupOfEmbryo, method) {
  groups <- setdiff(unique(groupOfEmbryo), "control")
  if (length(groups) == 0L) return(NULL)
  testFun <- switch(method,
                    t = function(x, y) stats::t.test(x, y)$p.value,
                    wilcox = function(x, y) stats::wilcox.test(x, y)$p.value)
  pMat <- matrix(NA_real_, nrow = nrow(rawMat), ncol = length(groups),
                 dimnames = list(rownames(rawMat), groups))
  ctrl <- groupOfEmbryo == "control"
  for (g in groups) {
    for (i in seq_len(nrow(rawMat))) {
      x <- rawMat[i, ctrl]; y <- rawMat[i, groupOfEmbryo == g]
      x <- x[!is.na(x)]; y <- y[!is.na(y)]
      if (length(x) < 2L || length(y) < 2L) next   # need n>=2 per side
      p <- suppressWarnings(tryCatch(testFun(x, y), error = function(e) NA_real_))
      pMat[i, g] <- if (is.finite(p)) p else NA_real_  # e.g. all-identical wilcox
    }
  }
  pMat
}

#' plotEmbryoHeatmap()
#' Faceted embryo-by-cell heatmap of any per-cell attribute
#'
#' @param data0 single data.frame or list of data.frames (one per embryo);
#'   the "control" group (per-cell control mean, cellSet reference).
#' @param data1 optional: single data.frame, OR list of data.frames (group
#'   "experiment"), OR a named list of such lists (one group per entry).
#' @param attribute single column name present in every data.frame.
#' @param agg aggregation function (see constraints above the helpers);
#'   default: max-min for attribute == "time", else mean.
#' @param aggArgs list of further arguments for agg.
#' @param tail "drop" (default) | "keep" | character vector of boundary
#'   cells to keep; applies to cells whose max(time) == max(time) of their
#'   data.frame.
#' @param head same, for cells whose min(time) == min(time) of their
#'   data.frame.
#' @param cellSet what cell to display?
#'   "control" (default) only from control embryos
#'   "intersection" cells occuring in all groups
#'   "union" union of cells from all groups, cells not found in control will be normalized to mean of the rest and name will be red
#' @param cellOccurrence single whole number >= 0; 0 == "union" (warning).
#' @param lineages lineage spec(s) for grepCells(); subsets all data.frames
#'   BEFORE cellSet selection (times = "ALL"); no override of cellSet.
#' @param cells exact cell name(s); also passed to grepCells(); always
#'   plotted, even if not found or filtered out by cellSet/cellOccurrence.
#' @param colLimits colour limits on the ratio scale, straddling 1.
#' @param legendBreaks numeric legend breaks.
#' @param legendName legend title (default: "<attribute>\n divided by\n control mean").
#' @param showEmbryoNames show embryo names under the facets (default FALSE,
#'   matching axis.text.x = element_blank()).
#' @param na_col colour for NA cells (default "black").
#' @param return_heatmap whether to return a ComplexHeatmap object that can be re-styled & drawn
#'   or directly draw the plot to grid graphics system (standard ::draw() behavior)
#' @param newpage whether to open a new grid graphics page (default TRUE) or plot on top of a current one,
#'   should use default TRUE. Use FALSE only if you know what your are doing
#' @param test `"wilcox"` or `"t"`.
#'   Perform the indicated test on the averaged plotting value for each cell for each data1 group,
#'   p_value below 0.05 will be indicated with left annotation (brighter for lower p)
#' @return the drawn ComplexHeatmap (re-style/re-draw as needed) or plotted out with nothing useful returned
#' @import ComplexHeatmap
#' @import circlize
plotEmbryoHeatmap <- function(
    data0,
    data1 = NULL,
    attribute,
    agg = NULL,
    aggArgs = list(),
    tail = "drop",
    head = "drop",
    cellSet = c("control", "intersection", "union"),
    cellOccurrence = 3L,
    lineages = NULL,
    cells = NULL,
    showEmbryoNames = FALSE,
    colLimits = NULL,
    legendBreaks = NULL,
    legendName = NULL,
    na_col = "black",
    test = NULL,
    return_heatmap = TRUE, newpage = TRUE
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("This function needs the bioconductor packages ComplexHeatmap")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("This function needs the package circlize")
  }

  ## -- validate arguments -----------------------------------------------
  cellSet <- match.arg(cellSet)
  if (missing(attribute) || !is.character(attribute) || length(attribute) != 1L) {
    stop("`attribute` must be a single column name.", call. = FALSE)
  }
  if (is.null(agg)) {
    agg <- if (identical(attribute, "time")) function(x) max(x) - min(x) else mean
  }
  if (!is.function(agg)) stop("`agg` must be a function or NULL.", call. = FALSE)
  if (is.null(aggArgs)) aggArgs <- list()
  else if (!is.list(aggArgs)) stop("`aggArgs` must be a list.", call. = FALSE)
  if (!identical(tail, "drop") && !identical(tail, "keep") && !is.character(tail)) {
    stop("`tail` must be 'drop', 'keep', or a character vector of cell names ",
         "kept at the last timepoint.", call. = FALSE)
  }
  if (!identical(head, "drop") && !identical(head, "keep") && !is.character(head)) {
    stop("`head` must be 'drop', 'keep', or a character vector of cell names ",
         "kept at the first timepoint.", call. = FALSE)
  }
  if (length(cellOccurrence) != 1L || is.na(cellOccurrence) ||
      cellOccurrence %% 1 != 0 || cellOccurrence < 0) {
    stop("`cellOccurrence` must be a single whole number >= 0.", call. = FALSE)
  }
  cellOccurrence <- as.integer(cellOccurrence)
  if (is.null(data1) && cellSet != "control") {
    message("`data1` is NULL: enforcing cellSet = 'control'.")
    cellSet <- "control"
  } else if (cellOccurrence == 0L && cellSet != "union") {
    # with a threshold of 0, "occurs in >= 0 data frames" is trivially true
    # for every observed cell, so all three methods collapse to the union
    warning("cellOccurrence = 0 makes cellSet = '", cellSet,
            "' equivalent to cellSet = 'union'; using union cells.", call. = FALSE)
    cellSet <- "union"
  }
  if (!is.null(legendBreaks)) {
    if (!is.numeric(legendBreaks) || anyNA(legendBreaks) || any(legendBreaks <= 0)) {
      stop("`legendBreaks` must be positive numbers (ratio scale).", call. = FALSE)
    }
    legendBreaks <- sort(unique(legendBreaks))
  }
  if (!is.null(colLimits)) {
    colLimits <- sort(colLimits)
    if (length(colLimits) != 2L || !is.numeric(colLimits) || anyNA(colLimits) ||
        colLimits[1] >= 1 || colLimits[2] <= 1) {
      stop("`colLimits` must be two numbers straddling 1.", call. = FALSE)
    }
  }

  ### -- `test` ----
  if (!is.null(test)) {
    if (!is.list(test)) stop("`test` must be NULL or a list.", call. = FALSE)
    defaults <- list(method = "wilcox", cutoff = 0.05,
                     palette = circlize::colorRamp2(
                       seq(0, 1, length.out = 5),
                       c("#440153", "#3e4989", "#21918c", "#5ec962", "#fde725")
                     ),
                     na.color = "gray44")
    unknown <- setdiff(names(test), names(defaults))
    if (length(unknown) > 0L) {
      stop("unknown field(s) in `test`: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    test <- utils::modifyList(defaults, test)   # partial lists keep defaults
    test$method <- match.arg(test$method, c("t", "wilcox"))
    if (!is.numeric(test$cutoff) || length(test$cutoff) != 1L ||
        is.na(test$cutoff) || test$cutoff <= 0 || test$cutoff >= 1) {
      stop("`test$cutoff` must be a single number strictly between 0 and 1.",
           call. = FALSE)
    }
    if (!is.function(test$palette)) {
      stop("`test$palette` must be a colorRamp2() function on the ",
           "significance scale s = 1 - p/cutoff (1 = most significant).",
           call. = FALSE)
    }
  }

  ## -- parse data0/data1 into groups of embryo data.frames ---------------
  groups <- .parseGroups(data0, data1)
  groupNames <- make.unique(names(groups))
  names(groups) <- groupNames

  ## -- optional lineages/cells subsetting BEFORE cellSet selection -------
  if (!is.null(lineages) || !is.null(cells)) {
    groups <- groups |> lapply(function(group) {
      group |> lapply(function(embdf) {
        grepCells(CDData = embdf, cells = cells, lineages = lineages, times = "ALL")
      })
    })
  }

  ## -- one named vector of per-cell values per embryo ---------------------
  embryoValues <- seq_along(groups) |> lapply(
    function(groupId) {
      group <- groups[[groupId]]
      seq_along(group) |> lapply(function(embID) {
        embdf <- group[[embID]]
        tryCatch(
          .summariseEmbryo(embdf, attribute = attribute, agg = agg,
                           aggArgs = aggArgs, tail = tail, head = head),
          error = function(e) stop(
            "group '", groupNames[gi], "', embryo '", names(embryos)[ei], "': ",
            conditionMessage(e), call. = FALSE
          )
        )
    })
  }) |>
    unlist(recursive = FALSE)
  embryoNames <- make.unique(unlist(lapply(groups, names)))
  groupOfEmbryo <- rep(groupNames, times = vapply(groups, length, integer(1)))
  names(embryoValues) <- embryoNames

  ## -- cells x embryos matrix unfiltered ------------------------------------------
  allCells <- sort(unique(unlist(lapply(embryoValues, names))))
  rawMat <- matrix(NA_real_, nrow = length(allCells),
                   ncol = length(embryoValues), dimnames = list(allCells, embryoNames))
  for (i in seq_along(embryoValues)) {
    v <- embryoValues[[i]]
    rawMat[names(v), i] <- v
  }

  ## -- cell set selection (shared y axis, enforced on ALL groups) ----------
  selectedCells <- .selectCellSet(rawMat, groupOfEmbryo, cellSet, cellOccurrence)
  forcedCells <- if (is.null(cells)) character(0) else as.character(cells)
  plotCells <- sort(unique(c(selectedCells, forcedCells)))
  if (length(plotCells) == 0L) {
    warning("no cells selected (cellSet = '", cellSet, "', cellOccurrence = ",
            cellOccurrence, "); nothing to plot.", call. = FALSE)
    return(invisible(NULL))
  }

  rawMat <- rawMat[plotCells, , drop = FALSE]

  ## -- normalise against the per-cell data0 (control) mean -----------------
  controlMean <- rowMeans(rawMat[, groupOfEmbryo == "control", drop = FALSE],
                          na.rm = TRUE)
  controlMean[is.nan(controlMean)] <- NA_real_
  normMat <- if (nrow(rawMat) == 0L) rawMat else rawMat / controlMean
  missing <- is.na(controlMean)
  if (any(missing)) {
    fallback <- rowMeans(rawMat[missing, , drop = FALSE], na.rm = TRUE)
    warning("no data0 embryo has cell(s): ", paste(names(fallback), collapse = ", "),
            "; normalised by the mean of the ", colSums(!is.na(rawMat[missing, , drop = FALSE])),
            " other embryo(s) that do.", call. = FALSE)
    controlMean[missing] <- fallback
  }

  ## -- `cells` argument mention unobserved cells: plot as NA ----------
  extraCells <- setdiff(plotCells, rownames(normMat))
  if (length(extraCells) > 0L) {
    normMat <- rbind(
      normMat,
      matrix(NA_real_, nrow = length(extraCells), ncol = ncol(normMat),
             dimnames = list(extraCells, colnames(normMat)))
    )
  }
  normMat <- normMat[plotCells, , drop = FALSE]
  unnormalisable <- plotCells[is.na(controlMean[plotCells])]
  if (length(unnormalisable) > 0L) {
    warning("no non-dropped values in any data0 embryo for cell(s): ",
            paste(unnormalisable, collapse = ", "),
            "; they cannot be normalised and are shown as NA.", call. = FALSE)
  }

  ## -- "Overall" row: per-embryo mean of the normalised values -------------
  overall <- suppressWarnings(colMeans(normMat, na.rm = TRUE))
  overall[is.nan(overall)] <- NA_real_
  plotMat <- rbind(Overall = overall, normMat)

  ### auto colLimits / legendBreaks: log2-symmetric, always centred on 1 ----
  colLimitsAuto <- is.null(colLimits)
  vals <- plotMat[is.finite(plotMat) & plotMat > 0]   # ratio scale assumed
  if (colLimitsAuto) {
    if (length(vals) == 0L) {
      warning("no finite positive values; falling back to colLimits = c(0.5, 4).",
              call. = FALSE)
      colLimits <- c(0.5, 4)
    } else {
      r <- log2(vals)
      k <- stats::quantile(abs(r), 0.95, names = FALSE)  # outlier-robust; use max(abs(r)) for literal range
      k <- max(0.5, ceiling(k * 4) / 4)                  # floor at 0.5-2, snap up to 0.25 steps
      colLimits <- c(2^-k, 2^k)                          # e.g. k=1.5 -> 0.354-2.83
    }
  }
  if (is.null(legendBreaks)) {
    brks <- 2^pretty(c(-log2(colLimits[2]), log2(colLimits[2])), n = 4)
    brks <- unique(signif(brks[brks >= colLimits[1] & brks <= colLimits[2]], 3))
    if (!any(abs(brks - 1) < 1e-12)) brks <- sort(c(brks, 1))
    legendBreaks <- brks
  } else if (colLimitsAuto &&
             any(legendBreaks < colLimits[1] | legendBreaks > colLimits[2])) {
    ## user breaks stick out of auto limits: widen symmetrically so nothing clips
    k <- max(abs(log2(c(colLimits, legendBreaks))))
    colLimits <- c(2^-k, 2^k)
    warning("`legendBreaks` outside the data-driven colLimits; widened to ",
            paste(signif(colLimits, 3), collapse = "-"), ".", call. = FALSE)
  }

  ## -- row label colour: red for cells absent from every data0 embryo -----
  controlMissing <- rowSums(!is.na(rawMat[, groupOfEmbryo == "control", drop = FALSE])) == 0L
  miss <- controlMissing[plotCells]
  miss[is.na(miss)] <- TRUE               # never-observed cells = missing too
  labelCol <- c("black", ifelse(miss, "red", "black"))  # 1st entry = "Overall"

  ## -- colour scale: diverging scale, linear in log2 space, centred on 1 ---
  # values outside colLimits are squished onto the extreme colours)
  rampBreaks <- log2(c(colLimits[1], 1, colLimits[2]))
  colRamp <- circlize::colorRamp2(rampBreaks,
                                  c(rgb(0, 0, 1), rgb(1, 1, 1), rgb(1, 0, 0)))
  colFun <- function(x) {
    out <- rep(na_col, length(x))
    keep <- !is.na(x)
    if (any(keep)) {
      lg <- log2(x[keep]) |>
        pmax(rampBreaks[1]) |>
        pmin(rampBreaks[3])
      out[keep] <- colRamp(lg)
    }
    out
  }
  attr(colFun, "breaks") <- legendBreaks

  ### -- significance annotation, one layer per data1 group -----
  rightAnn <- NULL
  if (!is.null(test)) {
    if (is.null(data1)) {
      warning("`test` ignored: no data1 group to compare against data0.", call. = FALSE)
    } else {
      rawOverall <- colMeans(rawMat, na.rm = TRUE) # test the mean of each embryo as well
      rawOverall[is.nan(rawOverall)] <- NA_real_   # embryos with no values at all
      testMat <- rbind(Overall = rawOverall, rawMat)

      pMat <- .pValues(testMat, groupOfEmbryo, test$method)
      # significance score s in (0,1]: s -> 1 as p -> 0 (high palette end,
      # i.e. yellow for the default); p > cutoff or no test -> NA -> na.color
      sMat <- pMat * NA_real_
      below <- !is.na(pMat) & pMat <= test$cutoff
      sMat[below] <- 1 - pMat[below] / test$cutoff
      palBrks <- if (is.null(attr(test$palette, "breaks"))) {
        seq(0, 1, length.out = 5)
      } else attr(test$palette, "breaks")
      # legend ticks shown as p-values, reversed from the s scale
      annGroups <- colnames(pMat)
      legParams <- setNames(lapply(annGroups, function(g) {
        list(title = paste0(g, " ", test$method, " p"),
             at = palBrks,
             labels = as.character(round(pmax(0, 1 - palBrks) * test$cutoff, 3)),
             show_legend = g == annGroups[1])           # one shared legend
      }), annGroups)

      rightAnn <- do.call(ComplexHeatmap::rowAnnotation, c(
        setNames(lapply(annGroups, function(g) sMat[rownames(plotMat), g]), annGroups),
        list(
          col = setNames(rep(list(test$palette), length(annGroups)), annGroups),
          na_col = test$na.color,
          border = TRUE,
          simple_anno_size = grid::unit(4, "mm"),
          annotation_legend_param = legParams
        )
      ))

    }
  }

  ## -- heatmap ----
  colSplit <- factor(groupOfEmbryo, levels = groupNames)
  rowSplit <- factor(c("Overall", rep("Cells", nrow(normMat))),
                     levels = c("Overall", "Cells"))

  ht <- ComplexHeatmap::Heatmap(
    plotMat,
    name = "value",
    col = colFun,
    na_col = na_col,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, row_names_side = "left", row_split = rowSplit,
    row_title = NULL, row_gap = grid::unit(3, "mm"),
    row_names_gp = grid::gpar(col = labelCol),
    show_column_names = showEmbryoNames, column_split = colSplit,
    column_title = groupNames, column_gap = grid::unit(3, "mm"),
    border = TRUE,
    right_annotation = rightAnn,
    heatmap_legend_param = list(
      title = legendName,
      at = legendBreaks,
      labels = as.character(legendBreaks))
  )

  if (return_heatmap) return(invisible(ht))
  ComplexHeatmap::draw(ht, newpage = newpage)
  invisible(NULL)
}
