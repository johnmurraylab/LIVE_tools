#' drawEmb
#' This file contains several functions relevant to drawing a 3D (interactive) embryo,
#' able to represent expression levels and highlight lineages
#' if deliver multiple "lineages" must also deliver a list "lineageColors" with elements corresponding to "lineages"


#' drawEmbVal draw the nucleus positions of the embryo in 3D scatter plots, colored by expression
#'
#' @param embryoCD the embryo dataframe with cell, time, and blot(expression) column, mandatory
#' @param time which time to plot the embryo, mandatory
#' @param valCol column that will color the nucleus by
#' @param lineages a list of lineages to highlight (plot with the given symbol in "shapes" parameter)
#' @param shapes a list of shapes to plot each lineage given in "lineages"
#' @param ReporterForAll color the nucleus points by expression or not (default TRUE)
#' @param colorScheme a list of at least 2 vectors made of a number between 0 and 1 to specify expression fraction, and a color value to specify color fro that fraction
#' @param maxBlot upper limit for expression coloring
#' @param minBlot lower limit for expresion coloring
#' @param xSize size of each x unit in embryoCD dataframe
#' @param ySize size of each y unit in embryoCD dataframe
#' @param zSize size of each z unit in embryoCD dataframe
#' @param center which point in the 3d space does the camera aim at
#' @param viewPoint where to put the imaginative camera initially
#' @param showAxis whether to display the x,y,z axis or not
#'
#' @return a list made of: the plotly figure, a dataframe of the highlighted cells, and a dataframe of cells not highlighted
#' @export
#' @import plotly
#' @import viridis
#'
#' @examples drawEmbVal(Embryo, lineages = c("MS", "E", "C"), shapes = c("circle","x","square"), time = 139)
drawEmbVal<-function(embryoCD, time, valCol = "blot",
                     lineages=c("P0"),
                     shapes = NULL,
                     ReporterForAll=F, colorScheme = NULL,
                     maxBlot = NULL, minBlot = NULL,
                     xSize=1, ySize=1, zSize=1, #whiteSpace = c(0.05,2),
                     center = list(x=0,y=0,z=0), viewPoint = list(x=0,y=0,z=1.8),
                     showAxis = F){
  embDat <- embryoCD|>grepCells(lineages = "ALL", times = time)
  row.names(embDat)<- NULL
  embDat$x <- embDat$x*xSize
  embDat$y <- embDat$y*ySize
  embDat$z <- embDat$z*zSize
  cell_value <- aggregate(embryoCD[,valCol], by = list(embryoCD[,"cell"]), FUN = "mean")
  names(cell_value)<-c("cell", "value")
  embDat[,valCol]<-NULL
  embDat<- left_join(embDat, cell_value, by = "cell")
  if(is.null(maxBlot)){maxBlot <- max(embDat$value)}
  if(is.null(minBlot)){minBlot <- min(embDat$value)}
  if(identical(lineages, NULL)){lineages <- list("ABa", "ABp","MS", "E", "C", "D", "P4")} #default lineages
  if(length(shapes) != length(lineages)){
    print("\'shapes\' parameters not properly specified, using circles for all")
    # shapesSet <- schema(F)$traces$scatter$attributes$marker$symbol$values
    # shapesSet <- grep("[a-z]", shapesSet, value = TRUE)
    # shapesSet <- shapesSet[-grep("dot", shapesSet)]
    # shapesSet <- shapesSet[-grep("open", shapesSet)]
    # shapes <- seq(from = 1, by = 1, length.out = length(lineages))
    # shapes = shapesSet[shapes]
    shapes <- rep("circle", times = length(lineages))
  }
  if(is.null(colorScheme)){colorScheme =  "Viridis"}
  else if(colorScheme=="blackBlueOrange"){
    colorScheme<-list(c(0,rgb(0,0,0)), c(0.3,rgb(0,0.2,1)), c(0.7,rgb(1,0,0.2)), c(1,rgb(1,0.7,0.3)))}
  else if(colorScheme=="blackOrange"){colorScheme<-list(c(0,rgb(0,0,0)), c(0.6,rgb(1,0.1,0.1)), c(1,rgb(1,0.7,0.3)))}

  # rangeDim <- function(dimVals){
  #   span <- max(dimVals)-min(dimVals)
  #   return(span*c(-whiteSpace[1],whiteSpace[1]) + c(-whiteSpace[2],whiteSpace[2]) + range(dimVals))
  # }
  # xrange <- rangeDim(embryoCD[,"x"])
  # yrange <- rangeDim(embryoCD[,"y"])
  # zrange <- rangeDim(embryoCD[,"z"])
  if(showAxis){
    xtitle <- "x"
    ytitle <- "y"
    ztitle <- "z"
  }
  else{
    xtitle <- ""
    ytitle <- ""
    ztitle <- ""
  }
  fig<- plotly::plot_ly() |> #initiate a plotly figure object
    plotly::layout(
      scene = list(
        aspectmode = "data",
        xaxis = list(title = xtitle, showgrid = showAxis, showticklabels=showAxis),
        yaxis = list(title = ytitle, showgrid = showAxis, showticklabels=showAxis),
        zaxis = list(title = ztitle, showgrid = showAxis, showticklabels=showAxis),
        camera = list(eye = viewPoint, center=center, up = list(x = 0, y = 1, z = 0)
                      #,projection=list(type="orthographic")
        )),
      paper_bgcolor = rgb(1,1,1)
    )

  selectCells <- NULL
  for (i in seq_along(lineages)) { #plot highlighted cells
    lineage <- lineages[[i]]
    thisSymbol <- shapes[[i]]
    thisCells <- grepCells(CDData = embDat, lineages = lineage, dataReturn = F)
    fig <- fig%>%AddGroupExp(groupName = lineage, data = embDat, selectedCells = thisCells, symbol = thisSymbol,
                          colorScheme = colorScheme, colorMin = minBlot, colorMax = maxBlot)
    selectCells <- union(selectCells, thisCells)
  }

  otherCells <- embDat|>rownames()|>as.integer()
  otherCells <- otherCells[!otherCells%in%selectCells]
  if(is.null(lineages)){alpha <- 1}
  else{alpha <- 0.4}
  fig<-fig%>%plotly::add_trace(name="other cells",
                       x=embDat[otherCells,"x"],
                       y=embDat[otherCells,"y"],
                       z=embDat[otherCells,"z"], #cells that are not selected will be plotted transparent
                       type = "scatter3d", mode="markers",
                       marker = list(color=embDat[otherCells,"value"], size = 6, opacity = alpha,
                                     symbol = "circle-open",
                                     #line=list(color = "grey", width=3),
                                     colorscale=colorScheme, cmin=minBlot, cmax=maxBlot)
  )

  return(list(fig, embDat[selectCells,], embDat[otherCells,]))
}

#' AddGroupExp
#' @description
#' auxiliary function for `drawEmbVal`
#'
#' @import plotly
#' @return modified plotly figure
AddGroupExp <-function(plotlyFig, groupName, data, selectedCells, symbol,
                   colorScheme, colorMin, colorMax, opacity=1){
  plotlyFig<-plotlyFig|>plotly::add_trace(name=groupName,
             x=data[selectedCells,"x"],
             y=data[selectedCells,"y"],
             z=data[selectedCells,"z"], #cells that are not selected will be plotted transparent
             type = "scatter3d", mode="markers",
             marker = list(color=data[selectedCells,"value"], size = 7.5, opacity = 1,
                           symbol = symbol,
                           colorscale=colorScheme, cmin=colorMin, cmax=colorMax,
                           colorbar = list(title="expression", x = 0)
             )
  )
  return(plotlyFig)
}

#' drawEmbLine draw the nucleus positions of the embryo in 3D scatter plots, colored by lineage
#'
#' @param embryoCD the embryo dataframe with cell, time column, mandatory
#' @param time which time to plot the embryo, mandatory
#' @param lineages a list of lineages to highlight (plot with the given color and opacity in "colors" and "opacity_s" parameter)
#' @param colors a list of colors to plot each lineage given in "lineages"
#' @param opacity_s a list of opacity for each lineage given in "lineages"
#' @param xSize size of each x unit in embryoCD dataframe
#' @param ySize size of each y unit in embryoCD dataframe
#' @param zSize size of each z unit in embryoCD dataframe
#' @param otherOpacity opacity for cells not specified
#' @param cellSize single numeric, size of the data points
#' @param center which point in the 3d space does the camera aim at
#' @param viewPoint where to put the imaginative camera initially
#' @param showAxis wherther to diaply axis and axis labels
#' @param time which time to plot the embryo, mandatory
#'
#' @return a list made of: the plotly figure, a dataframe of the highlighted cells, and a dataframe of cells not highlighted
#' @export
#' @import plotly
#' @import viridis
#'
#' @examples drawEmbVal(Embryo, lineages = c("MS", "E", "C"), colors = c(rgb(1,0,0), rgb(1,0,1), rgb(0,0,1)), time = 139)
drawEmbLine <- function(embryoCD, time, lineages=NULL, colors = NULL, opacity_s = NULL,
                        xSize=1, ySize=1, zSize=1,
                        otherOpacity = 0.2, cellSize = 7.5,
                        center = list(x=0,y=0,z=0), viewPoint = list(x=0,y=0,z=1.8),
                        showAxis = F){
  embDat <- embryoCD|>grepCells(lineages = "ALL", times = time)
  row.names(embDat)<- NULL
  embDat$x <- embDat$x*xSize
  embDat$y <- embDat$y*ySize
  embDat$z <- embDat$z*zSize
  if(identical(lineages, NULL)){lineages <- list("ABa", "ABp","MS", "E", "C", "D", "P4")} #default lineages
  traceCount <- length(lineages)
  if(length(colors) != traceCount){
    print("\'colors\' argument not properly specified")
    colors <- viridis::viridis_pal(option = "H")(traceCount)
  }
  if(length(opacity_s) != traceCount){
    print("\'opacity_s\' argument not properly specified")
    opacity_s <- rep(1, traceCount)
  }
  if(showAxis){
    xtitle <- "x"
    ytitle <- "y"
    ztitle <- "z"
  }
  else{
    xtitle <- ""
    ytitle <- ""
    ztitle <- ""
  }
  fig<-plotly::plot_ly() |> #initiate a plotly figure object
    plotly::layout(
      scene = list(
        aspectmode = "data",
        xaxis = list(title = xtitle, showgrid = showAxis, showticklabels=showAxis),
        yaxis = list(title = ytitle, showgrid = showAxis, showticklabels=showAxis),
        zaxis = list(title = ztitle, showgrid = showAxis, showticklabels=showAxis),
        camera = list(eye = viewPoint, center=center, up = list(x = 0, y = 1, z = 0)
                      #,projection=list(type="orthographic")
                      )),
      paper_bgcolor = rgb(1,1,1)
    )
  selectCells <- NULL
  for (i in seq_along(lineages)) { #add each lineage as a trace
    lineage <- lineages[[i]]
    color <- colors[[i]]
    opacity <- opacity_s[[i]]
    thisCells <- grepCells(CDData = embDat, lineages = lineage, dataReturn = F)
    fig <- fig|>AddGroupLine(groupName = lineage,
                             data = embDat, selectedCells = thisCells,
                             color=color, opacity = opacity, cellSize=cellSize)
    selectCells <- union(selectCells, thisCells)
  }
  otherCells <- embDat|>rownames()|>as.integer()
  otherCells <- otherCells[!otherCells%in%selectCells]
  fig<-fig%>%add_trace(
    name="other cells",
    x=embDat[otherCells,"x"],
    y=embDat[otherCells,"y"],
    z=embDat[otherCells,"z"], #cells that are not selected will be plotted transparent
    type = "scatter3d", mode="markers",
    marker = list(color="black", size = cellSize, opacity = otherOpacity, symbol = "circle")
  )
  return(list(fig, embDat[selectCells,], embDat[-selectCells,]))
}

#' AddGroupLine
#' @description
#' auxiliary function for `drawEmbLine`
#'
#' @import plotly
#' @return modified plotly figure
AddGroupLine <- function(plotlyFig, groupName, data, selectedCells, color, opacity, cellSize){
  out <- plotlyFig|>plotly::add_trace(
    name = groupName,
    x=data[selectedCells,"x"],
    y=data[selectedCells,"y"],
    z=data[selectedCells,"z"],
    type = "scatter3d", mode="markers",
    marker = list(color=color, size = cellSize, opacity = opacity, symbol = 'circle')
  )
 out
}

#' printEmbImg Generates a .png image from drawEmbExp or drawEmbLine output figure
#' Require Reticulate that connects to a python distribution with "kaleido" and "plotly" installed
#'  recommend setting up a python venv for this and specify the python kernel by:
#'  reticulate::use_virtualenv("the/directory/holding/python/executable")
#'
#' @param fig plotly figure
#' @param output_file output file name
#' @param center coordinate the camera face
#' @param viewPoint from where does the camera look at the embryo
#' @param up which direction (as unit vector) the upper side of the camer points to
#' @param width in pixels, width of output image
#' @param height in pixels, height of output image
#' @param keepLabels Boolean, to keep the color bar and legends or not
#' @export
#' @import plotly
saveEmbImg <- function(
    fig, output_file = "embryo_DV_view.png",
    center = list(0,0,0), viewPoint = list(x=0,y=0,z=1.8), up = list(x = 0, y = 1, z = 0),
    width = 600, height = 450, keepLabels = F
    ){
  fig <- fig|>
    plotly::layout(
      margin = list(l = 0, r = 0, t = 0, b = 0),
      scene = list(
        domain = list(x = c(0, 1), y = c(0, 1)),
        camera = list(eye = viewPoint, center=center, up = up)
      )
    )
  if(!keepLabels){fig<-fig|> plotly::hide_colorbar() |> plotly::hide_legend()}
  plotly::save_image(fig, file = output_file, width = width, height = height, scale = 2)
}
