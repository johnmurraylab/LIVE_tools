#' drawEmb
#' This file contains several functions relevant to drawing a 3D (interactive) embryo, 
#' able to represent expression levels and highlight lineages
#' if deliver multiple "lineages" must also deliver a list "lineageColors" with elements corresponding to "lineages"
#' 
#' @import plotly
source('LineageProcessing.R')
if(!require('plotly')){stop("Package 'plotly' could not be loaded")}
if (!requireNamespace("webshot", quietly = TRUE)) {
  stop("Package 'webshot' is required. Install with install.packages('webshot')")
}
if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required. Install with install.packages('magick')")
}
if (!requireNamespace("animation", quietly = TRUE)) {
  stop("Package 'animation' is required. Install with install.packages('animation')")
}

drawEmbExp<-function(embryoCD, time, lineages=NULL, symbols = NULL, 
                     ReporterForAll=F, colorScheme = "blackBlueOrange", 
                     maxBlot = NULL, minBlot = NULL, 
                     xSize=0.08651, ySize=0.08651, zSize=0.5, 
                     center = list(x=0,y=0,z=0), viewPoint = list(x=0,y=0,z=2.5)
                     ){
  embDat <- embryoCD|>grepCells(lineages = "ALL", times = time)
  embDat$x <- embDat$x*xSize
  embDat$y <- embDat$y*ySize
  embDat$z <- embDat$z*zSize
  if(is.null(maxBlot)){maxBlot <- max(embDat$blot)}
  if(is.null(minBlot)){minBlot <- min(embDat$blot)}
  if(length(symbols) != length(lineages)){stop("symbols argument not properly specified")}
  if(colorScheme=="blackBlueOrange"){#default color scheme for reporter intensity (blue red white)
    colorScheme<-list(c(0,rgb(0,0,0)), c(0.3,rgb(0,0.2,1)), c(0.7,rgb(1,0,0.2)), c(1,rgb(1,0.7,0.3)))}
  else if(colorScheme=="blackOrange"){colorScheme<-list(c(0,rgb(0,0,0)), c(0.6,rgb(1,0.1,0.1)), c(1,rgb(1,0.7,0.3)))}
  
  fig<-plot_ly() |> #initiate a plotly figure object
    layout(
      scene = list(aspectmode = "data", 
        xaxis = list(title = '', showgrid = FALSE, showticklabels=F),
        yaxis = list(title = '', showgrid = FALSE, showticklabels=F),
        zaxis = list(title = '', showgrid = FALSE, showticklabels=F),
        camera = list(eye = viewPoint, center=center, up = list(x = 0, y = 1, z = 0),
                      projection=list(type="orthographic"))), 
      paper_bgcolor=rgb(1,1,1))
  
  selectCells <- NULL
  for (i in seq_along(lineages)) { #plot highlighted cells
    lineage <- lineages[[i]]
    thisSymbol <- symbols[[i]]
    thisCellRE <- lineageRE(lineage)
    thisCells <- grep(thisCellRE, embDat$cell)
    fig <- fig%>%AddGroupExp(groupName = lineage, data = embDat, selectedCells = thisCells, symbol = thisSymbol, 
                          colorScheme = colorScheme, colorMin = minBlot, colorMax = maxBlot)
    selectCells <- union(selectCells, thisCells)
  }
  
  if(is.null(lineages)){remainCells<-embDat}
  else{remainCells <- embDat[-selectCells,]}
  fig<-fig%>%add_trace(name="other cells",
                       x=remainCells[,"x"], 
                       y=remainCells[,"y"], 
                       z=remainCells[,"z"], #cells that are not selected will be plotted transparent
                       type = "scatter3d", mode="markers",
                       marker = list(color=remainCells[,"blot"], size = 10, opacity = 1, 
                                     symbol = "circle-open",
                                     #line=list(color = "grey", width=3), 
                                     colorscale=colorScheme, cmin=minBlot, cmax=maxBlot)
  )
  
  return(list(fig, embDat[selectCells,], remainCells[,]))
}

AddGroupExp <-function(plotlyFig, groupName, data, selectedCells, symbol, 
                   colorScheme, colorMin, colorMax, opacity=1){
  plotlyFig<-plotlyFig|>add_trace(name=groupName,
             x=data[selectedCells,"x"], 
             y=data[selectedCells,"y"], 
             z=data[selectedCells,"z"], #cells that are not selected will be plotted transparent
             type = "scatter3d", mode="markers",
             marker = list(color=data[selectedCells,"blot"], size = 10, opacity = 1, 
                           symbol = symbol, 
                           colorscale=colorScheme, cmin=colorMin, cmax=colorMax, 
                           colorbar = list(title="expression", x = 0)
             )
  )
  return(plotlyFig)
}

drawEmbLine <- function(embryoCD, time, lineages, colors, opacitys,
                        xSize=0.08651, ySize=0.08651, zSize=0.5, 
                        otherOpacity = 0.2, cellSize = 10,
                        center = list(x=0,y=0,z=0), viewPoint = list(x=0,y=0,z=2.5), 
                        showAxis = F){
  embDat <- embryoCD|>grepCells(lineages = "ALL", times = time)
  embDat$x <- embDat$x*xSize
  embDat$y <- embDat$y*ySize
  embDat$z <- embDat$z*zSize
  if(length(colors) != length(lineages)){stop("\'colors\' argument not properly specified")}
  if(length(opacitys) != length(lineages)){stop("\'opacitys\' argument not properly specified")}
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
  fig<-plot_ly() |> #initiate a plotly figure object
    layout(
      scene = list(
        aspectmode = "data", 
        xaxis = list(title = xtitle, showgrid = showAxis, showticklabels=showAxis),
        yaxis = list(title = ytitle, showgrid = showAxis, showticklabels=showAxis),
        zaxis = list(title = ztitle, showgrid = showAxis, showticklabels=showAxis), 
        camera = list(eye = viewPoint, center=center, up = list(x = 0, y = 1, z = 0),
                      projection=list(type="orthographic"))), 
      paper_bgcolor = rgb(1,1,1)
    )
  selectCells <- NULL
  for (i in seq_along(lineages)) { #add each lineage as a trace
    lineage <- lineages[[i]]
    color <- colors[[i]]
    opacity <- opacitys[[i]]
    thisCellRE <- lineageRE(lineage)
    thisCells <- grep(thisCellRE, embDat$cell)
    fig <- fig|>AddGroupLine(groupName = lineage, 
                             data = embDat, selectedCells = thisCells,
                             color=color, opacity = opacity, cellSize=cellSize)
    selectCells <- union(selectCells, thisCells)
  }
  fig<-fig%>%add_trace(
    name="other cells",
    x=embDat[-selectCells,"x"], 
    y=embDat[-selectCells,"y"], 
    z=embDat[-selectCells,"z"], #cells that are not selected will be plotted transparent
    type = "scatter3d", mode="markers",
    marker = list(color="black", size = cellSize, opacity = otherOpacity, symbol = "circle")
  )
  return(list(fig, embDat[selectCells,], embDat[-selectCells,]))
}

AddGroupLine <- function(plotlyFig, groupName, data, selectedCells, color, opacity, cellSize){
  out <- plotlyFig|>add_trace(
    name = groupName, 
    x=data[selectedCells,"x"], 
    y=data[selectedCells,"y"], 
    z=data[selectedCells,"z"], 
    type = "scatter3d", mode="markers",
    marker = list(color=color, size = cellSize, opacity = opacity, symbol = 'circle')
  )
 out
}


#' createRotationGIF
#' NOT IMPLEMENTED!!
#' Generates a GIF showing rotation around the x-axis.
#' need Reticulate that connects to a python distribution with "kaleido" and "plotly" installed
#' @param fig The plotly figure object generated by drawEmbExp or drawEmbLine.
#' @param output_file Path to save the output GIF. Default: "animation.gif".
#' @param radius Camera distance from the center. Default: 2.5.
#' @param num_frames Number of frames in the animation. Default: 36.
#' @param duration Total duration of the GIF in seconds. Default: 5.
#' @param center the coordinate that camera rotates around, default (0,0,0)
#' @param width Width of the output GIF in pixels. Default: 800.
#' @param height Height of the output GIF in pixels. Default: 800.
#' @param dpi DPI for the output images. Default: 150.
#' 
#' @note Requires 'webshot', 'magick', and 'animation' packages. Ensure ImageMagick and PhantomJS are installed and accessible in the system PATH.
#'
# RotatingGIF <- function(fig, output_file = "animation.gif", 
#                               radius = 2.5, num_frames = 36, 
#                               duration = 5, center = list(0,0,0),
#                               width = 800, height = 800, dpi = 150) {
#   # Create temporary directory for frames
#   temp_dir <- tempdir()
#   dir.create(temp_dir, showWarnings = FALSE)
#   # Generate angles for rotation
#   angles <- seq(0, 2 * pi, length.out = num_frames)
#   frame_files <- vector(mode='list', length=length(angles))
#   
#   # Generate each frame
#   for (i in seq_along(angles)) {
#     theta <- angles[i]
#     eye_x <- center["x"]
#     eye_y <- radius * sin(theta)
#     eye_z <- radius * cos(theta)
#     fig_upd <- fig %>% layout(
#       scene = list
#       (camera = list(
#         center = center,
#         eye = list(x = eye_x, y = eye_y, z = eye_z),
#         up = list(x = 0, y = 0, z = 1)
#         )
#       )
#     )
#     # Save frame as PNG
#     if(i<10){fileNO<-paste0("00",i)}
#     else if(i<100){fileNO<-paste0("0",i)}
#     else{fileNO<-paste0("0",i)}
#     frame_files[i]<-paste0(temp_dir,"/img",fileNO,".png")
#     plotly::save_image(
#       fig_upd, 
#       file = frame_files[i], 
#       width = width, 
#       height = height, 
#       scale = dpi / 72
#     )
#   }
#   # Use magick to create GIF
    #saw memory overflow here
#   imgs <- list.files(temp_dir, full.names = TRUE)
#   img_list <- lapply(imgs, magick::image_read)
#   img_joined <- magick::image_join(img_list)
#   img_animated <- magick::image_animate(img_joined, fps = num_frames/duration)
#   magick::image_write(image = img_animated, path = output_file)
#   message("GIF saved to ", output_file)
#   file.remove(frame_files)
#   return(img_animated)
# }