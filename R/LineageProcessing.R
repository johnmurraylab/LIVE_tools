#' LineageProcessing
#' This file contains several utility functions that helps processing data in CD-like dataFrames

#' readEmbryoTable
#' @description
#' load a single embryo cell tracking file into data.framne
#'
#' @param file file directory of the embryo cell tracking table
#' @param time.file embryo imaging time record files
#'
#' @return a list of: `[[1]]` the embryo cell tracking table `[[2]]` Boolean, if a time file is loaded
#'
#' @import data.table
#' @export
readEmbryoTable <- function(file, time.file = NULL){
  temp_dir <- tempdir()
  if(grepl(pattern = ".zip$", file)){ #acetree output zip file
    zip_file <- unzip(file, list = T) #get valid timepoint records
    timepoint_files <- grep("^nuclei/t\\d{3}", zip_file[zip_file$Length>1,"Name"], value = TRUE)
    unzip(file, list = F, files = timepoint_files, exdir = temp_dir)
    tables <- lapply(timepoint_files, function(file_path) {
      timepoint <- as.integer(sub("^.*/t(\\d{3})-nuclei", "\\1", file_path))
      file_content <- data.table::fread(paste0(temp_dir,"/",file_path), header = F, sep = ",")
      file_content[, time := timepoint]
      timepoint <- file_content[
        !V10%in%c("polar1","polar2","") & !grepl("^Nuc",V10),
        .(
          x = V6,y = V7,z = V8, size = V9,
          cell = V10, time = time,
          gweight = V11, cross = V16, global = V17,
          local = V16-V18, blot = V16-V19
        )
      ]
      timepoint})
    CD_df <- data.table::rbindlist(tables) |> data.frame()
  }
  else{CD_df <- read.csv(file)}#assume .csv table
  timed <- F
  if(!is.null(time.file)){#modify the time attributes
    if(file.exists(time.file)){
      print(paste0("getting time data from: ", time.file))
      timed <- T
    }
    TIME <- read.table(time.file, header = F)[,2:3]
    colnames(TIME)<-c("time","realTime")
    CD_df <- merge(CD_df, TIME, by = "time", all.x = TRUE)
    CD_df$time <- CD_df$realTime/60
    CD_df <- CD_df[,!names(CD_df)%in%"realTime"]
  }
  return(list(CD_df, timed))
}

#' timeAvg
#' @description
#'  averages the selected attribute for all cells of the same time, can add corresponding standard deviation attribute
#' @param datList list of CD-like dataFrames
#' @param attribute the attribute to be averaged
#' @return dataFrame,
#' @export
timeAvg <- function(datList, attribute){
  datList |> dplyr::group_by(time) |>
      dplyr::summarise(stdev:=sd(!!sym(attribute)), !!attribute:=mean(!!sym(attribute)))
}

#' alignTime
#' @description
#' 'align' the cell time of a list of CD-like dataFrame by a specific cell stage
#' so that there is a "common" reference point (like common start time of lineage)
#' @param datList list of CD-like dataFrames
#' @param alignCell the cell used to find common time point to align
#' @param alignPoint a 0 to 1 ratio number showing what ratio of the alignCell life to align together (1 means the end, 0 means the start)
#' @param align_t the time of the `alignCell` at `alignPoint` after alignment, if the input is `"mean"` the align time will be the mean value of `alignPoint` across all embryos
#' @param alignBlot whether and method to align the blot values together default `FALSE` for not align, "mean" for aligning the mean values together (by addition/subtraction)
#' @return modified list of CD-like dataFrames
#' @export
alignTime <- function(datList, alignCell, alignPoint = 1, align_t = 0, alignBlot = F){
  rawTimes <- datList|>lapply(function(df){
    df[with(df, cell==alignCell),'time']|>quantile(alignPoint)})
  if(align_t=='mean'){align_t<-mean(unlist(rawTimes))}
  else if(align_t=='max'){align_t<-max(unlist(rawTimes))}
  embryos <- names(datList)
  times<-embryos |> lapply(function(embName){
    displace <- align_t-rawTimes[[embName]] #move all times by the same amount
    return(datList[[embName]][,'time']+displace)
  })
  if(alignBlot=="mean"){
    rawBlots <- datList|>lapply(function(df){df[with(df, expr = cell==alignCell),'blot']|>mean()})
    blotMean <- mean(unlist(rawBlots))
    blotDiff <- rawBlots|>lapply(function(val){blotMean-val})
    names(blotDiff)<-embryos
    datList <- embryos|>lapply(function(embName){
      df <- datList[[embName]]
      df$blot <- df$blot+blotDiff[[embName]]
      df
    })
    names(datList) <- embryos
  }
  for(i in 1:length(datList)){datList[[i]][,'time']<-times[[i]]}
  return(datList)
}

#' grepCells
#' @description
#' extract the data of given cells from the given CD table
#' @param CDData a CD-like dataFrame
#' @param cells cells to retrieve
#' @param lineages mother cell of lineages to retrieve (inclusive to mother)
#' @param times scope of time to select from
#' @param timesDiffThreshold define how much difference between the input `time` and retrieved `time` is allowed, if there are no exact `time` match, the function will look for closest `time` value
#'
#' @return if `dataReturn==TRUE`: dataFrame of selected cells, if `dataReturn==FALSE`: return the row index in the original `CDData`
#'
#' @export
grepCells <- function(CDData, cells=NULL, lineages=NULL,
                      times = "ALL", timesDiffThreshold = 1.5, dataReturn = TRUE){
  selectCells <- integer(0)
  CellREs <- NULL
  if((is.null(cells) && is.null(lineages)) || is.null(times)){
    if (dataReturn) {return(CDData[0, ])}
    return(integer(0))
  }
  else if("ALL"%in%cells || "ALL"%in%lineages){
    selectCells<-CDData|>rownames()|>as.integer()}
  else{
    CellREs <- cells|> lapply(function(line){paste0("^", gsub("x", "[a-z]", line),"$")})
    CellREs <- c(CellREs, lineageRE(lineages))
  }
  for (RE in CellREs){selectCells <- union(selectCells, grep(RE, CDData$cell))}

  if(!"ALL"%in%times){#by default, no filter by times
    timeMatches <- NULL
    for(t in times){
      if (t %in% CDData$time) { # Exact match found
        tMatches <- which(CDData$time == t)
      } else {# Find nearest time
        timeDiffs <- abs(CDData$time - t)
        minDiff <- min(timeDiffs)
        if(minDiff>timesDiffThreshold){tMatches=NULL}
        else{tMatches <- which(timeDiffs-minDiff < 0.001)}
      }
      timeMatches <- union(timeMatches, tMatches)
    }
    selectCells <- intersect(selectCells,timeMatches)
  }
  if(dataReturn){return(CDData[selectCells, ])}
  else{return(selectCells)}
}

#' RePosition
#' @description
#'  orientate one embryo at a given time point so that
#'  the embryo is centered at (0,0,0)
#'  coordinate the AP axis is aligned to the x-axis,
#'  DV to z-axis and LR to y-axis.
#'  Note that other time points are not adjusted - use totalRePosition to align all time points
#'
#'
#' @param CDFrame a CD-like dataFrame
#' @param time the time to align
#' @param indicatorP lineage on the relative posterior of embryo center
#' @param indicatorD lineage on the relative dorsal of embryo center, use NULL not NaN if not considering
#' @param indicatorV lineage on the relative ventral of embryo center, use NULL not NaN if not considering
#' @param indicatorL lineage on the relative left of embryo center, use NULL not NaN if not considering
#' @param indicatorR lineage on the relative right of embryo center, use NULL not NaN if not considering
#' @param xSize size of each unit on x-axis in original CDFrame
#' @param ySize size of each unit on y-axis in original CDFrame
#' @param zSize size of each unit on z-axis in original CDFrame
#'
#' @return the slice of the provided CDFrame at time==time, rotated
#' @export
RePosition <- function(CDFrame, time, indicatorP = "C", indicatorD = "Cxa",
                       indicatorV = "MSxxp", indicatorL = NULL, indicatorR = NULL,
                       xSize=0.08651, ySize=0.08651, zSize=0.5){
  CDFrame <- CDFrame |> grepCells(lineages = "ALL", times = time)
  CDFrame$x <- CDFrame$x*xSize
  CDFrame$y <- CDFrame$y*ySize
  CDFrame$z <- CDFrame$z*zSize
  #move the embryo to 0,0,0 by averaging the top 95% and bottom 5% values
  CDFrame$x <- CDFrame$x-mean(quantile(CDFrame$x, probs = c(0.05, 0.95)))
  CDFrame$y <- CDFrame$y-mean(quantile(CDFrame$y, probs = c(0.05, 0.95)))
  CDFrame$z <- CDFrame$z-mean(quantile(CDFrame$z, probs = c(0.05, 0.95)))

  transformMatrix <- rotationVec(CDFrame, indicatorP, indicatorD, indicatorV, indicatorL, indicatorR)

  coords <- as.matrix(CDFrame[, c("x", "y", "z")])
  transformed_coords <- coords %*% transformMatrix
  CDFrame$x <- transformed_coords[, 1]
  CDFrame$y <- transformed_coords[, 2]
  CDFrame$z <- transformed_coords[, 3]
  row.names(CDFrame)<- NULL
  return(CDFrame)
}

#' totalRePosition
#' @description
#'  center the embryo at 0,0,0 coordinate and rotate the embryo across all time points
#'  this uses the same strategy as RePosition and defines the orientation based on the specified
#'  time point, then applies the same rotation to all time points
#'
#' @param CDFrame a CD-like dataFrame
#' @param time the time to generate rotation parameters, nucleus across all time points will be rotated with the same parameter
#' @param indicatorP lineage on the relative posterior of embryo center
#' @param indicatorD lineage on the relative dorsal of embryo center, use NULL not NaN if not considering
#' @param indicatorV lineage on the relative ventral of embryo center, use NULL not NaN if not considering
#' @param indicatorL lineage on the relative left of embryo center, use NULL not NaN if not considering
#' @param indicatorR lineage on the relative right of embryo center, use NULL not NaN if not considering
#' @param xSize size of each unit on x-axis in original CDFrame
#' @param ySize size of each unit on y-axis in original CDFrame
#' @param zSize size of each unit on z-axis in original CDFrame
#'
#' @return the provided CDFrame rotated
#' @export
totalRePosition <- function(CDFrame, time, indicatorP = "C", indicatorD = "Cxa",
                            indicatorV = "MSxxp", indicatorL = NULL, indicatorR = NULL,
                            xSize=0.08651, ySize=0.08651, zSize=0.5){
  CDFrame$x <- CDFrame$x*xSize
  CDFrame$y <- CDFrame$y*ySize
  CDFrame$z <- CDFrame$z*zSize

  #align the embryo center to (0,0,0)
  positionAgg_FUN <- function(X){mean(quantile(X, probs = c(0.05, 0.95)))}
  centerCord <- aggregate(cbind(x,y,z)~time, data = CDFrame, FUN = positionAgg_FUN)
  cellCount <- count(CDFrame, time)
  t_start <- min(cellCount[cellCount[,"n"]>=50,"time"]) #embryos with >50 cells have more stable center coordination calculated by this method
  coord_start <- centerCord[centerCord[,"time"]==t_start,c("x","y","z")]
  centerCord[centerCord[,"time"]<t_start,c("x","y","z")] <- coord_start
  CDFrame <- left_join(CDFrame, centerCord, by = "time", suffix = c("", "_c"))
  CDFrame[,c("x","y","z")] <- CDFrame[,c("x","y","z")] - CDFrame[,c("x_c","y_c","z_c")]
  CDFrame <- CDFrame[,!names(CDFrame)%in%c("x_c","y_c","z_c")]
  CD_timeRef <- CDFrame |> grepCells(lineages = "ALL", times = time)

  transformMatrix <- rotationVec(CD_timeRef, indicatorP, indicatorD, indicatorV, indicatorL, indicatorR)
  coords <- as.matrix(CDFrame[, c("x", "y", "z")])
  transformed_coords <- coords %*% transformMatrix
  CDFrame$x <- transformed_coords[, 1]
  CDFrame$y <- transformed_coords[, 2]
  CDFrame$z <- transformed_coords[, 3]
  return(CDFrame)
}


rotationVec <- function(CDFrame, indicatorP, indicatorD, indicatorV, indicatorL, indicatorR){
  # 1) find the AB axis for this embryo
  pca1 <- prcomp(CDFrame[,c('x','y','z')]) #prcomp() always returns unit vectors
  APVect <- pca1$rotation[c(1,2,3)] #get the PC1 transformation UNIT vector
  PCells <- grepCells(CDFrame, lineages = indicatorP) #make sure that positive value for principle component 1 transformation vector is toward the A side
  PCellsVect <- c(mean(PCells$x), mean(PCells$y), mean(PCells$z))
  if(sum(PCellsVect*APVect)>0){APVect <- -APVect}
  # 2) find the unit transformation vector toward D side
  DVVect <- c(0,0,0)
  if(!is.null(indicatorD)){ #directional contribution from each indicator lineage
    DCells <- grepCells(CDFrame, lineages = indicatorD)
    DVVect_D <- c(mean(DCells$x), mean(DCells$y), mean(DCells$z)) |> OrthgonalProj(APVect)
    DVVect <- DVVect + DVVect_D
  }
  if(!is.null(indicatorV)){
    VCells <- grepCells(CDFrame, lineages = indicatorV)
    VDVect_V <- c(mean(VCells$x), mean(VCells$y), mean(VCells$z)) |> OrthgonalProj(APVect)
    DVVect <- DVVect - VDVect_V
  }
  if(!is.null(indicatorR)){
    RCells <- grepCells(CDFrame, lineages = indicatorR)
    DVVect_R <- c(mean(RCells$x), mean(RCells$y), mean(RCells$z)) |> cross_product(APVect)
    DVVect <- DVVect - DVVect_R
  }
  if(!is.null(indicatorL)){
    LCells <- grepCells(CDFrame, lineages = indicatorL)
    VDVect_L <- c(mean(LCells$x), mean(LCells$y), mean(LCells$z)) |> cross_product(APVect)
    DVVect <- DVVect + VDVect_L
  }
  #do not scale before adding together so that lineages further from AP make more contributions
  DVVect <- DVVect/sqrt(sum(DVVect^2))#scale into a unit vector
  # 3) use vector cross product to get RL transformation vector (+ side is the right side)
  RLVect <- cross_product(APVect, DVVect)
  transformMatrix <- cbind(APVect,RLVect,DVVect)
  return(transformMatrix)
}

#'OrthgonalProj
#'@description
#'Gram-Schmidt Process gives the orthogonal component to AP axis
OrthgonalProj <- function(v,u){
  vOrth <- v - (sum(v*u)*u)
}

#' cross_product
#' @description
#' vector cross-product for orthogonal vector
cross_product <- function(v, u) {
  c(v[2] * u[3] - v[3] * u[2],
    v[3] * u[1] - v[1] * u[3],
    v[1] * u[2] - v[2] * u[1])
}

#' blotRange return the range (min, max) of blot values in a given CD-like dataframe
#' @param datList a list of embryo cell tracking tables
#' @return a named list with a value 'max' and a value 'min'
#' @export
blotRange <- function(datList){
  maxBlot <- datalist$blot|>max()
  minBlot <- datalist$blot|>min()
  list(max = maxBlot, min = minBlot)
}

#' aggEmbryos merge a list of CD-like dataFrame into one dataFrame, with options to keep embryo identifier
#' @param datList list of CD-like dataFrames
#' @param labelEmbryo Boolean, whether to add a column to identify the embryo each cell come from
#' @return a dataframe combining all CD-like dataframes in the first dimension
#' @export
#' @import dplyr
aggEmbryos <- function(datList, labelEmbryo = TRUE){
  if(!is.data.frame(datList[[1]])){ #will do aggregation if a LargeList of DataFrame(s) is delivered
    return(datList)
  }
  if(labelEmbryo){
    out <- lapply(names(datList),
                  function(name) {
      datList[[name]] |> dplyr::mutate(embryo = name)
    })
  }
  else{out<-datList}
  out|>dplyr::bind_rows()
}

#' blotAvgFun
#' @description
#' get the average blot value for each cell, can choose a cutoff to select for top blot values only
#' @param datList a list of embryo cell tracking tables
#' @param cutoff the minimum percentile of blot values to consider (account for the underestimations during cell division)
blotAvgFun <- function(datList, cutoff = 0.5){
  CD_merged <- aggEmbryos(datList, labelEmbryo = TRUE)
  output <- CD_merged |>
    dplyr::group_by(cell, embryo) |>
    dplyr::gsummarise(
      avgBlot = mean(blot[blot >= quantile(blot, cutoff)], na.rm = TRUE)
    )
  output
}


#' lineageNames
#' @description
#' returns all possible starts of daughter lineages (including the mother cell itself) under a given lineage/mother cell
#' @param lineage a string or list of lineage name / mother cell name
#'
#' @return a list of regular expressions for cell
#' @export
#'
#' @examples lineageNames("EMS")
lineageNames <- function(lineage){
  P4 <- c("P4", "Z2", "Z3")
  P3 <- c("P3", "D", P4)
  P2 <- c("P2", "C", P3)
  EMS <- c("EMS", "MS", "E")
  E <- c("E$", "Ex")
  P1 <- c("P1", EMS, P2)
  P0 <- c("AB", P1)
  lineageDict <- list("P4" = P4, "P3" = P3, "P2" = P2, "P1" = P1, "P0" = P0, "EMS" = EMS, "E" = E)
  out<- lineage
  if(lineage %in% names(lineageDict)){out<-unlist(lineageDict[lineage])}
  else{out <- lineage}
  out
}

#' lineageRE returns all regular expressions needed to grep cells under a given lineage/mother cell
#'
#' @param lineage a string or list of lineage name / mother cell name, x as single-letter wildcard
#'
#' @return a list of regular expressions for cell
#' @export
#' @examples lineageRE(c("EMS", "C"))
lineageRE <- function(lineage){
  lineage<- lineage |> lapply(lineageNames)|>unlist()
  REs <- lineage |> lapply(function(line){paste0("^", gsub("x", "[a-z]", line))})
  return(REs)
}
