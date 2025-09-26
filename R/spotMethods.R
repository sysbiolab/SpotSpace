
#' @title  Constructor of PathwaySpace-class Objects
#' 
#' @description \code{buildSpotSpace} is a constructor of
#' PathwaySpace-class objects.
#' 
#' @param spot_coord A data frame with spot coordinates. It must include 
#' row names and coordinates assigned to `x` and `y` column names.
#' @param raster_image A raster image serving as the background on which 
#' spot coordinates are mapped.
#' @param nrc A single positive integer indicating the number of rows and 
#' columns of a square image matrix onto which spot signals will be projected.
#' @param mar A single numeric value (in [0,1]) indicating the size of the 
#' outer margins as a fraction of the image matrix defined by `nrc`.
#' @param rotate.xy Logical; whether to rotate the `raster_image`.
#' @param flip.y Logical; whether to flip the `raster_image` along the y-axis.
#' @param flip.x Logical; whether to flip the `raster_image` along the x-axis.
#' @param verbose Logical; whether to display detailed messages.
#' @author Sysbiolab Team
#' @seealso \code{\link[PathwaySpace]{buildPathwaySpace}}
#' @examples
#' # See examples in the SpotSpace's vignette:
#' vignette("SpotSpace")
#' 
#' @importFrom PathwaySpace buildPathwaySpace
#' @importFrom RGraphSpace GraphSpace getGraphSpace
#' @importFrom igraph make_empty_graph V 'V<-'
#' @importFrom grDevices is.raster as.raster col2rgb
#' @importFrom methods is
#' @importFrom RANN nn2
#' @importFrom Seurat SCTransform
#' @importFrom SeuratObject GetTissueCoordinates GetImage GetAssayData
#' @importFrom patchwork wrap_plots
#' @aliases buildSpotSpace
#' @export
#' 
buildSpotSpace <- function(spot_coord, raster_image, mar = 0.1, nrc = 500, 
  rotate.xy = TRUE, flip.y = TRUE, flip.x = FALSE,
  verbose = TRUE) {
  
  if(verbose) message("Validating arguments...")
  
  #--- validate argument types
  .validate.spot.args("singleNumber", "mar", mar)
  .validate.spot.args("singleNumber", "nrc", nrc)
  .validate.spot.args("singleLogical", "rotate.xy", rotate.xy)
  .validate.spot.args("singleLogical", "flip.y", flip.y)
  .validate.spot.args("singleLogical", "flip.x", flip.x)
  .validate.spot.args("singleLogical", "verbose", verbose)
  
  #--- build GraphSpace-class
  g_lt <- .graphFromImageCoordinates(coord = spot_coord, 
    image = raster_image, rotate.xy = rotate.xy, flip.y = flip.y, 
    flip.x = flip.x)
  gs <- GraphSpace(g = g_lt$g, image = g_lt$image, mar = mar, 
    verbose = verbose)
  ps <- buildPathwaySpace(gs, nrc = nrc, verbose = verbose)
  return(ps)
}


#-------------------------------------------------------------------------------
#' @title  getNearestSpot
#' 
#' @description \code{getNearestSpot} retrieves the nearest neighbor for 
#' each spot based on Euclidean distances.
#' 
#' @param ps A \code{\link[PathwaySpace]{PathwaySpace}} class object.
#' @seealso \code{\link[RANN]{nn2}}
#' @examples
#' # See examples in the SpotSpace's vignette:
#' vignette("SpotSpace")
#' 
#' @importFrom RANN nn2
#' @aliases getNearestSpot
#' @export
#' 
getNearestSpot <- function(ps){
  if(!is(ps, "GraphSpace")){
    stop("'ps' should be either a 'GraphSpace' or 'PathwaySpace' class object.")
  }
  gxy <- getGraphSpace(ps, "nodes")
  nnpg <- nn2(gxy[,c("x","y")], gxy[,c("x","y")], k=2)
  nn.idx <- nnpg$nn.idx[,2]
  nn.dists <- nnpg$nn.dists[,2]
  nn <- data.frame(from=gxy$name, to=gxy$name[nn.idx], dist=nn.dists)
  return(nn)
}
