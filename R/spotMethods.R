
#-------------------------------------------------------------------------------
#' Constructor for GraphSpace-class objects from spot coordinates
#'
#' \code{buildSpotSpace} constructs a \code{GraphSpace-class} object 
#' by integrating spatial spot coordinates with a raster background. This 
#' wrapper streamlines data preparation for graph and spatial workflows, 
#' such as in the \code{PathwaySpace} package.
#' 
#' @param spot_coord A data frame with spot coordinates. It must include 
#' row names and coordinates assigned to `x` and `y` column names.
#' @param raster_image A raster image serving as the background on which 
#' spot coordinates are mapped.
#' @param mar A single numeric value in \code{[0, 0.5]} for the graph margins.
#' @param flip.v Logical; whether to vertically flip the background image  
#' matrix (top-to-bottom) to align with the graph coordinate system.
#' @param flip.h Logical; whether to horizontally flip the background image  
#' matrix (left-to-right) to align with the graph coordinate system.
#' @param flip.y Logical; whether to flip the node coordinates along the y-axis.
#' @param flip.x Logical; whether to flip the node coordinates along the x-axis.
#' @param rotate.xy Logical; whether to rotate graph coordinates.
#' @param verbose Logical; whether to display detailed messages.
#' @param crop_coord `r lifecycle::badge("deprecated")` Deprecated since
#' SpotSpace 0.0.7; use \link[RGraphSpace]{cropGraphSpace} instead.
#' @param nrc `r lifecycle::badge("deprecated")` Deprecated since
#' SpotSpace 0.0.7; use \link[PathwaySpace]{buildPathwaySpace} instead.
#' 
#' @details
#' Spatial alignment between spot coordinates and background images 
#' often requires orientation adjustments due to differing coordinate 
#' conventions (e.g., top-left vs. bottom-left origins). \code{buildSpotSpace} 
#' automates the initial object creation and calls \code{normalizeGraphSpace} 
#' internally to handle these adjustments via \code{rotate.xy}, \code{flip.x}, 
#' and \code{flip.y}.
#' 
#' @return A \code{\link[RGraphSpace:GraphSpace-class]{GraphSpace}} object with 
#' integrated spatial and image data.
#' 
#' @author Sysbiolab Team
#' 
#' @seealso \link[RGraphSpace]{GraphSpace}, 
#' \link[RGraphSpace]{normalizeGraphSpace}, 
#' \link[PathwaySpace]{buildPathwaySpace}
#' @examples
#' # See examples in the PathwaySpace's online tutorials:
#' # https://sysbiolab.github.io/PathwaySpace/
#' 
#' @importFrom PathwaySpace buildPathwaySpace 
#' @importFrom RGraphSpace normalizeGraphSpace
#' @importFrom RGraphSpace GraphSpace getGraphSpace
#' @importFrom igraph make_empty_graph V 'V<-'
#' @importFrom grDevices is.raster as.raster col2rgb
#' @importFrom methods is
#' @importFrom RANN nn2
#' @importFrom Seurat SCTransform
#' @importFrom SeuratObject GetTissueCoordinates GetImage GetAssayData
#' @importFrom patchwork wrap_plots
#' @importFrom scales rescale
#' @importFrom lifecycle deprecated is_present deprecate_soft
#' @aliases buildSpotSpace
#' @export
buildSpotSpace <- function(spot_coord, raster_image, mar = 0.1, 
  flip.v = FALSE, flip.h = FALSE, flip.x = FALSE, flip.y = FALSE,
  rotate.xy = FALSE, verbose = TRUE, 
  crop_coord = deprecated(), 
  nrc = deprecated()){
  
  ### deprecate
  if (lifecycle::is_present(crop_coord)) {
    deprecate_soft("0.0.7", "buildSpotSpace(crop_coord)", 
      "cropGraphSpace()")
  }
  if (lifecycle::is_present(nrc)) {
    deprecate_soft("0.0.7", "buildSpotSpace(nrc)", 
      "buildPathwaySpace(nrc)")
  }
  
  if(!is.data.frame(spot_coord)){
    stop("'spot_coord' must be a 'data.frame'.",
      call. = FALSE)
  }
  if(!all(c("x","y") %in% colnames(spot_coord))){
    stop("'spot_coord' is missing 'x' and 'y' coordinates.",
      call. = FALSE)
  }
 
  #--- build GraphSpace-class
  gs <- GraphSpace(g = spot_coord, verbose = verbose)
  
  if(missing(raster_image)){
    gs <- normalizeGraphSpace(gs, mar = mar,  
      flip.x = flip.x, flip.y = flip.y, 
      rotate.xy = rotate.xy,
      verbose = verbose)
  } else {
    gs <- normalizeGraphSpace(gs, raster_image, mar = mar, 
      flip.x = flip.x, flip.y = flip.y, rotate.xy = rotate.xy, 
      flip.v = flip.v, flip.h = flip.h, 
      verbose = verbose)
  }
  
}
