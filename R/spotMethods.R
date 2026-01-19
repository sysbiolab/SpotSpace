
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
#' @param crop_coord An optional numeric vector of length four specifying a  
#' cropping region (xmin, xmax, ymin, ymax), with values in normalized 
#' coordinates [0,1].
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
#' @importFrom scales rescale
#' @aliases buildSpotSpace
#' @export
#' 
buildSpotSpace <- function(spot_coord, raster_image, mar = 0.1, 
  nrc = 500, rotate.xy = TRUE, flip.y = TRUE, flip.x = FALSE, 
  crop_coord = c(0, 1, 0, 1), verbose = TRUE) {
  
  if(verbose) message("Validating arguments...")
  
  #--- validate argument types
  .validate.spot.args("singleNumber", "mar", mar)
  .validate.spot.args("singleNumber", "nrc", nrc)
  .validate.spot.args("singleLogical", "rotate.xy", rotate.xy)
  .validate.spot.args("singleLogical", "flip.y", flip.y)
  .validate.spot.args("singleLogical", "flip.x", flip.x)
  .validate.spot.args("numeric_vec", "crop_coord", crop_coord)
  .validate.spot.args("singleLogical", "verbose", verbose)
  if(missing(raster_image) || is.null(raster_image)){
    g <- .graphFromCoordinates(coord = spot_coord, 
      rotate.xy = rotate.xy, flip.y = flip.y, 
      flip.x = flip.x)
    g_lt <- list(g=g, image=NULL)
  } else {
    if(!is.matrix(raster_image) && !is.raster(raster_image)){
      stop("'raster_image' should a raster image or matrix." )
    }
    g_lt <- .graphFromImageCoordinates(coord = spot_coord, 
      image = raster_image, rotate.xy = rotate.xy, flip.y = flip.y, 
      flip.x = flip.x)
  }
  if(length(crop_coord)!=4){
    stop("'crop_coord' should be a numeric vector of length = 4.")
  }
  if(any(crop_coord < 0) || any(crop_coord > 1)){
    stop("'crop_coord' should be in [0,1].")
  }
  #--- build GraphSpace-class
  gs <- GraphSpace(g = g_lt$g, image = g_lt$image, mar = mar, 
    verbose = verbose)
  
  #-------------------------------------------------------------------------------
  if( !missing(crop_coord) ){
    gs <- .crop_gspace(gs, crop_coord, mar)
  }

  ps <- buildPathwaySpace(gs, nrc = nrc, verbose = verbose)
  
  return(ps)
}

.crop_gspace <- function(gs, crop_coord, mar){
  
  # remove nodes
  nodes <- gs@nodes
  cx <- nodes$x >= crop_coord[1] & nodes$x <= crop_coord[2]
  cy <- nodes$y >= crop_coord[3] & nodes$y <= crop_coord[4]
  nodes <- nodes[ which(cx & cy), ]
  
  # remove edges
  idx <- (gs@edges$name1 %in% gs@nodes$name) & (gs@edges$name2 %in% gs@nodes$name) 
  gs@edges <- gs@edges[idx,]
  
  # center nodes
  nodes$x <- nodes$x - mean(range(nodes$x))
  nodes$y <- nodes$y - mean(range(nodes$y))
  from <- range(c(nodes$x, nodes$y))
  to <- c(mar, 1-mar)
  nodes$x <- scales::rescale(nodes$x, from = from, to=c(0,1))
  nodes$y <- scales::rescale(nodes$y, from = from, to=c(0,1))
  
  # update graph
  idx <- V(gs@graph)$name %in% rownames(nodes)
  gs@graph <- igraph::delete_vertices(gs@graph, which(!idx))
  idx <- match(rownames(nodes), V(gs@graph)$name)
  V(gs@graph)$x[idx] <- nodes$x
  V(gs@graph)$y[idx] <- nodes$y
  gs@nodes <- nodes
  
  # crop image
  if(gs@pars$image.layer){
    nrow_mat <- nrow(gs@image)
    ncol_mat <- ncol(gs@image)
    xmin <- max(1, floor(crop_coord[1] * ncol_mat) + 1)
    xmax <- min(ncol_mat, ceiling(crop_coord[2] * ncol_mat))
    ymin <- max(1, floor(crop_coord[3] * nrow_mat) + 1)
    ymax <- min(nrow_mat, ceiling(crop_coord[4] * nrow_mat))
    gs@image <- gs@image[ymin:ymax, xmin:xmax, drop = FALSE]
  }
  
  return(gs)
}

#-------------------------------------------------------------------------------
#' @title  getNearestNode
#' 
#' @description \code{getNearestNode} retrieves the nearest neighbor for 
#' each spot based on Euclidean distances.
#' 
#' @param ps A \code{\link[PathwaySpace]{PathwaySpace}} class object.
#' @seealso \code{\link[RANN]{nn2}}
#' @examples
#' # See examples in the SpotSpace's vignette:
#' vignette("SpotSpace")
#' 
#' @importFrom RANN nn2
#' @aliases getNearestNode
#' @export
#' 
getNearestNode <- function(ps){
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
