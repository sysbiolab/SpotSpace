
#-------------------------------------------------------------------------------
.nearest.nodes <- function(gs){
  gxy <- getGraphSpace(gs, "nodes")
  nnpg <- nn2(gxy[,c("x","y")], gxy[,c("x","y")], k=2)
  nn.idx <- nnpg$nn.idx[,2]
  nn.dists <- nnpg$nn.dists[,2]
  nn <- data.frame(from=gxy$name, to=gxy$name[nn.idx], dist=nn.dists)
  return(nn)
}

#-------------------------------------------------------------------------------
.graphFromImageCoordinates <- function(coord, image, rotate.xy = TRUE,
  flip.y = TRUE, flip.x = FALSE){
  
  # Check attributes
  if(!is.data.frame(coord) || is.null(rownames(coord)) ){
    stop("'coord' should be a rownamed data frame.")
  }
  if(is.null(colnames(coord)) || is.null(rownames(coord)) ){
    stop("'coord' should be a row- and col-named data frame.")
  }
  
  # Check attributes
  attr <- colnames(coord)
  if(!all(c("x","y") %in% attr)){
    stop("'coord' is missing 'x' and 'y' coordinates.")
  }
  attr <- attr[!attr%in%c("x","y")]
  
  # Set image as raster
  if(!is.raster(image)) image <- as.raster(image)
  
  # Rotated coordinates
  if(rotate.xy){
    coord$x2 <- coord$y
    coord$y2 <- coord$x
  } else {
    coord$x2 <- coord$x
    coord$y2 <- coord$y
  }
  
  # Flip y-coordinates over image axis
  if(flip.y){
    y <- coord$y2
    coord$y2 <- -(y - max(y)) + nrow(image) - max(y) + 1
  }
  
  # Flip x-coordinates over image axis
  if(flip.x){
    x <- coord$x2
    coord$x2 <- -(x - max(x)) + ncol(image) - max(x) + 1
  }
  
  # Initialize a graph using 'spots' as vertices, with no edges
  g <- make_empty_graph(n = nrow(coord), directed = FALSE)
  V(g)$name <- rownames(coord)
  
  # Add coordinates as vertex attributes
  V(g)$x <- coord$x2
  V(g)$y <- coord$y2
  V(g)$nodeSize <- 1.2
  
  if(length(attr)>0){
    for(name in attr){
      igraph::vertex_attr(g, name) <- coord[[name]]
    }
  }
  
  res <- list(g=g, image=image)
  
  return(res)
  
}

#-------------------------------------------------------------------------------
.graphFromCoordinates <- function(coord, rotate.xy = TRUE,
  flip.y = TRUE, flip.x = FALSE){
  
  # Rotated coordinates
  if(rotate.xy){
    coord$x2 <- coord$y
    coord$y2 <- coord$x
  } else {
    coord$x2 <- coord$x
    coord$y2 <- coord$y
  }
  
  # Flip y-coordinates
  if(flip.y){
    y <- coord$y2
    coord$y2 <- -(y - max(y)) - max(y) + 1
  }
  
  # Flip x-coordinates
  if(flip.x){
    x <- coord$x2
    coord$x2 <- -(x - max(x)) - max(x) + 1
  }
  
  # Initialize a graph using 'spots' as vertices, with no edges
  g <- make_empty_graph(n = nrow(coord), directed = FALSE)
  V(g)$name <- rownames(coord)
  
  # Add coordinates as vertex attributes
  V(g)$x <- coord$x2
  V(g)$y <- coord$y2
  V(g)$nodeSize <- 0.5
  
  return(g)
  
}
