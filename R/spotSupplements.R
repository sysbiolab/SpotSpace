
#-------------------------------------------------------------------------------
.nearest_nodes <- function(gs){
  gxy <- getGraphSpace(gs, "nodes")
  nnpg <- RANN::nn2(gxy[,c("x","y")], gxy[,c("x","y")], k=2)
  nn.idx <- nnpg$nn.idx[,2]
  nn.dists <- nnpg$nn.dists[,2]
  nn <- data.frame(from=gxy$name, to=gxy$name[nn.idx], dist=nn.dists)
  return(nn)
}
.edge_distances <- function(ps) {
  nodes <- getGraphSpace(ps, "nodes")
  edges <- getGraphSpace(ps, "edges")
  dx <- nodes[edges$vertex1, "x"] - nodes[edges$vertex2, "x"]
  dy <- nodes[edges$vertex1, "y"] - nodes[edges$vertex2, "y"]
  dist <- sqrt(dx^2 + dy^2)
  nn <- data.frame(from=edges$name1, to=edges$name2, dist=dist)
  return(nn)
}

#-------------------------------------------------------------------------------
.graphFromImageCoordinates <- function(coord, image, rotate.xy = TRUE,
  flip.y = TRUE, flip.x = FALSE, verbose = TRUE){
  
  # Check attributes
  if(!is.data.frame(coord) || is.null(rownames(coord)) ){
    stop("'coord' should be a rownamed data frame.", call. = FALSE)
  }
  if(is.null(colnames(coord)) || is.null(rownames(coord)) ){
    stop("'coord' should be a row- and col-named data frame.", call. = FALSE)
  }
  
  # Check attributes
  attr <- colnames(coord)
  if(!all(c("x","y") %in% attr)){
    stop("'coord' is missing 'x' and 'y' coordinates.", call. = FALSE)
  }
  attr <- attr[!attr%in%c("x","y")]
  
  # Set image as raster
  if(!is.raster(image)) image <- as.raster(image)
  
  # Rotated coordinates
  if(rotate.xy){
    if(verbose) message("Rotating xy-coordinates...")
    coord$x2 <- coord$y
    coord$y2 <- coord$x
  } else {
    coord$x2 <- coord$x
    coord$y2 <- coord$y
  }
  
  # Flip y-coordinates over image axis
  if(flip.y){
    if(verbose) message("Flipping y-coordinates...")
    y <- coord$y2
    y <- -(y - max(y)) + nrow(image) - max(y) + 1
    rg <- range(y)
    if(min(rg)<1 || max(rg)>nrow(image)){
      ms2 <- "Revise buildSpotSpace() arguments."
      if(rotate.xy){
        ms1 <- "Coordinate rotation/flip not compatible with the input image dimensions. "
      } else {
        ms1 <- "Coordinate flip not compatible with the input image dimensions. "
      }
      stop(paste0(ms1, ms2), call. = FALSE)
    }
    coord$y2 <- y
  }
  
  # Flip x-coordinates over image axis
  if(flip.x){
    if(verbose) message("Flipping x-coordinates...")
    x <- coord$x2
    x <- -(x - max(x)) + ncol(image) - max(x) + 1
    rg <- range(x)
    if(min(rg)<1 || max(rg)>ncol(image)){
      ms2 <- "Revise buildSpotSpace() arguments."
      if(rotate.xy){
        ms1 <- "Coordinate rotation/flip not compatible with the input image dimensions. "
      } else {
        ms1 <- "Coordinate flip not compatible with the input image dimensions. "
      }
      stop(paste0(ms1, ms2), call. = FALSE)
    }
    coord$x2 <- x
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
  flip.y = TRUE, flip.x = FALSE, verbose = TRUE){
  
  # Check attributes
  attr <- unique(colnames(coord))
  attr <- attr[!is.na(attr)]
  if(!all(c("x","y") %in% attr)){
    stop("'coord' is missing 'x' and 'y' coordinates.")
  }
  attr <- attr[!attr%in%c("x","y")]
  
  # Rotated coordinates
  if(rotate.xy){
    if(verbose) message("Rotating xy-coordinates...")
    coord$x2 <- coord$y
    coord$y2 <- coord$x
  } else {
    coord$x2 <- coord$x
    coord$y2 <- coord$y
  }
  
  # Flip y-coordinates
  if(flip.y){
    if(verbose) message("Flipping y-coordinates...")
    y <- coord$y2
    coord$y2 <- -(y - max(y)) - max(y) + 1
  }
  
  # Flip x-coordinates
  if(flip.x){
    if(verbose) message("Flipping x-coordinates...")
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
  
  if(length(attr)>0){
    for(name in attr){
      igraph::vertex_attr(g, name) <- coord[[name]]
    }
  }
  
  return(g)
  
}
