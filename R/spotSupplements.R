
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


