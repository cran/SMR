#############################################################
# Function to calculate nodes and weights of Gauss-Legendre #
# quadrature, where size is the number of points. It uses the  #
# Golub-Welsh algorithm                                     #
#############################################################
GaussLegendre_aux <- function(size)
{
  size <- as.integer(size)
  if (size < 0)
    stop("Must be a non-negative number of nodes!")
  if (size == 0)
     return(list(x = numeric(0), w = numeric(0)))
  savepoints <- c(16, 32, 64, 100, 150, 200)
  if (any(size == savepoints)) {
    points <- .points_GaussLegendre(size)
    return(list(nodes = points$nodes, weights = points$weights))
  } else {
    res <- gaussLegendre(size)
    return(list(nodes = as.vector(res[1,]), weights = as.vector(res[2,])))
  }
}



