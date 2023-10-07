#############################################################
# Function to calculate nodes and weights of Gauss-Legendre #
# quadrature, where size is the number of points. It uses the  #
# Golub-Welsh algorithm                                     #
#############################################################
#' #' @importFrom Rcpp sourceCpp
#' #' @export
#' cGaussLegendre <- function(size)
#' {
#'   size <- as.integer(size)
#'   if (size < 0)
#'     stop("Must be a non-negative number of nodes!")
#'   if (size == 0)
#'      return(list(x = numeric(0), w = numeric(0)))
#'   aux <- aux_gausslegendre(size)
#'   w <- as.vector(aux$aux2)
#'   w <- 2 * w^2
#'   x <- aux$aux1
#'   return(list(nodes = x, weights = w))
#' }

