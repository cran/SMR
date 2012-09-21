#############################################################
# Function to calculate nodes and weights of Gauss-Legendre #
# quadrature, where n is the number of points. It uses the  #
# Golub-Welsh algorithm                                     #
#############################################################
GaussLegendre <- function(n) 
{
  n <- as.integer(n)
  if (n < 0) 
    stop("Must be a non-negative number of nodes!")
  if (n == 0) 
     return(list(x = numeric(0), w = numeric(0)))
  i  <- 1:n
  j   <- 1:(n-1)
  mu0 <- 2
  b <- j / (4 * j^2 - 1)^0.5
  A <- rep(0, n * n)
  A[(n + 1) * (j - 1) + 2] <- b
  A[(n + 1) * j] <- b
  dim(A) <- c(n, n)
  sd <- eigen(A, symmetric = TRUE)
  w <- rev(as.vector(sd$vectors[1, ]))
  w <- mu0 * w^2
  x <- rev(sd$values)
  return(list(nodes = x, weights = w))
}

