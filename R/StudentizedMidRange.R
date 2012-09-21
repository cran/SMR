# computes the CDF of the  normal midrange
# input: the quantile q in R and  the number of means n >= 2. 
pNMR <- function(q, n, np = 32)
{
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # interval 1: (-infty; q - 8]
    y  <- (q - 8) + (1 + x) / (x - 1) # from -1, 1 to -infty, q - 8
    aux <- pnorm(2*q-y)-pnorm(y)
    aux[aux <= 0] <- 1.0e-300
    fy <- log(dnorm(y)) + (n-1) * log(aux)
    fy <- exp(fy)
    fy <- n * (2 / (x - 1)^2) * fy 
    I  <- sum(fy * w) 
    # interval 2: [q-8; q]
    a  <- q - 8
    b  <- q
    y  <- (b - a) / 2 * x + (a + b) / 2 # from -1, 1 to (q - 8),  q
    aux <- pnorm(2*q-y)-pnorm(y)
    aux[aux <= 0] <- 1.0e-300
    fy <- log(dnorm(y)) + (n-1) * log(aux)
    fy <- exp(fy)
    fy <- (b - a) / 2 * n * fy 
    I  <- I + sum(fy * w) 
    return(I) 
}

# auxiliar function 2
pMR_aux2 <- function(s, q, n, nu, np = 32)
{
   
    Ii <- pNMR(s * q, n, np)    
    fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + (nu - 1) * log(s) - nu * s * s / 2
    fx <- exp(fx)
    I  <- fx * Ii
    return(I)  
}


# computes the CDF of the  externally studentized normal midrange
# input: the quantile q in R and  the number of means n >= 2,
# degrees of freedom nu > 0, and number of points of the 
# gaussian quadrature (np)
pMR <- function(q, n, nu, np = 32)
{
    if (nu == Inf) return(pNMR(q, n, np))
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # Integration interval from 0 to 1
    y  <- 0.5 * x + 0.5 # from [-1; 1] to [0; 1]
    y  <- matrix(y, np, 1)
    fy <-  0.5 * apply(y, 1, pMR_aux2, q, n, nu, np) 
    I  <- sum(fy * w)
    # integration interval from 1 to +infty
    y <- 1+(1+x)/(1-x) # from [-1; 1] to [1; +infty)
    y <- matrix(y, np, 1)
    fy.aux <- apply(y, 1, pMR_aux2, q, n, nu, np) 
    fy <- log(fy.aux)+log(2)-2*log(1-x)
    fy <- exp(fy)
    I <- I+sum(fy * w) 
    return(I)
}

# computes the PDF of the normal midrange
# input: the quantile q in R and  the number of means n >= 2,
# degrees of freedom nu > 0 and number of points of the 
# gaussian quadrature (np). We use two divisions of the 
# integration limits. The first division was between inf and q-8,
# and the second division was between q-8 and q.
dNMR <- function(q, n, np = 32)
{
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # interval 1: (-infty; q - 8]
    y  <- (q - 8) + (1 + x) / (x - 1) # from [-1; 1] to (-infty; q - 8]
    aux <- pnorm(2 * q - y) - pnorm(y)
    aux[aux <= 0] <- 1.0e-300
    fy <- log(dnorm(y)) + (n - 2) * log(aux) + log(dnorm(2 * q - y)) + log(2) +
          log(n) + log(n-1)
    fy <- exp(fy)
    fy <- (2 / (x - 1)^2) * fy 
    I  <- sum(fy * w) 
    # interval 2: [q-8; q]
    a  <- q - 8
    b  <- q
    y  <- (b - a) / 2 * x + (a + b) / 2 #from [-1; 1] to [(q - 8);  q]
    aux <- pnorm(2 * q - y) - pnorm(y)
    aux[aux <= 0] <- 1.0e-300
    fy <- log(dnorm(y)) + (n - 2) * log(aux) + log(dnorm(2 * q - y)) + log(2) +
          log(n) + log(n-1)
    fy <- exp(fy)
    fy <- (b - a) / 2 * fy 
    I  <- I + sum(fy * w) 
    return(I) 
}

# auxiliar function 3
pMR_aux3 <- function(s, q, n, nu, np = 32)
{
   
    Ii <- dNMR(s * q, n, np)    
    fx <- nu/2 * log(nu) - lgamma(nu/2) - (nu/2-1) * log(2) + 
          (nu - 1) * log(s) - nu * s * s / 2
    fx <- exp(fx) * s
    I  <- fx * Ii
    return(I)  
}

# computes the PDF of the  externally studentized normal midrange
# input: the quantile q in R and the number of means n >= 2,
# degrees of freedom nu > 0,  degrees of freedom nu > 0 
# and number of points of the gaussian quadrature (np). 
# We use two divisions of the integration limits. 
# The first division was between inf and q-8, 
# and the second division was between q-8 and q.
dMR <- function(q, n, nu, np = 32)
{
    if (nu == Inf) return(dNMR(q, n, np))
    xx <- GaussLegendre(np)
    x  <- xx$nodes
    w  <- xx$weights
    # integration limit between 0 and 1
    y  <- 0.5 * x + 0.5 # changing the interval [-1,1] to [0,1]
    y  <- matrix(y, np, 1)
    fy <-  0.5 * apply(y, 1, pMR_aux3, q, n, nu, np) 
    I  <- sum(fy * w) 
    # integration limit [1,Infty] - transformation y <- -ln(x), 0 < x < exp(-1)
    b <- exp(-1)
    a <- 0
    x <- (b - a) / 2 * x + (b + a) / 2 # changing the interval [-1,1] to [0,exp(0,-1)]        
    y  <- -log(x) # # changing the interval [0,exp(-1)] to [1,infty]
    y  <- matrix(y, np, 1)
    fy <- (b - a) / 2 * apply(y, 1, pMR_aux3, q, n, nu, np) / x 
    I  <- I + sum(fy * w) 
    return(I)  
}


# computes the quantiles of the normal midrange
# inputs: the cumulative probability (0 < p < 1), 
# the number of means (n >= 2) n and number 
# of points of the gaussian quadrature (np)
qNMR <- function(p, n, np = 32)
{
   if (p < 0.5) q0 <- -0.5 else
   if (p > 0.5) q0 <-  0.5 else q0 <- 0 # arbitrary initial value
   found <- F
   maxIt <- 5000
   it <- 0
   eps <- 1e-13
   while ((found == F) & (it <= maxIt))
   {
      q1 <- q0 - (pNMR(q0, n, np) - p) / dNMR(q0, n, np)
      if (abs(q1-q0) <= eps) found = T
      it <- it + 1
      q0 <- q1
   }
   return(q1)
}

# computes the quantiles of the  externally studentized normal midrange
# inputs: the cumulative probability (0 < p < 1), 
# the number of means (n >= 2) n, degrees of freedom nu > 0 
# and number of points of the gaussian quadrature (np)
qMR <- function(p, n, nu, np = 32)
{
   if (nu == Inf) return(qNMR(p, n, np))  
   if (p < 0.5) q0 <- -0.5 else
   if (p > 0.5) q0 <-  0.5 else q0 <- 0 # arbitrary initial value
   found <- F
   maxIt <- 20
   it <- 0
   eps <- 1e-13
   while ((found == F) & (it <= maxIt))
   {
      q1 <- q0 - (pMR(q0, n, nu, np) - p) / dMR(q0, n, nu, np)
      if (abs(q1-q0) <= eps) found = T
      it <- it + 1
      q0 <- q1
   }
   return(q1)
}

# computes the vector of random numbers of the normal midrange or
#  externally studentized normal midrange.
# Inputs -  N: size of the vector to be simulated, n: sample size
# nu in the interval [0, Inf], is the degrees of freedom.
rMR <- function(N, n, nu = Inf)
{
    if (nu == Inf) X <- (matrix(rnorm(N * n), N, n)) else
    X <- (matrix(rnorm(N * n), N, n)) / (rchisq(N, nu) / nu)^0.5 
    midrange <- function(x) 
    {
       return((max(x)+ min(x)) / 2)
    }
    res <- apply(X, 1, midrange)
    return(res)
}

