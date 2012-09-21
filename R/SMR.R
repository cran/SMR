dSMR <- function(x, n, nu, np = 32)
{
  nn <- length(x) 
  if (nn == length(n))   xx <- cbind(x, n) else
    if (length(n) == 1) xx <- cbind(x, rep(n, times = nn))
  if (nn == length(nu))   xx <- cbind(xx, nu) else
    if (length(nu) == 1) xx <- cbind(xx, rep(nu, times = nn))
  if (nn == length(np))   xx <- cbind(xx, np) else
    if (length(np) == 1) xx <- cbind(xx, rep(np, times = nn))
  dtched <- function(xx) return(dMR(xx[1], xx[2], xx[3], xx[4]))
  d <- apply(xx, 1, dtched)
  return(d)
}

pSMR <- function(x, n, nu, np = 32)
{
    nn <- length(x) 
    if (nn == length(n))   xx <- cbind(x, n) else
      if (length(n) == 1) xx <- cbind(x, rep(n, times = nn))
    if (nn == length(nu))   xx <- cbind(xx, nu) else
      if (length(nu) == 1) xx <- cbind(xx, rep(nu, times = nn))
    if (nn == length(np))   xx <- cbind(xx, np) else
      if (length(np) == 1) xx <- cbind(xx, rep(np, times = nn))
    dtched <- function(xx) return(pMR(xx[1], xx[2], xx[3], xx[4]))
    p <- apply(xx, 1, dtched)
    return(p)
}

qSMR <- function(p, n, nu, np = 32)
{
  nn <- length(p) 
  if (any(p>=1) | any(p<=0)) stop("Warning: probabilities must be between 0 and 1!")     
  if (nn == length(n))    xx <- cbind(p, n) else
    if (length(n) == 1)   xx <- cbind(p, rep(n, times = nn))
  if (nn == length(nu))   xx <- cbind(xx, nu) else
    if (length(nu) == 1)  xx <- cbind(xx, rep(nu, times = nn))
  if (nn == length(np))   xx <- cbind(xx, np) else
    if (length(np) == 1)  xx <- cbind(xx, rep(np, times = nn))
  dtched <- function(xx) return(qMR(xx[1], xx[2], xx[3], xx[4]))
  q <- apply(xx, 1, dtched)  
  return(q)
}

rSMR <- function(N, n,  nu = Inf)
{
  x <- rMR(N, n, nu)
  return(x)
}