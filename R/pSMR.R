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
