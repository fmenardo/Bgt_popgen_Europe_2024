sfs_TajD <- function (sfs1) 
{
  n <- length(sfs1) + 1
  if (n < 4) {
    warning("Tajima test requires at least 4 sequences")
    D <- NA
    return(D)
  }
  sfs_numbers <- seq(along=sfs1)
  khat <- sum(sfs_numbers*(n-sfs_numbers)*sfs1)/choose(n,2)
  S <- sum(sfs1)
  if (S==0) {
    warning("no segregating sites")
    D <- NA
    return(D)
  }
  tmp <- 1:(n - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (n + 1)/(3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  D <- (khat - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
  return(D)
}
