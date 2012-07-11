knap <-
function(rho, k, n, q)
  {
    N <- sum(n)
    c1 <- (1-k)/rho * (n-q)
    c2 <- k/(N-rho) * q

    obj <- c1-c2
    o <- order(obj/n)

    keep1 <- n[o] <= rho
    ii <- which(keep1)

    s <- cumsum(n[o][ii])
    keep <- which(s <= rho)

    if (s[max(keep)] < rho) {
      new.rho <- s[max(keep)+1]
      return(knap(new.rho, k, n, q))
    }

    keep <- ii[keep]
    ids <- o[keep]

    if (length(ids) == length(n)) {
      ids <- ids[-length(ids)]
    }

    res <- rep(0, length(n))
    res[ids] <- 1
    res
  }

