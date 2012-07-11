getgrad <-
function(t, rho, k, n, q)
  {
    s <- sum(n[t==1])
    sN <- sum(n)-s

    sq <- sum(q[t==1])
    sQ <- sum(q)-sq

    g1 <- (1-k)/s^2 * (sq*n - s*q)
    g2 <- k/sN^2 * (sQ*n - sN*q)
    g1+g2
  }

