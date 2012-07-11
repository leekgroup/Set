get.obj <-
function(t, k, n, q)
  {
    res <- (1-k) * sum((n-q)[t==1])
    res1 <- res / sum(n[t==1])
    res <- k * sum(q[t==0])
    res2 <- res / sum(n[t==0])
    res1 + res2
  }

