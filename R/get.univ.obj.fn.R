get.univ.obj.fn <-
function(t, p, k ,n ,q)
  {
    a1 <- p%*%(n-q)
    b1 <- t%*%(n-q)
    c1 <- p%*%n
    d1 <- t%*%n

    a2 <- -p%*%q
    b2 <- sum(q) - t%*%q
    c2 <- -p%*%n
    d2 <- sum(n) - t%*%n

    function(sigma)
      {
        res <- (1-k) * (sigma*a1 + b1) / (sigma*c1 + d1)
        res <- res + k * (sigma*a2 + b2) / (sigma*c2 + d2)
        res
      }

  }

