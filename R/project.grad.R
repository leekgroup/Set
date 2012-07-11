project.grad <-
function(t, g, k, n, q, verbose=TRUE)
  {
    sigma.bar <- rep(Inf, length(t))
    sigma.bar[g<0] <- ((t-1)/g)[g<0]
    sigma.bar[g>0] <- (t/g)[g>0]

    ssigma.bar <- sort(unique(c(0,sigma.bar)))

    if (length(ssigma.bar) == 1)
      return(t)

    obj0 <- get.obj(t, k, n, q)
    OK <- FALSE

    if (verbose)
      message("gradient projection: ", length(ssigma.bar)-1, " intervals to check")

    p <- -g
    for (j in 2:length(ssigma.bar)) {
      p <- ifelse(ssigma.bar[j-1] < sigma.bar, -g, 0)
      univ.obj <- get.univ.obj.fn(t, p, k, n, q)
      res <- optimize(f=univ.obj,
                      interval=ssigma.bar[(j-1):j])
      if (verbose)
        message("j=", j, " obj0=", obj0, " obj=", res$objective)

      if (res$objective < obj0) {
        obj0 <- res$objective
        OK <- TRUE
      } else if (OK) {
        break
      }
    }
    if (!OK) return(t)
    return(t + res$minimum * p)
  }

