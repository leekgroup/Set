postProb <-
    function(obs.stat,
             ##method = "Storey",
             null,
             B=NULL,pos=FALSE,
             p0=1,K=length(obs.stat),dim.basis=10,...)
{
    ##if the null parameter is a character, then assume it's a random generation mechanism and get null statistics
    if(is.character(null))
    {
        rand.gen <- get(null, mode = "function")

        if(is.null(B))
        {
            stop("Need B = nr. of simulations")
        }

        null <- rand.gen(B*length(obs.stat),...)
    }

    stat <- obs.stat
    stat0 <- null

    ## This function calculates a posterior probability
    ## that the alternative model is true based on
    ## observed statistics (stat) and null statistics
    ## (stat0), based on window sizes of (K)
    m <- length(stat)
    m0 <- length(stat0)
    if(is.null(B)){B <- m0/m}

    ##labels whether number comes from observed or null
    z <- c(rep(1,m),rep(0,m0))

    ##if statistics are all positive, we take logs
    if(pos) {xx <- c(log(stat),log(stat0))}
    if(!pos) {xx <- c(stat,stat0)}

    ##order labels and statistics
    z <- z[order(xx)] ##ordered labels
    xx <- xx[order(xx)] ##ordered statistics

    ##splits statistics into K intervals
    ##finds mean number of observed stat in each interval
    ##finds median statistic value in each interval
    d <- rep(NA,K)
    b <- rep(NA,K)
    nK <- floor( (m + m0) / K)

    for(i in 1:K) {
        d[i] <- mean(z[(1+(i-1)*(nK)):(i*(nK))])
        ##d[i] now has the mean number of observed statistics in interval i
        b[i] <- median(xx[(1+(i-1)*(nK)):(i*(nK))])
        ##b[i] now has the median statistic in interval i
    }

    ##print(max(b))
    fit <- gam(d~s(b,
                   k=dim.basis,
                   fx=FALSE,bs="tp"),weights=c(rep(nK,K)),family=binomial())
    ##par(mfrow=c(1,2))
    ##plot(d~b)
    ##lines(fitted.values(fit)~b,col="red",lwd=2)
    ##model the mean number of observed statistics in each interval as a smooth function of the median statistics in each interval
    ##"cr" means cubic spline regression
    if(!pos) {
        lr <- -exp(-predict(object=fit,newdata=data.frame(b=stat)))/B
    }
    if(pos) {
        lr <- -exp(-predict(object=fit,newdata=data.frame(b=log(stat))))/B
    }
    ##plot(lr)
    postprob <- ifelse((1 + p0*lr) > 0, (1 + p0*lr), 0)

    p1.y <- postprob

    return(p1.y)
}

