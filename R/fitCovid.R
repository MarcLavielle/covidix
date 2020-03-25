
errpred <- function(p, d, a, tau=NULL, model=NULL){
  y.pred=model(p, d$t, tau=tau)
  e=sum((log(y.pred+a) - log(d$y+a))^2)
  #  e=sum((y.pred - d$y)^2)
  return(e)
}

linear <- function(x, t, tau, deriv=FALSE) {
  M <- length(tau)-1
  a   <- x[1]
  b   <- x[2]
  c   <- x[3]
  d   <- 0
  if (M>1) 
    d <- c(d, x[4:(M+2)])
  f <- exp( (a/100)*(t^2) + b*t + c )
  k <- 2*(a/100)*t + b
  for (m in (1:M)) {
    jm <- which(t >= tau[m] )
    f[jm] <- f[jm] * exp( d[m]*(t[jm]-tau[m])^2 )
    k[jm] <- k[jm] + 2*d[m]*(t[jm]-tau[m])
  }
  if (deriv) 
    return(list(f=f, k=k))
  else
    return(f)
}

quadratic <- function(x, t, tau, deriv=FALSE) {
  M <- length(tau)-1
  a   <- x[1]
  b   <- x[2]
  c   <- x[3]
  d   <- x[4]
  e   <- 0
  if (M>1) 
    e <- c(e, x[5:(M+3)])
  lf <-  (a/1000)*(t^3) +(b/100)*(t^2) + c*t + d 
  k <- 3*(a/1000)*(t^2) +2*(b/100)*t + c
  for (m in (1:M)) {
    jm <- which(t >= tau[m] )
    lf[jm] <- lf[jm]  +  e[m]*(t[jm]-tau[m])^3 
    k[jm] <- k[jm] + 3*e[m]*(t[jm]-tau[m])^2
  }
  f <- exp(lf)
  if (deriv) 
    return(list(f=f, k=k))
  else
    return(f)
}

linearlogit <- function(x, t, tau, deriv=FALSE) {
  M <- length(tau)-1
  a   <- x[1]
  b   <- x[2]
  c   <- x[3]
  d   <- 0
  if (M>1) 
    d <- c(d, x[4:(M+2)])
  gamma <- x[length(x)]
  f <- exp( (a/100)*(t^2) + b*t + c )
  fM <- exp( (a/100)*(tau[M]^2) + b*tau[M] + c )
  k <- 2*(a/100)*t + b
  bM <- 2*(a/100)
  aM <- b
  for (m in (1:(M-1))) {
    jm <- which(t >= tau[m] )
    f[jm] <- f[jm] * exp( d[m]*(t[jm]-tau[m])^2 )
    k[jm] <- k[jm] + 2*d[m]*(t[jm]-tau[m])
    fM <- fM * exp( d[m]*(tau[M]-tau[m])^2 )
    aM <- aM - 2*d[m]*tau[m]
    bM <- bM + 2*d[m]
  }
  jM <- which(t >= tau[M] )
  betaM <- (aM+bM*tau[M])*exp(gamma*tau[M])/(gamma-aM-bM*tau[M])
  A <- fM*(1 + betaM*exp(-gamma*tau[M]))
  f[jM] <- A/(1 + betaM*exp(-gamma*t[jM]))
  k[jM] <- betaM*gamma/(betaM+exp(gamma*t[jM]))
  if (deriv) 
    return(list(f=f, k=k))
  else
    return(f)
}

fitCovid <- function(data="covid19.csv", country="France", model="linearlogit", p=0.005, ra=1, 
                     M.max=4, d.tau=4, nb.day=2, type=c("confirmed", "deaths")) {
  # model:  "linear", "quadradic", "linearlogit"
  
  model.in <- model
  country.in <- country
  Dc <- Do <- NULL
  for (type.in in type) {
    d <- subset(read.csv(data), country==country.in & type==type.in)
    d <- d[,c('date', 'day', 'value')]
  #  if (is.null(a.in)) 
      a <- sd(d$value)*ra
    S <- p*max(d$value)
    day0 <- min(which(d$value>S))
    dy0 <- ifelse(day0>1, d$value[day0-1], 0)
    d <- subset(d, value>S)
    names(d) <- c("date", "t", "y")
    
    t.ori <- d$t
    t.min <- min(t.ori)
    t.max <- max(t.ori) - t.min
    d$t <- d$t - t.min
    tau.ini <- c(-0.5, t.max+0.5)
    if (model.in=="linear") {
      model <- linear
      nls1 <- nls(y ~ exp((a/100)*t^2 + b*t + c), start=c(a=0, b=0.2, c=5), data=d)
      r <- optim(coef(nls1), errpred, model=model, d=d, a=a, tau=tau.ini)
      dM <- 1
    } else if (model.in=="quadratic") {
      model <- quadratic
      nls1 <- nls(y ~ exp((a/1000)*t^3 + (b/100)*t^2 + c*t + d), start=c(a=0, b=0, c=0.2, d=5), data=d)
      r <- optim(coef(nls1), errpred, model=model, d=d, a=a, tau=tau.ini)
      dM <- 2
    } else {
      model <- linearlogit
      nls1 <- nls(y ~ exp((a/100)*t^2 + b*t + c), start=c(a=0, b=0.2, c=5), data=d)
      r <- optim(coef(nls1), errpred, model=linear, d=d, a=a, tau=tau.ini)
      r$value <- Inf
      dM <- 1
    }
    p0 <- r$par
    R <-  r$value
    Tau <- c(tau.ini, rep(0,M.max-1))
    P <- c(p0, rep(0,M.max-1))
    for (M in (2:M.max)) {
      tau.ini <- (t.max+1)*(0:M)/M - 0.5
      p.ini <- c(p0, rep(0, M-1))
      if (model.in=="linearlogit")
        p.ini[length(p.ini)] <- 2
      r <- optim(p.ini, errpred, model=model, d=d, a=a, tau=tau.ini)
      tau.est <- tau.ini
      r.min <- r$value
      p.est <- r$par
      
      test <- TRUE
      while (test) {
        r.min0 <- r.min
        for (m in (2:M)) {
          tau.c <- tau.est
          taum.min <- tau.est[m]
          if (tau.est[m+1]-d.tau > tau.est[m-1]+d.tau )
            for ( taum in seq((tau.est[m-1]+d.tau),(tau.est[m+1]-d.tau), by=0.5) ) {
              tau.c[m] <- taum
              r.c <- optim(p.ini, errpred, model=model, d=d, a=a, tau=tau.c)
              if (r.c$value < r.min) {
                r.min <- r.c$value 
                p.est <- r.c$par
                taum.min <- taum
              }
            }
          tau.est[m] <- taum.min
        }
        if (r.min>=r.min0)
          test <- FALSE
      }
      R <- c(R, r.min)
      Tau <- rbind(Tau, c(tau.est, rep(0,M.max-M)))
      P <- rbind(P, c(p.est, rep(0,M.max-M)))
    }
    
    BIC <- nrow(d)*log(R/nrow(d)) + log(nrow(d))*seq(1,2*M.max,by=2)
    M.est <- which.min(BIC)
    tau.est <- Tau[M.est, seq(1, M.est+1)]
    p.est <- P[M.est, seq(1, dM+1+M.est)]
    
    d$dy <- diff(c(dy0, d$y))
    tc <- seq(min(d$t), max(d$t)+nb.day, by=0.1)
    fc <- model(p.est, tc, tau.est, deriv=TRUE)
    dc <- data.frame(day=tc + t.min, f=fc$f, df=fc$k*fc$f, type=type.in)
    dc$date <- as.Date(dc$day, "2020-01-21")
    d$date <- as.Date(d$t+t.min, "2020-01-21")
    d$type <- type.in
    Do <- rbind(Do, d)
    Dc <- rbind(Dc, dc)
  }
  Dmo <- melt(Do, id=list("date", "t", "type"))
  Dmo$group <- as.factor(with(Dmo, paste(type, variable)))
  levels(Dmo$group) <- c("Daily number of confirmed cases", "Total number of confirmed cases", 
                         "Daily number of deaths", "Total number of deaths")
  Dmo$group = factor(Dmo$group,levels(Dmo$group)[c(2, 1, 4, 3)])
  Dmc <- melt(Dc, id=list("date", "day", "type"))
  Dmc$group <- as.factor(with(Dmc, paste(type, variable)))
  levels(Dmc$group) <- c("Daily number of confirmed cases", "Total number of confirmed cases", 
                         "Daily number of deaths", "Total number of deaths")
  Dmc$group = factor(Dmc$group,levels(Dmc$group)[c(2, 1, 4, 3)])
  
  Dc$k <- Dc$df/Dc$f
  
  pl1 <- ggplot(Dc, aes(date, k)) + geom_line(color="red", size=0.75) + 
    xlab("") + ylab("") +  expand_limits(y=0) + facet_wrap( ~type , scales = "free_y", ncol=2) 
  
  
  pl2 <- ggplot(Dmo, aes(date, value)) + geom_point(color="blue") + 
    geom_line(data=Dmc, aes(date,value), color="red", size=0.75) + 
    xlab("") + ylab("") + expand_limits(y=0) + facet_wrap( ~group , scales = "free_y", dir="v", ncol=2)
  
  
  return(list(pl.fit=pl2, pl.rate=pl1, Dc=Dc, Do=Do))
}

