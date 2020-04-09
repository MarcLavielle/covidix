errpred2 <- function(pp,tau,y,dataIn,ilog,cw){
  pp[ilog] <- exp(pp[ilog])
  dataIn$individual_parameters$value=matrix(c(pp,tau), nrow=1)
  res = simulx(data=dataIn)
  y.pred1=res[[1]][,2]
  dy.pred1=diff(y.pred1)
  y.pred2=res[[2]][,2]
  dy.pred2=diff(y.pred2)
  s1=mean((y.pred1-y[[1]])^2)
  s2=mean((y.pred2-y[[2]])^2)
  sd1=mean((dy.pred1-y[[3]])^2)
  sd2=mean((dy.pred2-y[[4]])^2)
  n1 <- length(y.pred1)
  n2 <- length(y.pred2)
  e <- cw[1]*n1*log(s1) + cw[2]*(n1-1)*log(sd1) + cw[3]*n2*log(s2) + cw[4]*(n2-1)*log(s2) 
  return(e)
}

mlxoptim2 <- function(model=NULL, output=NULL, data=NULL, initial=NULL, ilog=NULL, cw=rep(1,4), tau=NULL, optim=T) {
  y.obs <- list(data[[1]]$value, data[[2]]$value, diff(data[[1]]$value), diff(data[[2]]$value))
  out1 <- list(name=output[1], time=data[[1]]$day)
  out2 <- list(name=output[2], time=data[[2]]$day)
  dd <- simulx(model=model,
               parameter=c(initial, tau),
               output=list(out1, out2),
               settings=list(data.in=TRUE,load.design=TRUE))
  l0 <- initial
  l0[ilog] <- log(l0[ilog])
  sse.ini <- errpred2(l0,tau,y.obs,dd,ilog=ilog,cw=cw)
  if (optim){ 
    r <- optim(l0, errpred2, y=y.obs, tau=tau, dataIn=dd, ilog=ilog, cw=cw) #, control=list(trace=1))
    pest=r$par
    pest[ilog]=exp(pest[ilog])
    sse.est <- r$value
  } else {
    sse.est <- sse.ini
    pest <- initial
  }
  final = list(parameter.ini=initial, sse.ini=sse.ini, par=pest, value=sse.est)
  return(final)
}

#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------

fitCovid <- function(data=NULL, country="France", nb.day=14,
                     M.max=3, d.tau=6, l.tau=6, cw=rep(1,4), dir=".", K1=10, ndt=5) {
  cw.ini <- c(1,0.5,0.1,0.1)
  out.name <- c("W", "D")
  file.out <- gsub(" ","",paste0(dir,"/",country,".RData"))
  # ilog <- c(1, 3, 4, 5, 6, 7, 8, 9)
  ilog <- 1:9
  delta <- 1
  country.in <- country
  Dc <- Do <- NULL
  Tau.est <- list()
  d0 <- subset(data, country==country.in)
  d0$country <- d0$percentage <- NULL
  d <- d0[,c('day', 'value','type')]
  d$day <- d$day-min(d$day)+1
  
  d1 <- subset(d, type=="confirmed")
  d2 <- subset(d, type=="deaths")
  
  dy01 <- d1$value[1]
  dy02 <- d2$value[1]
  n <- table(d0$type)
  d0$dy <- c( c(dy01, diff(d0$value[1:n[1]])), c(dy02, diff(d0$value[(n[1]+1):(n[1]+n[2])]) ))
  
  d <- list(d1, d2)
  out <- list(list(name="W", time=d1$day), list(name="D", time=d2$day) , list(name="ks", time=d1$day))
  
  model <- sir0
  W0=d1$value[1]
  D0=d2$value[1]
  
  p.ini0 <- c(k0=2, a1=0.01, kd=0.06, kr=0.4, d=2, I0=W0/5, W0=W0, D0=D0, L0=D0*2)
  # if (file.exists(file.out)) {
  #   load(file=file.out)
  #   p.ini[1:9] <- P[1, 1:9]
  # }  
  ilogM <- ilog[1:9]
  
  r <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini0, ilog=ilogM, cw=cw.ini, tau=NULL)
  r.min <- mlxoptim2(model=model, output=out.name, data=d, initial=r$par, ilog=ilogM, cw=cw, tau=NULL)
  v.min <- r.min$value
  for (k in (1:K1)) {
    pk <- p.ini0*(1+rnorm(length(p.ini0))*0.2)
    r <- mlxoptim2(model=model, output=out.name, data=d, initial=pk, ilog=ilogM, cw=cw.ini, tau=NULL)
    r <- mlxoptim2(model=model, output=out.name, data=d, initial=r$par, ilog=ilogM, cw=cw, tau=NULL)
    if (r$value < v.min) {
      v.min <- r$value
      p.est <- r$par
      r.min <- r
    }
  }
  r <- r.min
  # model <- sir1
  # p.ini1 <- c(r.min$par, AW=0, AD=0, phi=phi)
  # ilogM <- ilog
  r <- mlxoptim2(model=model, output=out.name, data=d, initial=r$par, ilog=ilogM, cw=cw, tau=NULL)
  print(c(v.min, r$value))
  p.est <- r$par
  res <- simulx(model=model,
                parameter=p.est,
                output=out,
                settings=list(load.design=TRUE))
  y1 <- d1$value
  y2 <- d2$value
  y.pred1=res[[1]][,2]
  
  s1=mean((y.pred1-y1)^2)
  y.pred2=res[[2]][,2]
  s2=mean((y.pred2-y2)^2)
  print(c(log(s1), log(s2), r$value))
  d[[1]]$pred <- y.pred1
  d[[2]]$pred <- y.pred2
  
  dd <- list()
  dd[[1]] <- data.frame(day=d[[1]]$day[2:nrow(d[[1]])], value=diff(d[[1]]$value), pred=diff(d[[1]]$pred))
  dd[[2]] <- data.frame(day=d[[2]]$day[2:nrow(d[[2]])], value=diff(d[[2]]$value), pred=diff(d[[2]]$pred))
  
  pl1 <- ggplot(data=d[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl2 <- ggplot(data=d[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl3 <- ggplot(data=dd[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl4 <- ggplot(data=dd[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl5 <- ggplot(data=res[['ks']]) + geom_line(aes(time,ks), color="red")
  grid.arrange(pl1, pl2, pl5, pl3, pl4, ncol=3)
  
  t.ori <- d[[1]]$day
  t.min <- min(t.ori)
  t.max <- max(t.ori) - t.min + 1
  tau.ini <- c(t.min-0.5, t.max+0.5)
  
  M.max <- min(M.max, floor((t.max-t.min)/l.tau))
  p0 <- r$par
  W0=d1$value[1]
  D0=d2$value[1]
  p0['W0']=W0
  p0['D0']=D0
  p0['L0']=max(p0['L0'], 1)
  p0['I0']=max(p0['I0'], 1)
  
  #r <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini1, ilog=ilogM, cw=c(0,1,0,1), tau=NULL, optim=F)
  
  R1 <-  r$value
  
  p.iniM <- c(p0, AW=0, AD=0, phi=0)
  q.iniM <- NULL
  r <- mlxoptim2(model=sir1, output=out.name, data=d, initial=p.iniM, ilog=ilogM, cw=cw, tau=q.iniM)
  R2 <-  r$value
  
  Tauc <- matrix(c(tau.ini, rep(0,M.max-1)), nrow=1)
  Pc <- matrix(c(r$par, rep(0,M.max-1)), nrow=1)
  if (M.max >1)
  colnames(Pc) <- c(names(r$par), paste0("a",2:M.max))
  else
    colnames(Pc) <- names(r$par)
  
  if (M.max>1) {
    for (M in (2:M.max)) {
      print(M)
      
      tau.ini <- (t.max-t.min+1)*(0:M)/M + t.min - 0.5
      p.ini <- c(p0, rep(0,(M-1)))
      #      p.ini <- c(p0, rep(0.1,(M-1)))
      ilogM <- ilog
      
      names(tau.ini) <- paste0("tau", 0:M)
      names(p.ini)[(length(p0)+1):length(p.ini)] <- paste0("a",2:M)
      p.ini['a2'] <- -p0['k0']*(1-exp(-p0['a1']*tau.ini[M]))/tau.ini[M]
      q.ini <- c(tau.ini[-c(1,M+1)], AW=0, AD=0, phi=0)
      eval(parse(text=paste0("model <- sir",M)))
      r <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini, ilog=ilogM, cw=cw, tau=q.ini)
      tau.est <- tau.ini
      r.min <- r$value
      p.est <- r$par
      r.c <- r
      
      test <- TRUE
      while (test) {
        r.min0 <- r.min
        for (m in (2:M)) {
          tau.c <- tau.est
          taum.min <- tau.est[m]
          for ( taum in seq((tau.est[m]-d.tau),(tau.est[m]+d.tau), by=delta) ) {
            if ( (taum-tau.est[m-1]>l.tau) & (tau.est[m+1]-taum>l.tau) ) {
              tau.c[m] <- c(tau1=taum)
              q.c <- c(tau.c[-c(1,M+1)], AW=0, AD=0, phi=0)
              #              r.c <- mlxoptim2(model=model, output=out.name, data=d, initial=r$par, ilog=ilogM, cw=cw.ini, tau=tau.c[-c(1,M+1)])
              r.c <- mlxoptim2(model=model, output=out.name, data=d, initial=r.c$par, ilog=ilogM, cw=cw, tau=q.c)
              if (r.c$value < r.min) {
                r.min <- r.c$value 
                p.est <- r.c$par
                taum.min <- taum
              }
            }
          }
          tau.est[m] <- taum.min
        }
        print(c(tau.est[2:M], r.min))
        if ( r.min0-r.min <= 1 )
          test <- FALSE
      }
      
      R1 <- c(R1, r.min)
      q.est <- c(tau.est[2:M], AW=0, AD=0, phi=0)
      res <- simulx(model=model,
                    parameter=c(p.est, q.est),
                    output=out,
                    settings=list(load.design=TRUE))
      y.pred1=res[[1]][,2]
      y.pred2=res[[2]][,2]
      dd1 <- data.frame(day=d[[1]]$day[2:nrow(d[[1]])], value=diff(d[[1]]$value), pred=diff(y.pred1))
      dd2 <- data.frame(day=d[[2]]$day[2:nrow(d[[2]])], value=diff(d[[2]]$value), pred=diff(y.pred2))
      dd1$e <- dd1$value-dd1$pred
      dd1$e <- dd1$e/sd(dd1$e)
      dd2$e <- dd2$value-dd2$pred
      dd2$e <- dd2$e/sd(dd2$e)
      dd <- rbind(dd1, dd2)
      omega <- 2*pi/7
      lmd1 <- lm(e ~ -1 + I(sin(omega*day)) + I(cos(omega*day)), data=dd)
      phi <- as.vector(atan(lmd1$coefficients[1]/lmd1$coefficients[2]))
      p.iniM <- c(p.est[1:9], AW=0, AD=0, phi=phi, p.est[10:length(p.est)])
      q.iniM <- tau.est[2:M]
      ilogM <- ilog
      r <- mlxoptim2(model=model, output=out.name, data=d, initial=p.iniM, ilog=ilogM, cw=cw, tau=q.iniM)
      
      
      R2 <- c(R2, r$value)
      Tauc <- rbind(Tauc, c(tau.est, rep(0,M.max-M)))
      p.est <- r$par
     # names(p.est)[(length(p0)+1):length(p.est)] <- paste0("a",2:M)
      Pc <- rbind(Pc, c(p.est, rep(0,M.max-M)))
      colnames(Pc)[1:length(p.est)] <- names(p.est)
      names(tau.est) <- paste0("tau",0:M)
      res <- simulx(model=model,
                    parameter=c(p.est, tau.est[2:M]),
                    output=out,
                    settings=list(load.design=TRUE))
      d[[1]]$pred <- res[[1]][,2]
      d[[2]]$pred <- res[[2]][,2]
      
      dd <- list()
      dd[[1]] <- data.frame(day=d[[1]]$day[2:nrow(d[[1]])], value=diff(d[[1]]$value), pred=diff(d[[1]]$pred))
      dd[[2]] <- data.frame(day=d[[2]]$day[2:nrow(d[[2]])], value=diff(d[[2]]$value), pred=diff(d[[2]]$pred))
      
      pl1 <- ggplot(data=d[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl2 <- ggplot(data=d[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl3 <- ggplot(data=dd[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl4 <- ggplot(data=dd[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl5 <- ggplot(data=res[['ks']]) + geom_line(aes(time,ks), color="red")
      grid.arrange(pl1, pl2, pl5, pl3, pl4, ncol=3)
      
    }
  }
  
  Tau <- Tauc
  P <- Pc
  BIC <- R2 + log(nrow(d[[1]]) + nrow(d[[2]]))*seq(1,2*M.max,by=2)
  M.est <- which.min(BIC)
  tau.est <- Tau[M.est, seq(1, M.est+1)]
  eval(parse(text=paste0("model <- sir",M.est)))
  
  np0 <- ncol(P)-nrow(P)+1
  p.est <- P[M.est, seq(1, np0+M.est-1)]
  #names(p.est) <- names(p0)
  # if (length(p.est)>np0)
  #   names(p.est)[(np0+1):length(p.est)] <- paste0("a",2:M.est)
  names(tau.est) <- paste0("tau",0:M.est)
  eval(parse(text=paste0("model <- sir",M.est)))
  #model <- sir1e
  
  delta.t <- 1/ndt
  tc <- seq(min(d1$day), max(d2$day) + nb.day, delta.t)
  outc <- list(name=c("W", "D", "ks"), time=tc)
  
  p.est0 <- p.est
  p.est0['AW'] <- p.est0['AD'] <- 0
  res0 <- simulx(model=model,
                 parameter=c(p.est0, tau.est[2:M.est]),
                 output=outc,
                 settings=list(load.design=TRUE))
  res <- simulx(model=model,
                parameter=c(p.est, tau.est[2:M.est]),
                output=outc)
  ks <- res[["ks"]][,2]
  t.min <- min(d0$day) - 1
  dc1 <- data.frame(day=res[["W"]][,1]+t.min, pred0=res0[["W"]][,2], pred=res[["W"]][,2], variable="y", type="confirmed")
  dc2 <- data.frame(day=res[["D"]][,1]+t.min, pred0=res0[["D"]][,2], pred=res[["D"]][,2], variable="y", type="deaths")
  dW0 <- dc1$pred0[(ndt+1):nrow(dc1)] - dc1$pred0[1:(nrow(dc1)-ndt)]
  dD0 <- dc2$pred0[(ndt+1):nrow(dc2)] - dc2$pred0[1:(nrow(dc2)-ndt)]
  dW <- dc1$pred[(ndt+1):nrow(dc1)] - dc1$pred[1:(nrow(dc1)-ndt)]
  dD <- dc2$pred[(ndt+1):nrow(dc2)] - dc2$pred[1:(nrow(dc2)-ndt)]
  dc3 <- data.frame(day=dc1$day[(ndt+1):nrow(dc1)], pred0=dW0, pred=dW, variable="dy", type="confirmed")
  dc4 <- data.frame(day=dc2$day[(ndt+1):nrow(dc1)], pred0=dD0, pred=dD, variable="dy", type="deaths")
  dc <- rbind(dc1, dc2, dc3, dc4)
  dc$date <- as.Date(dc$day, "2020-01-21")
  
  res0 <- simulx(model=model,
                 parameter=c(p.est0, tau.est[2:M.est]),
                 output=out,
                 settings=list(load.design=TRUE))
  res <- simulx(model=model,
                parameter=c(p.est, tau.est[2:M.est]),
                output=out)
  d1$variable="y"
  d1$pred0 <- res0[["W"]][,2]
  d1$pred <- res[["W"]][,2]
  d2$variable="y"
  d2$pred0 <- res0[["D"]][,2]
  d2$pred <- res[["D"]][,2]
  d <- rbind(d1, d2)[, c(1, 2, 5, 6, 4,3)]
  dd1 <- data.frame(day=d1$day[2:nrow(d1)], value=diff(d1$value), pred0=diff(d1$pred0), 
                    pred=diff(d1$pred), variable="dy", type="confirmed")
  dd2 <- data.frame(day=d2$day[2:nrow(d2)], value=diff(d2$value), pred0=diff(d2$pred0), 
                    pred=diff(d2$pred), variable="dy", type="deaths")
  dd <- rbind(dd1, dd2)
  d <- rbind(d, dd)
  d$day <- d$day+ t.min
  d$date <- as.Date(d$day, "2020-01-21")
  d$variable <- as.factor(d$variable)
  d$e <- d$value-d$pred
  d <- d %>% group_by(variable, type) %>% mutate(sd=sd(e))
  
  R0 <- data.frame(day=dc1$day, date=as.Date(dc1$day, "2020-01-21"), r0=ks/(p.est['kr']+p.est['kd']))
  
  Tau <- Tau + t.min 
  tau.est <- tau.est + t.min 
  print(BIC)
  save(Tau, P, R1, R2, BIC, cw, R0, d, dc, tau.est, file=file.out, version=2)
  
  return(list(file=file.out, p.est=p.est, tau.est=tau.est))
}
