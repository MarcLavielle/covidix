errpred2 <- function(pp,args,y,dataIn,ilog,cw){
  pp.sir <- pp[!(names(pp) %in% c('AW', 'AD', 'phi'))]
  if (length(pp)>length(pp.sir)) {
    AW <- exp(pp['AW'])
    AD <- exp(pp['AD'])
    phi <- pp['phi']
  } else {
    AW = AD = phi = 0
  }
  omega <- 2*pi/7
  pp.sir[ilog] <- exp(pp.sir[ilog])
  if (pp.sir['d']>7)
    e <- Inf
  else {
    dataIn$individual_parameters$value=matrix(c(pp.sir,args), nrow=1)
    res = simulx(data=dataIn)
    y.pred1=res[[1]][,2]
    t1 <- res[[1]][-1,1]
    dy.pred1=diff(y.pred1)*(1 + AW*cos(omega*t1 + phi))
    y.pred1 <-  cumsum(c(y.pred1[1], dy.pred1))
    y.pred2=res[[2]][,2]
    t2 <- res[[2]][-1,1]
    dy.pred2=diff(y.pred2)*(1 + AD*cos(omega*t2 + phi))
    y.pred2 <- cumsum(c(y.pred2[1], dy.pred2))
    s1=mean((y.pred1-y[[1]])^2)
    s2=mean((y.pred2-y[[2]])^2)
    sd1=mean((dy.pred1-y[[3]])^2)
    sd2=mean((dy.pred2-y[[4]])^2)
    n1 <- length(y.pred1)
    n2 <- length(y.pred2)
    e <- cw[1]*n1*log(s1) + cw[2]*(n1-1)*log(sd1) + cw[3]*n2*log(s2) + cw[4]*(n2-1)*log(sd2)
  }
  return(e)
}

mlxoptim2 <- function(model=NULL, output=NULL, data=NULL, initial=NULL, ilog=NULL, cw=rep(1,4), args=NULL, optim=T) {
  y.obs <- list(data[[1]]$value, data[[2]]$value, diff(data[[1]]$value), diff(data[[2]]$value))
  out1 <- list(name=output[1], time=data[[1]]$day)
  out2 <- list(name=output[2], time=data[[2]]$day)
  ini.sir <- initial[!(names(initial) %in% c('AW', 'AD', 'phi'))]
  dd <- simulx(model=model,
               parameter=c(ini.sir, args),
               output=list(out1, out2),
               settings=list(data.in=TRUE,load.design=TRUE))
  l0 <- initial
  l0[ilog] <- log(l0[ilog])
  if ('AW' %in% names(l0)) {
    l0['AW'] <- log(initial['AW']+0.001)
    l0['AD'] <- log(initial['AD']+0.001)
  }
  sse.ini <- errpred2(l0,args,y.obs,dd,ilog=ilog,cw=cw)
  if (optim){ 
    r <- tryCatch( {
      optim(l0, errpred2, y=y.obs, args=args, dataIn=dd, ilog=ilog, cw=cw) #, control=list(trace=1))
    },
    error=function(cond) {
      return(list(value=Inf, par=NULL))
    })
    pest=r$par
    if (!is.null(pest))
      pest[ilog]=exp(pest[ilog])
    # sse.est <- r$value
    # r <- optim(l0, errpred2, y=y.obs, args=args, dataIn=dd, ilog=ilog, cw=cw) #, control=list(trace=1))
    # pest=r$par
    # pest[ilog]=exp(pest[ilog])
    if ('AW' %in% names(l0)) {
      pest['AW'] <- exp(pest['AW'])
      pest['AD'] <- exp(pest['AD'])
    }
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

fitCovid <- function(data=NULL, country=NULL, nb.day=14, estim.tau=T, estim.p0=T, M=NULL,
                     M.max=3, d.tau=6, l.tau=6, cw=c(1,1,1,1), dir=".", K1=5, ndt=5, file.out=NULL) {
  
  if (is.null(country))
    country <- as.character(data$country[1])
  
  if (is.null(file.out))
    file.out <- gsub(" ","",paste0(dir,"/",country,".RData"))
  if (!estim.tau) {
    if (file.exists(file.out) ) {
      load(file.out)
      #P.ini <- P
      Tau.ini <- Tau[c(1,3,7),] - Tau[1,1] + 0.5
      #estim.p0 <- F
    } else 
      estim.tau <- T
  }
  if (!estim.p0) {
    if (file.exists(file.out) ) {
      load(file.out)
      p.ini.1B <- P[1,][1:12]
    } else 
      estim.p0 <- T
  }
  if (!is.null(M))
    M.est <- M.max <- M
  else
    M.est <- NULL
  
  period <- 7
  omega <- 2*pi/period
  # ilog <- c(1, 3, 4, 5, 6, 7, 8, 9)
  # cw.ini <- c(1,0.5,0.1,0.1)
  cw.ini <- cw
  out.name <- c("W", "D")
  ilog <- 1:6
  delta <- 2
  country.in <- country
  Dc <- Do <- NULL
  Tau.est <- list()
  d0 <- subset(data, country==country.in)
  d0$country <- d0$percentage <- NULL
  d <- d0[,c('day', 'value','type')]
  d$day <- d$day-min(d$day)+1
  day.max <- max(d$day)
  
  d1 <- subset(d, type=="confirmed")
  d2 <- subset(d, type=="deaths")
  
  # dy01 <- d1$value[1]
  # dy02 <- d2$value[1]
  # n <- table(d0$type)
  # d0$dy <- c( c(dy01, diff(d0$value[1:n[1]])), c(dy02, diff(d0$value[(n[1]+1):(n[1]+n[2])]) ))
  
  d <- list(d1, d2)
  out <- list(list(name="W", time=d1$day), list(name="D", time=d2$day) , list(name="ks", time=d1$day))
  
  model <- sir1
  
  W0=d1$value[1]
  D0=d2$value[1]
  args1 <- c(tmax=day.max, W0=W0, D0=D0, L0=D0/2)
  if (estim.p0) {
    
    #p.ini0 <- c(k0=2, a1=0.5, kd=0.4, kr=0.2, d=6, I0=W0/2)
    p.ini0 <- c(k0=2, a1=0.01, kd=0.06, kr=0.4, d=2, I0=W0/2)
    # if (file.exists(file.out)) {
    #   load(file=file.out)
    #   p.ini[1:9] <- P[1, 1:9]
    # }  
    
    r <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini0, ilog=ilog, cw=cw.ini, args=args1)
    r.min <- mlxoptim2(model=model, output=out.name, data=d, initial=r$par, ilog=ilog, cw=cw, args=args1)
    v.min <- r.min$value
    for (k in (1:K1)) {
      pk <- p.ini0*(1+rnorm(length(p.ini0))*0.2)
      r <- mlxoptim2(model=model, output=out.name, data=d, initial=pk, ilog=ilog, cw=cw.ini, args=args1)
      r <- mlxoptim2(model=model, output=out.name, data=d, initial=r$par, ilog=ilog, cw=cw, args=args1)
      if (r$value < v.min) {
        v.min <- r$value
        p.est <- r$par
        r.min <- r
      }
    }
    r <- r.min
    p.ini.1A <- r$par
    # p.ini.1A['W0']=d1$value[1]
    # p.ini.1A['D0']=d2$value[1]
    # p.ini.1A['L0']=max(p.ini.1A['L0'], 1)
    # p.ini.1A['I0']=max(p.ini.1A['I0'], 1)
    
  } else {
    p.ini.1A <- p.ini.1B[1:6]
  }
  
  dn.tot <- nrow(d[[1]]) + nrow(d[[2]]) 
  coef.bic <- log(dn.tot)
  
  t.ori <- d[[1]]$day
  t.min <- min(t.ori)
  t.max <- max(t.ori) - t.min + 1
  tau.ini <- c(t.min-0.5, t.max+0.5)
  M.max <- min(M.max, floor((t.max-t.min)/l.tau))
  
  rA <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini.1A, ilog=ilog, cw=cw, args=args1)
  p.est.1A <- rA$par
  #devA <- mlxoptim2(model=model, output=out.name, data=d, initial=p.est.1A, ilog=ilog, cw=c(0,1,0,1), args=args1, optim=F)$value
  devA <- rA$value
  
  res <- simulx(model=model,
                parameter=c(p.est.1A, args1),
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
  omega <- 2*pi/period
  lmd1 <- lm(e ~ -1 + I(sin(omega*day)) + I(cos(omega*day)), data=dd)
  phi <- as.vector(-atan(lmd1$coefficients[1]/lmd1$coefficients[2]))
  
  p.ini.1B <- c(p.est.1A, AW=0.2, AD=0.2, phi=phi)
  rB <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini.1B, ilog=ilog, cw=cw, args=args1)
  p.est.1B <- rB$par
  #  devB <- mlxoptim2(model=model, output=out.name, data=d, initial=p.est.1B, ilog=ilog, cw=c(0,1,0,1), args=args1, optim=F)$value
  devB <- rB$value
  DEV <- c(devA, devB)
  
  tau1A <- tau1B <- c(tau.ini, rep(0,M.max-1))
  Tau <- rbind(tau1A, tau1B)
  colnames(Tau) <- paste0("tau",0:M.max)
  pa1 <- c(p.est.1A, AW=0, AD=0, phi=0, a0=0, a2=0, args1)
  pb1 <- c(p.est.1B, a0=0, a2=0, args1) 
  P <- rbind(pa1, pb1)
  
  DIM <- apply(P, MARGIN=1, function(x) {sum(x !=0)}) + apply(Tau, MARGIN=1, function(x) {sum(x !=0)})
  BIC <- DEV + coef.bic*DIM + coef.bic*2*(P[,'AW'] !=0)
  print(rbind(DEV, BIC))
  
  i.model <- which.min(BIC)
  p.est <- P[i.model,]
  res <- simulx(model=model,
                parameter=c(p.est),
                output=out,
                settings=list(load.design=TRUE))
  
  y1 <- d1$value
  y2 <- d2$value
  y.pred1=res[[1]][,2]
  
  s1=mean((y.pred1-y1)^2)
  y.pred2=res[[2]][,2]
  s2=mean((y.pred2-y2)^2)
  d[[1]]$pred <- y.pred1
  d[[2]]$pred <- y.pred2
  
  dd <- list()
  dd[[1]] <- data.frame(day=d[[1]]$day[-1], value=diff(d[[1]]$value), pred=diff(d[[1]]$pred))
  dd[[2]] <- data.frame(day=d[[2]]$day[-1], value=diff(d[[2]]$value), pred=diff(d[[2]]$pred))
  dd[[2]] <- data.frame(day=d[[2]]$day[-1], value=diff(d[[2]]$value), pred=diff(d[[2]]$pred))
  dd[[1]]$pred <- dd[[1]]$pred*(1 + p.est['AW']*cos(omega*dd[[1]]$day + p.est['phi']))
  dd[[2]]$pred <- dd[[2]]$pred*(1 + p.est['AD']*cos(omega*dd[[2]]$day + p.est['phi']))
  d[[1]]$pred <- cumsum(c(d[[1]]$pred[1], dd[[1]]$pred))
  d[[2]]$pred <- cumsum(c(d[[2]]$pred[1], dd[[2]]$pred))
  res[['ks']]$Reff <- res[['ks']]$ks/(p.est['kr']+p.est['kd'])
  
  pl1 <- ggplot(data=d[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl2 <- ggplot(data=d[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl3 <- ggplot(data=dd[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl4 <- ggplot(data=dd[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
  pl5 <- ggplot(data=res[['ks']]) + geom_line(aes(time,Reff), color="red")
  grid.arrange(pl1, pl2, pl5, pl3, pl4, ncol=3)
  
  if (M.max>1) {
    for (M in (2:M.max)) {
      print(M)
      
      if (estim.tau) {
        tau.ini <- (t.max-t.min+1)*(0:M)/M + t.min - 0.5
        if (M==2) {
          tau.ini[2] <- tau.ini[3] - l.tau
        }  else {
          tau.ini[2] <- tau.ini[1] + l.tau
          tau.ini[3] <- tau.est[2] 
        }
      } else 
        tau.ini <- Tau.ini[M, 1:(M+1)]
      names(tau.ini) <- paste0("tau", 0:M)
      if (M==2) {
        p.ini <- p.est.1A
        argsM <- c(args1,tau.ini[-c(1,M+1)], a0=0, a2=0)
      }  else {
        p.ini <- p.est.A
        argsM <- c(args1, tau.ini[-c(1,M+1)], a0=0, a2=0)
      }
      
      eval(parse(text=paste0("model <- sir",M)))
      r <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini, ilog=ilog, cw=cw, args=argsM)
      tau.est <- tau.ini
      r.min <- r$value
      p.estM <- r$par
      r.c <- r
      #      r.c$par['D0'] <- max(r.c$par.c['D0'] ,1)
      #      r.c$par['L0'] <- max(r.c$par.c['L0'] ,1)
      
      test <- estim.tau
      #   test <- F
      while (test) {
        r.min0 <- r.min
        for (m in (2:M)) {
          tau.c <- tau.est
          taum.min <- tau.est[m]
          for ( taum in seq((tau.est[m]-d.tau), (tau.est[m]+d.tau), by=delta) ) {
            if ( (taum-tau.est[m-1]>l.tau) & (tau.est[m+1]-taum>l.tau) ) {
              tau.c[m] <- taum
              argsM <- c(args1,tau.c[-c(1,M+1)], a0=0, a2=0)
              r.c <- mlxoptim2(model=model, output=out.name, data=d, initial=r.c$par, ilog=ilog, cw=cw, args=argsM)
              #              r.c$par['D0'] <- max(r.c$par.c['D0'] ,1)
              #              r.c$par['L0'] <- max(r.c$par.c['L0'] ,1)
              if (r.c$value < r.min) {
                r.min <- r.c$value 
                p.estM <- r.c$par
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
      
      argsM <- c(args1, tau.est[-c(1,M+1)], a0=0, a2=0)
      rA <- mlxoptim2(model=model, output=out.name, data=d, initial=p.estM, ilog=ilog, cw=cw, args=argsM)
      p.est.A <- rA$par
      #devA <- mlxoptim2(model=model, output=out.name, data=d, initial=p.est.A, ilog=ilog, cw=c(0,1,0,1), args=argsM, optim=F)$value
      devA <- rA$value
      
      p.ini.B <- c(p.est.A, p.est.1B[c('AW', 'AD', 'phi')])
      rB <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini.B, ilog=ilog, cw=cw, args=argsM)
      p.est.B <- rB$par
      #devB <- mlxoptim2(model=model, output=out.name, data=d, initial=p.est.B, ilog=ilog, cw=c(0,1,0,1), args=argsM, optim=F)$value
      devB <- rB$value
      
      DEVM <- c(devA, devB)
      pA <- c(p.est.A, AW=0, AD=0, phi=0, a0=0, a2=0, args1)
      pB <- c(p.est.B, a0=0, a2=0, args1) 
      PM <- rbind(pA, pB)
      tauA <- tauB <- c(tau.est, rep(0,M.max-M))
      TauM <- rbind(tauA, tauB)
      colnames(TauM) <- paste0("tau",0:(M.max))
      
      if (M==2) {
        list.arg <- list(c(a0=0))
        list.est <- list(c(a2=0))
      } else {
        list.arg <- list(c(a0=0), c(a2=0), NULL)
        list.est <- list(c(a2=0), c(a0=0), c(a0=0, a2=0))
      }
      
      for (k in seq_len(length(list.arg))) {
        argsM <- c(list.arg[[k]], tau.est[-c(1,M+1)], args1)
        p.ini.A <- c(p.est.A, list.est[[k]])
        p.ini.B <- c(p.est.B, list.est[[k]])
        rA <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini.A, ilog=ilog, cw=cw, args=argsM)
        p.est.Ak <- rA$par
        # devA <- mlxoptim2(model=model, output=out.name, data=d, initial=p.est.Ak, ilog=ilog, cw=c(0,1,0,1), args=argsM, optim=F)$value
        devA <- rA$value
        
        rB <- mlxoptim2(model=model, output=out.name, data=d, initial=p.ini.B, ilog=ilog, cw=cw, args=argsM)
        p.est.Bk <- rB$par
        # devB <- mlxoptim2(model=model, output=out.name, data=d, initial=p.est.Bk, ilog=ilog, cw=c(0,1,0,1), args=argsM, optim=F)$value
        devB <- rB$value
        
        DEVM <- c(DEVM, devA, devB)
        pA <- c(p.est.Ak, AW=0, AD=0, phi=0, list.arg[[k]], args1)
        pB <- c(p.est.Bk, a0=0, list.arg[[k]], args1) 
        iA=match(colnames(P),names(pA))
        iB=match(colnames(P),names(pB))
        PM <- rbind(PM, pA[iA], pB[iB])
        
        #      save(list = ls(all.names = TRUE), file = "temp.RData", envir = environment())
        
        tauA <- tauB <- c(tau.est, rep(0,M.max-M))
        TauM <- rbind(TauM, tauA, tauB)
      }
      DIMM <- apply(PM, MARGIN=1, function(x) {sum(x !=0)}) +  apply(TauM, MARGIN=1, function(x) {sum(x !=0)})
      BICM <- DEVM + coef.bic*DIMM + coef.bic*2*(PM[,'AW'] !=0)
      i.model <- which.min(BICM)
      p.est <- PM[i.model,]
      tau.est <- TauM[i.model,]
      res <- simulx(model=model,
                    parameter=c(args1, p.est, tau.est[2:M]),
                    output=out,
                    settings=list(load.design=TRUE))
      d[[1]]$pred <- res[[1]][,2]
      d[[2]]$pred <- res[[2]][,2]
      
      dd <- list()
      dd[[1]] <- data.frame(day=d[[1]]$day[-1], value=diff(d[[1]]$value), pred=diff(d[[1]]$pred))
      dd[[2]] <- data.frame(day=d[[2]]$day[-1], value=diff(d[[2]]$value), pred=diff(d[[2]]$pred))
      dd[[1]]$pred <- dd[[1]]$pred*(1 + p.est['AW']*cos(omega*dd[[1]]$day + p.est['phi']))
      dd[[2]]$pred <- dd[[2]]$pred*(1 + p.est['AD']*cos(omega*dd[[2]]$day + p.est['phi']))
      d[[1]]$pred <- cumsum(c(d[[1]]$pred[1], dd[[1]]$pred))
      d[[2]]$pred <- cumsum(c(d[[2]]$pred[1], dd[[2]]$pred))
      res[['ks']]$Reff <- res[['ks']]$ks/(p.est['kr']+p.est['kd'])
      
      pl1 <- ggplot(data=d[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl2 <- ggplot(data=d[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl3 <- ggplot(data=dd[[1]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl4 <- ggplot(data=dd[[2]]) + geom_point(aes(day,value),color="blue") + geom_line(aes(day,pred), color="red")
      pl5 <- ggplot(data=res[['ks']]) + geom_line(aes(time,Reff), color="red")
      grid.arrange(pl1, pl2, pl5, pl3, pl4, ncol=3)
      
      print(tau.est[2:M])
      print(rbind(DEVM, BICM))
      
      DIM <- c(DIM, DIMM)
      DEV <- c(DEV, DEVM)
      BIC <- c(BIC, BICM)
      P <- rbind(P, PM)
      Tau <- rbind(Tau, TauM)
    }
  }
  
  if (!is.null(M.est)) {
    M <- apply(Tau, MARGIN=1, function(x) {sum(x !=0)})
    jM <- which(M!=(M.est+1))
    BICM <- BIC
    BICM[jM] = Inf
    i.model <- which.min(BICM)
  } else {
    i.model <- which.min(BIC)
  }
  p.est <- P[i.model,]
  tau.est <- Tau[i.model,]
  M.est <- sum(tau.est>0) - 1
  eval(parse(text=paste0("model <- sir",M.est)))
  # res <- simulx(model=model,
  #               parameter=c(args1, p.est, tau.est[2:M.est]),
  #               output=out,
  #               settings=list(load.design=TRUE))
  
  delta.t <- 1/ndt
  tc <- seq(min(d1$day), max(d2$day) + nb.day, delta.t)
  outc <- list(name=c("W", "D", "ks"), time=tc)
  
  resc <- simulx(model=model,
                 parameter=c(args1, p.est, tau.est[2:M.est]),
                 output=outc)
  r.eff <- resc[["ks"]][,2]/(p.est['kr']+p.est['kd'])
  t.min <- min(d0$day) - 1
  dc1 <- data.frame(day=resc[["W"]][,1]+t.min, pred0=resc[["W"]][,2], variable="y", type="confirmed")
  dc2 <- data.frame(day=resc[["D"]][,1]+t.min, pred0=resc[["D"]][,2], variable="y", type="deaths")
  dW0 <- dc1$pred0[(ndt+1):nrow(dc1)] - dc1$pred0[1:(nrow(dc1)-ndt)]
  dD0 <- dc2$pred0[(ndt+1):nrow(dc2)] - dc2$pred0[1:(nrow(dc2)-ndt)]
  # dW <- dc1$pred[(ndt+1):nrow(dc1)] - dc1$pred[1:(nrow(dc1)-ndt)]
  # dD <- dc2$pred[(ndt+1):nrow(dc2)] - dc2$pred[1:(nrow(dc2)-ndt)]
  dc3 <- data.frame(day=dc1$day[(ndt+1):nrow(dc1)], pred0=dW0, variable="dy", type="confirmed")
  dc4 <- data.frame(day=dc2$day[(ndt+1):nrow(dc2)], pred0=dD0, variable="dy", type="deaths")
  dc3$pred <- dc3$pred0*(1 + p.est['AW']*cos(omega*(tc[1:nrow(dc3)]+1) + p.est['phi']))
  dc4$pred <- dc4$pred0*(1 + p.est['AD']*cos(omega*(tc[1:nrow(dc3)]+1) + p.est['phi']))
  dc1$pred <- dc2$pred <- 0
  for (j in (1:ndt)) {
    dc1$pred[seq(j, nrow(dc1), by=ndt)] <- cumsum(c(dc1$pred[j], dc3$pred[seq(j, nrow(dc3), by=ndt)]))
    dc2$pred[seq(j, nrow(dc2), by=ndt)] <- cumsum(c(dc2$pred[j], dc4$pred[seq(j, nrow(dc4), by=ndt)]))
  }
  
  dc <- rbind(dc1, dc2, dc3, dc4)
  dc$date <- as.Date(dc$day, "2020-01-21")
  
  res <- simulx(model=model,
                parameter=c(args1, p.est, tau.est[2:M.est]),
                output=out)
  d1$variable="y"
  d1$pred0 <- res[["W"]][,2]
  d2$variable="y"
  d2$pred0 <- res[["D"]][,2]
  dd1 <- data.frame(day=d1$day[2:nrow(d1)], value=diff(d1$value), pred0=diff(d1$pred0), variable="dy", type="confirmed")
  dd2 <- data.frame(day=d2$day[2:nrow(d2)], value=diff(d2$value), pred0=diff(d2$pred0), variable="dy", type="deaths")
  dd1$pred <- dd1$pred0*(1 + p.est['AW']*cos(omega*((1:nrow(dd1))+1)+ p.est['phi']))
  dd2$pred <- dd2$pred0*(1 + p.est['AD']*cos(omega*((1:nrow(dd2))+1) + p.est['phi']))
  d1$pred <- cumsum(c(d1$pred0[1], dd1$pred))
  d2$pred <- cumsum(c(d2$pred0[1], dd2$pred))
  d <- rbind(d1, d2)
  dd <- rbind(dd1, dd2)[, c(1, 2, 5, 4, 3, 6)]
  d <- rbind(d, dd)
  d$day <- d$day+ t.min
  d$date <- as.Date(d$day, "2020-01-21")
  d$variable <- as.factor(d$variable)
  d$e <- d$value-d$pred
  d <- d %>% group_by(variable, type) %>% mutate(sd=sd(e))
  
  R0 <- data.frame(day=dc1$day, date=as.Date(dc1$day, "2020-01-21"), r0=r.eff)
  
  Tau <- Tau + t.min 
  tau.est <- tau.est + t.min 
  print(BIC)
  save(Tau, P, DEV, BIC, DIM, cw, R0, d, dc, tau.est, file=file.out, version=2)
  
  return(list(file=file.out, p.est=p.est, tau.est=tau.est))
}
