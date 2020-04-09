# setwd("C:/Users/Marc/OneDrive/Xpop/covidix/shiny/R")
# load("../data4/Italy.RData")
bootix <- function(file=NULL, nb.day=14, level=0.90, n.mc=200, ndt=5) {
  
  load(file)
  out.name <- c("W", "D")
  ilog <- 1:9
  cw <- c(1, 1, 1, 1)
  delta.t <- 1/ndt
  
  M.est <- length(tau.est)-1
  
  eval(parse(text=paste0("model <- sir",M.est)))
  
  np0 <- ncol(P)-nrow(P)+1
  p.est <- P[M.est, seq(1, np0+M.est-1)]
  if (length(p.est)>np0)
    names(p.est)[(np0+1):length(p.est)] <- paste0("a",2:M.est)
  names(tau.est) <- paste0("tau",0:M.est)
  eval(parse(text=paste0("model <- sir",M.est)))
  
  
  d.min <- min(d$day) 
  d.max <- max(d$day) 
  d$day <- d$day - d.min + 1
  tc <- seq(min(d$day), max(d$day) + nb.day, delta.t)
  outc <- list(name=c("W", "D", "ks"), time=tc)
  out1 <- list(name=c("W"), time=subset(d, type=="confirmed" & variable=="y")$day )
  out2 <- list(name=c("D"), time=subset(d, type=="deaths" & variable=="y")$day )
  res <- simulx(model=model,
                parameter=c(p.est, tau.est[2:M.est]),
                output=list(out1, out2))
  
  
  d1 <- res[["W"]][,2]
  d2 <- res[["D"]][,2]
  n1 <- length(d1)
  n2 <- length(d2)
  sy1 <- subset(d, type=="confirmed" & variable=="y")[['sd']][1]
  sdy1 <- subset(d, type=="confirmed" & variable=="dy")[['sd']][1]
  sy2 <- subset(d, type=="deaths" & variable=="y")[['sd']][1]
  sdy2 <- subset(d, type=="deaths" & variable=="dy")[['sd']][1]
  
  
  DM <- vector(mode = "list", length = 13)
  ks <- P <- NULL
  for (i.mc in 1:n.mc) {
    print(i.mc)
    dy1 <- c(d1[1] + rnorm(1, 0, sy1), diff(d1) + rnorm(n1-1, 0, sdy1))
    dy2 <- c(d2[1] + rnorm(1, 0, sy2), diff(d2) + rnorm(n2-1, 0, sdy2))
    y1 <- cumsum(dy1)
    y2 <- cumsum(dy2)
    
    d0 <- list()
    d[which(d$type=="confirmed" & d$variable=="y"),'value'] <- y1
    d[which(d$type=="deaths" & d$variable=="y"),'value'] <- y2
    d[which(d$type=="confirmed" & d$variable=="dy"),'value'] <- dy1[2:n1]
    d[which(d$type=="deaths" & d$variable=="dy"),'value'] <- dy2[2:n2]
    d0[[1]] <- subset(d, type=="confirmed" & variable=="y")
    d0[[2]] <- subset(d, type=="deaths" & variable=="y")
    d0[[3]] <- subset(d, type=="confirmed" & variable=="dy")
    d0[[4]] <- subset(d, type=="deaths" & variable=="dy")
    
    r <- mlxoptim2(model=model, output=out.name, data=d0, initial=p.est, ilog=ilog, cw=cw, tau=tau.est[-c(1,M.est+1)])
    
    p.est <- r$par
    p.est0 <- p.est
    p.est0['AW'] <- p.est0['AD'] <- 0
    res0 <- simulx(model=model,
                   parameter=c(p.est0, tau.est[2:M.est]),
                   output=outc,
                   settings=list(load.design=TRUE))
    res <- simulx(model=model,
                  parameter=c(p.est, tau.est[2:M.est]),
                  output=outc)
    P <- rbind(P, p.est)
    
    W0 <- res0[["W"]][,2]
    W <- res[["W"]][,2]
    D0 <- res0[["D"]][,2]
    D <- res[["D"]][,2]
    DM[[1]] <- cbind(DM[[1]], W0)
    DM[[2]] <- cbind(DM[[2]], W)
    DM[[3]] <- cbind(DM[[3]], D0)
    DM[[4]] <- cbind(DM[[4]], D)
    DM[[5]] <- cbind(DM[[5]], W + rnorm(1, 0, sy1))
    DM[[6]] <- cbind(DM[[6]], D + rnorm(1, 0, sy2))
    
    dW0 <- W0[(ndt+1):length(W0)] - W0[1:(length(W0)-ndt)]
    dD0 <- D0[(ndt+1):length(D0)] - D0[1:(length(D0)-ndt)]
    dW <- W[(ndt+1):length(W)] - W[1:(length(W)-ndt)]
    dD <- D[(ndt+1):length(D)] - D[1:(length(D)-ndt)]
    DM[[7]] <- cbind(DM[[7]], dW0)
    DM[[8]] <- cbind(DM[[8]], dW)
    DM[[9]] <- cbind(DM[[9]], dD0)
    DM[[10]] <- cbind(DM[[10]], dD)
    DM[[11]] <- cbind(DM[[11]], dW + rnorm(1, 0, sdy1))
    DM[[12]] <- cbind(DM[[12]], dD + rnorm(1, 0, sdy2))
    
    r0=res[["ks"]][,2]/(p.est['kr']+p.est['kd'])
    DM[[13]] <- cbind(DM[[13]], r0)
  }
  Q <- list()
  nc <- length(tc)
  for (q in (1: length(DM))) {
     nq <- nrow(DM[[q]])
     qDM <- apply(DM[[q]], MARGIN=1, quantile, probs=c((1-level)/2, 0.5, (1+level)/2),na.rm=TRUE)
     Q[[q]] <- data.frame(day = tc[(nc-nq+1):nc] + d.min - 1, t(qDM))
  }
  for (q in (1: 12)) {
    Q[[q]]$type <- "deaths"
    Q[[q]]$variable <- "y"
  }
  for (q in c(1,2,5,7,8,11)) 
    Q[[q]]$type <- "confirmed"
  for (q in c(7:12)) 
    Q[[q]]$variable <- "dy"
  
  Q.pred0.conf <- rbind(Q[[1]], Q[[3]], Q[[7]], Q[[9]])
  Q.pred.conf <- rbind(Q[[2]], Q[[4]], Q[[8]], Q[[10]])
  Q.pred.pred <- rbind(Q[[5]], Q[[6]], Q[[11]], Q[[12]])
  Qu <- list(pred0.conf=Q.pred0.conf, pred.conf=Q.pred.conf, 
             pred.pred=Q.pred.pred, Reff=Q[[13]])
  
  load(file=file)
  save(Q, Qu, P, DM, Tau, P, R1, R2, BIC, cw, R0, d, dc, tau.est, file=file, version=2)
  return(list(quantiles=Qu, parameter=P))
}
