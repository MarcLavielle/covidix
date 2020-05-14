fitPlot <- function(country="France", file=NULL, file.test=NULL, nb.day=2, plot.tau=F, plot.pred=T, plot.pred0=T,
                    plot.day=T, log.scale=F, level=0.9, pred.int=F, conf.int=F) {
  if (is.null(file))
    file <- gsub(" ","",paste0("data/",country,".RData"))
  
  
  # if (!is.null(file.train))
  #   load(file.train)
  # else   
    load(file)
  
  R0 <- subset(R0, day<=max(d$day)+nb.day)
  pl1 <- ggplot(R0)  + 
    xlab("date") + geom_hline(aes(yintercept=1), colour="blue", linetype="dashed") + 
    geom_line(aes(date, r0), color="red", size=1) +
 #   scale_y_continuous(name="Reff", breaks=c(1))
  scale_y_continuous(name="Reff")
  
  M.est <- length(tau.est)-1
  if (plot.day) {
    Dmo <- subset(d, variable=="dy")
    Dmc <- subset(dc, variable=="dy" & day<=max(d$day)+nb.day)
    pl2 <- ggplot(Dmo) + 
      xlab("") + ylab("") + expand_limits(y=0) +
      facet_wrap( ~type, scales = "free_y", ncol=1)
    
    if (pred.int) {
      sd <- Dmo %>% group_by(type) %>% summarize(sd=mean(sd))
      Dmc$lwr <- Dmc$pred + qnorm((1-level)/2)*subset(sd, type=="confirmed")$sd
      Dmc$upr <- Dmc$pred + qnorm((1+level)/2)*subset(sd, type=="confirmed")$sd
      id <- which(Dmc$type=="deaths")
      Dmc$lwr[id] <- Dmc$pred[id] + qnorm((1-level)/2)*subset(sd, type=="deaths")$sd
      Dmc$upr[id] <- Dmc$pred[id] + qnorm((1+level)/2)*subset(sd, type=="deaths")$sd
      # Dmc$lwr <- subset(Qu[['pred.pred']], variable=="dy" & day<=max(d$day)+nb.day)[,2]
      # Dmc$upr <- subset(Qu[['pred.pred']], variable=="dy" & day<=max(d$day)+nb.day)[,4]
      pl2 <- pl2 + geom_ribbon(data=Dmc, aes(x=date, ymin = lwr, ymax = upr), 
                               alpha=0.1, inherit.aes=F, fill="blue") 
    }
    # if (conf.int) {
    #   Dmc$lwr <- subset(Qu[['pred0.conf']], variable=="dy" & day<=max(d$day)+nb.day)[,2]
    #   Dmc$upr <- subset(Qu[['pred0.conf']], variable=="dy" & day<=max(d$day)+nb.day)[,4]
    #   pl2 <- pl2 + geom_ribbon(data=Dmc, aes(x=date, ymin = lwr, ymax = upr), 
    #                            alpha=0.4, inherit.aes=F, fill="#339900") 
    # }
    
    if (plot.pred0) {
      if (plot.pred)
        pl2 <- pl2 + geom_line(data=Dmc, aes(date,pred, color="#339900"), size=0.75)  
      else
        pl2 <- pl2 + geom_line(data=Dmc, aes(date,pred0, color="#339900"), size=0.75)  
    }
    if (plot.pred0)
      pl2 <- pl2 + geom_point( aes(date, value, color="#993399"), size=2) +
      scale_colour_manual(values = c( "#339900", "#993399"),
                          labels= c("prediction", "data"),
                          guide = guide_legend(override.aes = list(linetype = c( "solid", "blank"), shape = c(NA, 16)))) +
      theme(legend.title = element_blank(), legend.position=c(0.1, 0.9))
    else
      pl2 <- pl2 + geom_point( aes(date, value), color="#993399", size=2) 
    
  } else {
    Dmo <- subset(d, variable=="y")
    Dmc <- subset(dc, variable=="y" & day<=max(d$day)+nb.day)
    pl2 <- ggplot(Dmo) + 
      geom_line(data=Dmc, aes(date,pred), color="#339900", size=0.75) + 
      geom_point( aes(date, value), color="#993399", size=2) + 
      xlab("") + ylab("") + expand_limits(y=0) +
      facet_wrap( ~type, scales = "free_y")
    
  }
  if (plot.tau & M.est>1) {
    tau.est <- tau.est[2: (length(tau.est)-1)]
    Tau.est1 <- data.frame(day=tau.est, type="confirmed", variable="y")
    Tau.est1 <- rbind(Tau.est1, data.frame(day=tau.est, type="confirmed", variable="dy"))
    Tau.est2 <- data.frame(day=tau.est, type="deaths", variable="y")
    Tau.est2 <- rbind(Tau.est2, data.frame(day=tau.est, type="deaths", variable="dy"))
    Tau.est <- rbind(Tau.est1, Tau.est2)
    Tau.est$group <- as.factor(with(Tau.est, paste(type, variable)))
    levels(Tau.est$group) <- c("Daily number of confirmed cases", "Cumulated number of confirmed cases", 
                               "Daily number of deaths", "Cumulated number of deaths")
    Tau.est$group = factor(Tau.est$group,levels(Tau.est$group)[c(2, 1, 4, 3)])
    Tau.est$date <- as.Date(Tau.est$day, "2020-01-21")
    pl1 <- pl1 + geom_vline(data=Tau.est, aes(xintercept=date), color="red", linetype="dashed")
    pl2 <- pl2 + geom_vline(data=Tau.est, aes(xintercept=date), color="red", linetype="dashed")
  }
  if (log.scale)
    pl2 <- pl2 + scale_y_log10()
  
  if (!is.null(file.test)) {
    day.max <- max(d$day)
    load(file.test)
    if (plot.day) 
      pl2 <- pl2 + geom_point(data=subset(d, variable=='dy' & day>day.max), aes(date, value), color="red", size=3) 
    else
      pl2 <- pl2 + geom_point(data=subset(d, variable=='y' & day>day.max), aes(date, value), color="red", size=3) 
    
  }
  
  
  
  return(list(pl.R0=pl1, pl.fit=pl2))
}
