corfit <- function(data=NULL, country=NULL, shift=NULL, log.scale=FALSE, variable=c("total", "daily"),
                   n.max=7, p=0, point.size=1.5, line.size=0.75, predict=TRUE, interval=TRUE, level=0.95, nc.min=100, nd.min=5) {
  
  if (is.null(country)) {
    d0 <- data
    country <- levels(d0$country)
  } else {
    country.in <- country
    d0 <- subset(data,  country==country.in)
  }
  # d0 <- filterCovid(file.in="covid19.csv", country=country, nc.min=100, nd.min=2)
  # d <- subset(read.csv(data), country==country.in)
  d0 <- subset(d0, value>nc.min | type=="deaths")
  d0 <- subset(d0, value>nd.min | type=="confirmed")
  d <- d0[,c('day', 'value','type')]
  d1 <- subset(d, type=="confirmed")
  d2 <- subset(d, type=="deaths")
  
  S1 <- p*max(d1$value)
  day01 <- min(which(d1$value>S1))
  dy01 <- ifelse(day01>1, d1$value[day01-1], d1$value[day01])
  S2 <- p*max(d2$value)
  day02 <- min(which(d2$value>S2))
  dy02 <- ifelse(day02>1, d2$value[day02-1], d2$value[day02])
  n <- table(d0$type)
  d0$dy <- c( c(dy01, diff(d0$value[1:n[1]])), c(dy02, diff(d0$value[(n[1]+1):(n[1]+n[2])]) ))
  
  
  n1 <- nrow(d1)
  L <- 6
  X <- cbind(1, d1$day, d1$day^2)
  dL <- seq(L,n1, by=L)
  for (dl in dL) {
    X <- cbind(X, c(rep(0, dl-1), (1:(n1-dl+1))^2))
  }
  d1$data <- d1$value
  d1$value <- X%*%solve(t(X)%*%X)%*%t(X)%*%d1$value
  
  d1 <- subset(d1, value>S1)
  d2 <- subset(d2, value>S2)
  omega <- 2*pi/7
  if (is.null(shift)) {
    d1k <- d1
    r1 <- r2 <- NULL
    dx <- min(c(nrow(d1), nrow(d2), 7))
    for (k in 0:dx) {
      d1k$day <- d1$day+k
      dk <- merge(d1k,d2, by=c("day"), all=TRUE)
      ddk <- drop_na(dk)
      lmk <- summary(lm(value.y ~ -1 + value.x + I(cos(omega*day)) + I(sin(omega*day)) , data=ddk))
      r1 <- c(r1, lmk$r.squared)
      # r1 <- c(r1, cor(dk[,'value.x'], dk[,'value.y'], use="pairwise.complete.obs"))
      # sk <- sum(ddk[,'value.x']*ddk[,'value.y'])/sum(ddk[,'value.x']^2)
      # r2 <- c(r2, 1-sum((ddk[,'value.x']*sk-ddk[,'value.y'])^2)/var(ddk[,'value.y'])/(nrow(ddk)-1))
    }
    r1 <- signif(r1,digit=3)
    # r2 <- signif(r2,digit=2)
    shift <- which.max(r1)
    #shift <- which.max(r1)
  }
  
  d1k <- d1
  d1k$day <- d1$day+shift
  d.tot <- merge(d1k,d2, by=c("day"), all=TRUE)
  d.tot <- d.tot[shift:nrow(d.tot),]
  d.tot$shift.day <- d.tot$day-shift
  lmxy <- lm(value.y ~ -1 + value.x + I(cos(omega*day)) + I(sin(omega*day)), data=d.tot)
  scale <- lmxy$coefficients[1]
  pred.tot <- predict(lmxy, newdata = d.tot, interval = "prediction")
  conf.tot <- predict(lmxy, newdata = d.tot, interval = "confidence")
  # # colnames(pred.tot) <- paste0(colnames(pred.tot),".tot")
  d.tot <- cbind(d.tot, pred.tot)
  d.tot$variable <- "Total numbers"
  # M <- model.matrix(lmxy)
  # V <- vcov(lmxy)
  # Md <- apply(M, MARGIN=2, diff)
  # 
  # sxy <- summary(lmxy)
  # se.beta <- sxy$coefficients[1,2]
  # sigma <- sxy$sigma
  # df <- sxy$df[2]
  
  d1k <-  data.frame(day=d1k$day[2:nrow(d1k)], value=diff(d1k$value), data=diff(d1k$data), type="confirmed")
  d2k <-  data.frame(day=d2$day[2:nrow(d2)], value=diff(d2$value), type="death")
  d.day <- merge(d1k,d2k, by=c("day"), all=TRUE)
  d.day <- d.day[shift:nrow(d.day),]
  d.day$shift.day <- d.day$day-shift
  # se.day <- se.beta*d.day$value.x + sqrt(2)*sigma
  # p.day <- d.day$value.x*scale
  # pred.day <- data.frame(fit=p.day, lwr=p.day+se.day*qt((1-level)/2,df), upr=p.day+se.day*qt((1+level)/2,df) )
  
  lmxy <- lm(value.y ~ -1 + value.x + I(cos(omega*day)) + I(sin(omega*day)), data=d.day)
  slmxy <- summary(lmxy)
  #  scale <- lmxy$coefficients
  pred.day <- predict(lmxy, newdata = d.day, interval = "prediction")
  conf.day <- predict(lmxy, newdata = d.day, interval = "confidence")
  
  d.day <- cbind(d.day, pred.day)
  d.day$variable <- "Daily numbers"
  
  # if (ra==1 | predict==FALSE) {
  Dmo <- NULL
  if ("total" %in% variable)
    Dmo <- rbind(Dmo, d.tot)
  if ("daily" %in% variable)
    Dmo <- rbind(Dmo, d.day)
  Dmo$variable = factor(Dmo$variable) 
  Dmo$variable = factor(Dmo$variable ,levels(Dmo$variable)[c(2, 1)])
  # } else {
  #   Dmo <- d.tot
  # }
  Dmo$lwr <- pmax(0, Dmo$lwr)
  Dmo[c('value.x', 'value.y', 'fit', 'lwr', 'upr')] <- (Dmo[c('value.x', 'value.y', 'fit', 'lwr', 'upr')])
  
  nl.max <- 6
  n1 <- min(d0$day)
  n2 <- max(d0$day)
  br <- round(c(seq(n1,n2, length.out = nl.max)))
  lb <- gsub("2020-0", "", as.Date(br, "2020-01-21"))
  
  if (predict==TRUE) {
    pl <- ggplot(Dmo) + geom_point(aes(day,value.y), color="#993399") +
      geom_line(aes(day, fit), color= "#339900") +
      facet_wrap(~variable, scales = "free_y", ncol=1) + theme(legend.position = c(0.1,0.85)) +  
      scale_x_continuous("date", breaks = br, label= lb) + ylab("number of deaths") 
    if (log.scale) 
      pl <- pl + scale_y_continuous(trans="log10")
    if (interval)
      pl <- pl + geom_ribbon(aes(x=day, ymin=lwr, ymax=upr), alpha=0.1, inherit.aes=F, fill="blue") 
    
  } else {
    colors <- c("deaths" = "#993399", "confirmed" = "#339900")
    pl <- ggplot(Dmo) + 
      geom_point(aes(day,value.y, color="deaths")) + geom_line(aes(day,value.y, color="deaths")) +
      geom_point(aes(day,data*scale, color="confirmed")) + geom_line(aes(day,data*scale, color="confirmed")) +
      scale_color_manual(values = colors) + labs(x="date", color = "Legend") +
      facet_wrap(~variable, scales = "free_y", ncol=1) + theme(legend.position = c(0.1,0.85)) +  
      scale_x_continuous( sec.axis = sec_axis(~ . - shift, breaks = br, label= lb), breaks = br, label= lb) 
    
    if (log.scale) 
      pl <- pl + scale_y_continuous("Number of deaths", sec.axis = sec_axis(~ . / scale, name="Number of confirmed"), trans="log10")
    else
      pl <- pl + scale_y_continuous("Number of deaths", sec.axis = sec_axis(~ . / scale, name="Number of confirmed"))
  }
  
  return(list(plot=pl, coeff=c(shift=shift, scale=signif(scale[[1]],3)), res=Dmo))
}
