rwCovid <- function(file.out=NULL, file.correction=NULL) {
  
  #us=read.csv(url("https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"))
  
  # https://github.com/CSSEGISandData/COVID-19/ 
  
  url.confirmed <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
  
  url.deaths <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
  
  url.recovered <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
  # deprecated:
  # url.recovered <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv"
  
  # download.file("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv",
  #               destfile = "Confirmed.csv")
  # download.file("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv",
  #               destfile = "Deaths.csv")
  
  type.in <- c("confirmed", "deaths")
  #  type.in <- c("confirmed", "deaths", "recovered")
  #  type.in <- type
  d <- NULL
  
  for (type in type.in) {
    
    # ---- read the data
    if (type=="confirmed")
      dk <- read.csv(url(url.confirmed))
    else if (type=="deaths")
      dk <- read.csv(url(url.deaths))
    else
      dk <- read.csv(url(url.recovered))
    
    # ---- remove unncessary columns
    dk$Lat <- dk$Long <- NULL
    
    # ---- reorder the rows
    dk <- dk[order(dk[,2], dk[,1]),]
    
    # ---- redefine country/region as single name
    i1 <- which(as.character(dk[,1]) == as.character(dk[,2]))
    dk[i1,1] <- ""
    i2 <- which(!(dk[,1]==""))
    fi2 <- paste(dk[i2,2],dk[i2,1])
    i2 <- c(i2, grep("Korea", dk[,2]))
    fi2 <- c(fi2, "South Korea")
    levels(dk[,2]) <- c(levels(dk[,2]), fi2)
    dk[i2,2] <- fi2
    dk[,1]  <- NULL
    names(dk)[1] <- "country"
    
    # ---- remove some rows 
    i0 <- unique(c(grep("Princess",dk[,1]), grep(",",dk[,1]), grep("Virgin Islands",dk[,1])))
    dk <- dk[-i0,]
    
    # ---- reformat the data
    n.day <- ncol(dk)-1
    n.country <- nrow(dk)
    dk <- melt(dk, id=list("country"), variable.name="date")
    dk$country <- droplevels(dk$country)
    
    # ---- reformat the date and add the day
    nk <- sub("X","",dk[['date']])
    dk <- dk %>% 
      mutate(date=as.Date(nk,format = "%m.%d.%y")) %>%
      arrange(country,date) %>%
      mutate(day=rep(1:n.day, n.country)) %>% drop_na() %>%
      mutate(type=type) %>%
      select(country, date, day, value, type)
    
    d <- rbind(d, dk)
  }
  
  d.corr <- d
  
  idg <- which(d$country=="Germany" & d$day==81 & d$type=="deaths")
  d.corr[idg,]$value <- 2871
  
  
  if (!is.null(file.correction)) {
    for (j in (1:nrow(correction))) {
      cj <-correction[j,]
      ij <- which(d.corr$country==as.character(cj$country) & d.corr$day==cj$day & d.corr$type=="deaths")
      d.corr[ij,]$value <- d.corr[ij-1,]$value + cj$deaths
      ij <- which(d.corr$country==as.character(cj$country) & d.corr$day==cj$day & d.corr$type=="confirmed")
      d.corr[ij,]$value <- d.corr[ij-1,]$value + cj$confirmed
    }
  }
  
  # idf <- which(d$country=="France" & d$type=="deaths" & d$day>=71 )
  # nf <- diff(c(0,d[idf,]$value))
  # nf.corr <- nf*0.7
  # nf.corr.d <- c(nf[1], 471, 588, 441, 357, 595, 607, 541, 412, 554, 345, 310, 335, 
  #                541, 514, 417, 418, 364, 227, 444, 387, 336, 311, 305, 198)
  # nf.corr[1:length(nf.corr.d)] <- nf.corr.d
  # d.corr[idf,]$value <- cumsum(nf.corr)
  # 
  # idf <- which(d$country=="France" & d$type=="confirmed" & d$day>=73)
  # nf <- diff(c(0,d[idf,]$value))
  # nf.corr <- nf*0.75
  # nf.corr.d <- c(nf[1],  4060,  2886,  3116,  3777,  3881,  4286, 4332, 3114, 1613, 2673, 5497, 2633, 
  #                2641, 405, 2569, 785, 2051, 2667, 1827, 1653, 1773, 1537)
  # nf.corr[1:length(nf.corr.d)] <- nf.corr.d
  # d.corr[idf,]$value <- cumsum(nf.corr)
  # 
  # idf <- which(d$country=="France" & d$type=="deaths" & d$day>=71 )
  # d.new <- d.corr[idf[-1],]
  # d.new$deaths <- diff(d.corr[idf,]$value)
  # d.new$value <- d.new$type <- NULL
  # idf <- which(d$country=="France" & d$type=="confirmed" & d$day>=71 )
  # d.new$confirmed <- diff(d.corr[idf,]$value)
  # 
  # idf <- which(d$country=="Spain" & d$type=="deaths" & d$day>59 )
  # nf.corr <- diff(c(0,d[idf,]$value))
  # nf.corr.d <- c(nf.corr[1], 394, 462, 514, 738, 655, 769, 832, 838, 812, 849, 
  #                864, 950, 932, 809, 674, 637, 743, 757, 683, 605, 
  #                510, 619, 517, 567, 523, 551, 585, 565, 410, 399, 430, 435, 440, 367, 378 )
  # nf.corr[1:length(nf.corr.d)] <- nf.corr.d
  # d.corr[idf,]$value <- cumsum(nf.corr)
  # idf <- which(d$country=="Spain" & d$type=="confirmed" & d$day>59 )
  # nf.corr <- diff(c(0,d[idf,]$value))
  # nf.corr.d <- c(nf.corr[1], 3646, 4517, 6584, 7937, 8578, 7871, 8189, 6549, 6398, 9222,
  #                7719, 8102, 7472, 7026, 6023, 4273, 5478, 6180, 5756, 4576,
  #                4830, 4167, 3477, 3045, c(5092, 5183, 5252, 4499, 4218)-1680, 2585, 3352, 2460, 
  #                2881, 2796, 2944 )
  # nf.corr[1:length(nf.corr.d)] <- nf.corr.d
  # d.corr[idf,]$value <- cumsum(nf.corr)
  # 
  
  # idf <- which(d$country=="Spain" & d$type=="deaths" & d$day>=59 )
  # d.new2 <- d.corr[idf[-1],]
  # d.new2$deaths <- diff(d.corr[idf,]$value)
  # d.new2$value <- d.new2$type <- NULL
  # idf <- which(d$country=="Spain" & d$type=="confirmed" & d$day>=59 )
  # d.new2$confirmed <- diff(d.corr[idf,]$value)
  # 
  # d.new <- rbind(d.new, d.new2)
  
  # idf <- which(d$country=="France" & d$type=="recovered" & d$day>=43 & d$day <=62)
  # d.corr[idf,]$value <- exp(seq(log(12),log(2200), length.out=20))
  # 
  # idf <- which(d$country=="Spain" & d$type=="recovered" & d$day>85 )
  # nf.corr <- diff(c(0,d[idf,]$value))
  # nf.corr.d <- c(nf.corr[1], 3502, 3166, 2695, 3530 )
  # nf.corr[1:length(nf.corr.d)] <- nf.corr.d
  # d.corr[idf,]$value <- cumsum(nf.corr)
  
  if (!is.null(file.out))
    write.csv(d.corr, file=file.out, quote=F, row.names = F)
  return(d.corr)
}
