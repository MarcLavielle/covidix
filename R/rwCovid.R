rwCovid <- function(file.out=NULL) {
  
  # https://github.com/CSSEGISandData/COVID-19/ 
  
  url.confirmed <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
  
  url.deaths <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
  
  # deprecated:
  # url.recovered <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv"
  
  # download.file("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv",
  #               destfile = "Confirmed.csv")
  # download.file("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv",
  #               destfile = "Deaths.csv")
  
  type.in <- c("confirmed", "deaths")
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
  idf <- which(d$country=="France" & d$type=="deaths" & d$day>=71 )
  nf <- diff(c(0,d[idf,]$value))
  nf.corr <- nf*0.7
  nf.corr.d <- c(nf[1], 471, 588, 441, 357, 595, 607, 541, 412)
  nf.corr[1:length(nf.corr.d)] <- nf.corr.d
  d.corr[idf,]$value <- cumsum(nf.corr)
  
  idf <- which(d$country=="France" & d$type=="confirmed" & d$day>=73)
  nf <- diff(c(0,d[idf,]$value))
  nf.corr <- nf*0.75
  nf.corr.d <- c(nf[1],  4060,  2886,  3116,  3777,  3881,  4286)
  nf.corr[1:length(nf.corr.d)] <- nf.corr.d
  d.corr[idf,]$value <- cumsum(nf.corr)
  
  #   idf <- which(d$country=="France" & d$type=="confirmed" & d$day>=74)
  # d.corr[idf,]$value <- d[idf,]$value - 21555
  # idf <- which(d$country=="France" & d$type=="confirmed" & d$day>=76)
  # nf.d <- c(76455,3777, 3881, 2231)
  # d.corr[idf,]$value <- cumsum(nf.d)

  if (!is.null(file.out))
    write.csv(d.corr, file=file.out, quote=F, row.names = F)
  return(d.corr)
}
