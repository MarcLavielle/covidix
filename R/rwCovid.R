rwCovid <- function(file.out=NULL, data.correction=NULL) {
  
  url.confirmed <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
  
  url.deaths <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
  
  url.recovered <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
  
  type.in <- c("confirmed", "deaths")
  #  type.in <- c("confirmed", "deaths", "recovered")
  d <- NULL
  for (type in type.in) {
    
    # ---- read the data
    if (type=="confirmed")
      dk <- read.csv(url(url.confirmed), stringsAsFactors=T)
    else if (type=="deaths")
      dk <- read.csv(url(url.deaths), stringsAsFactors=T)
    else
      dk <- read.csv(url(url.recovered), stringsAsFactors=T)
    
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
    #dk$country <- as.factor(dk$country)
    
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
  
  
  if (!is.null(data.correction)) {
    data.correction$date <- as.Date(data.correction$date)
    for (ctr in levels(data.correction$country)) {
      dc <- subset(data.correction, country==ctr)
      for (tp in c("confirmed", "deaths")) {
        id <- which(d.corr$country==ctr & d.corr$type==tp)
        df <- d.corr[id,]
        ifc <- match(dc$date,df$date)
        dfv <- c(df$value[1], diff(df$value))
        dfv[ifc] <- dc[[tp]]
        d.corr[id,]$value <- cumsum(dfv)
      }
    }
    # for (j in (1:nrow(data.correction))) {
    #   cj <-data.correction[j,]
    #   ij <- which(d.corr$country==as.character(cj$country) & d.corr$day==cj$day & d.corr$type=="deaths")
    #   d.corr[ij,]$value <- d.corr[ij-1,]$value + cj$deaths
    #   ij <- which(d.corr$country==as.character(cj$country) & d.corr$day==cj$day & d.corr$type=="confirmed")
    #   d.corr[ij,]$value <- d.corr[ij-1,]$value + cj$confirmed
    # }
  }
  
  if (!is.null(file.out))
    write.csv(d.corr, file=file.out, quote=F, row.names = F)
  return(d.corr)
}
