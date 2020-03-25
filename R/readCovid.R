writeCovid <- function(file.out="covid19.csv") {
  
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
  

   # dk$type <- type
    
    d <- rbind(d, dk)
  }
  
  write.csv(d, file=file.out, quote=F, row.names = F)
}
