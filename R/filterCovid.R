filterCovid <- function(type=c("confirmed", "deaths"), country=NULL,
                        file.in="covid19.csv", file.out=NULL, data=NULL,
                        nc.min=0, nd.min=0, fc.min=0, fd.min=0,
                        day.max=Inf, day.min=0) {
  
  if (is.null(data))
    d <- read.csv(file.in)
  else
    d <- data
  if (is.null(country)) {
    country <- levels(d$country)
  } else if (length(country)==1 && grepl("\\*", country)) {
    country <- sub("*","", country)
    country <- unique (grep(paste(country,collapse="|"), levels(d$country), value=TRUE))
  }
  
  if (is.null(type))  type <- levels(d$type)
  
  type.in <- type
  country.in <- country
  d <- subset(d, type %in% type.in & country %in% country.in)
  
  d <- d %>%  filter(day>day.min & day<=day.max) #%>%
  # group_by(type) %>% 
  # filter(value >= ifelse(type=="confirmed",nc.min, nd.min)) %>%
  # ungroup()
  d <- droplevels.data.frame(d)
  
  d1 <- subset(d, type=="confirmed")
  d2 <- subset(d, type=="deaths")
  for (country.in in levels(d$country)) {
    d1c <- subset(d1, country==country.in  & value>=nc.min)$day
    d2c <- subset(d2, country==country.in &  value>=nd.min)$day
    d1min <-ifelse(length(d1c)>0, min(d1c), 0)
    d2min <-ifelse(length(d2c)>0, min(d2c), 0)
    dc.min <- max(d1min, d2min)
    if (sum(d$day<dc.min)>0) {
      if (dc.min>0)
        d <- d[-which(d$country==country.in & d$day<dc.min), ]
      else
        d <- d[-which(d$country==country.in), ]
    }
#    print(c(country.in, nrow(d)))
  }
  
  d <- d %>%  group_by(country, type) %>%
    filter(n() >= ifelse(type=="confirmed",fc.min,fd.min)) %>%
    ungroup()
  
  d <- droplevels.data.frame(d)
  
  iv <- 1
  while (length(iv)>0) {
    dv <- diff(d$value)
    nd <- nrow(d)
    
    rdv <- dv[2:(nd-2)] /pmax((dv[1:(nd-3)] + dv[3:(nd-1)])/2, 1)
    pv <- 0.05
    # iv <- which(dv[1:(nd-2)]>0 & dv[2:(nd-1)]==0 & dv[3:(nd)]>0 & 
    iv <- which( rdv < pv  & dv[3:(nd-1)]>dv[1:(nd-3)] & dv[3:(nd-1)]>dv[2:(nd-2)] &
                   d$country[1:(nd-3)]==d$country[2:(nd-2)] & 
                   d$country[2:(nd-2)]==d$country[3:(nd-1)] & 
                   d$type[1:(nd-3)]==d$type[2:(nd-2)] & 
                   d$type[2:(nd-2)]==d$type[3:(nd-1)] )
    
    d$value[iv+2] <-   (d$value[iv+2] + d$value[iv+3])/2
  }
  d['percentage'] <- d['value']
  i0 <- NULL
  for (country in levels(d[['country']])) {
    for (type in levels(d[['type']])) {
      ic <- which(d[['country']]==country & d[['type']]==type)
      if (length(ic)==0)
        i0 <- c(i0, country)
      else {
        max.v1 <- max(d[ic,'value'])
        d[ic,'percentage'] <- d[ic,'value']/max.v1*100
      }
    }
  }
  d <- subset(d, !(country %in% i0))
  d <- droplevels.data.frame(d)
  
  if (!is.null(file.out))
    write.csv(d, file=file.out, quote=F, row.names = F)
  
  return(d)
}
