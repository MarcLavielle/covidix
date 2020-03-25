filterCovid <- function(type=c("confirmed", "deaths"), country=NULL,
                        file.in="covid19.csv", file.out=NULL,
                        nc.min=0, nd.min=0, fc.min=0, fd.min=0,
                        day.max=Inf, day.min=0) {

  d <- read.csv(file.in)

  if (is.null(country)) {
    country <- levels(d$country)
  } else if (length(country)==1 && grep("*", country)) {
    country <- sub("*","", country)
    country <- unique (grep(paste(country,collapse="|"), levels(d$country), value=TRUE))
  }

  if (is.null(type))  type <- levels(d$type)

  type.in <- type
  country.in <- country
  d <- subset(d, type %in% type.in & country %in% country.in)

  d <- d %>%  filter(day>=day.min & day<=day.max) %>%
    group_by(type) %>%
    filter(value >= ifelse(type=="confirmed",nc.min,nd.min)) %>%
    ungroup() %>%
    group_by(country, type) %>%
    filter(n() >= ifelse(type=="confirmed",fc.min,fd.min)) %>%
    ungroup()

  d <- droplevels.data.frame(d)

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
