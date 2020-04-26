plotCovid <- function(data=d, type=c("confirmed", "deaths"), nc.min=200, nd.min=20,
                      p.min=0.005, point.size=1, line.size=0.75, log.scale=TRUE, plot=FALSE) {

  type.in <- type
  d <- data
  a.min <- aggregate(d$value, by=list(d$type, d$country), min)
  cmin <- aggregate(a.min$x, by=list(a.min$Group.1), max)[,2] + 1
  cmin <- pmax(cmin, c(100, 10))

  d.min1 <- d  %>%
    group_by(country, type) %>%
    slice(max(which(value <= ifelse(type=="confirmed",cmin[1],cmin[2])))) %>%
    mutate(value1=value, day1=day) %>% select(-c(day, value))
  d.min2 <- d  %>%
    group_by(country, type) %>%
    filter(value >= ifelse(type=="confirmed",nc.min,nd.min)) %>%
    slice(min(which(value >= ifelse(type=="confirmed",cmin[1],cmin[2])))) %>%
    mutate(value2=value, day2=day) %>% select(-c(day, value))
  d.min <- full_join(d.min1,d.min2, by=c("country","type")) %>% drop_na() %>%
    mutate(value=ifelse(type=="confirmed",cmin[1],cmin[2])) %>%
    mutate(day=ifelse(day1==day2,day1,day1+(day2-day1)*(value-value1)/(value2-value1)))

  d.min <- full_join(d,d.min, by=c("country","type")) %>%
    slice(which(value.x >= ifelse(type=="confirmed",cmin[1],cmin[2]))) %>%
    drop_na()

  d <- d %>% group_by(country, type) %>%
    mutate(prop=value/max(value)) %>%
    filter(prop>p.min)
  d$country <- droplevels(d$country)
  d$type <- droplevels(d$type)

  if ("confirmed" %in% levels(d$type))
    x.lab <- paste0("Number of days after the ", cmin[1], "th confirmed case  ")
  if ("death" %in% levels(d$type))
    x.lab <- c(x.lab, paste0(" ;  Number of days after the ", cmin[2], "th death "))

    if (length(type.in)>1) {
      dd <- merge(subset(d, type=="confirmed" & value>nc.min)[,c(1,3,4)],
                  subset(d, type=="deaths")[,c(1,3,4)],
                  by=c("country", "day"))
      names(dd)[3:4] <- c("confirmed", "death")
      dd$rate <- dd[["death"]]/dd[["confirmed"]]
      dd$country <- droplevels(dd$country)
  }

  j <- 1
  pl <- list()

  nl.max <- 15
  n1 <- min(d$day)
  n2 <- max(d$day)
  br <- round(c(seq(n1,n2, length.out = nl.max)))
  lb <- gsub("2020-0", "", as.Date(br, "2020-01-21"))

  if ("confirmed" %in% type) {
    pl[[j]] <- ggplot(data=subset(d, type=="confirmed"), aes(day, value)) +
      geom_point(size=point.size, color="red") + geom_line(size=line.size, color="blue") +
      scale_x_continuous(breaks = br, label= lb ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.7, size=7)) +
      facet_wrap(~ country, scales="free")  + ylab("Total number of confirmed cases")
  }
  if ("deaths" %in% type) {
    j <- j+1
    pl[[j]] <- ggplot(data=subset(d, type=="deaths"), aes(day, value)) +
      geom_point(size=point.size, color="red") + geom_line(size=line.size, color="blue") +
      scale_x_continuous(breaks = br, label= lb ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.7, size=7)) +
      facet_wrap(~ country, scales="free")  + ylab("Total number of deaths")
  }
  if ("recovered" %in% type) {
    j <- j+1
    pl[[j]] <- ggplot(data=subset(d, type=="deaths"), aes(day, value)) +
      geom_point(size=point.size, color="red") + geom_line(size=line.size, color="blue") +
      scale_x_continuous(breaks = br, label= lb ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.7, size=7)) +
      facet_wrap(~ country, scales="free")  + ylab("Total number of deaths")
  }
  
  j <- j+1
  pl[[j]] <- ggplot(data=d.min, aes(day.x-day.y, value.x, color=country, shape=country)) +
    geom_line(size=line.size)  +  geom_point(size=point.size) +
   scale_shape_manual(values=1:nlevels(d$country)) +
    facet_wrap(~ type, scales="free", ncol=2)  + ylab("#") +
    scale_x_continuous(name=x.lab)
  if (log.scale)
    pl[[j]] <-   pl[[j]] + scale_y_log10()

  j <- j+1
  if (length(type.in)>1) {
    pl[[j]] <- ggplot(data=dd, aes(confirmed, rate, color=country, shape=country)) +
      geom_line(size=line.size) + geom_point(size=point.size) + scale_shape_manual(values=1:nlevels(d$country)) +
      xlab("Total number of confirmed cases") + ylab("mortality rate")
    if (log.scale)
      pl[[j]] <-  pl[[j]] + scale_x_log10()
  }

  if (plot)
    for (j in 1:length(pl))
      print(pl[[j]])

  return(pl)
}
