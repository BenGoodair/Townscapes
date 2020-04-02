###load the packages###
library(splines)
library(grid)
library(lubridate)
library(plyr)
library(rgdal)
library(maptools)
library(raster)
library(dplyr)
library(sp)
library(tmaptools)
library(base)
library(lattice)
library(spatial)
library(classInt)
library(spatialEco)
library(RColorBrewer)
library(ggplot2)
library(tibble)
library(broom)
library(margins)
library(Ecdat)
library(readr)
library(rgeos)
library(colorspace)
library(ggplot2)
library(latticeExtra)
library(rMaps)
library(leaflet)
library(spatstat)
library(tidyverse)
library(ggrepel)
library(coefplot)
library(gcookbook)
library(ggthemes)
library(vcd)
library(graphics)
library(corrplot)
library(corrgram)
library(foreign)
library(plm)
library(clubSandwich)
library(ggeffects)
library(ggiraph)
library(ggiraphExtra)
library(stargazer)
library(ggmap)
library(feisr)
library(pglm)
library(survival)
library(sampling)
library("psycho")
library(reshape2)
library(ggpubr)
library(bife)
library(alpaca)

### download the data###


buas <- readShapeSpatial("bua")
buanames <- buas[-c(1,4,5,6,7,8,9)]
names(buanames)[names(buanames)=="bua11cd"] <- "BUA11CD"
names(buas)[names(buas)=="bua11cd"] <- "BUA11CD"
perseq <- seq(-100, 100, by=20)
ukcrs <- CRS("+init=epsg:27700")
popbua <- read_csv("ewpop16.csv")
popbua <- merge(popbua, buas, by="bua11nm")
popbua <- popbua[-c(1,3,5,6,7,8,9, 10)]
scotpop <- read_csv("scotpop.csv")
popbua <- rbind(popbua, scotpop)
popbua$BUA11CD <- as.character(popbua$BUA11CD)
proj4string(buas) <- crs(ukcrs)
collar <- (brewer.pal( 7, "RdBu"))
citylayer <- list( "sp.polygons",citybua, fill= "gray77", col="transparent")
theme_set(theme_minimal())

###UK geography###

scotbua <- readShapeSpatial("localities")
names(scotbua)[names(scotbua)=="code"] <- "BUA11CD"
names(scotbua)[names(scotbua)=="Shape_Area"] <- "st_areasha"
names(scotbua)[names(scotbua)=="Shape_Leng"] <- "st_lengths"
names(scotbua)[names(scotbua)=="OBJECTID"] <- "objectid"
names(scotbua)[names(scotbua)=="name"] <- "bua11nm"
buauk <- buas[-c(4:7)]
proj4string(scotbua) <- crs(ukcrs)
uktowns <- rbind(scotbua, buauk)
scotbuaname <- scotbua[-c(1,4,5)]
scotbuaname <- as.data.frame(scotbuaname)
buanames <- as.data.frame(buanames)
buanames <- rbind(buanames, scotbuaname)
citybua <- merge(uktowns, popbua, by="BUA11CD")
citybua <- citybua[which(citybua$poplsoa>175000), ]

buasize <- buas[-c(1,3,4,5,6,7,9)]
buasize$size <- (buasize$st_areasha/1000)/1000
buasize <- buasize[-c(2)]
###distance code###
citycent <- SpatialPointsDataFrame(gCentroid(citybua, byid=TRUE), citybua@data, match.ID=FALSE)
ucityz <- as.ppp(citycent)
uktownz <- merge(uktowns, popbua, by="BUA11CD")
uktownz <- uktownz[which(uktownz$poplsoa>10000 & uktownz$poplsoa<175000), ]
uktowncent <- SpatialPointsDataFrame(gCentroid(uktownz, byid=TRUE), uktownz@data, match.ID=FALSE)
uktownz <- as.ppp(uktowncent)
uktowncent$distance <- nncross(uktownz, ucityz, what="dist", k=1)

###newdistancecode###


##index of decline##

scotpop01 <- read_csv("scot01pop2.csv")
scotpop11 <- read_csv("scot11pops.csv")
scotyoung01 <- read_csv("scot01young.csv")
scotyoung11 <- read_csv("scot11young.csv")
scotbus10 <- read_csv("scotbus10.csv")
scotbus18 <- read_csv("scotbus18.csv")
scotdegree01 <- read_csv("scot01degree.csv")
scotdegree11 <- read_csv("scot11degree.csv")
scotemp01 <- read_csv("scot01emp.csv")
scotemp11 <- read_csv("scot11emp.csv")

engpop01 <- read_csv("eng01pop.csv")
engpop11 <- read_csv("eng11pop.csv")
engyoung01 <- read_csv("engyoung01.csv")
engyoung11 <- read_csv("eng11young.csv")
engbus10 <- read_csv("engbus10.csv")
engbus18 <- read_csv("eng18bus.csv")
engdegree01 <- read_csv("eng01degree.csv")
engdegree11 <- read_csv("eng11degree.csv")
engemp01 <- read_csv("eng01emp.csv")
engemp11 <- read_csv("eng11emp.csv")

engpop01 <- aggregate(engpop01[-1], engpop01["BUA11CD"], sum)
engpop11 <- aggregate(engpop11[-1], engpop11["BUA11CD"], sum)
engbus10 <- aggregate(engbus10[-1], engbus10["BUA11CD"], sum)
engbus18 <- aggregate(engbus18[-1], engbus18["BUA11CD"], sum)
engyoung01 <- aggregate(engyoung01[-1], engyoung01["BUA11CD"], sum)
engyoung11 <- aggregate(engyoung11[-1], engyoung11["BUA11CD"], sum)
engdegree01 <- aggregate(engdegree01[-1], engdegree01["BUA11CD"], sum)
engdegree11 <- aggregate(engdegree11[-1], engdegree11["BUA11CD"], sum)
engemp01 <- aggregate(engemp01[-1], engemp01["BUA11CD"], sum)
engemp11 <- aggregate(engemp11[-1], engemp11["BUA11CD"], sum)

scotbus10 <- aggregate(scotbus10[-1], scotbus10["BUA11CD"], sum)
scotbus18 <- aggregate(scotbus18[-1], scotbus18["BUA11CD"], sum)

engyoung01$young01 <- (engyoung01$fivenine01/engyoung01$pop01)*100
engyoung11$young11 <- (engyoung11$fivenine11/engyoung11$pop11)*100

engyoung <- merge(engyoung01, engyoung11, by="BUA11CD")
engyoung <- engyoung[-c(2,3,5,6)]
scotyoung <- merge(scotyoung01, scotyoung11, by="bua11nm")
engpop <- merge(engpop01, engpop11, by="BUA11CD")
scotpop <- merge(scotpop01, scotpop11, by="bua11nm")
scotdegree <- merge(scotdegree01, scotdegree11, by="bua11nm")
engdegree01$degree01 <- (engdegree01$degree/engdegree01$sixsevenpop)*100
sixseven11 <- engemp11[-c(3)]
engdegree11 <- merge(engdegree11, sixseven11, by="BUA11CD")
engdegree11$degree11 <- (engdegree11$degree/engdegree11$sixseven11)*100
engdegree <- merge(engdegree01, engdegree11, by="BUA11CD")
engdegree <- engdegree[-c(2,3,5,6)]

engbus <- merge(engbus10, engbus18, by="BUA11CD")
scotbus <- merge(scotbus10, scotbus18, by="BUA11CD")

engemp01$emp01 <- (engemp01$employed/engemp01$sixsevenpop)*100
engemp11$emp11 <- (engemp11$employ/engemp11$sixseven11)*100
engemp <- merge(engemp01, engemp11, by="BUA11CD")
engemp <- engemp[-c(2,3,5,6)]

scotemp <- merge(scotemp01, scotemp11, by="bua11nm")


scotemp$empchange <- scotemp$emp11-scotemp$emp01
engemp$empchange <- engemp$emp11-engemp$emp01
engdegree$degreechange <- engdegree$degree11-engdegree$degree01
scotdegree$degreechange <- scotdegree$degree11-scotdegree$degree01
engpop$popchange <- engpop$pop11-engpop$pop01
scotpop$popchange <- scotpop$pop11-scotpop$pop01
engyoung$youngchange <- engyoung$young11-engyoung$young01
scotyoung$youngchange <- scotyoung$young11-scotyoung$young01
engbus$buschange <- engbus$bus18-engbus$bus10
scotbus$buschange <- scotbus$bus18-scotbus$bus10

scotemp <- scotemp[-c(2,3)]
engemp <- engemp[-c(2,3)]
engdegree <- engdegree[-c(2,3)]
scotdegree <- scotdegree[-c(2,3)]
engpop <- engpop[-c(2,3)]
engyoung <- engyoung[-c(2,3)]
scotyoung <- scotyoung[-c(2,3)]
engbus <- engbus[-c(2,3)]
scotpop <- scotpop[-c(2,3)]
scotbus <- scotbus[-c(2,3)]

scotemp <- merge(scotemp, scotbuaname, by="bua11nm")
scotdegree <- merge(scotdegree, scotbuaname, by="bua11nm")
scotpop <- merge(scotpop, scotbuaname, by="bua11nm")
scotyoung <- merge(scotyoung, scotbuaname, by="bua11nm")

scotemp <- scotemp[-c(1)]
scotdegree <- scotdegree[-c(1)]
scotpop <- scotpop[-c(1)]
scotyoung <- scotyoung[-c(1)]

emp <- rbind(scotemp, engemp)
young <- rbind(engyoung, scotyoung)
bus <- rbind(engbus,scotbus)
degree <- rbind(engdegree, scotdegree)
pop <- rbind(scotpop, engpop)

indexdecline <- merge(emp, young, by="BUA11CD")
indexdecline <- merge(indexdecline, bus, by="BUA11CD")
indexdecline <- merge(indexdecline, pop, by="BUA11CD")
indexdecline <- merge(indexdecline, degree, by="BUA11CD")
indexdecline <- merge(indexdecline, popbua, by="BUA11CD")
indexdecline <- indexdecline[which(indexdecline$poplsoa<175000&indexdecline$poplsoa>10000), ]

indexdecline<- indexdecline %>% 
  psycho::standardize() 

indexdecline$growth <- indexdecline$empchange + indexdecline$youngchange + indexdecline$buschange + indexdecline$popchange +indexdecline$degreechange
indexdecline <- merge(indexdecline, buanames, by="BUA11CD")
indexdecline <- indexdecline[-c(2,3,4,5,6,7,9)]



towns <- as.data.frame(uktowncent)
towns <- towns[-c(2,4,5,8,9)]
towns <- merge(towns, indexdecline, by="BUA11CD")
towns$distance <- towns$distance/1000

towns <- towns %>% mutate(close=ifelse(towns$distance<20, 1, 0))
towns <- towns %>% mutate(middle=ifelse(towns$distance>20 & towns$distance<30, 1, 0))
towns <- towns %>% mutate(far=ifelse(towns$distance>30, 1, 0))

towns$close <- factor(towns$close)
towns$middle <- factor(towns$middle)
towns$far <- factor(towns$far)

towns <- towns %>% mutate(distancegroup=ifelse(towns$close==1, "Neighbouring", ifelse(towns$middle==1, "Middle", ifelse(towns$far==1, "Isolated", 0))))

towns$distancegroup <- factor(towns$distancegroup, levels = c("Isolated", "Middle", "Neighbouring"))

towns <- towns[-c(6,7,8)]

towns <- towns %>% mutate(declining=ifelse(towns$growth< -1, 1, 0))
towns <- towns %>% mutate(stagnant=ifelse(towns$growth> -1 & towns$growth<1, 1, 0))
towns <- towns %>% mutate(increase=ifelse(towns$growth>1, 1, 0))

towns$declining <- factor(towns$declining)
towns$stagnant <- factor(towns$stagnant)
towns$increase <- factor(towns$increase)

towns <- towns %>% mutate(declinegroup=ifelse(towns$declining==1, "Declining", ifelse(towns$stagnant==1, "Stagnant", ifelse(towns$increase==1, "Improving", 0))))

towns$declinegroup <- factor(towns$declinegroup, levels = c("Declining", "Stagnant", "Improving"))
towns <- towns[-c(7,8, 9)]

country <- read_csv("country.csv")
country <- country[-c(1)]
towns <- merge(towns, country, by="BUA11CD", all.x=TRUE)
towns$country <- ifelse(is.na(towns$country), 
                        'England', towns$country)

townnames <- read_csv("townnames.csv")
towns <- merge(towns, townnames, by="BUA11CD")



region <- read_csv("region2.csv")
region <- region %>% distinct(BUA11CD, .keep_all = TRUE)
towns <- merge(towns, region, by="BUA11CD", all.x=TRUE)
towns$Region <- ifelse(is.na(towns$Region), 
                       'Scotland', towns$Region)

towns <- merge(towns, buasize, by="BUA11CD", all.x=TRUE)
###jobdensity###

scotjobs <- read_csv("jobsscotinter.csv")
scotlookup <- read_csv("interset.csv")
scotjobs <- merge(scotjobs, scotlookup, by="inter")
scotjobs <- scotjobs[-c(1)]
scotjobs <- aggregate(scotjobs[-2], scotjobs["BUA11CD"], sum)
scotsixpop <- read_csv("sixsixsixpopscot.csv")

scotjobs <- merge(scotjobs, scotsixpop, by="BUA11CD")
scotjobs$jobdensity <- scotjobs$jobs/scotjobs$sixsixsixpop
scotjobs <- scotjobs[-c(2,3)]

engjobs <- read_csv("jobsengwa.csv")
engjobs <- aggregate(engjobs[-1], engjobs["BUA11CD"], sum)
engsixpop <- read_csv("sixsixsixpopengwa.csv")
engsixpop <- merge(engsixpop, buanames, by="bua11nm")
engsixpop <- engsixpop[-c(1)]
engjobs <- merge(engjobs, engsixpop, by="BUA11CD")
engjobs$jobdensity <- engjobs$jobs/engjobs$sixsixsixpop
engjobs <- engjobs[-c(2,3)]

jobdense <- rbind(engjobs, scotjobs)
towns <- merge(towns, jobdense, by="BUA11CD")

towns <- towns %>% mutate(residential=ifelse(towns$jobdensity<0.5, 1, 0))
towns <- towns %>% mutate(versatile=ifelse(towns$jobdensity>0.5 & towns$jobdensity<0.7, 1, 0))
towns <- towns %>% mutate(working=ifelse(towns$jobdensity>0.7, 1, 0))

towns$residential <- factor(towns$residential)
towns$versatile <- factor(towns$versatile)
towns$working <- factor(towns$working)

towns <- towns %>% mutate(workinggroup=ifelse(towns$working==1, "Working", ifelse(towns$versatile==1, "Partially Residential", ifelse(towns$residential==1, "Residential", 0))))

towns$workinggroup <- factor(towns$workinggroup, levels = c("Working", "Partially Residential", "Residential"))



##counts##


###map tries###

northeast <- towns[which(towns$Region=='North East'), ]
northw <- towns[which(towns$Region=='North West'), ]
southe <- towns[which(towns$Region=='South East'), ]
southw <- towns[which(towns$Region=='South West'), ]
york <- towns[which(towns$Region=='Yorkshire and The Humber'), ]
scot <- towns[which(towns$Region=='Scotland'), ]
wal <- towns[which(towns$Region=='Wales'), ]
westm <- towns[which(towns$Region=='West Midlands'), ]
eastm <- towns[which(towns$Region=='East Midlands'), ]
easteng <- towns[which(towns$Region=='East of England'), ]


uktownsname <- merge(uktowns, northeast, by="BUA11CD")
uktownsnw <- merge(uktowns, northw, by="BUA11CD")
uktownssc <- merge(uktowns, scot, by="BUA11CD")
uktownswa <- merge(uktowns, wal, by="BUA11CD")
uktownssw <- merge(uktowns, southw, by="BUA11CD")
uktownsse <- merge(uktowns, southe, by="BUA11CD")
uktownsyo <- merge(uktowns, york, by="BUA11CD")
uktownsem <- merge(uktowns, eastm, by="BUA11CD")
uktownswm <- merge(uktowns, westm, by="BUA11CD")
uktownsee <- merge(uktowns, easteng, by="BUA11CD")

uktownssc <- uktownssc[-c(15)]
uktownsname <- sp.na.omit(uktownsname)
uktownsnw <- sp.na.omit(uktownsnw)
uktownssc <- sp.na.omit(uktownssc)
uktownswa <- sp.na.omit(uktownswa)
uktownsyo <- sp.na.omit(uktownsyo)
uktownssw <- sp.na.omit(uktownssw)
uktownsse <- sp.na.omit(uktownsse)
uktownsem <- sp.na.omit(uktownsem)
uktownswm <- sp.na.omit(uktownswm)
uktownsee <- sp.na.omit(uktownsee)

health_2018 <- read_delim("health_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                          col_types = cols('PointX Classification Code' = col_number()))
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Easting`), ]
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Northing`), ]
names(health_2018)[names(health_2018)=="PointX Classification Code"] <- "code"

schools_2018 <- read_delim("schools_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                           col_types = cols('PointX Classification Code' = col_number()))
schools_2018 <- schools_2018 [complete.cases(schools_2018$`Feature Easting`), ]
schools_2018 <- schools_2018 [complete.cases(schools_2018$`Feature Northing`), ]
names(schools_2018)[names(schools_2018)=="PointX Classification Code"] <- "code"

gov_2018 <- read_delim("gov_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"

station_2018 <- read_delim("station_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                           col_types = cols('PointX Classification Code' = col_number()))
station_2018 <- station_2018 [complete.cases(station_2018$`Feature Easting`), ]
station_2018 <- station_2018 [complete.cases(station_2018$`Feature Northing`), ]
names(station_2018)[names(station_2018)=="PointX Classification Code"] <- "code"


mental18 <- health_2018[which(health_2018$code==05280372), ]
mental18 <- subset(mental18, "Feature Easting" != "" | "Feature Northing" != "")
mental18$poi_ID <- 1:nrow(mental18)
coords <- cbind(Easting = as.numeric(as.character(mental18$`Feature Easting`)),
                Northing = as.numeric(as.character(mental18$`Feature Northing`)))
mental18sp <- SpatialPointsDataFrame(coords, data = data.frame(mental18$Name,
                                                               mental18$poi_ID), proj4string = CRS("+init=epsg:27700"))
mental18sp <- mental18sp[!is.na(over(mental18sp, geometry(uktownsname))), ]

station18 <- station_2018[which(station_2018$code==010570738), ]
station18 <- subset(station18, "Feature Easting" != "" | "Feature Northing" != "")
station18$poi_ID <- 1:nrow(station18)
coords <- cbind(Easting = as.numeric(as.character(station18$`Feature Easting`)),
                Northing = as.numeric(as.character(station18$`Feature Northing`)))
station18sp <- SpatialPointsDataFrame(coords, data = data.frame(station18$Name,
                                                                station18$poi_ID), proj4string = CRS("+init=epsg:27700"))
station18sp <- station18sp[!is.na(over(station18sp, geometry(uktownsname))), ]



jobcent18 <- gov_2018[which(gov_2018$code==06330418), ]
jobcent18 <- subset(jobcent18, "Feature Easting" != "" | "Feature Northing" != "")
jobcent18$poi_ID <- 1:nrow(jobcent18)
coords <- cbind(Easting = as.numeric(as.character(jobcent18$`Feature Easting`)),
                Northing = as.numeric(as.character(jobcent18$`Feature Northing`)))
jobcent18sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent18$Name,
                                                                jobcent18$poi_ID), proj4string = CRS("+init=epsg:27700"))
jobcent18sp <- jobcent18sp[!is.na(over(jobcent18sp, geometry(uktownsname))), ]


furthered18 <- schools_2018[which(schools_2018$code==05310376), ]
furthered18 <- subset(furthered18, "Feature Easting" != "" | "Feature Northing" != "")
furthered18$poi_ID <- 1:nrow(furthered18)
coords <- cbind(Easting = as.numeric(as.character(furthered18$`Feature Easting`)),
                Northing = as.numeric(as.character(furthered18$`Feature Northing`)))
furthered18sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered18$Name,
                                                                  furthered18$poi_ID), proj4string = CRS("+init=epsg:27700"))
furthered18sp <- furthered18sp[!is.na(over(furthered18sp, geometry(uktownsname))), ]


GPs18 <- health_2018[which(health_2018$code==05280354), ]
GPs18 <- subset(GPs18, "Feature Easting" != "" | "Feature Northing" != "")
GPs18$poi_ID <- 1:nrow(GPs18)
coords <- cbind(Easting = as.numeric(as.character(GPs18$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs18$`Feature Northing`)))
GPs18sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs18$Name,
                                                            GPs18$poi_ID), proj4string = CRS("+init=epsg:27700"))
GPs18sp <- GPs18sp[!is.na(over(GPs18sp, geometry(uktownsname))), ]



hospitals18 <- health_2018[which(health_2018$code==05280371), ]
hospitals18 <- subset(hospitals18, "Feature Easting" != "" | "Feature Northing" != "")
hospitals18$poi_ID <- 1:nrow(hospitals18)
coords <- cbind(Easting = as.numeric(as.character(hospitals18$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals18$`Feature Northing`)))
hospitals18sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals18$Name,
                                                                  hospitals18$poi_ID), proj4string = CRS("+init=epsg:27700"))
hospitals18sp <- hospitals18sp[!is.na(over(hospitals18sp, geometry(uktownsname))), ]


hmrc18 <- gov_2018[which(gov_2018$code==06330417), ]
hmrc18 <- subset(hmrc18, "Feature Easting" != "" | "Feature Northing" != "")
hmrc18$poi_ID <- 1:nrow(hmrc18)
coords <- cbind(Easting = as.numeric(as.character(hmrc18$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc18$`Feature Northing`)))
hmrc18sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc18$Name,
                                                             hmrc18$poi_ID), proj4string = CRS("+init=epsg:27700"))
hmrc18sp <- hmrc18sp[!is.na(over(hmrc18sp, geometry(uktownsname))), ]




school18 <- read_delim("schools_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
school18 <- school18[complete.cases(school18$`Feature Easting`), ]
school18 <- school18[complete.cases(school18$`Feature Northing`), ]
names(school18)[names(school18)=="PointX Classification Code"] <- "code"
school18 <- school18[which(school18$code==05310379|school18$code==05310375), ]
school18 <- subset(school18, "Feature Easting" != "" | "Feature Northing" != "")
school18$poi_ID <- 1:nrow(school18)
coords <- cbind(Easting = as.numeric(as.character(school18$"Feature Easting")),
                Northing = as.numeric(as.character(school18$"Feature Northing")))
school18sp <- SpatialPointsDataFrame(coords, data = data.frame(school18$Name,
                                                               school18$poi_ID), proj4string = CRS("+init=epsg:27700"))
school18sp <- school18sp[!is.na(over(school18sp, geometry(uktownsname))), ]


health18 <- read_delim("health_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
health18 <- health18[complete.cases(health18$`Feature Easting`), ]
health18 <- health18[complete.cases(health18$`Feature Northing`), ]
names(health18)[names(health18)=="PointX Classification Code"] <- "code"
health18 <- health18[which(health18$code>05280000 &health18$code<05290000), ]
health18 <- subset(health18, "Feature Easting" != "" | "Feature Northing" != "")
health18$poi_ID <- 1:nrow(health18)
coords <- cbind(Easting = as.numeric(as.character(health18$"Feature Easting")),
                Northing = as.numeric(as.character(health18$"Feature Northing")))
health18sp <- SpatialPointsDataFrame(coords, data = data.frame(health18$Name,
                                                               health18$poi_ID), proj4string = CRS("+init=epsg:27700"))
health18sp <- health18sp[!is.na(over(health18sp, geometry(uktownsname))), ]



hall18 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
hall18 <- hall18[complete.cases(hall18$`Feature Easting`), ]
hall18 <- hall18[complete.cases(hall18$`Feature Northing`), ]
names(hall18)[names(hall18)=="PointX Classification Code"] <- "code"
hall18 <- hall18[which(hall18$code==06340456), ]
hall18 <- subset(hall18, "Feature Easting" != "" | "Feature Northing" != "")
hall18$poi_ID <- 1:nrow(hall18)
coords <- cbind(Easting = as.numeric(as.character(hall18$"Feature Easting")),
                Northing = as.numeric(as.character(hall18$"Feature Northing")))
hall18sp <- SpatialPointsDataFrame(coords, data = data.frame(hall18$Name,
                                                             hall18$poi_ID), proj4string = CRS("+init=epsg:27700"))
hall18sp <- hall18sp[!is.na(over(hall18sp, geometry(uktownsname))), ]



post18 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
post18 <- post18[complete.cases(post18$`Feature Easting`), ]
post18 <- post18[complete.cases(post18$`Feature Northing`), ]
names(post18)[names(post18)=="PointX Classification Code"] <- "code"
post18 <- post18[which(post18$code==06340457), ]
post18 <- subset(post18, "Feature Easting" != "" | "Feature Northing" != "")
post18$poi_ID <- 1:nrow(post18)
coords <- cbind(Easting = as.numeric(as.character(post18$"Feature Easting")),
                Northing = as.numeric(as.character(post18$"Feature Northing")))
post18sp <- SpatialPointsDataFrame(coords, data = data.frame(post18$Name,
                                                             post18$poi_ID), proj4string = CRS("+init=epsg:27700"))
post18sp <- post18sp[!is.na(over(post18sp, geometry(uktownsname))), ]



bus18 <- read_delim("bust18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
bus18 <- bus18[complete.cases(bus18$`Feature Easting`), ]
bus18 <- bus18[complete.cases(bus18$`Feature Northing`), ]
names(bus18)[names(bus18)=="PointX Classification Code"] <- "code"
bus18 <- bus18[which(bus18$code==10590732), ]
bus18 <- subset(bus18, "Feature Easting" != "" | "Feature Northing" != "")
bus18$poi_ID <- 1:nrow(bus18)
coords <- cbind(Easting = as.numeric(as.character(bus18$"Feature Easting")),
                Northing = as.numeric(as.character(bus18$"Feature Northing")))
bus18sp <- SpatialPointsDataFrame(coords, data = data.frame(bus18$Name,
                                                            bus18$poi_ID), proj4string = CRS("+init=epsg:27700"))
bus18sp <- bus18sp[!is.na(over(bus18sp, geometry(uktownsname))), ]



wifi18 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
wifi18 <- wifi18[complete.cases(wifi18$`Feature Easting`), ]
wifi18 <- wifi18[complete.cases(wifi18$`Feature Northing`), ]
names(wifi18)[names(wifi18)=="PointX Classification Code"] <- "code"
wifi18 <- wifi18[which(wifi18$code==06340802), ]
wifi18 <- subset(wifi18, "Feature Easting" != "" | "Feature Northing" != "")
wifi18$poi_ID <- 1:nrow(wifi18)
coords <- cbind(Easting = as.numeric(as.character(wifi18$"Feature Easting")),
                Northing = as.numeric(as.character(wifi18$"Feature Northing")))
wifi18sp <- SpatialPointsDataFrame(coords, data = data.frame(wifi18$Name,
                                                             wifi18$poi_ID), proj4string = CRS("+init=epsg:27700"))
wifi18sp <- wifi18sp[!is.na(over(wifi18sp, geometry(uktownsname))), ]



nursery18 <- read_delim("educatjun18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
nursery18 <- nursery18[complete.cases(nursery18$`Feature Easting`), ]
nursery18 <- nursery18[complete.cases(nursery18$`Feature Northing`), ]
names(nursery18)[names(nursery18)=="PointX Classification Code"] <- "code"
nursery18 <- nursery18[which(nursery18$code==05320397), ]
nursery18 <- subset(nursery18, "Feature Easting" != "" | "Feature Northing" != "")
nursery18$poi_ID <- 1:nrow(nursery18)
coords <- cbind(Easting = as.numeric(as.character(nursery18$"Feature Easting")),
                Northing = as.numeric(as.character(nursery18$"Feature Northing")))
nursery18sp <- SpatialPointsDataFrame(coords, data = data.frame(nursery18$Name,
                                                                nursery18$poi_ID), proj4string = CRS("+init=epsg:27700"))


nursery18sp <- nursery18sp[!is.na(over(nursery18sp, geometry(uktownsname))), ]

nursery18sp <- spTransform(nursery18sp, CRS("+init=epsg:4326"))
health18sp <- spTransform(health18sp, CRS("+init=epsg:4326"))

bus18sp <- spTransform(bus18sp, CRS("+init=epsg:4326"))
mymap <- leaflet(data=bus18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5)
addMarkers(mymap, data=bus18sp)%>%
  addScaleBar(mymap, position="bottomright")

icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = "red")

mymap <- leaflet(data=nursery18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addMarkers(mymap, data=nursery18sp, label= as.character(nursery18sp$nursery18.Name), labelOptions = labelOptions(noHide = T))

mymap <- leaflet(data=health18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addMarkers(mymap, data=health18sp, label= as.character(health18sp$health18.Name), labelOptions = labelOptions(noHide = T))


mymap <- leaflet(data=health18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addMarkers(mymap, data=health18sp)


uktownsname <- spTransform(uktownsname, CRS("+init=epsg:4326"))
uktownsname <- SpatialPointsDataFrame(gCentroid(uktownsname, byid=TRUE), uktownsname@data, match.ID=FALSE)

getColorp <- function(uktownsname) {
  sapply(uktownsname$growth, function(growth) {
    if(growth > 1) {
      "lightblue"
    } else if(growth > -1) {
      "blue"
    } else {
      "red"
    } })
}

icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownsname))

mymap <- leaflet(data=uktownsname) %>%
  addTiles() %>%  
  setView(lng=-1.5774599, lat=55.1, zoom=9) 
mymap <- addAwesomeMarkers(mymap, icon = icons, data=uktownsname, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("lightblue", "blue", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Improvement")
mapshot(mymap, file="figure 5.png", zoom=20, vwidth=680, vheight=980)  

uktownsname <- spTransform(uktownsname, CRS("+init=epsg:4326"))
uktownsname <- SpatialPointsDataFrame(gCentroid(uktownsname, byid=TRUE), uktownsname@data, match.ID=FALSE)

getColorp <- function(uktownsname) {
  sapply(uktownsname$jobdensity, function(jobdensity) {
    if(jobdensity > 0.7) {
      "lightblue"
    } else if(jobdensity > 0.5) {
      "blue"
    } else {
      "red"
    } })
}

icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownsname))

mymap <- leaflet(data=uktownsname) %>%
  addTiles() %>%  
  setView(lng=-1.5774599, lat=55.1, zoom=9) 
mymap <- addAwesomeMarkers(mymap, icon = icons, data=uktownsname, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("lightblue", "blue", "red"),
    labels = c("Working", "Partially Residential", "Residential"), opacity = 1,
    title = "Town Type")
mapshot(mymap, file="figure 8.png", zoom=20, vwidth=680, vheight=980)  





GPs18 <- health_2018[which(health_2018$code==05280369), ]
GPs18 <- subset(GPs18, "Feature Easting" != "" | "Feature Northing" != "")
GPs18$poi_ID <- 1:nrow(GPs18)
coords <- cbind(Easting = as.numeric(as.character(GPs18$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs18$`Feature Northing`)))
GPs18sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs18$Name,
                                                            GPs18$poi_ID), proj4string = CRS("+init=epsg:27700"))
GPs18sp <- GPs18sp[!is.na(over(GPs18sp, geometry(uktownsnw))), ]




bus18 <- read_delim("bust18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
bus18 <- bus18[complete.cases(bus18$`Feature Easting`), ]
bus18 <- bus18[complete.cases(bus18$`Feature Northing`), ]
names(bus18)[names(bus18)=="PointX Classification Code"] <- "code"
bus18 <- bus18[which(bus18$code==10590732), ]
bus18 <- subset(bus18, "Feature Easting" != "" | "Feature Northing" != "")
bus18$poi_ID <- 1:nrow(bus18)
coords <- cbind(Easting = as.numeric(as.character(bus18$"Feature Easting")),
                Northing = as.numeric(as.character(bus18$"Feature Northing")))
bus18sp <- SpatialPointsDataFrame(coords, data = data.frame(bus18$Name,
                                                            bus18$poi_ID), proj4string = CRS("+init=epsg:27700"))
bus18sp <- bus18sp[!is.na(over(bus18sp, geometry(uktownsnw))), ]




uktownsnw <- SpatialPointsDataFrame(gCentroid(uktownsnw, byid=TRUE), uktownsnw@data, match.ID=FALSE)
uktownsnw <- spTransform(uktownsnw, CRS("+init=epsg:4326"))



getColorp <- function(uktownsnw) {
  sapply(uktownsnw$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })
}

icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownsnw))

mymap <- leaflet(data=uktownsnw) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownsnw, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")



GPs18sp <- spTransform(GPs18sp, CRS("+init=epsg:4326"))
mymap <- leaflet(data=GPs18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5)
addMarkers(mymap, data=GPs18sp)%>%
  addScaleBar(mymap, position="bottomright")


bus18sp <- spTransform(bus18sp, CRS("+init=epsg:4326"))
mymap <- leaflet(data=bus18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5)
addMarkers(mymap, data=bus18sp)%>%
  addScaleBar(mymap, position="bottomright")


uktownssc <- SpatialPointsDataFrame(gCentroid(uktownssc, byid=TRUE), uktownssc@data, match.ID=FALSE)
uktownssc <- spTransform(uktownssc, CRS("+init=epsg:4326"))




getColorp <- function(uktownssc) {
  sapply(uktownssc$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })
}

icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownssc))

mymap <- leaflet(data=uktownssc) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownssc, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")




gov_2018 <- read_delim("gov_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"

station_2018 <- read_delim("station_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                           col_types = cols('PointX Classification Code' = col_number()))
station_2018 <- station_2018 [complete.cases(station_2018$`Feature Easting`), ]
station_2018 <- station_2018 [complete.cases(station_2018$`Feature Northing`), ]
names(station_2018)[names(station_2018)=="PointX Classification Code"] <- "code"

health_2018 <- read_delim("health_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                          col_types = cols('PointX Classification Code' = col_number()))
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Easting`), ]
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Northing`), ]
names(health_2018)[names(health_2018)=="PointX Classification Code"] <- "code"

infra2018 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Easting`), ]
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Northing`), ]
names(infra2018)[names(infra2018)=="PointX Classification Code"] <- "code"



library18 <- infra2018[which(infra2018$code==06340458), ]
library18 <- subset(library18, "Feature Easting" != "" | "Feature Northing" != "")
library18$poi_ID <- 1:nrow(library18)
coords <- cbind(Easting = as.numeric(as.character(library18$`Feature Easting`)),
                Northing = as.numeric(as.character(library18$`Feature Northing`)))
library18sp <- SpatialPointsDataFrame(coords, data = data.frame(library18$Name,
                                                                library18$poi_ID), proj4string = CRS("+init=epsg:27700"))

library18sp <- library18sp[!is.na(over(library18sp, geometry(uktownssc))), ]



station18 <- station_2018[which(station_2018$code==010570738), ]
station18 <- subset(station18, "Feature Easting" != "" | "Feature Northing" != "")
station18$poi_ID <- 1:nrow(station18)
coords <- cbind(Easting = as.numeric(as.character(station18$`Feature Easting`)),
                Northing = as.numeric(as.character(station18$`Feature Northing`)))
station18sp <- SpatialPointsDataFrame(coords, data = data.frame(station18$Name,
                                                                station18$poi_ID), proj4string = CRS("+init=epsg:27700"))
station18sp <- station18sp[!is.na(over(station18sp, geometry(uktownssc))), ]


hospitals18 <- health_2018[which(health_2018$code==05280371), ]
hospitals18 <- subset(hospitals18, "Feature Easting" != "" | "Feature Northing" != "")
hospitals18$poi_ID <- 1:nrow(hospitals18)
coords <- cbind(Easting = as.numeric(as.character(hospitals18$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals18$`Feature Northing`)))
hospitals18sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals18$Name,
                                                                  hospitals18$poi_ID), proj4string = CRS("+init=epsg:27700"))
hospitals18sp <- hospitals18sp[!is.na(over(hospitals18sp, geometry(uktownssc))), ]


fire18 <- gov_2018[which(gov_2018$code==06330414), ]
fire18 <- subset(fire18, "Feature Easting" != "" | "Feature Northing" != "")
fire18$poi_ID <- 1:nrow(fire18)
coords <- cbind(Easting = as.numeric(as.character(fire18$`Feature Easting`)),
                Northing = as.numeric(as.character(fire18$`Feature Northing`)))
fire18sp <- SpatialPointsDataFrame(coords, data = data.frame(fire18$Name,
                                                                fire18$poi_ID), proj4string = CRS("+init=epsg:27700"))
fire18sp <- fire18sp[!is.na(over(fire18sp, geometry(uktownssc))), ]

police18 <- gov_2018[which(gov_2018$code==06330422), ]
police18 <- subset(police18, "Feature Easting" != "" | "Feature Northing" != "")
police18$poi_ID <- 1:nrow(police18)
coords <- cbind(Easting = as.numeric(as.character(police18$`Feature Easting`)),
                Northing = as.numeric(as.character(police18$`Feature Northing`)))
police18sp <- SpatialPointsDataFrame(coords, data = data.frame(police18$Name,
                                                                police18$poi_ID), proj4string = CRS("+init=epsg:27700"))
police18sp <- police18sp[!is.na(over(police18sp, geometry(uktownssc))), ]

GPs18 <- health_2018[which(health_2018$code==05280369), ]
GPs18 <- subset(GPs18, "Feature Easting" != "" | "Feature Northing" != "")
GPs18$poi_ID <- 1:nrow(GPs18)
coords <- cbind(Easting = as.numeric(as.character(GPs18$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs18$`Feature Northing`)))
GPs18sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs18$Name,
                                                            GPs18$poi_ID), proj4string = CRS("+init=epsg:27700"))
GPs18sp <- GPs18sp[!is.na(over(GPs18sp, geometry(uktownssc))), ]

police18sp <- spTransform(police18sp, CRS("+init=epsg:4326"))
fire18sp <- spTransform(fire18sp, CRS("+init=epsg:4326"))
library18sp <- spTransform(library18sp, CRS("+init=epsg:4326"))
hospitals18sp <- spTransform(hospitals18sp, CRS("+init=epsg:4326"))
station18sp <- spTransform(station18sp, CRS("+init=epsg:4326"))
GPs18sp <- spTransform(GPs18sp, CRS("+init=epsg:4326"))

icons <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'green')
iconz <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'pink')
icona <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'blue')
iconq <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'orange')
iconw <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'red')
icont <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'purple')


mymap <- leaflet(data=police18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, data=police18sp, icon = icons, label = "Police Station", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
  addAwesomeMarkers( data=hospitals18sp,icon = iconq, label = "Hospital", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
  addAwesomeMarkers( data=library18sp,icon = iconw, label = "Library", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
  addAwesomeMarkers( data=fire18sp,icon = icona, label = "Fire Station", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
  addAwesomeMarkers( data=station18sp,icon = iconz, label = "Train Station",labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
  addScaleBar( position="bottomright")
    # addLegend(
    # position = "bottomright",
    # colors = c("green", "orange", "red", "blue", "pink", "purple"),
    # labels = c("Police Stations", "Hospitals", "Libraries", "Fire Stations", "Train Stations", "Doctor's Surgeries"), opacity = 1,
    # title = "Public Service"  ) 
    # 





uktownswa <- spTransform(uktownswa, CRS("+init=epsg:4326"))
uktownswa <- SpatialPointsDataFrame(gCentroid(uktownswa, byid=TRUE), uktownswa@data, match.ID=FALSE)

getColorp <- function(uktownswa) {
  sapply(uktownswa$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })}


icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownswa))

mymap <- leaflet(data=uktownswa) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownswa, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")

uktownswa <- merge(uktownswa, wales,by="BUA11CD")
uktownswa <- uktownswa%>% mutate(hosp =ifelse(uktownswa$hospitals==0,1,0)) 
uktownswa <- uktownswa%>% mutate(lib =ifelse(uktownswa$librarys==0,1,0)) 
uktownswa <- uktownswa%>% mutate(fir =ifelse(uktownswa$fires==0,1,0)) 
uktownswa <- uktownswa%>% mutate(pol =ifelse(uktownswa$polices==0,1,0)) 
uktownswa <- uktownswa%>% mutate(stat =ifelse(uktownswa$stations==0,1,0)) 
uktownswa <- uktownswa%>% mutate(fur =ifelse(uktownswa$furthereds==0,1,0)) 

uktownswa$missing <- uktownswa$hosp+uktownswa$fir+uktownswa$fur+uktownswa$pol+uktownswa$lib+uktownswa$stat
uktownswa$missing <- as.numeric(uktownswa$missing)

getColorp <- function(uktownswa) {
  sapply(uktownswa$missing, function(missing) {
    if(missing==5 ) {
      "Brown"
    } else if(missing==4) {
      "purple"
    }else if(missing==3) {
      "pink"
    }else if(missing==2)  {
      "blue"
    }else if(missing==1) {
      "black"
    } else {
      "green"
    } })
}

icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownswa))

mymap <- leaflet(data=uktownswa) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownswa, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("#33FF33", "#CCFF99", "#FFFF99", "#FFFF00", "#FFCC00", "#FF9900", "#FF3300"),
    labels = c("0", "1", "2", "3","4" ,"5", "6"), opacity = 1,
    title = "Town Decline")


uktownswm <- spTransform(uktownswm, CRS("+init=epsg:4326"))
uktownswm <- SpatialPointsDataFrame(gCentroid(uktownswm, byid=TRUE), uktownswm@data, match.ID=FALSE)

getColorp <- function(uktownswm) {
  sapply(uktownswm$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })}


icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownswm))

mymap <- leaflet(data=uktownswm) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownswm, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")

uktownsem <- spTransform(uktownsem, CRS("+init=epsg:4326"))
uktownsem <- SpatialPointsDataFrame(gCentroid(uktownsem, byid=TRUE), uktownsem@data, match.ID=FALSE)

getColorp <- function(uktownsem) {
  sapply(uktownsem$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })}


icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownsem))

mymap <- leaflet(data=uktownsem) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownsem, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")

uktownsee <- spTransform(uktownsee, CRS("+init=epsg:4326"))
uktownsee <- SpatialPointsDataFrame(gCentroid(uktownsee, byid=TRUE), uktownsee@data, match.ID=FALSE)

getColorp <- function(uktownsee) {
  sapply(uktownsee$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })}


icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownsee))

mymap <- leaflet(data=uktownsee) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownsee, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")

uktownsyo <- spTransform(uktownsyo, CRS("+init=epsg:4326"))
uktownsyo <- SpatialPointsDataFrame(gCentroid(uktownsyo, byid=TRUE), uktownsyo@data, match.ID=FALSE)

getColorp <- function(uktownsyo) {
  sapply(uktownsyo$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })}


icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownsyo))

mymap <- leaflet(data=uktownsyo) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownsyo, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")

uktownsse <- spTransform(uktownsse, CRS("+init=epsg:4326"))
uktownsse <- SpatialPointsDataFrame(gCentroid(uktownsse, byid=TRUE), uktownsse@data, match.ID=FALSE)

getColorp <- function(uktownsse) {
  sapply(uktownsse$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })}


icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownsse))

mymap <- leaflet(data=uktownsse) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownsse, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")

uktownssw <- spTransform(uktownssw, CRS("+init=epsg:4326"))
uktownssw <- SpatialPointsDataFrame(gCentroid(uktownssw, byid=TRUE), uktownssw@data, match.ID=FALSE)

getColorp <- function(uktownssw) {
  sapply(uktownssw$growth, function(growth) {
    if(growth > 1) {
      "green"
    } else if(growth > -1) {
      "orange"
    } else {
      "red"
    } })}


icons <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(uktownssw))

mymap <- leaflet(data=uktownssw) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownssw, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")




gov_2018 <- read_delim("gov_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"

station_2018 <- read_delim("station_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                           col_types = cols('PointX Classification Code' = col_number()))
station_2018 <- station_2018 [complete.cases(station_2018$`Feature Easting`), ]
station_2018 <- station_2018 [complete.cases(station_2018$`Feature Northing`), ]
names(station_2018)[names(station_2018)=="PointX Classification Code"] <- "code"

health_2018 <- read_delim("health_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                          col_types = cols('PointX Classification Code' = col_number()))
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Easting`), ]
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Northing`), ]
names(health_2018)[names(health_2018)=="PointX Classification Code"] <- "code"

infra2018 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Easting`), ]
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Northing`), ]
names(infra2018)[names(infra2018)=="PointX Classification Code"] <- "code"



library18 <- infra2018[which(infra2018$code==06340458), ]
library18 <- subset(library18, "Feature Easting" != "" | "Feature Northing" != "")
library18$poi_ID <- 1:nrow(library18)
coords <- cbind(Easting = as.numeric(as.character(library18$`Feature Easting`)),
                Northing = as.numeric(as.character(library18$`Feature Northing`)))
library18sp <- SpatialPointsDataFrame(coords, data = data.frame(library18$Name,
                                                                library18$poi_ID), proj4string = CRS("+init=epsg:27700"))

library18sp <- library18sp[!is.na(over(library18sp, geometry(uktownswm))), ]



station18 <- station_2018[which(station_2018$code==010570738), ]
station18 <- subset(station18, "Feature Easting" != "" | "Feature Northing" != "")
station18$poi_ID <- 1:nrow(station18)
coords <- cbind(Easting = as.numeric(as.character(station18$`Feature Easting`)),
                Northing = as.numeric(as.character(station18$`Feature Northing`)))
station18sp <- SpatialPointsDataFrame(coords, data = data.frame(station18$Name,
                                                                station18$poi_ID), proj4string = CRS("+init=epsg:27700"))
station18sp <- station18sp[!is.na(over(station18sp, geometry(uktownswm))), ]


hospitals18 <- health_2018[which(health_2018$code==05280371), ]
hospitals18 <- subset(hospitals18, "Feature Easting" != "" | "Feature Northing" != "")
hospitals18$poi_ID <- 1:nrow(hospitals18)
coords <- cbind(Easting = as.numeric(as.character(hospitals18$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals18$`Feature Northing`)))
hospitals18sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals18$Name,
                                                                  hospitals18$poi_ID), proj4string = CRS("+init=epsg:27700"))
hospitals18sp <- hospitals18sp[!is.na(over(hospitals18sp, geometry(uktownswm))), ]


fire18 <- gov_2018[which(gov_2018$code==06330414), ]
fire18 <- subset(fire18, "Feature Easting" != "" | "Feature Northing" != "")
fire18$poi_ID <- 1:nrow(fire18)
coords <- cbind(Easting = as.numeric(as.character(fire18$`Feature Easting`)),
                Northing = as.numeric(as.character(fire18$`Feature Northing`)))
fire18sp <- SpatialPointsDataFrame(coords, data = data.frame(fire18$Name,
                                                             fire18$poi_ID), proj4string = CRS("+init=epsg:27700"))
fire18sp <- fire18sp[!is.na(over(fire18sp, geometry(uktownswm))), ]

police18 <- gov_2018[which(gov_2018$code==06330422), ]
police18 <- subset(police18, "Feature Easting" != "" | "Feature Northing" != "")
police18$poi_ID <- 1:nrow(police18)
coords <- cbind(Easting = as.numeric(as.character(police18$`Feature Easting`)),
                Northing = as.numeric(as.character(police18$`Feature Northing`)))
police18sp <- SpatialPointsDataFrame(coords, data = data.frame(police18$Name,
                                                               police18$poi_ID), proj4string = CRS("+init=epsg:27700"))
police18sp <- police18sp[!is.na(over(police18sp, geometry(uktownswm))), ]

GPs18 <- health_2018[which(health_2018$code==05280369), ]
GPs18 <- subset(GPs18, "Feature Easting" != "" | "Feature Northing" != "")
GPs18$poi_ID <- 1:nrow(GPs18)
coords <- cbind(Easting = as.numeric(as.character(GPs18$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs18$`Feature Northing`)))
GPs18sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs18$Name,
                                                            GPs18$poi_ID), proj4string = CRS("+init=epsg:27700"))
GPs18sp <- GPs18sp[!is.na(over(GPs18sp, geometry(uktownswm))), ]

police18sp <- spTransform(police18sp, CRS("+init=epsg:4326"))
fire18sp <- spTransform(fire18sp, CRS("+init=epsg:4326"))
library18sp <- spTransform(library18sp, CRS("+init=epsg:4326"))
hospitals18sp <- spTransform(hospitals18sp, CRS("+init=epsg:4326"))
station18sp <- spTransform(station18sp, CRS("+init=epsg:4326"))
GPs18sp <- spTransform(GPs18sp, CRS("+init=epsg:4326"))

icons <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'green')
iconz <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'pink')
icona <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'blue')
iconq <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'orange')
iconw <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'red')
icont <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'purple')

mymap <- leaflet(data=police18sp) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
#addAwesomeMarkers(mymap, data=police18sp, icon = icons) %>%
  addAwesomeMarkers(mymap, data=hospitals18sp,icon = iconq) %>%
  addAwesomeMarkers( data=library18sp,icon = iconw) %>%
  addAwesomeMarkers( data=fire18sp,icon = icona) %>%
  addAwesomeMarkers( data=station18sp,icon = iconz) %>%
  addAwesomeMarkers( data=GPs18sp,icon = icont) %>%
  addScaleBar( position="bottomright")%>%
  addLegend(
    position = "bottomright",
    colors = c( "orange", "red", "blue", "pink", "purple"),
    labels = c("Hospitals", "Libraries", "Fire Stations", "Train Stations", "Doctor's Surgeries"), opacity = 1,
    title = "Public Service"  ) 


###map of ublic index#####
  
  northeast <- data18[which(data18$Region=='North East'), ]
  northw <- data18[which(data18$Region=='North West'), ]
  southe <- data18[which(data18$Region=='South East'), ]
  southw <- data18[which(data18$Region=='South West'), ]
  york <- data18[which(data18$Region=='Yorkshire and The Humber'), ]
  scot <- data18[which(data18$Region=='Scotland'), ]
  wal <- data18[which(data18$Region=='Wales'), ]
  westm <- data18[which(data18$Region=='West Midlands'), ]
  eastm <- data18[which(data18$Region=='East Midlands'), ]
  easteng <- data18[which(data18$Region=='East of England'), ]
  # 
  # growthsouthwest <- ggplot(data=york, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  #   geom_point()+
  #   labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in South West Towns and Distance from City")+
  #   scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  #   theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  #   geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
  #   geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  #   geom_text(aes(x=0, y=100, label="British Town Average"))+
  #   stat_smooth(method = 'lm', se=FALSE,aes(group = 1))
  # growthsouthwest <- ggplot(data=southw, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  #   geom_point()+
  #   labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in South West Towns and Distance from City")+
  #   scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  #   theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  #   geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
  #   geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  #   geom_text(aes(x=0, y=100, label="British Town Average"))+
  #   stat_smooth(method = 'lm', se=FALSE,aes(group = 1))
  
  
  
  uktownsname <- merge(uktowns, northeast, by="BUA11CD")
  uktownsnw <- merge(uktowns, northw, by="BUA11CD")
  uktownssc <- merge(uktowns, scot, by="BUA11CD")
  uktownswa <- merge(uktowns, wal, by="BUA11CD")
  uktownssw <- merge(uktowns, southw, by="BUA11CD")
  uktownsse <- merge(uktowns, southe, by="BUA11CD")
  uktownsyo <- merge(uktowns, york, by="BUA11CD")
  uktownsem <- merge(uktowns, eastm, by="BUA11CD")
  uktownswm <- merge(uktowns, westm, by="BUA11CD")
  uktownsee <- merge(uktowns, easteng, by="BUA11CD")
  
  uktownsname <- sp.na.omit(uktownsname)
  uktownsnw <- sp.na.omit(uktownsnw)
  uktownssc <- sp.na.omit(uktownssc)
  uktownswa <- sp.na.omit(uktownswa)
  uktownsyo <- sp.na.omit(uktownsyo)
  uktownssw <- sp.na.omit(uktownssw)
  uktownsse <- sp.na.omit(uktownsse)
  uktownsem <- sp.na.omit(uktownsem)
  uktownswm <- sp.na.omit(uktownswm)
  uktownsee <- sp.na.omit(uktownsee)
  
  
  uktownsname <- spTransform(uktownsname, CRS("+init=epsg:4326"))
  uktownsname <- SpatialPointsDataFrame(gCentroid(uktownsname, byid=TRUE), uktownsname@data, match.ID=FALSE)
  
  getColorp <- function(uktownsname) {
    sapply(uktownsname$services, function(services) {
      if(services > 1) {
        "lightblue"
      } else if(services > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownsname))
  
  mymap <- leaflet(data=uktownsname) %>%
    addTiles() %>%  
    setView(lng=-1.5774599, lat=55.1, zoom=9) 
 mymap <-  addAwesomeMarkers(mymap, icon = icons, data=uktownsname, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Higher", "Middling", "Lower"), opacity = 1,
      title = "Overall Level of Public Service Provision")
  
mapshot(mymap, file="figure 6.png", zoom=20, vwidth=680, vheight=980)  
  
  getColorp <- function(uktownsname) {
    sapply(uktownsname$publicchange, function(publicchange) {
      if(publicchange > 1) {
        "lightblue"
      } else if(publicchange > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownsname))
  
  mymap <- leaflet(data=uktownsname) %>%
    addTiles() %>%  
    setView(lng=-1.5774599, lat=55.1, zoom=9) 
  mymap <- addAwesomeMarkers(mymap, icon = icons, data=uktownsname, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Increasing", "Similar", "Decreasing"), opacity = 1,
      title = "Change to Public Service Provision 2011-18")
  
  mapshot(mymap, file="figure 7.png", zoom=20, vwidth=680, vheight=980)  
  
  
  
  ###scotland####
  
  
  uktownssc <- spTransform(uktownssc, CRS("+init=epsg:4326"))
  uktownssc <- SpatialPointsDataFrame(gCentroid(uktownssc, byid=TRUE), uktownssc@data, match.ID=FALSE)
  
  getColorp <- function(uktownssc) {
    sapply(uktownssc$publicchange, function(publicchange) {
      if(publicchange > 1) {
        "lightblue"
      } else if(publicchange > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssc))
  
  mymap <- leaflet(data=uktownssc) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=56.3975, zoom=8) 
mymap <-   addAwesomeMarkers(mymap, icon = icons, data=uktownssc, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Higher", "Middling", "Lower"), opacity = 1,
      title = "Change to Public Service Provision 2011-18")
  mapshot(mymap, file="figure6scotversiontwo.png", zoom=20, vwidth=750, vheight=999)  
  
  getColorp <- function(uktownssc) {
    sapply(uktownssc$services, function(services) {
      if(services > 1) {
        "lightblue"
      } else if(services > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssc))
  
  mymap <- leaflet(data=uktownssc) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=56.3975, zoom=8) 
  mymap <-  addAwesomeMarkers(mymap, icon = icons, data=uktownssc, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Higher", "Middling", "Lower"), opacity = 1,
      title = "Overall Level of Public Service Provision")
  mapshot(mymap, file="figure5scotversiontwo.png", zoom=20, vwidth=750, vheight=999)  
  
  getColorp <- function(uktownssc) {
    sapply(uktownssc$growth, function(growth) {
      if(growth > 1) {
        "lightblue"
      } else if(growth > -1) {
        "blue"
      } else {
        "red"
      } })
  }
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssc))
  
  mymap <- leaflet(data=uktownssc) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=56.3975, zoom=8) 
  mymap <-  addAwesomeMarkers(mymap, icon = icons, data=uktownssc, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
      title = "Town Improvement")
  mapshot(mymap, file="figure4scotversiontwo.png", zoom=20, vwidth=750, vheight=999)  
  
  getColorp <- function(uktownssc) {
    sapply(uktownssc$jobdensity, function(jobdensity) {
      if(jobdensity > 0.7) {
        "lightblue"
      } else if(jobdensity > 0.5) {
        "blue"
      } else {
        "red"
      } })
  }
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssc))
  
  mymap <- leaflet(data=uktownssc) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=56.3975, zoom=8) 
  mymap <-  addAwesomeMarkers(mymap, icon = icons, data=uktownssc, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Working", "Partially Residential", "Residential"), opacity = 1,
      title = "Town Type")
  mapshot(mymap, file="figure7scotversiontwo.png", zoom=20, vwidth=750, vheight=999)  
  


  gov_2018 <- read_delim("gov_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                         col_types = cols('PointX Classification Code' = col_number()))
  gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
  gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
  names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"

  station_2018 <- read_delim("station_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                             col_types = cols('PointX Classification Code' = col_number()))
  station_2018 <- station_2018 [complete.cases(station_2018$`Feature Easting`), ]
  station_2018 <- station_2018 [complete.cases(station_2018$`Feature Northing`), ]
  names(station_2018)[names(station_2018)=="PointX Classification Code"] <- "code"

  health_2018 <- read_delim("health_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                            col_types = cols('PointX Classification Code' = col_number()))
  health_2018 <- health_2018 [complete.cases(health_2018$`Feature Easting`), ]
  health_2018 <- health_2018 [complete.cases(health_2018$`Feature Northing`), ]
  names(health_2018)[names(health_2018)=="PointX Classification Code"] <- "code"

  infra2018 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                          col_types = cols('PointX Classification Code' = col_number()))
  infra2018 <- infra2018 [complete.cases(infra2018$`Feature Easting`), ]
  infra2018 <- infra2018 [complete.cases(infra2018$`Feature Northing`), ]
  names(infra2018)[names(infra2018)=="PointX Classification Code"] <- "code"



  library18 <- infra2018[which(infra2018$code==06340458), ]
  library18 <- subset(library18, "Feature Easting" != "" | "Feature Northing" != "")
  library18$poi_ID <- 1:nrow(library18)
  coords <- cbind(Easting = as.numeric(as.character(library18$`Feature Easting`)),
                  Northing = as.numeric(as.character(library18$`Feature Northing`)))
  library18sp <- SpatialPointsDataFrame(coords, data = data.frame(library18$Name,
                                                                  library18$poi_ID), proj4string = CRS("+init=epsg:27700"))

  library18sp <- library18sp[!is.na(over(library18sp, geometry(uktownssc))), ]



  station18 <- station_2018[which(station_2018$code==010570738), ]
  station18 <- subset(station18, "Feature Easting" != "" | "Feature Northing" != "")
  station18$poi_ID <- 1:nrow(station18)
  coords <- cbind(Easting = as.numeric(as.character(station18$`Feature Easting`)),
                  Northing = as.numeric(as.character(station18$`Feature Northing`)))
  station18sp <- SpatialPointsDataFrame(coords, data = data.frame(station18$Name,
                                                                  station18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  station18sp <- station18sp[!is.na(over(station18sp, geometry(uktownssc))), ]


  hospitals18 <- health_2018[which(health_2018$code==05280371), ]
  hospitals18 <- subset(hospitals18, "Feature Easting" != "" | "Feature Northing" != "")
  hospitals18$poi_ID <- 1:nrow(hospitals18)
  coords <- cbind(Easting = as.numeric(as.character(hospitals18$`Feature Easting`)),
                  Northing = as.numeric(as.character(hospitals18$`Feature Northing`)))
  hospitals18sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals18$Name,
                                                                    hospitals18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  hospitals18sp <- hospitals18sp[!is.na(over(hospitals18sp, geometry(uktownssc))), ]


  fire18 <- gov_2018[which(gov_2018$code==06330414), ]
  fire18 <- subset(fire18, "Feature Easting" != "" | "Feature Northing" != "")
  fire18$poi_ID <- 1:nrow(fire18)
  coords <- cbind(Easting = as.numeric(as.character(fire18$`Feature Easting`)),
                  Northing = as.numeric(as.character(fire18$`Feature Northing`)))
  fire18sp <- SpatialPointsDataFrame(coords, data = data.frame(fire18$Name,
                                                               fire18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  fire18sp <- fire18sp[!is.na(over(fire18sp, geometry(uktownssc))), ]

  police18 <- gov_2018[which(gov_2018$code==06330422), ]
  police18 <- subset(police18, "Feature Easting" != "" | "Feature Northing" != "")
  police18$poi_ID <- 1:nrow(police18)
  coords <- cbind(Easting = as.numeric(as.character(police18$`Feature Easting`)),
                  Northing = as.numeric(as.character(police18$`Feature Northing`)))
  police18sp <- SpatialPointsDataFrame(coords, data = data.frame(police18$Name,
                                                                 police18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  police18sp <- police18sp[!is.na(over(police18sp, geometry(uktownssc))), ]

  GPs18 <- health_2018[which(health_2018$code==05280369), ]
  GPs18 <- subset(GPs18, "Feature Easting" != "" | "Feature Northing" != "")
  GPs18$poi_ID <- 1:nrow(GPs18)
  coords <- cbind(Easting = as.numeric(as.character(GPs18$`Feature Easting`)),
                  Northing = as.numeric(as.character(GPs18$`Feature Northing`)))
  GPs18sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs18$Name,
                                                              GPs18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  GPs18sp <- GPs18sp[!is.na(over(GPs18sp, geometry(uktownssc))), ]

  police18sp <- spTransform(police18sp, CRS("+init=epsg:4326"))
  fire18sp <- spTransform(fire18sp, CRS("+init=epsg:4326"))
  library18sp <- spTransform(library18sp, CRS("+init=epsg:4326"))
  hospitals18sp <- spTransform(hospitals18sp, CRS("+init=epsg:4326"))
  station18sp <- spTransform(station18sp, CRS("+init=epsg:4326"))
  GPs18sp <- spTransform(GPs18sp, CRS("+init=epsg:4326"))
  # 
  # icons <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'green')
  # iconz <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'pink')
  # icona <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'blue')
  # iconq <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'orange')
  # iconw <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'red')
  # icont <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'purple')
   icons <- makeIcon("train.png")
   iconz <- makeIcon("library.png")
   icona <- makeIcon("fire.png")
   iconq <- makeIcon("police.png")
   iconw <- makeIcon("hospital.png")
   icont <- makeIcon("train.png", "train.png", 18,18)
   
  

  mymap <- leaflet(data=police18sp) %>%
    addTiles() %>%
    setView(lng=-1.3823, lat=53.0975, zoom=5)
  addAwesomeMarkers(mymap, data=police18sp, icon = iconq) %>%
    addMarkers( data=hospitals18sp,icon = iconw) %>%
    addMarkers( data=library18sp,icon = iconz) %>%
    addMarkers( data=fire18sp,icon = icona) %>%
    addMarkers( data=station18sp,icon = icont) %>%
    addScaleBar( position="bottomright")
  #   addLegend(
  #     position = "bottomright",
  #     colors = c("green", "orange", "red", "blue", "pink", "purple"),
  #     labels = c("Police Stations", "Hospitals", "Libraries", "Fire Stations", "Train Stations", "Doctor's Surgeries"), opacity = 1,
  #     title = "Public Service"  ) 
  # 
  # 
  
  ###wales####
  
  
  uktownswa <- spTransform(uktownswa, CRS("+init=epsg:4326"))
  uktownswa <- SpatialPointsDataFrame(gCentroid(uktownswa, byid=TRUE), uktownswa@data, match.ID=FALSE)
  
  getColorp <- function(uktownswa) {
    sapply(uktownswa$publicchange, function(publicchange) {
      if(publicchange > 1) {
        "lightblue"
      } else if(publicchange > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownswa))
  
  mymap <- leaflet(data=uktownswa) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=52.2675, zoom=8) 
  mymap <-   addAwesomeMarkers(mymap, icon = icons, data=uktownswa, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Higher", "Middling", "Lower"), opacity = 1,
      title = "Change to Public Service Provision 2011-18")
  mapshot(mymap, file="figure6wales.png", zoom=20, vwidth=550, vheight=750)  
  
  getColorp <- function(uktownswa) {
    sapply(uktownswa$services, function(services) {
      if(services > 1) {
        "lightblue"
      } else if(services > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownswa))
  
  mymap <- leaflet(data=uktownswa) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=52.2675, zoom=8) 
  mymap <-  addAwesomeMarkers(mymap, icon = icons, data=uktownswa, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Higher", "Middling", "Lower"), opacity = 1,
      title = "Overall Level of Public Service Provision")
  mapshot(mymap, file="figure7wales.png", zoom=20, vwidth=550, vheight=750)  
  
  getColorp <- function(uktownswa) {
    sapply(uktownswa$growth, function(growth) {
      if(growth > 1) {
        "lightblue"
      } else if(growth > -1) {
        "blue"
      } else {
        "red"
      } })
  }
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownswa))
  
  mymap <- leaflet(data=uktownswa) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=52.2675, zoom=8) 
  mymap <-  addAwesomeMarkers(mymap, icon = icons, data=uktownswa, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
      title = "Town Improvement")
  mapshot(mymap, file="figure8wales.png", zoom=20, vwidth=550, vheight=750)  
  
  getColorp <- function(uktownswa) {
    sapply(uktownswa$jobdensity, function(jobdensity) {
      if(jobdensity > 0.7) {
        "lightblue"
      } else if(jobdensity > 0.5) {
        "blue"
      } else {
        "red"
      } })
  }
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownswa))
  
  mymap <- leaflet(data=uktownswa) %>%
    addTiles() %>%  
    setView(lng=-3.6923, lat=52.2675, zoom=8) 
  mymap <-  addAwesomeMarkers(mymap, icon = icons, data=uktownswa, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Working", "Partially Residential", "Residential"), opacity = 1,
      title = "Town Type")
  mapshot(mymap, file="figure9wales.png", zoom=20, vwidth=550, vheight=750)  
  
  
  
  
  ####South West####
  
  
  uktownssw <- spTransform(uktownssw, CRS("+init=epsg:4326"))
  uktownssw <- SpatialPointsDataFrame(gCentroid(uktownssw, byid=TRUE), uktownssw@data, match.ID=FALSE)
  
  getColorp <- function(uktownssw) {
    sapply(uktownssw$services, function(services) {
      if(services > 1) {
        "lightblue"
      } else if(services > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssw))
  
  mymap <- leaflet(data=uktownssw) %>%
    addTiles() %>%  
    setView(lng=-1.3823, lat=53.0975, zoom=5) 
  addAwesomeMarkers(mymap, icon = icons, data=uktownssw, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Higher", "Middling", "Lower"), opacity = 1,
      title = "Overall Level of Public Service Provision")
  
  getColorp <- function(uktownssw) {
    sapply(uktownssw$growth, function(growth) {
      if(growth > 1) {
        "lightblue"
      } else if(growth > -1) {
        "blue"
      } else {
        "red"
      } })
  }
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssw))
  
  mymap <- leaflet(data=uktownssw) %>%
    addTiles() %>%  
    setView(lng=-1.3823, lat=53.0975, zoom=5) 
  addAwesomeMarkers(mymap, icon = icons, data=uktownssw, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
      title = "Town Improvement")
  
  getColorp <- function(uktownssw) {
    sapply(uktownssw$jobdensity, function(jobdensity) {
      if(jobdensity > 0.7) {
        "lightblue"
      } else if(jobdensity > 0.5) {
        "blue"
      } else {
        "red"
      } })
  }
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssw))
  
  mymap <- leaflet(data=uktownssw) %>%
    addTiles() %>%  
    setView(lng=-1.3823, lat=53.0975, zoom=5) 
  addAwesomeMarkers(mymap, icon = icons, data=uktownssw, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Working", "Partially Residential", "Residential"), opacity = 1,
      title = "Town Type")
  
  getColorp <- function(uktownssw) {
    sapply(uktownssw$publicchange, function(publicchange) {
      if(publicchange > 1) {
        "lightblue"
      } else if(publicchange > -1) {
        "blue"
      } else {
        "red"
      } })}
  
  
  icons <- awesomeIcons(
    icon = 'ios-close',
    iconColor = 'black',
    library = 'ion',
    markerColor = getColorp(uktownssw))
  
  mymap <- leaflet(data=uktownssw) %>%
    addTiles() %>%  
    setView(lng=-1.3823, lat=53.0975, zoom=5) 
  addAwesomeMarkers(mymap, icon = icons, data=uktownssw, labelOptions = labelOptions(noHide = T))%>%
    addLegend(
      position = "bottomright",
      colors = c("lightblue", "blue", "red"),
      labels = c("Higher", "Middling", "Lower"), opacity = 1,
      title = "Change to Public Service Provision 2011-18")
  
  
  gov_2018 <- read_delim("gov_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                         col_types = cols('PointX Classification Code' = col_number()))
  gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
  gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
  names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"
  
  station_2018 <- read_delim("station_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                             col_types = cols('PointX Classification Code' = col_number()))
  station_2018 <- station_2018 [complete.cases(station_2018$`Feature Easting`), ]
  station_2018 <- station_2018 [complete.cases(station_2018$`Feature Northing`), ]
  names(station_2018)[names(station_2018)=="PointX Classification Code"] <- "code"
  
  health_2018 <- read_delim("health_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                            col_types = cols('PointX Classification Code' = col_number()))
  health_2018 <- health_2018 [complete.cases(health_2018$`Feature Easting`), ]
  health_2018 <- health_2018 [complete.cases(health_2018$`Feature Northing`), ]
  names(health_2018)[names(health_2018)=="PointX Classification Code"] <- "code"
  
  infra2018 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,
                          col_types = cols('PointX Classification Code' = col_number()))
  infra2018 <- infra2018 [complete.cases(infra2018$`Feature Easting`), ]
  infra2018 <- infra2018 [complete.cases(infra2018$`Feature Northing`), ]
  names(infra2018)[names(infra2018)=="PointX Classification Code"] <- "code"
  
  
  
  library18 <- infra2018[which(infra2018$code==06340458), ]
  library18 <- subset(library18, "Feature Easting" != "" | "Feature Northing" != "")
  library18$poi_ID <- 1:nrow(library18)
  coords <- cbind(Easting = as.numeric(as.character(library18$`Feature Easting`)),
                  Northing = as.numeric(as.character(library18$`Feature Northing`)))
  library18sp <- SpatialPointsDataFrame(coords, data = data.frame(library18$Name,
                                                                  library18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  
  library18sp <- library18sp[!is.na(over(library18sp, geometry(uktownssw))), ]
  
  
  
  station18 <- station_2018[which(station_2018$code==010570738), ]
  station18 <- subset(station18, "Feature Easting" != "" | "Feature Northing" != "")
  station18$poi_ID <- 1:nrow(station18)
  coords <- cbind(Easting = as.numeric(as.character(station18$`Feature Easting`)),
                  Northing = as.numeric(as.character(station18$`Feature Northing`)))
  station18sp <- SpatialPointsDataFrame(coords, data = data.frame(station18$Name,
                                                                  station18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  station18sp <- station18sp[!is.na(over(station18sp, geometry(uktownssw))), ]
  
  
  hospitals18 <- health_2018[which(health_2018$code==05280371), ]
  hospitals18 <- subset(hospitals18, "Feature Easting" != "" | "Feature Northing" != "")
  hospitals18$poi_ID <- 1:nrow(hospitals18)
  coords <- cbind(Easting = as.numeric(as.character(hospitals18$`Feature Easting`)),
                  Northing = as.numeric(as.character(hospitals18$`Feature Northing`)))
  hospitals18sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals18$Name,
                                                                    hospitals18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  hospitals18sp <- hospitals18sp[!is.na(over(hospitals18sp, geometry(uktownssw))), ]
  
  
  fire18 <- gov_2018[which(gov_2018$code==06330414), ]
  fire18 <- subset(fire18, "Feature Easting" != "" | "Feature Northing" != "")
  fire18$poi_ID <- 1:nrow(fire18)
  coords <- cbind(Easting = as.numeric(as.character(fire18$`Feature Easting`)),
                  Northing = as.numeric(as.character(fire18$`Feature Northing`)))
  fire18sp <- SpatialPointsDataFrame(coords, data = data.frame(fire18$Name,
                                                               fire18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  fire18sp <- fire18sp[!is.na(over(fire18sp, geometry(uktownssw))), ]
  
  police18 <- gov_2018[which(gov_2018$code==06330422), ]
  police18 <- subset(police18, "Feature Easting" != "" | "Feature Northing" != "")
  police18$poi_ID <- 1:nrow(police18)
  coords <- cbind(Easting = as.numeric(as.character(police18$`Feature Easting`)),
                  Northing = as.numeric(as.character(police18$`Feature Northing`)))
  police18sp <- SpatialPointsDataFrame(coords, data = data.frame(police18$Name,
                                                                 police18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  police18sp <- police18sp[!is.na(over(police18sp, geometry(uktownssw))), ]
  
  GPs18 <- health_2018[which(health_2018$code==05280369), ]
  GPs18 <- subset(GPs18, "Feature Easting" != "" | "Feature Northing" != "")
  GPs18$poi_ID <- 1:nrow(GPs18)
  coords <- cbind(Easting = as.numeric(as.character(GPs18$`Feature Easting`)),
                  Northing = as.numeric(as.character(GPs18$`Feature Northing`)))
  GPs18sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs18$Name,
                                                              GPs18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  GPs18sp <- GPs18sp[!is.na(over(GPs18sp, geometry(uktownssw))), ]
  
  mental18 <- health_2018[which(health_2018$code==05280372), ]
  mental18 <- subset(mental18, "Feature Easting" != "" | "Feature Northing" != "")
  mental18$poi_ID <- 1:nrow(mental18)
  coords <- cbind(Easting = as.numeric(as.character(mental18$`Feature Easting`)),
                  Northing = as.numeric(as.character(mental18$`Feature Northing`)))
  mental18sp <- SpatialPointsDataFrame(coords, data = data.frame(mental18$Name,
                                                                 mental18$poi_ID), proj4string = CRS("+init=epsg:27700"))
  mental18sp <- mental18sp[!is.na(over(mental18sp, geometry(uktownssw))), ]
  
  
  police18sp <- spTransform(police18sp, CRS("+init=epsg:4326"))
  fire18sp <- spTransform(fire18sp, CRS("+init=epsg:4326"))
  library18sp <- spTransform(library18sp, CRS("+init=epsg:4326"))
  hospitals18sp <- spTransform(hospitals18sp, CRS("+init=epsg:4326"))
  station18sp <- spTransform(station18sp, CRS("+init=epsg:4326"))
  GPs18sp <- spTransform(GPs18sp, CRS("+init=epsg:4326"))
  mental18sp <- spTransform(mental18sp, CRS("+init=epsg:4326"))
  
  icons <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'green')
  iconz <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'pink')
  icona <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'blue')
  iconq <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'orange')
  iconw <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'red')
  icont <- makeAwesomeIcon(icon = 'ios-close', iconColor = 'black', library = 'ion', markerColor = 'purple')
  
  
  mymap <- leaflet(data=police18sp) %>%
    addTiles() %>%  
    setView(lng=-1.3823, lat=53.0975, zoom=5) 
  addAwesomeMarkers(mymap, data=mental18sp, icon = icons, label = "Mental Health", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
     addAwesomeMarkers( data=hospitals18sp,icon = iconq, label = "Hospital", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
     addAwesomeMarkers( data=library18sp,icon = iconw, label = "Library", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
    # addAwesomeMarkers( data=fire18sp,icon = icona, label = "Fire Station", labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
     addAwesomeMarkers( data=station18sp,icon = iconz, label = "Train Station",labelOptions = labelOptions(noHide = T, direction = "auto", opacity = 0.7, style = list("font-style"="strong","font-family"="ariel", "border-color" = "rgba(0,0,0,0.5"  ,"font-size"="15px"))) %>%
    addScaleBar( position="bottomright")
  
  