###load the packages###

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
library(googleway)
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

scotpop01 <- read_csv("scot01pop.csv")
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



region <- read_csv("region.csv")
region <- region %>% distinct(BUA11CD, .keep_all = TRUE)
towns <- merge(towns, region, by="BUA11CD", all.x=TRUE)
towns$Region <- ifelse(is.na(towns$Region), 
                       'Scotland', towns$Region)

towns <- merge(towns, buasize, by="BUA11CD", all.x=TRUE)
##counts##

hmrc18 <- read_csv("hmrc18.csv")
jobcent18 <- read_csv("jobcent18.csv")
station18 <- read_csv("station18.csv")
gps18 <- read_csv("gps18.csv")
hospital18 <- read_csv("hospital18.csv")
furthered18 <- read_csv("furthered18.csv")
mental18 <- read_csv("mental18.csv")

hmrc17 <- read_csv("hmrc17.csv")
jobcent17 <- read_csv("jobcent17.csv")
station17 <- read_csv("station17.csv")
gps17 <- read_csv("gps17.csv")
hospital17 <- read_csv("hospital17.csv")
furthered17 <- read_csv("furthered17.csv")
mental17 <- read_csv("mental17.csv")

hmrc16 <- read_csv("hmrc16.csv")
jobcent16 <- read_csv("jobcent16.csv")
station16 <- read_csv("station16.csv")
gps16 <- read_csv("gps16.csv")
hospital16 <- read_csv("hospital16.csv")
furthered16 <- read_csv("furthered16.csv")
mental16 <- read_csv("mental16.csv")

hmrc15 <- read_csv("hmrc15.csv")
jobcent15 <- read_csv("jobcent15.csv")
station15 <- read_csv("station15.csv")
gps15 <- read_csv("gps15.csv")
hospital15 <- read_csv("hospital15.csv")
furthered15 <- read_csv("furthered15.csv")
mental15 <- read_csv("mental15.csv")

hmrc14 <- read_csv("hmrc14.csv")
jobcent14 <- read_csv("jobcent14.csv")
station14 <- read_csv("station14.csv")
gps14 <- read_csv("gps14.csv")
hospital14 <- read_csv("hospital14.csv")
furthered14 <- read_csv("furthered14.csv")
mental14 <- read_csv("mental14.csv")

hmrc13 <- read_csv("hmrc13.csv")
jobcent13 <- read_csv("jobcent13.csv")
station13 <- read_csv("station13.csv")
gps13 <- read_csv("gps13.csv")
hospital13 <- read_csv("hospital13.csv")
furthered13 <- read_csv("furthered13.csv")
mental13 <- read_csv("mental13.csv")

hmrc12 <- read_csv("hmrc12.csv")
jobcent12 <- read_csv("jobcent12.csv")
station12 <- read_csv("station12.csv")
gps12 <- read_csv("gps12.csv")
hospital12 <- read_csv("hospital12.csv")
furthered12 <- read_csv("furthered12.csv")
mental12 <- read_csv("mental12.csv")

hmrc11 <- read_csv("hmrc11.csv")
jobcent11 <- read_csv("jobcent11.csv")
station11 <- read_csv("station11.csv")
gps11 <- read_csv("gps11.csv")
hospital11 <- read_csv("hospital11.csv")
furthered11 <- read_csv("furthered11.csv")
mental11 <- read_csv("mental11.csv")

station10 <- read_csv("station10.csv")
gps10 <- read_csv("gps10.csv")
hospital10 <- read_csv("hospital10.csv")
furthered10 <- read_csv("furthered10.csv")
mental10 <- read_csv("mental10.csv")

station09 <- read_csv("station09.csv")
gps09 <- read_csv("gps09.csv")
hospital09 <- read_csv("hospital09.csv")
furthered09 <- read_csv("furthered09.csv")
mental09 <- read_csv("mental09.csv")

station08 <- read_csv("station08.csv")
gps08 <- read_csv("gps08.csv")
hospital08 <- read_csv("hospital08.csv")
furthered08 <- read_csv("furthered08.csv")
mental08 <- read_csv("mental08.csv")


###combine###


station <- merge(station18, station17, by="BUA11CD")
station <- merge(station, station16, by="BUA11CD")
station <- merge(station, station15, by="BUA11CD")
station <- merge(station, station14, by="BUA11CD")
station <- merge(station, station13, by="BUA11CD")
station <- merge(station, station12, by="BUA11CD")
station <- merge(station, station11, by="BUA11CD")
station <- merge(station, station10, by="BUA11CD")
station <- merge(station, station09, by="BUA11CD")
station <- merge(station, station08, by="BUA11CD")

hospital <- merge(hospital18, hospital17, by="BUA11CD")
hospital <- merge(hospital, hospital16, by="BUA11CD")
hospital <- merge(hospital, hospital15, by="BUA11CD")
hospital <- merge(hospital, hospital14, by="BUA11CD")
hospital <- merge(hospital, hospital13, by="BUA11CD")
hospital <- merge(hospital, hospital12, by="BUA11CD")
hospital <- merge(hospital, hospital11, by="BUA11CD")
hospital <- merge(hospital, hospital10, by="BUA11CD")
hospital <- merge(hospital, hospital09, by="BUA11CD")
hospital <- merge(hospital, hospital08, by="BUA11CD")

mental <- merge(mental18, mental17, by="BUA11CD")
mental <- merge(mental, mental16, by="BUA11CD")
mental <- merge(mental, mental15, by="BUA11CD")
mental <- merge(mental, mental14, by="BUA11CD")
mental <- merge(mental, mental13, by="BUA11CD")
mental <- merge(mental, mental12, by="BUA11CD")
mental <- merge(mental, mental11, by="BUA11CD")
mental <- merge(mental, mental10, by="BUA11CD")
mental <- merge(mental, mental09, by="BUA11CD")
mental <- merge(mental, mental08, by="BUA11CD")

furthered <- merge(furthered18, furthered17, by="BUA11CD")
furthered <- merge(furthered, furthered16, by="BUA11CD")
furthered <- merge(furthered, furthered15, by="BUA11CD")
furthered <- merge(furthered, furthered14, by="BUA11CD")
furthered <- merge(furthered, furthered13, by="BUA11CD")
furthered <- merge(furthered, furthered12, by="BUA11CD")
furthered <- merge(furthered, furthered11, by="BUA11CD")
furthered <- merge(furthered, furthered10, by="BUA11CD")
furthered <- merge(furthered, furthered09, by="BUA11CD")
furthered <- merge(furthered, furthered08, by="BUA11CD")

gps <- merge(gps18, gps17, by="BUA11CD")
gps <- merge(gps, gps16, by="BUA11CD")
gps <- merge(gps, gps15, by="BUA11CD")
gps <- merge(gps, gps14, by="BUA11CD")
gps <- merge(gps, gps13, by="BUA11CD")
gps <- merge(gps, gps12, by="BUA11CD")
gps <- merge(gps, gps11, by="BUA11CD")
gps <- merge(gps, gps10, by="BUA11CD")
gps <- merge(gps, gps09, by="BUA11CD")
gps <- merge(gps, gps08, by="BUA11CD")

hmrc <- merge(hmrc18, hmrc17, by="BUA11CD")
hmrc <- merge(hmrc, hmrc16, by="BUA11CD")
hmrc <- merge(hmrc, hmrc15, by="BUA11CD")
hmrc <- merge(hmrc, hmrc14, by="BUA11CD")
hmrc <- merge(hmrc, hmrc13, by="BUA11CD")
hmrc <- merge(hmrc, hmrc12, by="BUA11CD")
hmrc <- merge(hmrc, hmrc11, by="BUA11CD")

jobcent <- merge(jobcent18, jobcent17, by="BUA11CD")
jobcent <- merge(jobcent, jobcent16, by="BUA11CD")
jobcent <- merge(jobcent, jobcent15, by="BUA11CD")
jobcent <- merge(jobcent, jobcent14, by="BUA11CD")
jobcent <- merge(jobcent, jobcent13, by="BUA11CD")
jobcent <- merge(jobcent, jobcent12, by="BUA11CD")
jobcent <- merge(jobcent, jobcent11, by="BUA11CD")

station <- station[-c(2,4,6,8,10,12,14,16,18,20,22)]
hospital <- hospital[-c(2,4,6,8,10,12,14,16,18,20,22)]
gps <- gps[-c(2,4,6,8,10,12,14,16,18,20,22)]
furthered <- furthered[-c(2,4,6,8,10,12,14,16,18,20,22)]
mental <- mental[-c(2,4,6,8,10,12,14,16,18,20,22)]
hmrc <- hmrc[-c(2,4,6,8,10,12,14,16)]
jobcent <- jobcent[-c(2,4,6,8,10,12,14,16)]


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns",  "hospital", "hmrc", "jobcent", "mental", "furthered", "gps", "station", "towns" )])

publiclosses <- merge(hospital, furthered, by="BUA11CD")
publiclosses <- merge(publiclosses, gps, by="BUA11CD")
publiclosses <- merge(publiclosses, mental, by="BUA11CD")
publiclosses <- merge(publiclosses, station, by="BUA11CD")
publiclosses <- merge(publiclosses, jobcent, by="BUA11CD")

publiclosses <- merge(publiclosses, towns, by="BUA11CD")


losthosp <- publiclosses[which(publiclosses$hospitals18==0), ]
losthosp$nohosp <- losthosp$hospitals17+losthosp$hospitals16+losthosp$hospitals15+losthosp$hospitals14+losthosp$hospitals13+losthosp$hospitals12+losthosp$hospitals11+losthosp$hospitals10+losthosp$hospitals09+losthosp$hospitals08
lostgp <- publiclosses[which(publiclosses$GPs18==0), ]
lostgp$nogp <- lostgp$GPs17+lostgp$GPs16+lostgp$GPs15+lostgp$GPs14+lostgp$GPs13+lostgp$GPs12+lostgp$GPs11+lostgp$GPs10+lostgp$GPs09+lostgp$GPs08


hospital <- merge(hospital, popbua, by="BUA11CD")
furthered <- merge(furthered, popbua, by="BUA11CD")
hmrc <- merge(hmrc, popbua, by="BUA11CD")
gps <- merge(gps, popbua, by="BUA11CD")
mental <- merge(mental, popbua, by="BUA11CD")
station <- merge(station, popbua, by="BUA11CD")
jobcent <- merge(jobcent, popbua, by="BUA11CD")



hospital$"2018-06-01" <- hospital$hospitals18
hospital$"2017-06-01" <- hospital$hospitals17
hospital$"2016-06-01" <- hospital$hospitals16
hospital$"2015-06-01" <- hospital$hospitals15
hospital$"2014-09-01" <- hospital$hospitals14
hospital$"2013-06-01" <- hospital$hospitals13
hospital$"2012-06-01" <- hospital$hospitals12
hospital$"2011-06-01" <- hospital$hospitals11
hospital$"2010-09-01" <- hospital$hospitals10
hospital$"2009-06-01" <- hospital$hospitals09
hospital$"2008-03-01" <- hospital$hospitals08

jobcent$"2018-06-01" <- jobcent$jobcent18
jobcent$"2017-06-01" <- jobcent$jobcent17
jobcent$"2016-06-01" <- jobcent$jobcent16
jobcent$"2015-06-01" <- jobcent$jobcent15
jobcent$"2014-09-01" <- jobcent$jobcent14
jobcent$"2013-06-01" <- jobcent$jobcent13
jobcent$"2012-06-01" <- jobcent$jobcent12
jobcent$"2011-06-01" <- jobcent$jobcent11

hmrc$"2018-06-01" <- hmrc$hmrc18
hmrc$"2017-06-01" <- hmrc$hmrc17
hmrc$"2016-06-01" <- hmrc$hmrc16
hmrc$"2015-06-01" <- hmrc$hmrc15
hmrc$"2014-09-01" <- hmrc$hmrc14
hmrc$"2013-06-01" <- hmrc$hmrc13
hmrc$"2012-06-01" <- hmrc$hmrc12
hmrc$"2011-06-01" <- hmrc$hmrc11

mental$"2018-06-01" <- mental$mental18
mental$"2017-06-01" <- mental$mental17
mental$"2016-06-01" <- mental$mental16
mental$"2015-06-01" <- mental$mental15
mental$"2014-09-01" <- mental$mental14
mental$"2013-06-01" <- mental$mental13
mental$"2012-06-01" <- mental$mental12
mental$"2011-06-01" <- mental$mental11
mental$"2010-09-01" <- mental$mental10
mental$"2009-06-01" <- mental$mental09
mental$"2008-03-01" <- mental$mental08

gps$"2018-06-01" <- gps$GPs18
gps$"2017-06-01" <- gps$GPs17
gps$"2016-06-01" <- gps$GPs16
gps$"2015-06-01" <- gps$GPs15
gps$"2014-09-01" <- gps$GPs14
gps$"2013-06-01" <- gps$GPs13
gps$"2012-06-01" <- gps$GPs12
gps$"2011-06-01" <- gps$GPs11
gps$"2010-09-01" <- gps$GPs10
gps$"2009-06-01" <- gps$GPs09
gps$"2008-03-01" <- gps$GPs08

station$"2018-06-01" <- station$station18
station$"2017-06-01" <- station$station17
station$"2016-06-01" <- station$station16
station$"2015-06-01" <- station$station15
station$"2014-09-01" <- station$station14
station$"2013-06-01" <- station$station13
station$"2012-06-01" <- station$station12
station$"2011-06-01" <- station$station11
station$"2010-09-01" <- station$station10
station$"2009-06-01" <- station$station09
station$"2008-03-01" <- station$station08

furthered$"2018-06-01" <- furthered$furthered18
furthered$"2017-06-01" <- furthered$furthered17
furthered$"2016-06-01" <- furthered$furthered16
furthered$"2015-06-01" <- furthered$furthered15
furthered$"2014-09-01" <- furthered$furthered14
furthered$"2013-06-01" <- furthered$furthered13
furthered$"2012-06-01" <- furthered$furthered12
furthered$"2011-06-01" <- furthered$furthered11
furthered$"2010-09-01" <- furthered$furthered10
furthered$"2009-06-01" <- furthered$furthered09
furthered$"2008-03-01" <- furthered$furthered08

furthered <- furthered[-c(2,3,4,5,6,7,8,9,10,11,12,13)]
station <- station[-c(2,3,4,5,6,7,8,9,10,11,12,13)]
gps <- gps[-c(2,3,4,5,6,7,8,9,10,11,12,13)]
hospital <- hospital[-c(2,3,4,5,6,7,8,9,10,11,12,13)]
mental <- mental[-c(2,3,4,5,6,7,8,9,10,11,12,13)]
jobcent <- jobcent[-c(2,3,4,5,6,7,8,9,10)]
hmrc <- hmrc[-c(2,3,4,5,6,7,8,9,10)]

furthereds <- furthered[-c(1)] 

furthereds <- reshape(furthereds, idvar = "BUA11CD", ids = furthered$BUA11CD,
                      times = names(furthereds), timevar = "date",
                      varying = list(names(furthereds)),v.names="furthereds", new.row.names = 1:((dim(furthereds)[2])*(dim(furthereds)[1])),direction = "long")

furthereds$date <- as.Date(furthereds$date)


gpss <- gps[-c(1)] 

gpss <- reshape(gpss, idvar = "BUA11CD", ids = gps$BUA11CD,
                times = names(gpss), timevar = "date",
                varying = list(names(gpss)),v.names="gpss", new.row.names = 1:((dim(gpss)[2])*(dim(gpss)[1])),direction = "long")

gpss$date <- as.Date(gpss$date)


jobcents <- jobcent[-c(1)] 

jobcents <- reshape(jobcents, idvar = "BUA11CD", ids = jobcent$BUA11CD,
                    times = names(jobcents), timevar = "date",
                    varying = list(names(jobcents)),v.names="jobcents", new.row.names = 1:((dim(jobcents)[2])*(dim(jobcents)[1])),direction = "long")

jobcents$date <- as.Date(jobcents$date)


hmrcs <- hmrc[-c(1)] 

hmrcs <- reshape(hmrcs, idvar = "BUA11CD", ids = hmrc$BUA11CD,
                 times = names(hmrcs), timevar = "date",
                 varying = list(names(hmrcs)),v.names="hmrcs", new.row.names = 1:((dim(hmrcs)[2])*(dim(hmrcs)[1])),direction = "long")

hmrcs$date <- as.Date(hmrcs$date)


stations <- station[-c(1)] 

stations <- reshape(stations, idvar = "BUA11CD", ids = station$BUA11CD,
                    times = names(stations), timevar = "date",
                    varying = list(names(stations)),v.names="stations", new.row.names = 1:((dim(stations)[2])*(dim(stations)[1])),direction = "long")

stations$date <- as.Date(stations$date)


mentals <- mental[-c(1)] 

mentals <- reshape(mentals, idvar = "BUA11CD", ids = mental$BUA11CD,
                   times = names(mentals), timevar = "date",
                   varying = list(names(mentals)),v.names="mentals", new.row.names = 1:((dim(mentals)[2])*(dim(mentals)[1])),direction = "long")

mentals$date <- as.Date(mentals$date)


hospitals <- hospital[-c(1)] 

hospitals <- reshape(hospitals, idvar = "BUA11CD", ids = hospital$BUA11CD,
                     times = names(hospitals), timevar = "date",
                     varying = list(names(hospitals)),v.names="hospitals", new.row.names = 1:((dim(hospitals)[2])*(dim(hospitals)[1])),direction = "long")

hospitals$date <- as.Date(hospitals$date)


publicbine <- merge(mentals, hospitals, by=c("BUA11CD", "date"))
publicbine <- merge(publicbine, stations, by=c("BUA11CD", "date"))
publicbine <- merge(publicbine, furthereds, by=c("BUA11CD", "date"))
publicbine <- merge(publicbine, gpss, by=c("BUA11CD", "date"))
publicbine <- merge(publicbine, hmrcs, by=c("BUA11CD", "date"), all.x=TRUE)
publicbine <- merge(publicbine, jobcents, by=c("BUA11CD", "date"), all.x=TRUE)

##jsa##

#scotjsa <- read_csv("scotjsa.csv")
#jsa08 <- read_csv("jsa08.csv")
#jsa14 <- read_csv("jsa14.csv")#

#scotjsa <- scotjsa[complete.cases(scotjsa$BUA11CD), ]
#jsa08 <- jsa08[complete.cases(jsa08$BUA11CD), ]
#jsa14 <- jsa14[complete.cases(jsa14$BUA11CD), ]
#
#scotjsa <-  aggregate(scotjsa[-1], scotjsa["BUA11CD"], sum)
#jsa08 <-  aggregate(jsa08[-1], jsa08["BUA11CD"], sum)
#jsa14 <-  aggregate(jsa14[-1], jsa14["BUA11CD"], sum)
#
#
#scotjsa <- merge(scotjsa, popbua, by="BUA11CD")
#jsa08 <- merge(jsa08, popbua, by="BUA11CD")
#jsa14 <- merge(jsa14, popbua, by="BUA11CD")


#scotjsa <- scotjsa[which(scotjsa$poplsoa<175000&scotjsa$poplsoa>10000), ]
#jsa08 <- jsa08[which(jsa08$poplsoa<175000&jsa08$poplsoa>10000), ]
#jsa14 <- jsa14[which(jsa14$poplsoa<175000&jsa14$poplsoa>10000), ]

#scotjsas <- scotjsa[-c(1,13)]
#jsa08s <- jsa08[-c(1,8)]
jsa14s <- jsa14[-c(1,7)]


#jsa08s <- reshape(jsa08s, idvar = "BUA11CD", ids = jsa08$BUA11CD,
#                  times = names(jsa08s), timevar = "date",
 #                 varying = list(names(jsa08s)),v.names="jsa", new.row.names = 1:((dim(jsa08s)[2])*(dim(jsa08s)[1])),direction = "long")
#
#
#jsa14s <- reshape(jsa14s, idvar = "BUA11CD", ids = jsa14$BUA11CD,
 #                 times = names(jsa14s), timevar = "date",
  #                varying = list(names(jsa14s)),v.names="jsa", new.row.names = 1:((dim(jsa14s)[2])*(dim(jsa14s)[1])),direction = "long")
#
#scotjsas <- reshape(scotjsas, idvar = "BUA11CD", ids = scotjsa$BUA11CD,
 #                   times = names(scotjsas), timevar = "date",
  #                  varying = list(names(scotjsas)),v.names="jsa", new.row.names = 1:((dim(scotjsas)[2])*(dim(scotjsas)[1])),direction = "long")
#
#jobseekers <- merge(jsa08s, jsa14s, by=c("BUA11CD", "date", "jsa"), all=TRUE)
#jobseekers <- merge(jobseekers, scotjsas, by=c("BUA11CD", "date", "jsa"), all=TRUE)


#publicbine <- merge(publicbine, jobseekers, by=c("BUA11CD", "date"))


##counts##

schools148 <- read_csv("school148.csv")
pbox148 <- read_csv("post148.csv")
wifi148 <- read_csv("wifi148.csv")
hall148 <- read_csv("hall148.csv")
bus148 <- read_csv("bus148.csv")
nursery148 <- read_csv("nursery148.csv")
health148 <- read_csv("health148.csv")

schools13 <- read_csv("school13.csv")
pbox13 <- read_csv("pbox13.csv")
wifi13 <- read_csv("wifi13.csv")
hall13 <- read_csv("hall13.csv")
bus13 <- read_csv("bus13.csv")
nursery13 <- read_csv("nursery13.csv")
health13 <- read_csv("health13.csv")

schools12 <- read_csv("school12.csv")
pbox12 <- read_csv("pbox12.csv")
wifi12 <- read_csv("wifi12.csv")
hall12 <- read_csv("hall12.csv")
bus12 <- read_csv("bus12.csv")
nursery12 <- read_csv("nursery12.csv")
health12 <- read_csv("health12.csv")

schools11 <- read_csv("school11.csv")
pbox11 <- read_csv("pbox11.csv")
wifi11 <- read_csv("wifi11.csv")
hall11 <- read_csv("hall11.csv")
bus11 <- read_csv("bus11.csv")
nursery11 <- read_csv("nursery11.csv")
health11 <- read_csv("health11.csv")

schools10 <- read_csv("school10.csv")
pbox10 <- read_csv("pbox10.csv")
wifi10 <- read_csv("wifi10.csv")
hall10 <- read_csv("hall10.csv")
nursery10 <- read_csv("nursery10.csv")
health10 <- read_csv("health10.csv")

schools09 <- read_csv("school09.csv")
pbox09 <- read_csv("pbox09.csv")
wifi09 <- read_csv("wifi09.csv")
hall09 <- read_csv("hall09.csv")
nursery09 <- read_csv("nursery09.csv")
health09 <- read_csv("health09.csv")

schools08 <- read_csv("school08.csv")
pbox08 <- read_csv("pbox08.csv")
wifi08 <- read_csv("wifi08.csv")
hall08 <- read_csv("hall08.csv")
nursery08 <- read_csv("nursery08.csv")
health08 <- read_csv("health08.csv")



bus <- merge(bus148, bus11, by="BUA11CD")
bus <- merge(bus, bus12, by="BUA11CD")
bus <- merge(bus, bus13, by="BUA11CD")

hall <- merge(hall148, hall08, by="BUA11CD")
hall <- merge(hall, hall09, by="BUA11CD")
hall <- merge(hall, hall10, by="BUA11CD")
hall <- merge(hall, hall11, by="BUA11CD")
hall <- merge(hall, hall12, by="BUA11CD")
hall <- merge(hall, hall13, by="BUA11CD")

pbox <- merge(pbox148, pbox08, by="BUA11CD")
pbox <- merge(pbox, pbox09, by="BUA11CD")
pbox <- merge(pbox, pbox10, by="BUA11CD")
pbox <- merge(pbox, pbox11, by="BUA11CD")
pbox <- merge(pbox, pbox12, by="BUA11CD")
pbox <- merge(pbox, pbox13, by="BUA11CD")

health <- merge(health148, health08, by="BUA11CD")
health <- merge(health, health09, by="BUA11CD")
health <- merge(health, health10, by="BUA11CD")
health <- merge(health, health11, by="BUA11CD")
health <- merge(health, health12, by="BUA11CD")
health <- merge(health, health13, by="BUA11CD")

wifi <- merge(wifi148, wifi08, by="BUA11CD")
wifi <- merge(wifi, wifi09, by="BUA11CD")
wifi <- merge(wifi, wifi10, by="BUA11CD")
wifi <- merge(wifi, wifi11, by="BUA11CD")
wifi <- merge(wifi, wifi12, by="BUA11CD")
wifi <- merge(wifi, wifi13, by="BUA11CD")

nursery <- merge(nursery148, nursery08, by="BUA11CD")
nursery <- merge(nursery, nursery09, by="BUA11CD")
nursery <- merge(nursery, nursery10, by="BUA11CD")
nursery <- merge(nursery, nursery11, by="BUA11CD")
nursery <- merge(nursery, nursery12, by="BUA11CD")
nursery <- merge(nursery, nursery13, by="BUA11CD")

schools <- merge(schools148, schools08, by="BUA11CD")
schools <- merge(schools, schools09, by="BUA11CD")
schools <- merge(schools, schools10, by="BUA11CD")
schools <- merge(schools, schools11, by="BUA11CD")
schools <- merge(schools, schools12, by="BUA11CD")
schools <- merge(schools, schools13, by="BUA11CD")

health <- health[-c(2,8,10,12,14,16,18)]
schools <- schools[-c(2,8,10,12,14,16,18)]
wifi <- wifi[-c(2,8,10,12,14,16,18)]
pbox <- pbox[-c(2,8,10,12,14,16,18)]
bus <- bus[-c(2,8,10,12)]
hall <- hall[-c(2,8,10,12,14,16,18)]
nursery <- nursery[-c(2,8,10,12,14,16,18)]


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns",  "bus", "wifi", "nursery", "hall", "schools", "health", "pbox", "towns", "publicbine" )])

hall <- merge(hall, popbua, by="BUA11CD")
nursery <- merge(nursery, popbua, by="BUA11CD")
bus <- merge(bus, popbua, by="BUA11CD")
wifi <- merge(wifi, popbua, by="BUA11CD")
pbox <- merge(pbox, popbua, by="BUA11CD")
schools <- merge(schools, popbua, by="BUA11CD")
health <- merge(health, popbua, by="BUA11CD")

publicchange <- merge(nursery, wifi, by="BUA11CD")
publicchange <- merge(publicchange, hall, by="BUA11CD")
publicchange <- merge(publicchange, pbox, by="BUA11CD")
publicchange <- merge(publicchange, schools, by="BUA11CD")
publicchange <- merge(publicchange, health, by="BUA11CD")
publicchange <- merge(publicchange, bus, by="BUA11CD")


publicchange$buschange <- publicchange$bus18-publicchange$bus11
publicchange$healthchange <- publicchange$health18-publicchange$health09
publicchange$schoolchange <- publicchange$school18-publicchange$schools09
publicchange$nurserychange <- publicchange$nursery18-publicchange$nursery09

publicchange <- publicchange[c("BUA11CD", "nurserychange", "healthchange", "schoolchange", "buschange")]
publicchange <- merge(publicchange, buanames, by="BUA11CD")

nursery$"2018-06-01" <- nursery$nursery18
nursery$"2017-06-01" <- nursery$nursery17
nursery$"2016-06-01" <- nursery$nursery16
nursery$"2015-06-01" <- nursery$nursery15
nursery$"2014-09-01" <- nursery$nursery14
nursery$"2013-06-01" <- nursery$nursery13
nursery$"2012-06-01" <- nursery$nursery12
nursery$"2011-06-01" <- nursery$nursery11
nursery$"2010-09-01" <- nursery$nursery10
nursery$"2009-06-01" <- nursery$nursery09
nursery$"2008-03-01" <- nursery$nursery08

schools$"2018-06-01" <- schools$school18
schools$"2017-06-01" <- schools$school17
schools$"2016-06-01" <- schools$school16
schools$"2015-06-01" <- schools$school15
schools$"2014-09-01" <- schools$school14
schools$"2013-06-01" <- schools$schools13
schools$"2012-06-01" <- schools$schools12
schools$"2011-06-01" <- schools$schools11
schools$"2010-09-01" <- schools$schools10
schools$"2009-06-01" <- schools$schools09
schools$"2008-03-01" <- schools$schools08

health$"2018-06-01" <- health$health18
health$"2017-06-01" <- health$health17
health$"2016-06-01" <- health$health16
health$"2015-06-01" <- health$health15
health$"2014-09-01" <- health$health14
health$"2013-06-01" <- health$health13
health$"2012-06-01" <- health$health12
health$"2011-06-01" <- health$health11
health$"2010-09-01" <- health$health10
health$"2009-06-01" <- health$health09
health$"2008-03-01" <- health$health08

bus$"2018-06-01" <- bus$bus18
bus$"2017-06-01" <- bus$bus17
bus$"2016-06-01" <- bus$bus16
bus$"2015-06-01" <- bus$bus15
bus$"2014-09-01" <- bus$bus14
bus$"2013-06-01" <- bus$bus13
bus$"2012-06-01" <- bus$bus12
bus$"2011-06-01" <- bus$bus11

pbox$"2018-06-01" <- pbox$post18
pbox$"2017-06-01" <- pbox$post17
pbox$"2016-06-01" <- pbox$post16
pbox$"2015-06-01" <- pbox$post15
pbox$"2014-09-01" <- pbox$post14
pbox$"2013-06-01" <- pbox$pbox13
pbox$"2012-06-01" <- pbox$pbox12
pbox$"2011-06-01" <- pbox$pbox11
pbox$"2010-09-01" <- pbox$pbox10
pbox$"2009-06-01" <- pbox$pbox09
pbox$"2008-03-01" <- pbox$pbox08

wifi$"2018-06-01" <- wifi$wifi18
wifi$"2017-06-01" <- wifi$wifi17
wifi$"2016-06-01" <- wifi$wifi16
wifi$"2015-06-01" <- wifi$wifi15
wifi$"2014-09-01" <- wifi$wifi14
wifi$"2013-06-01" <- wifi$wifi13
wifi$"2012-06-01" <- wifi$wifi12
wifi$"2011-06-01" <- wifi$wifi11
wifi$"2010-09-01" <- wifi$wifi10
wifi$"2009-06-01" <- wifi$wifi09
wifi$"2008-03-01" <- wifi$wifi08

hall$"2018-06-01" <- hall$hall18
hall$"2017-06-01" <- hall$hall17
hall$"2016-06-01" <- hall$hall16
hall$"2015-06-01" <- hall$hall15
hall$"2014-09-01" <- hall$hall14
hall$"2013-06-01" <- hall$halls13
hall$"2012-06-01" <- hall$halls12
hall$"2011-06-01" <- hall$halls11
hall$"2010-09-01" <- hall$halls10
hall$"2009-06-01" <- hall$halls09
hall$"2008-03-01" <- hall$halls08


#publicdense <- merge(nursery, wifi, by=c("BUA11CD"))
#publicdense <- merge(publicdense, hall, by=c("BUA11CD"))
#publicdense <- merge(publicdense, pbox, by=c("BUA11CD"))
#publicdense <- merge(publicdense, schools, by=c("BUA11CD"))
#publicdense <- merge(publicdense, health, by=c("BUA11CD"))
#publicdense <- merge(publicdense, bus, by=c("BUA11CD"))

#publicdenses <- publicdenses[-c(13,14,15,16,17,18,19,20,21,22,23,24,36,37,38,39,40,41,42,43,44,45,46,47,59,60,61,62,63,64,65,66,67,68,69,70,82,83,84,85,86,87,88,89,90,91,92,93,105,106,107,108,109,110,111,112,113,114,115,116,128,129,130,131,132,133,134,135,136,137,138,139, 149,150,151,152,153,154,155,156)]
#publicdenses <-  publicdenses %>%
#  psycho::standardize()


hall <- hall[-c(2,3,4,5,6,7,8,9,10,11,12)]
wifi <- wifi[-c(2,3,4,5,6,7,8,9,10,11,12)]
pbox <- pbox[-c(2,3,4,5,6,7,8,9,10,11,12)]
bus <- bus[-c(2,3,4,5,6,7,8,9)]
health <- health[-c(2,3,4,5,6,7,8,9,10,11,12)]
schools <- schools[-c(2,3,4,5,6,7,8,9,10,11,12)]
nursery <- nursery[-c(2,3,4,5,6,7,8,9,10,11,12)]



schoolss <- schools[-c(1,2)] 

schoolss <- reshape(schoolss, idvar = "BUA11CD", ids = schools$BUA11CD,
                    times = names(schoolss), timevar = "date",
                    varying = list(names(schoolss)),v.names="schoolss", new.row.names = 1:((dim(schoolss)[2])*(dim(schoolss)[1])),direction = "long")

schoolss$date <- as.Date(schoolss$date)
healths <- health[-c(1,2)] 

healths <- reshape(healths, idvar = "BUA11CD", ids = health$BUA11CD,
                   times = names(healths), timevar = "date",
                   varying = list(names(healths)),v.names="healths", new.row.names = 1:((dim(healths)[2])*(dim(healths)[1])),direction = "long")

healths$date <- as.Date(healths$date)
pboxs <- pbox[-c(1,2)] 

pboxs <- reshape(pboxs, idvar = "BUA11CD", ids = pbox$BUA11CD,
                 times = names(pboxs), timevar = "date",
                 varying = list(names(pboxs)),v.names="pboxs", new.row.names = 1:((dim(pboxs)[2])*(dim(pboxs)[1])),direction = "long")

pboxs$date <- as.Date(pboxs$date)
buss <- bus[-c(1,2)] 

buss <- reshape(buss, idvar = "BUA11CD", ids = bus$BUA11CD,
                times = names(buss), timevar = "date",
                varying = list(names(buss)),v.names="buss", new.row.names = 1:((dim(buss)[2])*(dim(buss)[1])),direction = "long")

buss$date <- as.Date(buss$date)
halls <- hall[-c(1,2)] 

halls <- reshape(halls, idvar = "BUA11CD", ids = hall$BUA11CD,
                 times = names(halls), timevar = "date",
                 varying = list(names(halls)),v.names="halls", new.row.names = 1:((dim(halls)[2])*(dim(halls)[1])),direction = "long")

halls$date <- as.Date(halls$date)
wifis <- wifi[-c(1,2)] 

wifis <- reshape(wifis, idvar = "BUA11CD", ids = wifi$BUA11CD,
                 times = names(wifis), timevar = "date",
                 varying = list(names(wifis)),v.names="wifis", new.row.names = 1:((dim(wifis)[2])*(dim(wifis)[1])),direction = "long")

wifis$date <- as.Date(wifis$date)
nurserys <- nursery[-c(1,2)] 

nurserys <- reshape(nurserys, idvar = "BUA11CD", ids = nursery$BUA11CD,
                    times = names(nurserys), timevar = "date",
                    varying = list(names(nurserys)),v.names="nurserys", new.row.names = 1:((dim(nurserys)[2])*(dim(nurserys)[1])),direction = "long")

nurserys$date <- as.Date(nurserys$date)

publicdense <- merge(nurserys, wifis, by=c("BUA11CD", "date"))
publicdense <- merge(publicdense, halls, by=c("BUA11CD", "date"))
publicdense <- merge(publicdense, pboxs, by=c("BUA11CD", "date"))
publicdense <- merge(publicdense, schoolss, by=c("BUA11CD", "date"))
publicdense <- merge(publicdense, healths, by=c("BUA11CD", "date"))
publicdense <- merge(publicdense, buss, by=c("BUA11CD", "date"), all.x=TRUE)

paneldata <- merge(publicdense, publicbine, by=c("BUA11CD", "date"), all.x=TRUE)


##counts##

creative148 <- read_csv("creative148.csv")
sport148 <- read_csv("sport148.csv")
cine148 <- read_csv("cine148.csv")
night148 <- read_csv("night148.csv")
bet148 <- read_csv("bet148.csv")

fire18 <- read_csv("fire18.csv")
fire17 <- read_csv("fire17.csv")
fire16 <- read_csv("fire16.csv")
fire15 <- read_csv("fire15.csv")
fire14 <- read_csv("fire14.csv")
fire13 <- read_csv("fire13.csv")
fire12 <- read_csv("fire12.csv")
fire11 <- read_csv("fire11.csv")
fire10 <- read_csv("fire10.csv")
fire09 <- read_csv("fire09.csv")
fire08 <- read_csv("fire08.csv")

police18 <- read_csv("police18.csv")
police17 <- read_csv("police17.csv")
police16 <- read_csv("police16.csv")
police15 <- read_csv("police15.csv")
police14 <- read_csv("police14.csv")
police13 <- read_csv("police13.csv")
police12 <- read_csv("police12.csv")
police11 <- read_csv("police11.csv")
police10 <- read_csv("police10.csv")
police09 <- read_csv("police09.csv")
police08 <- read_csv("police08.csv")

library18 <- read_csv("library18.csv")
library17 <- read_csv("library17.csv")
library16 <- read_csv("library16.csv")
library15 <- read_csv("library15.csv")
library14 <- read_csv("library14.csv")
library13 <- read_csv("library13.csv")
library12 <- read_csv("library12.csv")
library11 <- read_csv("library11.csv")
library10 <- read_csv("library10.csv")
library09 <- read_csv("library09.csv")
library08 <- read_csv("library08.csv")

creative13 <- read_csv("creative13.csv")
sport13 <- read_csv("sport13.csv")
cine13 <- read_csv("cine13.csv")
night13 <- read_csv("night13.csv")
bet13 <- read_csv("bet13.csv")

creative12 <- read_csv("creative12.csv")
sport12 <- read_csv("sport12.csv")
cine12 <- read_csv("cine12.csv")
night12 <- read_csv("night12.csv")
bet12 <- read_csv("bet12.csv")

creative11 <- read_csv("creative11.csv")
sport11 <- read_csv("sport11.csv")
cine11 <- read_csv("cine11.csv")
night11 <- read_csv("night11.csv")
bet11 <- read_csv("bet11.csv")

creative10 <- read_csv("creative10.csv")
sport10 <- read_csv("sport10.csv")
cine10 <- read_csv("cine10.csv")
night10 <- read_csv("night10.csv")
bet10 <- read_csv("bet10.csv")

creative09 <- read_csv("creative09.csv")
sport09 <- read_csv("sport09.csv")
cine09 <- read_csv("cine09.csv")
night09 <- read_csv("night09.csv")
bet09 <- read_csv("bet09.csv")

creative08 <- read_csv("creative08.csv")
sport08 <- read_csv("sport08.csv")
cine08 <- read_csv("cine08.csv")
night08 <- read_csv("night08.csv")
bet08 <- read_csv("bet08.csv")


###combine##

creative <- merge(creative148, creative08, by="BUA11CD")
creative <- merge(creative, creative09, by="BUA11CD")
creative <- merge(creative, creative10, by="BUA11CD")
creative <- merge(creative, creative11, by="BUA11CD")
creative <- merge(creative, creative12, by="BUA11CD")
creative <- merge(creative, creative13, by="BUA11CD")

cine <- merge(cine148, cine08, by="BUA11CD")
cine <- merge(cine, cine09, by="BUA11CD")
cine <- merge(cine, cine10, by="BUA11CD")
cine <- merge(cine, cine11, by="BUA11CD")
cine <- merge(cine, cine12, by="BUA11CD")
cine <- merge(cine, cine13, by="BUA11CD")

bet <- merge(bet148, bet08, by="BUA11CD")
bet <- merge(bet, bet09, by="BUA11CD")
bet <- merge(bet, bet10, by="BUA11CD")
bet <- merge(bet, bet11, by="BUA11CD")
bet <- merge(bet, bet12, by="BUA11CD")
bet <- merge(bet, bet13, by="BUA11CD")

night <- merge(night148, night08, by="BUA11CD")
night <- merge(night, night09, by="BUA11CD")
night <- merge(night, night10, by="BUA11CD")
night <- merge(night, night11, by="BUA11CD")
night <- merge(night, night12, by="BUA11CD")
night <- merge(night, night13, by="BUA11CD")

sport <- merge(sport148, sport08, by="BUA11CD")
sport <- merge(sport, sport09, by="BUA11CD")
sport <- merge(sport, sport10, by="BUA11CD")
sport <- merge(sport, sport11, by="BUA11CD")
sport <- merge(sport, sport12, by="BUA11CD")
sport <- merge(sport, sport13, by="BUA11CD")

fire <- merge(fire18, fire17, by="BUA11CD")
fire <- merge(fire, fire16, by="BUA11CD")
fire <- merge(fire, fire15, by="BUA11CD")
fire <- merge(fire, fire14, by="BUA11CD")
fire <- merge(fire, fire13, by="BUA11CD")
fire <- merge(fire, fire12, by="BUA11CD")
fire <- merge(fire, fire11, by="BUA11CD")
fire <- merge(fire, fire10, by="BUA11CD")
fire <- merge(fire, fire09, by="BUA11CD")
fire <- merge(fire, fire08, by="BUA11CD")

police <- merge(police18, police17, by="BUA11CD")
police <- merge(police, police16, by="BUA11CD")
police <- merge(police, police15, by="BUA11CD")
police <- merge(police, police14, by="BUA11CD")
police <- merge(police, police13, by="BUA11CD")
police <- merge(police, police12, by="BUA11CD")
police <- merge(police, police11, by="BUA11CD")
police <- merge(police, police10, by="BUA11CD")
police <- merge(police, police09, by="BUA11CD")
police <- merge(police, police08, by="BUA11CD")

library <- merge(library18, library17, by="BUA11CD")
library <- merge(library, library16, by="BUA11CD")
library <- merge(library, library15, by="BUA11CD")
library <- merge(library, library14, by="BUA11CD")
library <- merge(library, library13, by="BUA11CD")
library <- merge(library, library12, by="BUA11CD")
library <- merge(library, library11, by="BUA11CD")
library <- merge(library, library10, by="BUA11CD")
library <- merge(library, library09, by="BUA11CD")
library <- merge(library, library08, by="BUA11CD")

sport <- sport[-c(2,8,10,12,14,16,18)]
creative <- creative[-c(2,8,10,12,14,16,18)]
night <- night[-c(2,8,10,12,14,16,18)]
cine <- cine[-c(2,8,10,12,14,16,18)]
bet <- bet[-c(2,8,10,12,14,16,18)]
police <- police[-c(2,8,10,12,14,16,18)]
fire <- fire[-c(2,8,10,12,14,16,18)]
library <- library[-c(2,8,10,12,14,16,18)]

bet <- merge(bet, popbua, by="BUA11CD")
creative <- merge(creative, popbua, by="BUA11CD")
sport <- merge(sport, popbua, by="BUA11CD")
night <- merge(night, popbua, by="BUA11CD")
cine <- merge(cine, popbua, by="BUA11CD")
fire <- merge(fire, popbua, by="BUA11CD")
police <- merge(police, popbua, by="BUA11CD")
library <- merge(library, popbua, by="BUA11CD")



bet$"2018-06-01" <- bet$bet18
bet$"2017-06-01" <- bet$bet17
bet$"2016-06-01" <- bet$bet16
bet$"2015-06-01" <- bet$bet15
bet$"2014-09-01" <- bet$bet14
bet$"2013-06-01" <- bet$bet13
bet$"2012-06-01" <- bet$bet12
bet$"2011-06-01" <- bet$bet11
bet$"2010-09-01" <- bet$bet10
bet$"2009-06-01" <- bet$bet09
bet$"2008-03-01" <- bet$bet08


fire$"2018-06-01" <- fire$fire18
fire$"2017-06-01" <- fire$fire17
fire$"2016-06-01" <- fire$fire16
fire$"2015-06-01" <- fire$fire15
fire$"2014-09-01" <- fire$fire14
fire$"2013-06-01" <- fire$fire13
fire$"2012-06-01" <- fire$fire12
fire$"2011-06-01" <- fire$fire11
fire$"2010-09-01" <- fire$fire10
fire$"2009-06-01" <- fire$fire09
fire$"2008-03-01" <- fire$fire08

police$"2018-06-01" <- police$police18
police$"2017-06-01" <- police$police17
police$"2016-06-01" <- police$police16
police$"2015-06-01" <- police$police15
police$"2014-09-01" <- police$police14
police$"2013-06-01" <- police$police13
police$"2012-06-01" <- police$police12
police$"2011-06-01" <- police$police11
police$"2010-09-01" <- police$police10
police$"2009-06-01" <- police$police09
police$"2008-03-01" <- police$police08
library$"2018-06-01" <- library$library18
library$"2017-06-01" <- library$library17
library$"2016-06-01" <- library$library16
library$"2015-06-01" <- library$library15
library$"2014-09-01" <- library$library14
library$"2013-06-01" <- library$library13
library$"2012-06-01" <- library$library12
library$"2011-06-01" <- library$library11
library$"2010-09-01" <- library$library10
library$"2009-06-01" <- library$library09
library$"2008-03-01" <- library$library08

creative$"2018-06-01" <- creative$creative18
creative$"2017-06-01" <- creative$creative17
creative$"2016-06-01" <- creative$creative16
creative$"2015-06-01" <- creative$creative15
creative$"2014-09-01" <- creative$creative14
creative$"2013-06-01" <- creative$creative13
creative$"2012-06-01" <- creative$creative12
creative$"2011-06-01" <- creative$creative11
creative$"2010-09-01" <- creative$creative10
creative$"2009-06-01" <- creative$creative09
creative$"2008-03-01" <- creative$creative08

sport$"2018-06-01" <- sport$sport18
sport$"2017-06-01" <- sport$sport17
sport$"2016-06-01" <- sport$sport16
sport$"2015-06-01" <- sport$sport15
sport$"2014-09-01" <- sport$sport14
sport$"2013-06-01" <- sport$sport13
sport$"2012-06-01" <- sport$sport12
sport$"2011-06-01" <- sport$sport11
sport$"2010-09-01" <- sport$sport10
sport$"2009-06-01" <- sport$sport09
sport$"2008-03-01" <- sport$sport08

night$"2018-06-01" <- night$night18
night$"2017-06-01" <- night$night17
night$"2016-06-01" <- night$night16
night$"2015-06-01" <- night$night15
night$"2014-09-01" <- night$night14
night$"2013-06-01" <- night$night13
night$"2012-06-01" <- night$night12
night$"2011-06-01" <- night$night11
night$"2010-09-01" <- night$night10
night$"2009-06-01" <- night$night09
night$"2008-03-01" <- night$night08

cine$"2018-06-01" <- cine$cine18
cine$"2017-06-01" <- cine$cine17
cine$"2016-06-01" <- cine$cine16
cine$"2015-06-01" <- cine$cine15
cine$"2014-09-01" <- cine$cine14
cine$"2013-06-01" <- cine$cine13
cine$"2012-06-01" <- cine$cine12
cine$"2011-06-01" <- cine$cine11
cine$"2010-09-01" <- cine$cine10
cine$"2009-06-01" <- cine$cine09
cine$"2008-03-01" <- cine$cine08



cine <- cine[-c(2,3,4,5,6,7,8,9,10,12,11,13)]
creative <- creative[-c(2,3,4,5,6,7,8,9,10,12,11,13)]
night <- night[-c(2,3,4,5,6,7,8,9,10,12,11,13)]
sport <- sport[-c(2,3,4,5,6,7,8,9,10,12,11,13)]
bet <- bet[-c(2,3,4,5,6,7,8,9,10,12,11,13)]
fire <- fire[-c(2,3,4,5,6,7,8,9,10,12,11,13,14,15,16,17)]
police <- police[-c(2,3,4,5,6,7,8,9,10,12,11,13,14,15,16,17)]
library <- library[-c(2,3,4,5,6,7,8,9,10,12,11,13,14,15,16,17)]


cines <- cine[-c(1)]

cines <- reshape(cines, idvar = "BUA11CD", ids = cine$BUA11CD,
                 times = names(cines), timevar = "date",
                 varying = list(names(cines)),v.names="cines", new.row.names = 1:((dim(cines)[2])*(dim(cines)[1])),direction = "long")

cines$date <- as.Date(cines$date)

polices <- police[-c(1)]

polices <- reshape(polices, idvar = "BUA11CD", ids = police$BUA11CD,
                 times = names(polices), timevar = "date",
                 varying = list(names(polices)),v.names="polices", new.row.names = 1:((dim(polices)[2])*(dim(polices)[1])),direction = "long")

polices$date <- as.Date(polices$date)

librarys <- library[-c(1)]

librarys <- reshape(librarys, idvar = "BUA11CD", ids = library$BUA11CD,
                   times = names(librarys), timevar = "date",
                   varying = list(names(librarys)),v.names="librarys", new.row.names = 1:((dim(librarys)[2])*(dim(librarys)[1])),direction = "long")

librarys$date <- as.Date(librarys$date)

fires <- fire[-c(1)]

fires <- reshape(fires, idvar = "BUA11CD", ids = fire$BUA11CD,
                 times = names(fires), timevar = "date",
                 varying = list(names(fires)),v.names="fires", new.row.names = 1:((dim(fires)[2])*(dim(fires)[1])),direction = "long")

fires$date <- as.Date(fires$date)

bets <- bet[-c(1)]

bets <- reshape(bets, idvar = "BUA11CD", ids = bet$BUA11CD,
                times = names(bets), timevar = "date",
                varying = list(names(bets)),v.names="bets", new.row.names = 1:((dim(bets)[2])*(dim(bets)[1])),direction = "long")

bets$date <- as.Date(bets$date)

creatives <- creative[-c(1)]

creatives <- reshape(creatives, idvar = "BUA11CD", ids = creative$BUA11CD,
                     times = names(creatives), timevar = "date",
                     varying = list(names(creatives)),v.names="creatives", new.row.names = 1:((dim(creatives)[2])*(dim(creatives)[1])),direction = "long")

creatives$date <- as.Date(creatives$date)

nights <- night[-c(1)]

nights <- reshape(nights, idvar = "BUA11CD", ids = night$BUA11CD,
                  times = names(nights), timevar = "date",
                  varying = list(names(nights)),v.names="nights", new.row.names = 1:((dim(nights)[2])*(dim(nights)[1])),direction = "long")

nights$date <- as.Date(nights$date)

sports <- sport[-c(1)]

sports <- reshape(sports, idvar = "BUA11CD", ids = sport$BUA11CD,
                  times = names(sports), timevar = "date",
                  varying = list(names(sports)),v.names="sports", new.row.names = 1:((dim(sports)[2])*(dim(sports)[1])),direction = "long")

sports$date <- as.Date(sports$date)




cinezz <- cinez[-c(1)]

cinezz <- reshape(cinezz, idvar = "BUA11CD", ids = cinez$BUA11CD,
                  times = names(cinezz), timevar = "date",
                  varying = list(names(cinezz)),v.names="cinezz", new.row.names = 1:((dim(cinezz)[2])*(dim(cinezz)[1])),direction = "long")

cinezz$date <- as.Date(cinezz$date)

betzz <- betz[-c(1)]

betzz <- reshape(betzz, idvar = "BUA11CD", ids = betz$BUA11CD,
                 times = names(betzz), timevar = "date",
                 varying = list(names(betzz)),v.names="betzz", new.row.names = 1:((dim(betzz)[2])*(dim(betzz)[1])),direction = "long")

betzz$date <- as.Date(betzz$date)

creativezz <- creativez[-c(1)]

creativezz <- reshape(creativezz, idvar = "BUA11CD", ids = creativez$BUA11CD,
                      times = names(creativezz), timevar = "date",
                      varying = list(names(creativezz)),v.names="creativezz", new.row.names = 1:((dim(creativezz)[2])*(dim(creativezz)[1])),direction = "long")

creativezz$date <- as.Date(creativezz$date)

nightzz <- nightz[-c(1)]

nightzz <- reshape(nightzz, idvar = "BUA11CD", ids = nightz$BUA11CD,
                   times = names(nightzz), timevar = "date",
                   varying = list(names(nightzz)),v.names="nightzz", new.row.names = 1:((dim(nightzz)[2])*(dim(nightzz)[1])),direction = "long")

nightzz$date <- as.Date(nightzz$date)

sportzz <- sportz[-c(1)]

sportzz <- reshape(sportzz, idvar = "BUA11CD", ids = sportz$BUA11CD,
                   times = names(sportzz), timevar = "date",
                   varying = list(names(sportzz)),v.names="sportzz", new.row.names = 1:((dim(sportzz)[2])*(dim(sportzz)[1])),direction = "long")

sportzz$date <- as.Date(sportzz$date)

entertainz <- merge(cinezz, nightzz, by=c("BUA11CD", "date"))
entertainz <- merge(entertainz, betzz, by=c("BUA11CD", "date"))
entertainz <- merge(entertainz, sportzz, by=c("BUA11CD", "date"))
entertainz <- merge(entertainz, creativezz, by=c("BUA11CD", "date"))





entertain <- merge(cines, nights, by=c("BUA11CD", "date"))
entertain <- merge(entertain, bets, by=c("BUA11CD", "date"))
entertain <- merge(entertain, sports, by=c("BUA11CD", "date"))
entertain <- merge(entertain, creatives, by=c("BUA11CD", "date"))
fires$date <- as.Date(fires$date)
entertain <- merge(entertain, fires, by=c("BUA11CD", "date"), all.x=TRUE)
entertain <- merge(entertain, polices, by=c("BUA11CD", "date"), all.x=TRUE)
entertain <- merge(entertain, librarys, by=c("BUA11CD", "date"), all.x=TRUE)


entertain <- merge(entertain, towns, by="BUA11CD", all.x=TRUE)

paneldata <- merge(paneldata, entertain, by=c("BUA11CD", "date"), all.x=TRUE)






















































































###load the packages###

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


###distance code###
citycent <- SpatialPointsDataFrame(gCentroid(citybua, byid=TRUE), citybua@data, match.ID=FALSE)
ucityz <- as.ppp(citycent)
uktownz <- merge(uktowns, popbua, by="BUA11CD")
uktownz <- uktownz[which(uktownz$poplsoa>10000 & uktownz$poplsoa<175000), ]
uktowncent <- SpatialPointsDataFrame(gCentroid(uktownz, byid=TRUE), uktownz@data, match.ID=FALSE)
uktownz <- as.ppp(uktowncent)
uktowncent$distance <- nncross(uktownz, ucityz, what="dist", k=1)

##index of decline##

scotpop01 <- read_csv("scot01pop.csv")
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

indexdecline <- merge(emp, young, by="BUA11CD", all.x=TRUE)
indexdecline <- merge(indexdecline, bus, by="BUA11CD", all.x=TRUE)
indexdecline <- merge(indexdecline, pop, by="BUA11CD", all.x=TRUE)
indexdecline <- merge(indexdecline, degree, by="BUA11CD", all.x=TRUE)

indexdecline<- indexdecline %>% 
  psycho::standardize() 

indexdecline$growth <- indexdecline$empchange + indexdecline$youngchange + indexdecline$buschange + indexdecline$popchange +indexdecline$degreechange
indexdecline <- merge(indexdecline, buanames, by="BUA11CD", all.x=TRUE)
indexdecline <- merge(indexdecline, popbua, by="BUA11CD", all.x=TRUE)
indexdecline <- indexdecline[which(indexdecline$poplsoa<175000&indexdecline$poplsoa>10000), ]
indexdecline <- indexdecline[-c(2,3,4,5,6,8,9)]



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

towns <- towns %>% mutate(declining=ifelse(towns$growth< -0.8, 1, 0))
towns <- towns %>% mutate(stagnant=ifelse(towns$growth> -0.8 & towns$growth<0.25, 1, 0))
towns <- towns %>% mutate(increase=ifelse(towns$growth>0.25, 1, 0))

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



region <- read_csv("region.csv")
region <- region %>% distinct(BUA11CD, .keep_all = TRUE)
towns <- merge(towns, region, by="BUA11CD", all.x=TRUE)
towns$Region <- ifelse(is.na(towns$Region), 
                       'Scotland', towns$Region)


##counts##

hmrc18 <- read_csv("hmrc18.csv")
jobcent18 <- read_csv("jobcent18.csv")
station18 <- read_csv("station18.csv")
gps18 <- read_csv("gps18.csv")
hospital18 <- read_csv("hospital18.csv")
furthered18 <- read_csv("furthered18.csv")
mental18 <- read_csv("mental18.csv")

hmrc17 <- read_csv("hmrc17.csv")
jobcent17 <- read_csv("jobcent17.csv")
station17 <- read_csv("station17.csv")
gps17 <- read_csv("gps17.csv")
hospital17 <- read_csv("hospital17.csv")
furthered17 <- read_csv("furthered17.csv")
mental17 <- read_csv("mental17.csv")

hmrc16 <- read_csv("hmrc16.csv")
jobcent16 <- read_csv("jobcent16.csv")
station16 <- read_csv("station16.csv")
gps16 <- read_csv("gps16.csv")
hospital16 <- read_csv("hospital16.csv")
furthered16 <- read_csv("furthered16.csv")
mental16 <- read_csv("mental16.csv")

hmrc15 <- read_csv("hmrc15.csv")
jobcent15 <- read_csv("jobcent15.csv")
station15 <- read_csv("station15.csv")
gps15 <- read_csv("gps15.csv")
hospital15 <- read_csv("hospital15.csv")
furthered15 <- read_csv("furthered15.csv")
mental15 <- read_csv("mental15.csv")

hmrc14 <- read_csv("hmrc14.csv")
jobcent14 <- read_csv("jobcent14.csv")
station14 <- read_csv("station14.csv")
gps14 <- read_csv("gps14.csv")
hospital14 <- read_csv("hospital14.csv")
furthered14 <- read_csv("furthered14.csv")
mental14 <- read_csv("mental14.csv")

hmrc13 <- read_csv("hmrc13.csv")
jobcent13 <- read_csv("jobcent13.csv")
station13 <- read_csv("station13.csv")
gps13 <- read_csv("gps13.csv")
hospital13 <- read_csv("hospital13.csv")
furthered13 <- read_csv("furthered13.csv")
mental13 <- read_csv("mental13.csv")

hmrc12 <- read_csv("hmrc12.csv")
jobcent12 <- read_csv("jobcent12.csv")
station12 <- read_csv("station12.csv")
gps12 <- read_csv("gps12.csv")
hospital12 <- read_csv("hospital12.csv")
furthered12 <- read_csv("furthered12.csv")
mental12 <- read_csv("mental12.csv")

hmrc11 <- read_csv("hmrc11.csv")
jobcent11 <- read_csv("jobcent11.csv")
station11 <- read_csv("station11.csv")
gps11 <- read_csv("gps11.csv")
hospital11 <- read_csv("hospital11.csv")
furthered11 <- read_csv("furthered11.csv")
mental11 <- read_csv("mental11.csv")

station10 <- read_csv("station10.csv")
gps10 <- read_csv("gps10.csv")
hospital10 <- read_csv("hospital10.csv")
furthered10 <- read_csv("furthered10.csv")
mental10 <- read_csv("mental10.csv")

station09 <- read_csv("station09.csv")
gps09 <- read_csv("gps09.csv")
hospital09 <- read_csv("hospital09.csv")
furthered09 <- read_csv("furthered09.csv")
mental09 <- read_csv("mental09.csv")

station08 <- read_csv("station08.csv")
gps08 <- read_csv("gps08.csv")
hospital08 <- read_csv("hospital08.csv")
furthered08 <- read_csv("furthered08.csv")
mental08 <- read_csv("mental08.csv")


###combine###


station <- merge(station18, station17, by="BUA11CD")
station <- merge(station, station16, by="BUA11CD")
station <- merge(station, station15, by="BUA11CD")
station <- merge(station, station14, by="BUA11CD")
station <- merge(station, station13, by="BUA11CD")
station <- merge(station, station12, by="BUA11CD")
station <- merge(station, station11, by="BUA11CD")
station <- merge(station, station10, by="BUA11CD")
station <- merge(station, station09, by="BUA11CD")
station <- merge(station, station08, by="BUA11CD")

hospital <- merge(hospital18, hospital17, by="BUA11CD")
hospital <- merge(hospital, hospital16, by="BUA11CD")
hospital <- merge(hospital, hospital15, by="BUA11CD")
hospital <- merge(hospital, hospital14, by="BUA11CD")
hospital <- merge(hospital, hospital13, by="BUA11CD")
hospital <- merge(hospital, hospital12, by="BUA11CD")
hospital <- merge(hospital, hospital11, by="BUA11CD")
hospital <- merge(hospital, hospital10, by="BUA11CD")
hospital <- merge(hospital, hospital09, by="BUA11CD")
hospital <- merge(hospital, hospital08, by="BUA11CD")

mental <- merge(mental18, mental17, by="BUA11CD")
mental <- merge(mental, mental16, by="BUA11CD")
mental <- merge(mental, mental15, by="BUA11CD")
mental <- merge(mental, mental14, by="BUA11CD")
mental <- merge(mental, mental13, by="BUA11CD")
mental <- merge(mental, mental12, by="BUA11CD")
mental <- merge(mental, mental11, by="BUA11CD")
mental <- merge(mental, mental10, by="BUA11CD")
mental <- merge(mental, mental09, by="BUA11CD")
mental <- merge(mental, mental08, by="BUA11CD")

furthered <- merge(furthered18, furthered17, by="BUA11CD")
furthered <- merge(furthered, furthered16, by="BUA11CD")
furthered <- merge(furthered, furthered15, by="BUA11CD")
furthered <- merge(furthered, furthered14, by="BUA11CD")
furthered <- merge(furthered, furthered13, by="BUA11CD")
furthered <- merge(furthered, furthered12, by="BUA11CD")
furthered <- merge(furthered, furthered11, by="BUA11CD")
furthered <- merge(furthered, furthered10, by="BUA11CD")
furthered <- merge(furthered, furthered09, by="BUA11CD")
furthered <- merge(furthered, furthered08, by="BUA11CD")

gps <- merge(gps18, gps17, by="BUA11CD")
gps <- merge(gps, gps16, by="BUA11CD")
gps <- merge(gps, gps15, by="BUA11CD")
gps <- merge(gps, gps14, by="BUA11CD")
gps <- merge(gps, gps13, by="BUA11CD")
gps <- merge(gps, gps12, by="BUA11CD")
gps <- merge(gps, gps11, by="BUA11CD")
gps <- merge(gps, gps10, by="BUA11CD")
gps <- merge(gps, gps09, by="BUA11CD")
gps <- merge(gps, gps08, by="BUA11CD")

hmrc <- merge(hmrc18, hmrc17, by="BUA11CD")
hmrc <- merge(hmrc, hmrc16, by="BUA11CD")
hmrc <- merge(hmrc, hmrc15, by="BUA11CD")
hmrc <- merge(hmrc, hmrc14, by="BUA11CD")
hmrc <- merge(hmrc, hmrc13, by="BUA11CD")
hmrc <- merge(hmrc, hmrc12, by="BUA11CD")
hmrc <- merge(hmrc, hmrc11, by="BUA11CD")

jobcent <- merge(jobcent18, jobcent17, by="BUA11CD")
jobcent <- merge(jobcent, jobcent16, by="BUA11CD")
jobcent <- merge(jobcent, jobcent15, by="BUA11CD")
jobcent <- merge(jobcent, jobcent14, by="BUA11CD")
jobcent <- merge(jobcent, jobcent13, by="BUA11CD")
jobcent <- merge(jobcent, jobcent12, by="BUA11CD")
jobcent <- merge(jobcent, jobcent11, by="BUA11CD")

station <- station[-c(2,4,6,8,10,12,14,16,18,20,22)]
hospital <- hospital[-c(2,4,6,8,10,12,14,16,18,20,22)]
gps <- gps[-c(2,4,6,8,10,12,14,16,18,20,22)]
furthered <- furthered[-c(2,4,6,8,10,12,14,16,18,20,22)]
mental <- mental[-c(2,4,6,8,10,12,14,16,18,20,22)]
hmrc <- hmrc[-c(2,4,6,8,10,12,14,16)]
jobcent <- jobcent[-c(2,4,6,8,10,12,14,16)]










publicbine <- merge(mental, hospital, by="BUA11CD")
publicbine <- merge(publicbine, station, by="BUA11CD")
publicbine <- merge(publicbine, furthered, by="BUA11CD")
publicbine <- merge(publicbine, gps, by="BUA11CD")
publicbine <- merge(publicbine, hmrc, by="BUA11CD")
publicbine <- merge(publicbine, jobcent, by="BUA11CD")



##counts##

schools148 <- read_csv("school148.csv")
pbox148 <- read_csv("post148.csv")
wifi148 <- read_csv("wifi148.csv")
hall148 <- read_csv("hall148.csv")
bus148 <- read_csv("bus148.csv")
nursery148 <- read_csv("nursery148.csv")
health148 <- read_csv("health148.csv")

schools13 <- read_csv("school13.csv")
pbox13 <- read_csv("pbox13.csv")
wifi13 <- read_csv("wifi13.csv")
hall13 <- read_csv("hall13.csv")
bus13 <- read_csv("bus13.csv")
nursery13 <- read_csv("nursery13.csv")
health13 <- read_csv("health13.csv")

schools12 <- read_csv("school12.csv")
pbox12 <- read_csv("pbox12.csv")
wifi12 <- read_csv("wifi12.csv")
hall12 <- read_csv("hall12.csv")
bus12 <- read_csv("bus12.csv")
nursery12 <- read_csv("nursery12.csv")
health12 <- read_csv("health12.csv")

schools11 <- read_csv("school11.csv")
pbox11 <- read_csv("pbox11.csv")
wifi11 <- read_csv("wifi11.csv")
hall11 <- read_csv("hall11.csv")
bus11 <- read_csv("bus11.csv")
nursery11 <- read_csv("nursery11.csv")
health11 <- read_csv("health11.csv")

schools10 <- read_csv("school10.csv")
pbox10 <- read_csv("pbox10.csv")
wifi10 <- read_csv("wifi10.csv")
hall10 <- read_csv("hall10.csv")
nursery10 <- read_csv("nursery10.csv")
health10 <- read_csv("health10.csv")

schools09 <- read_csv("school09.csv")
pbox09 <- read_csv("pbox09.csv")
wifi09 <- read_csv("wifi09.csv")
hall09 <- read_csv("hall09.csv")
nursery09 <- read_csv("nursery09.csv")
health09 <- read_csv("health09.csv")

schools08 <- read_csv("school08.csv")
pbox08 <- read_csv("pbox08.csv")
wifi08 <- read_csv("wifi08.csv")
hall08 <- read_csv("hall08.csv")
nursery08 <- read_csv("nursery08.csv")
health08 <- read_csv("health08.csv")



bus <- merge(bus148, bus11, by="BUA11CD")
bus <- merge(bus, bus12, by="BUA11CD")
bus <- merge(bus, bus13, by="BUA11CD")

hall <- merge(hall148, hall08, by="BUA11CD")
hall <- merge(hall, hall09, by="BUA11CD")
hall <- merge(hall, hall10, by="BUA11CD")
hall <- merge(hall, hall11, by="BUA11CD")
hall <- merge(hall, hall12, by="BUA11CD")
hall <- merge(hall, hall13, by="BUA11CD")

pbox <- merge(pbox148, pbox08, by="BUA11CD")
pbox <- merge(pbox, pbox09, by="BUA11CD")
pbox <- merge(pbox, pbox10, by="BUA11CD")
pbox <- merge(pbox, pbox11, by="BUA11CD")
pbox <- merge(pbox, pbox12, by="BUA11CD")
pbox <- merge(pbox, pbox13, by="BUA11CD")

health <- merge(health148, health08, by="BUA11CD")
health <- merge(health, health09, by="BUA11CD")
health <- merge(health, health10, by="BUA11CD")
health <- merge(health, health11, by="BUA11CD")
health <- merge(health, health12, by="BUA11CD")
health <- merge(health, health13, by="BUA11CD")

wifi <- merge(wifi148, wifi08, by="BUA11CD")
wifi <- merge(wifi, wifi09, by="BUA11CD")
wifi <- merge(wifi, wifi10, by="BUA11CD")
wifi <- merge(wifi, wifi11, by="BUA11CD")
wifi <- merge(wifi, wifi12, by="BUA11CD")
wifi <- merge(wifi, wifi13, by="BUA11CD")

nursery <- merge(nursery148, nursery08, by="BUA11CD")
nursery <- merge(nursery, nursery09, by="BUA11CD")
nursery <- merge(nursery, nursery10, by="BUA11CD")
nursery <- merge(nursery, nursery11, by="BUA11CD")
nursery <- merge(nursery, nursery12, by="BUA11CD")
nursery <- merge(nursery, nursery13, by="BUA11CD")

schools <- merge(schools148, schools08, by="BUA11CD")
schools <- merge(schools, schools09, by="BUA11CD")
schools <- merge(schools, schools10, by="BUA11CD")
schools <- merge(schools, schools11, by="BUA11CD")
schools <- merge(schools, schools12, by="BUA11CD")
schools <- merge(schools, schools13, by="BUA11CD")

health <- health[-c(2,8,10,12,14,16,18)]
schools <- schools[-c(2,8,10,12,14,16,18)]
wifi <- wifi[-c(2,8,10,12,14,16,18)]
pbox <- pbox[-c(2,8,10,12,14,16,18)]
bus <- bus[-c(2,8,10,12)]
hall <- hall[-c(2,8,10,12,14,16,18)]
nursery <- nursery[-c(2,8,10,12,14,16,18)]





publicdense <- merge(nursery, wifi, by="BUA11CD")
publicdense <- merge(publicdense, hall, by="BUA11CD")
publicdense <- merge(publicdense, pbox, by="BUA11CD")
publicdense <- merge(publicdense, schools, by="BUA11CD")
publicdense <- merge(publicdense, health, by="BUA11CD")
publicdense <- merge(publicdense, bus, by="BUA11CD")

listdata <- merge(publicdense, publicbine, by="BUA11CD")


##counts##

creative148 <- read_csv("creative148.csv")
sport148 <- read_csv("sport148.csv")
cine148 <- read_csv("cine148.csv")
night148 <- read_csv("night148.csv")
bet148 <- read_csv("bet148.csv")

fire18 <- read_csv("fire18.csv")
fire17 <- read_csv("fire17.csv")
fire16 <- read_csv("fire16.csv")
fire15 <- read_csv("fire15.csv")
fire14 <- read_csv("fire14.csv")
fire13 <- read_csv("fire13.csv")
fire12 <- read_csv("fire12.csv")
fire11 <- read_csv("fire11.csv")
fire10 <- read_csv("fire10.csv")
fire09 <- read_csv("fire09.csv")
fire08 <- read_csv("fire08.csv")

police18 <- read_csv("police18.csv")
police17 <- read_csv("police17.csv")
police16 <- read_csv("police16.csv")
police15 <- read_csv("police15.csv")
police14 <- read_csv("police14.csv")
police13 <- read_csv("police13.csv")
police12 <- read_csv("police12.csv")
police11 <- read_csv("police11.csv")
police10 <- read_csv("police10.csv")
police09 <- read_csv("police09.csv")
police08 <- read_csv("police08.csv")


library18 <- read_csv("library18.csv")
library17 <- read_csv("library17.csv")
library16 <- read_csv("library16.csv")
library15 <- read_csv("library15.csv")
library14 <- read_csv("library14.csv")
library13 <- read_csv("library13.csv")
library12 <- read_csv("library12.csv")
library11 <- read_csv("library11.csv")
library10 <- read_csv("library10.csv")
library09 <- read_csv("library09.csv")
library08 <- read_csv("library08.csv")
creative13 <- read_csv("creative13.csv")
sport13 <- read_csv("sport13.csv")
cine13 <- read_csv("cine13.csv")
night13 <- read_csv("night13.csv")
bet13 <- read_csv("bet13.csv")

creative12 <- read_csv("creative12.csv")
sport12 <- read_csv("sport12.csv")
cine12 <- read_csv("cine12.csv")
night12 <- read_csv("night12.csv")
bet12 <- read_csv("bet12.csv")

creative11 <- read_csv("creative11.csv")
sport11 <- read_csv("sport11.csv")
cine11 <- read_csv("cine11.csv")
night11 <- read_csv("night11.csv")
bet11 <- read_csv("bet11.csv")

creative10 <- read_csv("creative10.csv")
sport10 <- read_csv("sport10.csv")
cine10 <- read_csv("cine10.csv")
night10 <- read_csv("night10.csv")
bet10 <- read_csv("bet10.csv")

creative09 <- read_csv("creative09.csv")
sport09 <- read_csv("sport09.csv")
cine09 <- read_csv("cine09.csv")
night09 <- read_csv("night09.csv")
bet09 <- read_csv("bet09.csv")

creative08 <- read_csv("creative08.csv")
sport08 <- read_csv("sport08.csv")
cine08 <- read_csv("cine08.csv")
night08 <- read_csv("night08.csv")
bet08 <- read_csv("bet08.csv")


###combine##

creative <- merge(creative148, creative08, by="BUA11CD")
creative <- merge(creative, creative09, by="BUA11CD")
creative <- merge(creative, creative10, by="BUA11CD")
creative <- merge(creative, creative11, by="BUA11CD")
creative <- merge(creative, creative12, by="BUA11CD")
creative <- merge(creative, creative13, by="BUA11CD")

cine <- merge(cine148, cine08, by="BUA11CD")
cine <- merge(cine, cine09, by="BUA11CD")
cine <- merge(cine, cine10, by="BUA11CD")
cine <- merge(cine, cine11, by="BUA11CD")
cine <- merge(cine, cine12, by="BUA11CD")
cine <- merge(cine, cine13, by="BUA11CD")

bet <- merge(bet148, bet08, by="BUA11CD")
bet <- merge(bet, bet09, by="BUA11CD")
bet <- merge(bet, bet10, by="BUA11CD")
bet <- merge(bet, bet11, by="BUA11CD")
bet <- merge(bet, bet12, by="BUA11CD")
bet <- merge(bet, bet13, by="BUA11CD")

night <- merge(night148, night08, by="BUA11CD")
night <- merge(night, night09, by="BUA11CD")
night <- merge(night, night10, by="BUA11CD")
night <- merge(night, night11, by="BUA11CD")
night <- merge(night, night12, by="BUA11CD")
night <- merge(night, night13, by="BUA11CD")

sport <- merge(sport148, sport08, by="BUA11CD")
sport <- merge(sport, sport09, by="BUA11CD")
sport <- merge(sport, sport10, by="BUA11CD")
sport <- merge(sport, sport11, by="BUA11CD")
sport <- merge(sport, sport12, by="BUA11CD")
sport <- merge(sport, sport13, by="BUA11CD")

fire <- merge(fire18, fire17, by="BUA11CD")
fire <- merge(fire, fire16, by="BUA11CD")
fire <- merge(fire, fire15, by="BUA11CD")
fire <- merge(fire, fire14, by="BUA11CD")
fire <- merge(fire, fire13, by="BUA11CD")
fire <- merge(fire, fire12, by="BUA11CD")
fire <- merge(fire, fire11, by="BUA11CD")
fire <- merge(fire, fire10, by="BUA11CD")
fire <- merge(fire, fire09, by="BUA11CD")
fire <- merge(fire, fire08, by="BUA11CD")

police <- merge(police18, police17, by="BUA11CD")
police <- merge(police, police16, by="BUA11CD")
police <- merge(police, police15, by="BUA11CD")
police <- merge(police, police14, by="BUA11CD")
police <- merge(police, police13, by="BUA11CD")
police <- merge(police, police12, by="BUA11CD")
police <- merge(police, police11, by="BUA11CD")
police <- merge(police, police10, by="BUA11CD")
police <- merge(police, police09, by="BUA11CD")
police <- merge(police, police08, by="BUA11CD")

library <- merge(library18, library17, by="BUA11CD")
library <- merge(library, library16, by="BUA11CD")
library <- merge(library, library15, by="BUA11CD")
library <- merge(library, library14, by="BUA11CD")
library <- merge(library, library13, by="BUA11CD")
library <- merge(library, library12, by="BUA11CD")
library <- merge(library, library11, by="BUA11CD")
library <- merge(library, library10, by="BUA11CD")
library <- merge(library, library09, by="BUA11CD")
library <- merge(library, library08, by="BUA11CD")

sport <- sport[-c(2,8,10,12,14,16,18)]
creative <- creative[-c(2,8,10,12,14,16,18)]
night <- night[-c(2,8,10,12,14,16,18)]
cine <- cine[-c(2,8,10,12,14,16,18)]
bet <- bet[-c(2,8,10,12,14,16,18)]
police <- police[-c(2,8,10,12,14,16,18)]
fire <- fire[-c(2,8,10,12,14,16,18)]
library <- library[-c(2,8,10,12,14,16,18)]








entertain <- merge(cine, night, by="BUA11CD")
entertain <- merge(entertain, bet, by="BUA11CD")
entertain <- merge(entertain, sport, by="BUA11CD")
entertain <- merge(entertain, creative, by="BUA11CD")
entertain <- merge(entertain, fire, by="BUA11CD")
entertain <- merge(entertain, police, by="BUA11CD")
entertain <- merge(entertain, library, by="BUA11CD")


entertain <- merge(entertain, towns, by="BUA11CD")

listdata <- merge(listdata, entertain, by="BUA11CD")















































































###map tries###

northeast <- towns[which(towns$Region=='North East'), ]
uktownsname <- merge(uktowns, northeast, by="BUA11CD")

uktownsname <- sp.na.omit(uktownsname)

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
    if(growth > 0.25) {
      "green"
    } else if(growth > -0.8) {
      "orange"
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
  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, icon = icons, data=uktownsname, labelOptions = labelOptions(noHide = T))%>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("Improving", "Stagnant", "Declining"), opacity = 1,
    title = "Town Decline")
  


getColorp <- function(uktownsnw) {
  sapply(uktownsname$growth, function(growth) {
    if(growth > 0.25) {
      "green"
    } else if(growth > -0.8) {
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

