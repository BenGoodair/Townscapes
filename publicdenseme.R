###load the packages###
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
library(psycho)
library(reshape2)
library(ggpubr)
library(bife)
library(alpaca)
library(psych)
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

towns <- towns %>% mutate(declining=ifelse(towns$growth< -0.8, 1, 0))
towns <- towns %>% mutate(stagnant=ifelse(towns$growth> -0.8 & towns$growth<0.25, 1, 0))
towns <- towns %>% mutate(increase=ifelse(towns$growth>0.25, 1, 0))

towns$declining <- factor(towns$declining)
towns$stagnant <- factor(towns$stagnant)
towns$increase <- factor(towns$increase)

towns <- towns %>% mutate(declinegroup=ifelse(towns$declining==1, "Declining Towns", ifelse(towns$stagnant==1, "Stagnant Towns", ifelse(towns$increase==1, "Prospering Towns", 0))))

towns$declinegroup <- factor(towns$declinegroup, levels = c("Declining Towns", "Stagnant Towns", "Prospering Towns"))
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

##Density#
engdense <- read_csv("popdensity.csv")
scotdense <- read_csv("popdensityscot.csv")

popdense <- rbind(engdense, scotdense)
towns <-merge(towns, popdense, by="bua11nm")
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


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns",  "bus", "wifi", "nursery", "hall", "schools", "health", "pbox", "towns" )])

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

nursery$"2018-06-01" <- nursery$nursery18/(nursery$poplsoa/1000)
nursery$"2017-06-01" <- nursery$nursery17/(nursery$poplsoa/1000)
nursery$"2016-06-01" <- nursery$nursery16/(nursery$poplsoa/1000)
nursery$"2015-06-01" <- nursery$nursery15/(nursery$poplsoa/1000)
nursery$"2014-09-01" <- nursery$nursery14/(nursery$poplsoa/1000)
nursery$"2013-06-01" <- nursery$nursery13/(nursery$poplsoa/1000)
nursery$"2012-06-01" <- nursery$nursery12/(nursery$poplsoa/1000)
nursery$"2011-06-01" <- nursery$nursery11/(nursery$poplsoa/1000)
nursery$"2010-09-01" <- nursery$nursery10/(nursery$poplsoa/1000)
nursery$"2009-06-01" <- nursery$nursery09/(nursery$poplsoa/1000)
nursery$"2008-03-01" <- nursery$nursery08/(nursery$poplsoa/1000)

schools$"2018-06-01" <- schools$school18/(schools$poplsoa/1000)
schools$"2017-06-01" <- schools$school17/(schools$poplsoa/1000)
schools$"2016-06-01" <- schools$school16/(schools$poplsoa/1000)
schools$"2015-06-01" <- schools$school15/(schools$poplsoa/1000)
schools$"2014-09-01" <- schools$school14/(schools$poplsoa/1000)
schools$"2013-06-01" <- schools$schools13/(schools$poplsoa/1000)
schools$"2012-06-01" <- schools$schools12/(schools$poplsoa/1000)
schools$"2011-06-01" <- schools$schools11/(schools$poplsoa/1000)
schools$"2010-09-01" <- schools$schools10/(schools$poplsoa/1000)
schools$"2009-06-01" <- schools$schools09/(schools$poplsoa/1000)
schools$"2008-03-01" <- schools$schools08/(schools$poplsoa/1000)

health$"2018-06-01" <- health$health18/(health$poplsoa/1000)
health$"2017-06-01" <- health$health17/(health$poplsoa/1000)
health$"2016-06-01" <- health$health16/(health$poplsoa/1000)
health$"2015-06-01" <- health$health15/(health$poplsoa/1000)
health$"2014-09-01" <- health$health14/(health$poplsoa/1000)
health$"2013-06-01" <- health$health13/(health$poplsoa/1000)
health$"2012-06-01" <- health$health12/(health$poplsoa/1000)
health$"2011-06-01" <- health$health11/(health$poplsoa/1000)
health$"2010-09-01" <- health$health10/(health$poplsoa/1000)
health$"2009-06-01" <- health$health09/(health$poplsoa/1000)
health$"2008-03-01" <- health$health08/(health$poplsoa/1000)

bus$"2018-06-01" <- bus$bus18/(bus$poplsoa/1000)
bus$"2017-06-01" <- bus$bus17/(bus$poplsoa/1000)
bus$"2016-06-01" <- bus$bus16/(bus$poplsoa/1000)
bus$"2015-06-01" <- bus$bus15/(bus$poplsoa/1000)
bus$"2014-09-01" <- bus$bus14/(bus$poplsoa/1000)
bus$"2013-06-01" <- bus$bus13/(bus$poplsoa/1000)
bus$"2012-06-01" <- bus$bus12/(bus$poplsoa/1000)
bus$"2011-06-01" <- bus$bus11/(bus$poplsoa/1000)


pbox$"2018-06-01" <- pbox$post18/(pbox$poplsoa/1000)
pbox$"2017-06-01" <- pbox$post17/(pbox$poplsoa/1000)
pbox$"2016-06-01" <- pbox$post16/(pbox$poplsoa/1000)
pbox$"2015-06-01" <- pbox$post15/(pbox$poplsoa/1000)
pbox$"2014-09-01" <- pbox$post14/(pbox$poplsoa/1000)
pbox$"2013-06-01" <- pbox$pbox13/(pbox$poplsoa/1000)
pbox$"2012-06-01" <- pbox$pbox12/(pbox$poplsoa/1000)
pbox$"2011-06-01" <- pbox$pbox11/(pbox$poplsoa/1000)
pbox$"2010-09-01" <- pbox$pbox10/(pbox$poplsoa/1000)
pbox$"2009-06-01" <- pbox$pbox09/(pbox$poplsoa/1000)
pbox$"2008-03-01" <- pbox$pbox08/(pbox$poplsoa/1000)

wifi$"2018-06-01" <- wifi$wifi18/(wifi$poplsoa/1000)
wifi$"2017-06-01" <- wifi$wifi17/(wifi$poplsoa/1000)
wifi$"2016-06-01" <- wifi$wifi16/(wifi$poplsoa/1000)
wifi$"2015-06-01" <- wifi$wifi15/(wifi$poplsoa/1000)
wifi$"2014-09-01" <- wifi$wifi14/(wifi$poplsoa/1000)
wifi$"2013-06-01" <- wifi$wifi13/(wifi$poplsoa/1000)
wifi$"2012-06-01" <- wifi$wifi12/(wifi$poplsoa/1000)
wifi$"2011-06-01" <- wifi$wifi11/(wifi$poplsoa/1000)
wifi$"2010-09-01" <- wifi$wifi10/(wifi$poplsoa/1000)
wifi$"2009-06-01" <- wifi$wifi09/(wifi$poplsoa/1000)
wifi$"2008-03-01" <- wifi$wifi08/(wifi$poplsoa/1000)

hall$"2018-06-01" <- hall$hall18/(hall$poplsoa/1000)
hall$"2017-06-01" <- hall$hall17/(hall$poplsoa/1000)
hall$"2016-06-01" <- hall$hall16/(hall$poplsoa/1000)
hall$"2015-06-01" <- hall$hall15/(hall$poplsoa/1000)
hall$"2014-09-01" <- hall$hall14/(hall$poplsoa/1000)
hall$"2013-06-01" <- hall$halls13/(hall$poplsoa/1000)
hall$"2012-06-01" <- hall$halls12/(hall$poplsoa/1000)
hall$"2011-06-01" <- hall$halls11/(hall$poplsoa/1000)
hall$"2010-09-01" <- hall$halls10/(hall$poplsoa/1000)
hall$"2009-06-01" <- hall$halls09/(hall$poplsoa/1000)
hall$"2008-03-01" <- hall$halls08/(hall$poplsoa/1000)


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



nurseryz <-  nursery %>%
  psycho::standardize()
schoolsz <-  schools %>%
  psycho::standardize()
healthz <-  health %>%
  psycho::standardize()
busz <-  bus %>%
  psycho::standardize()
pboxz <-  pbox %>%
  psycho::standardize()
wifiz <-  wifi %>%
  psycho::standardize()
hallz <-  hall %>%
  psycho::standardize()

schoolszz <- schoolsz[-c(1,2)] 

schoolszz <- reshape(schoolszz, idvar = "BUA11CD", ids = schoolsz$BUA11CD,
                    times = names(schoolszz), timevar = "date",
                    varying = list(names(schoolszz)),v.names="schoolszz", new.row.names = 1:((dim(schoolszz)[2])*(dim(schoolszz)[1])),direction = "long")

schoolszz$date <- as.Date(schoolszz$date)
healthzz <- healthz[-c(1,2)] 

healthzz <- reshape(healthzz, idvar = "BUA11CD", ids = healthz$BUA11CD,
                   times = names(healthzz), timevar = "date",
                   varying = list(names(healthzz)),v.names="healthzz", new.row.names = 1:((dim(healthzz)[2])*(dim(healthzz)[1])),direction = "long")

healthzz$date <- as.Date(healthzz$date)
pboxzz <- pboxz[-c(1,2)] 

pboxzz <- reshape(pboxzz, idvar = "BUA11CD", ids = pboxz$BUA11CD,
                 times = names(pboxzz), timevar = "date",
                 varying = list(names(pboxzz)),v.names="pboxzz", new.row.names = 1:((dim(pboxzz)[2])*(dim(pboxzz)[1])),direction = "long")

pboxzz$date <- as.Date(pboxzz$date)
buszz <- busz[-c(1,2)] 

buszz <- reshape(buszz, idvar = "BUA11CD", ids = busz$BUA11CD,
                times = names(buszz), timevar = "date",
                varying = list(names(buszz)),v.names="buszz", new.row.names = 1:((dim(buszz)[2])*(dim(buszz)[1])),direction = "long")

buszz$date <- as.Date(buszz$date)
hallzz <- hallz[-c(1,2)] 

hallzz <- reshape(hallzz, idvar = "BUA11CD", ids = hallz$BUA11CD,
                 times = names(hallzz), timevar = "date",
                 varying = list(names(hallzz)),v.names="hallzz", new.row.names = 1:((dim(hallzz)[2])*(dim(hallzz)[1])),direction = "long")

hallzz$date <- as.Date(hallzz$date)
wifizz <- wifiz[-c(1,2)] 

wifizz <- reshape(wifizz, idvar = "BUA11CD", ids = wifiz$BUA11CD,
                 times = names(wifizz), timevar = "date",
                 varying = list(names(wifizz)),v.names="wifizz", new.row.names = 1:((dim(wifizz)[2])*(dim(wifizz)[1])),direction = "long")

wifizz$date <- as.Date(wifizz$date)
nurseryzz <- nurseryz[-c(1,2)] 

nurseryzz <- reshape(nurseryzz, idvar = "BUA11CD", ids = nurseryz$BUA11CD,
                    times = names(nurseryzz), timevar = "date",
                    varying = list(names(nurseryzz)),v.names="nurseryzz", new.row.names = 1:((dim(nurseryzz)[2])*(dim(nurseryzz)[1])),direction = "long")

nurseryzz$date <- as.Date(nurseryzz$date)


publicdensez <- merge(nurseryzz, wifizz, by=c("BUA11CD", "date"))
publicdensez <- merge(publicdensez, hallzz, by=c("BUA11CD", "date"))
publicdensez <- merge(publicdensez, pboxzz, by=c("BUA11CD", "date"))
publicdensez <- merge(publicdensez, schoolszz, by=c("BUA11CD", "date"))
publicdensez <- merge(publicdensez, healthzz, by=c("BUA11CD", "date"))
publicdensez <- merge(publicdensez, buszz, by=c("BUA11CD", "date"), all.x=TRUE)

##jsa##

scotjsa <- read_csv("scotjsa.csv")
jsa08 <- read_csv("jsa08.csv")
jsa14 <- read_csv("jsa14.csv")

scotjsa <- scotjsa[complete.cases(scotjsa$BUA11CD), ]
jsa08 <- jsa08[complete.cases(jsa08$BUA11CD), ]
jsa14 <- jsa14[complete.cases(jsa14$BUA11CD), ]

scotjsa <-  aggregate(scotjsa[-1], scotjsa["BUA11CD"], sum)
jsa08 <-  aggregate(jsa08[-1], jsa08["BUA11CD"], sum)
jsa14 <-  aggregate(jsa14[-1], jsa14["BUA11CD"], sum)


scotjsa <- merge(scotjsa, popbua, by="BUA11CD")
jsa08 <- merge(jsa08, popbua, by="BUA11CD")
jsa14 <- merge(jsa14, popbua, by="BUA11CD")


scotjsa <- scotjsa[which(scotjsa$poplsoa<175000&scotjsa$poplsoa>10000), ]
jsa08 <- jsa08[which(jsa08$poplsoa<175000&jsa08$poplsoa>10000), ]
jsa14 <- jsa14[which(jsa14$poplsoa<175000&jsa14$poplsoa>10000), ]

scotjsas <- scotjsa[-c(1,13)]
jsa08s <- jsa08[-c(1,8)]
jsa14s <- jsa14[-c(1,7)]


jsa08s <- reshape(jsa08s, idvar = "BUA11CD", ids = jsa08$BUA11CD,
                  times = names(jsa08s), timevar = "date",
                  varying = list(names(jsa08s)),v.names="jsa", new.row.names = 1:((dim(jsa08s)[2])*(dim(jsa08s)[1])),direction = "long")


jsa14s <- reshape(jsa14s, idvar = "BUA11CD", ids = jsa14$BUA11CD,
                  times = names(jsa14s), timevar = "date",
                  varying = list(names(jsa14s)),v.names="jsa", new.row.names = 1:((dim(jsa14s)[2])*(dim(jsa14s)[1])),direction = "long")

scotjsas <- reshape(scotjsas, idvar = "BUA11CD", ids = scotjsa$BUA11CD,
                    times = names(scotjsas), timevar = "date",
                    varying = list(names(scotjsas)),v.names="jsa", new.row.names = 1:((dim(scotjsas)[2])*(dim(scotjsas)[1])),direction = "long")

jobseekers <- merge(jsa08s, jsa14s, by=c("BUA11CD", "date", "jsa"), all=TRUE)
jobseekers <- merge(jobseekers, scotjsas, by=c("BUA11CD", "date", "jsa"), all=TRUE)


publicdense <- merge(publicdense, jobseekers, by=c("BUA11CD", "date"))
publicdense <- merge(publicdense, towns, by="BUA11CD", all.x=TRUE)


try <- ggplot(data = publicdense, aes(x = date, y = healths, group = BUA11CD))
 healthplot <- try + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                                                      geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "per 1000 population", title="Average Health Providers")+ facet_grid(. ~ declinegroup)
try <- ggplot(data = publicdense, aes(x = date, y = schoolss, group = BUA11CD))
schoolplot <-try + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                                                     geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "per 1000 population", title="Average Schools")+ facet_grid(. ~ declinegroup)

try <- ggplot(data = publicdense, aes(x = date, y = nurserys, group = BUA11CD))
nurseryplot <-try + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                                                  geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "per 1000 population", title="Average Nurseries")+ facet_grid(. ~ declinegroup)

try <- ggplot(data = publicdense, aes(x = date, y = pboxs, group = BUA11CD))
pboxplot <-try + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                                                   geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "per 1000 population", title="Average Post Box")+ facet_grid(. ~ declinegroup)

try <- ggplot(data = publicdense, aes(x = date, y = buss, group = BUA11CD))
busplot <-try + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                                                   geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "per 1000 population", title="Average Bus stops")+ facet_grid(. ~ declinegroup)

#try <- ggplot(data = publicdense, aes(x = date, y = halls, group = BUA11CD))
#hallplot <-try + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
#                                                                                                                   geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "per 1000 population", title="Average Halls")+ facet_grid(. ~ declinegroup)
#
#try <- ggplot(data = publicdense, aes(x = date, y = wifis, group = BUA11CD))
#wifiplot <-try + stat_smooth(aes(group = 1)) + stat_summary(aes(group = 1),
#                                                                                                                   geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "per 1000 population", title="Average Wifi")+ facet_grid(. ~ declinegroup)
#
#ggarrange(wifiplot, hallplot, busplot, pboxplot, nurseryplot, schoolplot, healthplot, ncol=4, nrow=2)
#
publicdensez <- merge(publicdensez, towns, by="BUA11CD", all.x=TRUE)
publicdensez <- publicdensez[complete.cases(publicdensez$declinegroup), ]




try <- ggplot(data = publicdensez, aes(x = date, y = healthzz, group = BUA11CD))
healthplot <- try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                               geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Health per person", title="Health Providers per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
try <- ggplot(data = publicdensez, aes(x = date, y = schoolszz, group = BUA11CD))
schoolplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                              geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Schools per person", title="Schools per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = nurseryzz, group = BUA11CD))
nurseryplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                               geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Nurseries per person", title="Children's Nurseries per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = pboxzz, group = BUA11CD))
pboxplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                            geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Post Boxes per person", title="Post Box per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = buszz, group = BUA11CD))
busplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                           geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Bus stops per person", title="Bus stops per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = hallzz, group = BUA11CD))
hallplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                            geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Community Centres/ Halls per person", title="Community Centres/ Halls per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = wifizz, group = BUA11CD))
wifiplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                            geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Wifi Hotspots per person", title="Wifi Hotspots per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

ggarrange(wifiplot, hallplot, busplot, pboxplot, nurseryplot, schoolplot, healthplot, ncol=4, nrow=2)


try <- ggplot(data = publicdensez, aes(x = date, y = healthzz, group = BUA11CD))
healthplot <- try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                               geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Health Providers per person", title="Health Providers per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
try <- ggplot(data = publicdensez, aes(x = date, y = schoolszz, group = BUA11CD))
schoolplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                              geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Schools per person", title="Schools per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = nurseryzz, group = BUA11CD))
nurseryplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                               geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Nurseries per person", title="Nurseries per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = pboxzz, group = BUA11CD))
pboxplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                            geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Post Boxes per person", title="Post Box per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = buszz, group = BUA11CD))
busplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                           geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Bus stops per person", title="Bus stops per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = hallzz, group = BUA11CD))
hallplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                            geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Community Centres/ Halls per person", title="Community Centres/ Halls per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = wifizz, group = BUA11CD))
wifiplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                            geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Wifi Hotspots per person", title="Wifi Hotspots per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

ggarrange(wifiplot, hallplot, busplot, pboxplot, nurseryplot, schoolplot, healthplot, ncol=4, nrow=2)

publicdensez$country <- factor(publicdensez$country)


try <- ggplot(data = publicdensez, aes(x = date, y = healthzz, group = BUA11CD))
healthplot <- try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                          geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Health Providers per person", title="Average Health Providers")+ facet_grid(. ~ country)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
try <- ggplot(data = publicdensez, aes(x = date, y = schoolszz, group = BUA11CD))
schoolplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                         geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Schools per person", title="Average Schools")+ facet_grid(. ~ country)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = nurseryzz, group = BUA11CD))
nurseryplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                          geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Nurseries per person", title="Average Nurseries")+ facet_grid(. ~ country)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = pboxzz, group = BUA11CD))
pboxplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Post Boxes per person", title="Average Post Box")+ facet_grid(. ~ country)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = buszz, group = BUA11CD))
busplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                      geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Bus stops per person", title="Average Bus stops")+ facet_grid(. ~ country)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = hallzz, group = BUA11CD))
hallplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Community Centres/ Halls per person", title="Average Community Centres/ Halls")+ facet_grid(. ~ country)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = publicdensez, aes(x = date, y = wifizz, group = BUA11CD))
wifiplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Wifi Hotspots per person", title="Average Wifi Hotspots")+ facet_grid(. ~ country)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))


###18 analysis###

public18 <- publicdense[which(publicdense$date=='2018-06-01'), ]

school <- ggplot(data=public18, aes(x=reorder(BUA11CD,-schoolss), y=schoolss)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Schools per Person") +
  xlab("Town") + ylab("Schools (per 1000 population) ") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

nurserys <- ggplot(data=public18, aes(x=reorder(BUA11CD,-nurserys), y=nurserys, label=name)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Number of Children's Nurseries per Capita") +
  xlab("Town") + ylab("Nurseries (per 1000 population) ") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  geom_label_repel(data         = subset(public18, nurserys > 0.61),
                   nudge_x       = 4 - subset(public18, nurserys > 0.61)$nurserys,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")+
  geom_label_repel(data         = subset(public18, nurserys < 0.05),
                   nudge_x       = 4 - subset(public18, nurserys < 0.05)$nurserys,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")

halls <- ggplot(data=public18, aes(x=reorder(BUA11CD,-halls), y=halls, label=name) ) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Number of Community Centres/ Halls per Capita") +
  xlab("Town") + ylab("Halls (per 1000 population) ") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  geom_label_repel(data         = subset(public18, halls > 0.64),
                   nudge_x       = 4 - subset(public18, halls > 0.64)$halls,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")+
  geom_label_repel(data         = subset(public18, halls < 0.04),
                   nudge_x       = 4 - subset(public18, halls < 0.04)$halls,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")

pboxs <- ggplot(data=public18, aes(x=reorder(BUA11CD,-pboxs), y=pboxs, label=name)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Number of Post Boxes per Capita") +
  xlab("Town") + ylab("Post boxes (per 1000 population) ") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  geom_label_repel(data         = subset(public18, pboxs > 2.53),
                   nudge_x       = 4 - subset(public18, pboxs > 2.53)$pboxs,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")+
  geom_label_repel(data         = subset(public18, pboxs < 0.55),
                   nudge_x       = 4 - subset(public18, pboxs < 0.55)$pboxs,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")

wifis <- ggplot(data=public18, aes(x=reorder(BUA11CD,-wifis), y=wifis, label=name)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Wifi Hotspots per Person") +
  xlab("Town") + ylab("Wifi hotspots (per 1000 population) ") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  geom_label_repel(data         = subset(public18, wifis > 1.5),
                   nudge_x       = 4 - subset(public18, wifis > 1.5)$wifis,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 5,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")




buss <- ggplot(data=public18, aes(x=reorder(BUA11CD,-buss), y=buss, label=name)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Number of Bus Stops per Capita") +
  xlab("Town") + ylab("Bus stops (per 1000 population) ") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  geom_label_repel(data         = subset(public18, buss > 10.5),
                   
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")+
  geom_label_repel(data         = subset(public18, buss < 2),
                   
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")

health <- ggplot(data=public18, aes(x=reorder(BUA11CD,-healths), y=healths, label=name)) +
  geom_col()+
  theme_minimal()+ ggtitle("Number of Health Providers per Capita") +
  xlab("Town") + ylab("Healths providers (per 1000 population) ") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  geom_label_repel(data         = subset(public18, healths > 3.8),
                    nudge_x       = 4 - subset(public18, healths > 3.8)$healths,
                    size          = 2,
                    box.padding   = 1,
                    point.padding = 0.5,
                    force         = 3,
                    segment.size  = 0.2,
                    segment.color = "grey50",
                    direction     = "y")+
  geom_label_repel(data         = subset(public18, healths < 0.7),
                   nudge_x       = 4 - subset(public18, healths < 0.7)$healths,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 3,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")





heal <- ggplot(data = public18, aes(x = distance, y = healths, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Health providers (per 1000 population)", title="Health Providers and Distance from City")
bu <- ggplot(data = public18, aes(x = distance, y = buss, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Distance from nearest City (km)", y= "Bus Stops (per 1000 population)", title="Bus Stops and Distance from City")
wif <- ggplot(data = public18, aes(x = distance, y = wifis, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Distance from nearest City (km)", y= "Wifi Hotspots (per 1000 population)", title="Wifi Hotspots and Distance from City")
nurs <- ggplot(data = public18, aes(x = distance, y = nurserys, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Nurseries (per 1000 population)", title="Nurseries and Distance from City")
po <- ggplot(data = public18, aes(x = distance, y = pboxs, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Post Boxes (per 1000 population)", title="Post Boxes and Distance from City")
scho <- ggplot(data = public18, aes(x = distance, y = schoolss, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Schools (per 1000 population)", title="Schools and Distance from City")
hal <- ggplot(data = public18, aes(x = distance, y = halls, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Community Halls (per 1000 population)", title="Community Halls and Distance from City")
ggarrange(heal, bu, wif, nurs, scho, po, hal, ncol=2, nrow=4)

heal <- ggplot(data = public18, aes(x = growth, y = healths, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Health providers (per 1000 population)", title="Health Providers and Economic Growth") + scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
bu <- ggplot(data = public18, aes(x = growth, y = buss, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Growth Index", y= "Bus Stops (per 1000 population)", title="Bus Stops and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
wif <- ggplot(data = public18, aes(x = growth, y = wifis, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Growth Index", y= "Wifi Hotspots (per 1000 population)", title="Wifi Hotspots and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
nurs <- ggplot(data = public18, aes(x = growth, y = nurserys, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Nurseries (per 1000 population)", title="Nurseries and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
po <- ggplot(data = public18, aes(x = growth, y = pboxs, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Post Boxes (per 1000 population)", title="Post Boxes and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
scho <- ggplot(data = public18, aes(x = growth, y = schoolss, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Schools (per 1000 population)", title="Schools and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
hal <- ggplot(data = public18, aes(x = growth, y = halls, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Community Halls (per 1000 population)", title="Community Halls and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
ggarrange(heal, bu, wif, nurs, scho, po, hal, ncol=2, nrow=4)

public18$poplsoa <- public18$poplsoa/1000
public18$distance <- public18$distance/10



h <- lm(schoolss~growth+distance+poplsoa+jsa, data=public18)


f <- lm(healths~growth+distance+poplsoa+jsa, data=public18)


d <- lm(wifis~growth+distance+poplsoa+jsa, data=public18)


s <- lm(buss~growth+distance+poplsoa+jsa, data=public18)


a <- lm(pboxs~growth+distance+poplsoa+jsa, data=public18)


q <- lm(nurserys~growth+distance+poplsoa+jsa, data=public18)

t <- lm(halls~growth+distance+poplsoa+jsa, data=public18)


stargazer(q,t, a, s, d, f, h, style= "apsr", title="Results", align=TRUE, type="text")

###losing maps###

nurserychange <- ggplot(data=publicchange, aes(x=reorder(BUA11CD,-nurserychange), y=nurserychange)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Changes in number of Children's Nurseries") +
  xlab("Town") + ylab("Nursery change") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

healthchange <- ggplot(data=publicchange, aes(x=reorder(BUA11CD,-healthchange), y=healthchange)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Changes in number of Health Providers") +
  xlab("Town") + ylab("Health provider change") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

buschange <- ggplot(data=publicchange, aes(x=reorder(BUA11CD,-buschange), y=buschange)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Changes in number of Bus Stops") +
  xlab("Town") + ylab("Bus stop change") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

schoolchange <- ggplot(data=publicchange, aes(x=reorder(BUA11CD,-schoolchange), y=schoolchange)) +
  geom_col()+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()+ ggtitle("Changes in number of Schools") +
  xlab("Town") + ylab("School change") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())


  



publicchange <- publicchange %>% mutate(morenursery=ifelse(publicchange$nurserychange>0, 1,0))
publicchange <- publicchange %>% mutate(morehealth=ifelse(publicchange$healthchange>0, 1,0))
publicchange <- publicchange %>% mutate(moreschool=ifelse(publicchange$schoolchange>0, 1,0))
publicchange <- publicchange %>% mutate(morebus=ifelse(publicchange$buschange>0, 1,0))

publicchangez <- merge(uktowns, publicchange,  by="BUA11CD")
publicchangez <- sp.na.omit(publicchangez)
publicchangez <- SpatialPointsDataFrame(gCentroid(publicchangez, byid=TRUE), publicchangez@data, match.ID=FALSE)



spplot(publicchangez, "morenursery")
spplot(publicchangez, "moreschool")
spplot(publicchangez, "morebus")
spplot(publicchangez, "morehealth")


distecon <- ggplot(data = public18, aes(x = growth, y = distance, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Distance from nearest City (km)", title="Economic Growth and Distance from a City") + scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

econdense <- ggplot(data = towns, aes(x = growth, y = popdensity, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Population Density (people per hectare, 2011)", title="Economic Growth and Population Density") + scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
distancedense <- ggplot(data = towns, aes(x = distance, y = popdensity, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Population Density (people per hectare, 2011)", title="Distance from Nearest City and Population Density") 


growthplot <- ggplot(towns, aes(x = growth, y = country, label = name)) +
  geom_point(alpha = .10) + 
  labs(x = "Growth Index", y= "Country", title="Economic Growth") +scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))+
  geom_label_repel(data         = subset(towns, growth > 3),
                   nudge_x       = 4 - subset(towns, growth > 3)$growth,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 0.5,
                   force         = 5,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y") +
  geom_label_repel(data         = subset(towns, growth < -3),
                   nudge_x       = -3 - subset(towns, growth < -3)$growth,
                   size          = 2,
                   box.padding   = 1,
                   point.padding = 1,
                   force         = 5,
                   segment.size  = 0.2,
                   segment.color = "grey50",
                   direction     = "y")


publicdense <- publicdense[complete.cases(publicdense),]
try <- ggplot(data = publicdense, aes(x = date, y = pboxs, group = BUA11CD))
pboxplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Post Boxes per person", title="Post Box per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)


