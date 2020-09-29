###load the packages###
library(GGally)
library(ggstance)
library(jtools)
library(lme4)
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
library(pscyho)
library(reshape2)
library(ggpubr)
library(bife)
library(psych)
library(grid)
library(lubridate)
library(snakecase)
library(sjPlot)
library(mapview)

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
collar <- (brewer.pal( 10, "RdBu"))
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
theme_set(theme_minimal())

###downloadmastersets###

masterpanel <- read_csv("paneldata2406.csv")
masterpanel <- masterpanel[-c(34)]
##ffs###
region <- read_csv("region2.csv")
region <- region %>% distinct(BUA11CD, .keep_all = TRUE)
masterpanel <- merge(masterpanel, region, by = "BUA11CD", all.x=TRUE)
masterpanel$Region <- ifelse(is.na(masterpanel$Region), 
                             'Scotland', masterpanel$Region)

##index of decline####

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
###fffsss###
engdegree01try2 <- read_csv("eng01degree2ff.csv")
engdegree01try2 <- engdegree01try2[-c(1,2,4)]
engdegree01try2 <- aggregate(engdegree01try2[-1], engdegree01try2["BUA11CD"], sum)
engdegree01 <- engdegree01[-c(3)]
engdegree01 <- merge(engdegree01, engdegree01try2, by="BUA11CD")

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


masterpanel <- masterpanel[-c(29)]

masterpanel <- merge(masterpanel, indexdecline, by="BUA11CD")


masterpanel <- masterpanel[-c(2)]
masterpanel$date <- as.Date(masterpanel$date)

buasize <- buas[-c(1,3,4,5,6,7,9)]
buasize$size <- (buasize$st_areasha/1000)/1000
scotbuasize <- read_csv("scotbuasize.csv")
scotbuasize <- scotbuasize[-c(1)]
buasize <- buasize[-c(2)]
buasize <- as.data.frame(buasize)
buasize <- rbind(buasize, scotbuasize)

masterpanel <- masterpanel[-c(32)]
masterpanel <- merge(masterpanel, buasize, by="BUA11CD")
data18 <-masterpanel[which(masterpanel$date=='2018-06-01'),] 

masterlist <- read_csv("listdata2406.csv")
masterlist <- masterlist[-c(1)]
data18$busperpop <- data18$buss/(data18$poplsoa/1000)
data18$healthsper <- data18$healths/(data18$poplsoa/1000)
data18$pboxsper <- data18$pboxs/(data18$poplsoa/1000)
data18$nurseryper <- data18$nurserys/(data18$poplsoa/1000)
data18$schoolper <- data18$schoolss/(data18$poplsoa/1000)
data18$wifiper <- data18$wifis/(data18$poplsoa/1000)
data18$hallper <- data18$halls/(data18$poplsoa/1000)
data18$gpsper <- data18$gpss/(data18$poplsoa/1000)

masterpanel$busperpop <- masterpanel$buss/(masterpanel$poplsoa/1000)
masterpanel$healthsper <- masterpanel$healths/(masterpanel$poplsoa/1000)
masterpanel$pboxsper <- masterpanel$pboxs/(masterpanel$poplsoa/1000)
masterpanel$nurseryper <- masterpanel$nurserys/(masterpanel$poplsoa/1000)
masterpanel$schoolper <- masterpanel$schoolss/(masterpanel$poplsoa/1000)
masterpanel$wifiper <- masterpanel$wifis/(masterpanel$poplsoa/1000)
masterpanel$hallper <- masterpanel$halls/(masterpanel$poplsoa/1000)
masterpanel$gpsper <- masterpanel$gpss/(masterpanel$poplsoa/1000)

###public service index####

publicindex <- subset(data18, select = c(1,35,36,37,38,39,41,42))
publicindex<- publicindex %>% 
  psycho::standardize() 
publicindex$services <- publicindex$busperpop+publicindex$pboxsper+publicindex$healthsper+publicindex$nurseryper+publicindex$schoolper+publicindex$hallper+publicindex$gpsper
#publicindex <- merge(publicindex, buanames, by="BUA11CD")
publicindex <- publicindex[-c(2,3,4,5,6,7,8)]

data18 <- merge(data18, publicindex, by="BUA11CD")

###publicchangeindex###
publicchangeindex <- subset(masterlist, select = c(1,6,10,28,32,39,43,50,54,61,65,72,73, 120,127, 217,226))
publicchangeindex <- publicchangeindex[-c(6,7)]
publicchangeindex$nuserychange <- publicchangeindex$nursery18-publicchangeindex$nursery11
publicchangeindex$gpchange <- publicchangeindex$GPs18-publicchangeindex$GPs11
publicchangeindex$buschange <- publicchangeindex$bus18-publicchangeindex$bus11
publicchangeindex$schoolchange <- publicchangeindex$school18-publicchangeindex$schools11
publicchangeindex$policechange <- publicchangeindex$police18-publicchangeindex$police11
publicchangeindex$healthchange <- publicchangeindex$health18-publicchangeindex$health11
publicchangeindex$hallchange <- publicchangeindex$hall18-publicchangeindex$halls11

publicchangeindex <- publicchangeindex[-c(2:15)]

publicchangeindex<- publicchangeindex %>% 
  psycho::standardize() 

publicchangeindex$publicchange <- publicchangeindex$nuserychange+publicchangeindex$gpchange+publicchangeindex$buschange+publicchangeindex$schoolchange+publicchangeindex$policechange+publicchangeindex$healthchange+publicchangeindex$hallchange

publicchangeindex <- publicchangeindex[-c(2:8)]

data18 <- merge(data18, publicchangeindex, by="BUA11CD")


###controls###
scotdep <- read_csv("scothousedep.csv")
engdep <- read_csv("enghousedep.csv")

scotdep$housedep <- 100-((scotdep$nondep/scotdep$housetot)*100)
engdep$housedep <- 100-((engdep$nondep/engdep$housetot)*100)

scotdep <- scotdep[-c(2,3)]
engdep <- engdep[-c(2,3)]

housedep <- rbind(scotdep, engdep)

data18 <- merge(data18, housedep, by="bua11nm")

scotsevdep <- read_csv("severedepscot.csv")
engsevdep <- read_csv("severedepengwa.csv")

severedep <- rbind(engsevdep, scotsevdep)

data18 <- merge(data18, severedep, by="bua11nm")


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
data18 <- merge(data18, jobdense, by="BUA11CD")

data18 <- data18 %>% mutate(residential=ifelse(data18$jobdensity<0.5, 1, 0))
data18 <- data18 %>% mutate(versatile=ifelse(data18$jobdensity>0.5 & data18$jobdensity<0.7, 1, 0))
data18 <- data18 %>% mutate(working=ifelse(data18$jobdensity>0.7, 1, 0))

data18$residential <- factor(data18$residential)
data18$versatile <- factor(data18$versatile)
data18$working <- factor(data18$working)

data18 <- data18 %>% mutate(workinggroup=ifelse(data18$working==1, "Working", ifelse(data18$versatile==1, "Partially Residential", ifelse(data18$residential==1, "Residential", 0))))

data18$workinggroup <- factor(data18$workinggroup, levels = c("Working", "Partially Residential", "Residential"))




##Westmids####

nuts1 <- shapefile("nuts1")
nuts1 <- nuts1[which(nuts1$objectid!=12), ]

nuts1 <- spTransform(nuts1, CRS("+init=epsg:27700"))
nuts1s <- as.data.frame(nuts1)
nuts1s <- nuts1s %>% mutate(WM = ifelse(nuts1s$objectid==5, 1, 0))
nuts1ss <- merge(nuts1, nuts1s, by="objectid")
# nuts1ss$NE <- factor(nuts1ss$NE)

map <- spplot(nuts1ss, "WM", col.regions = c("white", "navyblue"), colorkey=FALSE)
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "wmmap.png", plot=map, width=6, height=8, dpi=600)


westmids <- data18[which(data18$Region=='West Midlands'), ]
westmidspanel <- masterpanel[which(masterpanel$Region=='West Midlands'), ]
westmidslist <- masterlist[which(masterlist$Region=='West Midlands'), ]
####wmdepgraph####
lsoabualookup <- read_csv("lsoabualookup.csv")
lsoapop <- read_csv("lsoapop18.csv")

imd <- read_csv("imd2019.csv")
imd <- imd[which(imd$Measurement=="Rank"), ]
imd <- imd[which(imd$`Indices of Deprivation`=="a. Index of Multiple Deprivation (IMD)"), ]
names(imd)[names(imd)=="FeatureCode"] <- "LSOA11CD"
imdbua <- merge(imd, lsoabualookup, by="LSOA11CD", all.x=TRUE)
imdbua <- merge(imdbua, lsoapop, by="LSOA11CD", all.x=TRUE)

imdpop <- imdbua[-c(1:6, 8,9,10,11)]
imdpop <- aggregate(imdpop[-c(1)], imdpop["BUA11CD"], sum)
names(imdpop)[names(imdpop)=="mid2018pop"] <- "lsoa18popbestfit"

imdbua <- merge(imdbua, imdpop, by="BUA11CD", all.x=TRUE)
imdbua$imdrank2019 <- imdbua$Value*(imdbua$mid2018pop/imdbua$lsoa18popbestfit)
imdbua <- imdbua[-c(2:13)]
imdbua <- aggregate(imdbua[-c(1)], imdbua["BUA11CD"], sum)
imdbua <- imdbua[order(imdbua$imdrank2019),]
imdbua$imdrank2019 <- 1:nrow(imdbua)

westmids <- merge(westmids, imdbua, by="BUA11CD", all.x=TRUE)
westmids <- westmids[order(westmids$imdrank2019),]
westmids$imdrank2019 <- 1:nrow(westmids)


imd2010 <- read_csv("imd2010.csv")
names(imd2010)[names(imd2010)=="LSOA CODE"] <- "LSOA11CD"

imd2010bua <- merge(imd2010, lsoabualookup, by="LSOA11CD", all.x=TRUE)
imd2010bua <- merge(imd2010bua, lsoapop, by="LSOA11CD", all.x=TRUE)

imd2010pop <- imd2010bua[-c(1:7, 9,10,11,12)]
imd2010pop <- aggregate(imd2010pop[-c(1)], imd2010pop["BUA11CD"], sum)
names(imd2010pop)[names(imd2010pop)=="mid2018pop"] <- "lsoa18popbestfit"

imd2010bua <- merge(imd2010bua, imdpop, by="BUA11CD", all.x=TRUE)
imd2010bua$imdrank2010 <- imd2010bua$`RANK OF IMD SCORE (where 1 is most deprived)`*(imd2010bua$mid2018pop/imd2010bua$lsoa18popbestfit)
imd2010bua <- imd2010bua[-c(2:14)]
imd2010bua <- aggregate(imd2010bua[-c(1)], imd2010bua["BUA11CD"], sum)
imd2010bua <- imd2010bua[order(imd2010bua$imdrank2010),]
imd2010bua$imdrank2010 <- 1:nrow(imd2010bua)

westmids <- merge(westmids, imd2010bua, by="BUA11CD", all.x=TRUE)
westmids <- westmids[order(westmids$imdrank2010),]
westmids$imdrank2010 <- 1:nrow(westmids)

imd2015 <- read_csv("imd2015lsoa.csv")
imd2015 <- imd2015[which(imd2015$Measurement=="Rank"), ]
imd2015 <- imd2015[which(imd2015$`Indices of Deprivation`=="a. Index of Multiple Deprivation (IMD)"), ]
names(imd2015)[names(imd2015)=="FeatureCode"] <- "LSOA11CD"
imd2015bua <- merge(imd2015, lsoabualookup, by="LSOA11CD", all.x=TRUE)
imd2015bua <- merge(imd2015bua, lsoapop, by="LSOA11CD", all.x=TRUE)

imd2015pop <- imd2015bua[-c(1:6, 8,9,10,11)]
imd2015pop <- aggregate(imd2015pop[-c(1)], imd2015pop["BUA11CD"], sum)
names(imd2015pop)[names(imd2015pop)=="mid2018pop"] <- "lsoa18popbestfit"

imd2015bua <- merge(imd2015bua, imd2015pop, by="BUA11CD", all.x=TRUE)
imd2015bua$imdrank2015 <- imd2015bua$Value*(imd2015bua$mid2018pop/imd2015bua$lsoa18popbestfit)
imd2015bua <- imd2015bua[-c(2:13)]
imd2015bua <- aggregate(imd2015bua[-c(1)], imd2015bua["BUA11CD"], sum)
imd2015bua <- imd2015bua[order(imd2015bua$imdrank2015),]
imd2015bua$imdrank2015 <- 1:nrow(imd2015bua)

westmids <- merge(westmids, imd2015bua, by="BUA11CD", all.x=TRUE)
westmids <- westmids[order(westmids$imdrank2015),]
westmids$imdrank2015 <- 1:nrow(westmids)

westmids <- westmids[order(westmids$imdrank2010),]

my_labs <- westmids$name
my_labs <- c(my_labs, rep("", 2*length(my_labs)))
my_labs2 <- westmids$name
my_labs2 <- c(rep("", 1*length(my_labs2)), my_labs2, rep("", 1*length(my_labs2)))
my_labs3 <- westmids$name
my_labs3 <- c(rep("", 2*length(my_labs3)), my_labs3)


effort1 <- ggparcoord(westmids,
                      columns = c(53, 54, 52), scale = "globalminmax")+
  theme(panel.grid = element_blank(), axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  scale_y_continuous(breaks=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,74), labels = c("Most Deprived Town", "5", "10","15","20", "25", "30", "35", "40", "45", "50",  "55", "60", "65", "70",  "Least Deprived Town"))+
  # scale_y_continuous(breaks=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,74), labels = c("Most Deprived Town","2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14","15", "16", "17", "18", "19", "20", "21","22", "23", "24", "25", "26", "27", "28","29", "30", "31", "32", "33", "34", "35","36", "37", "38", "39", "40", "41", "42", "43", "44", "45","46", "47", "48", "49", "50", "51", "52", "53", "54", "55","56", "57", "58", "59", "60", "61", "62", "63", "64", "65","66", "67", "68", "69", "70", "71", "72", "73", "Least Deprived Town"))+
  labs(x = "Year", y = "Multiple Deprivation Rank")+
  coord_cartesian(xlim = c(0.9, 3.2))+
  geom_text( aes(label=(my_labs)),  color = "navyblue", fontface = "bold", nudge_x= -0.27, size=2)+
  geom_text( aes(label=(my_labs2)),  color = "navyblue", fontface = "bold", nudge_x= 0 , size=2)+
  geom_text( aes(label=(my_labs3)),  color = "navyblue", fontface = "bold", nudge_x= 0.27 , size=2)+
  geom_line(size=0.4, color="grey", linetype="dotted", arrow = arrow(ends = "last", length = unit(0.1 , "inches")))+
  theme(axis.text.y = element_text(size = 10))
#ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "depplotwestmids.png", plot=effort1, width=6, height=13, dpi=600)

effort2 <- ggparcoord(westmids,
                      columns = c(53, 54, 52), scale = "globalminmax")+
  theme(panel.grid = element_blank(), axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  scale_y_continuous( breaks=c(1,5,10,15,20,25,30,35,38), labels = c("Most\nDeprived\nTown", "5", "10","15","20", "25", "30", "35", "Least\nDeprived\nTown"))+
  # scale_y_continuous(breaks=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,74), labels = c("Most Deprived Town","2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14","15", "16", "17", "18", "19", "20", "21","22", "23", "24", "25", "26", "27", "28","29", "30", "31", "32", "33", "34", "35","36", "37", "38", "39", "40", "41", "42", "43", "44", "45","46", "47", "48", "49", "50", "51", "52", "53", "54", "55","56", "57", "58", "59", "60", "61", "62", "63", "64", "65","66", "67", "68", "69", "70", "71", "72", "73", "Least Deprived Town"))+
  labs(x = "Year", y = "Multiple Deprivation Rank")+
  scale_x_discrete(labels = c("2010", "2015", "2019"))+
  coord_cartesian(xlim = c(0.9, 3.2))+
  geom_text( aes(label=(my_labs)),  color = "navyblue", fontface = "bold", nudge_x= -0.27, size=3)+
  geom_text( aes(label=(my_labs2)),  color = "navyblue", fontface = "bold", nudge_y= -0.0 , size=3)+
  geom_text( aes(label=(my_labs3)),  color = "navyblue", fontface = "bold", nudge_x= 0.27 , size=3)+
  geom_line(size=0.4, color="grey", linetype="solid", arrow = arrow(ends = "last", length = unit(0.1 , "inches")))+
  theme(axis.text.y = element_text(size = 10))
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "depplotwestmids.png", plot=effort2, width=6, height=5, dpi=600)

# effort3<- ggparcoord(Southwest,
#                      columns = c(53, 54, 52), scale = "globalminmax")+
#   theme(panel.grid = element_blank(), axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
#   scale_y_continuous(limits = c(55,74), breaks=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,74), labels = c("Most\nDeprived\nTown", "5", "10","15","20", "25", "30", "35", "40", "45", "50",  "55", "60", "65", "70",  "Least\nDeprived\nTown"))+
#   # scale_y_continuous(breaks=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,74), labels = c("Most Deprived Town","2", "3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13", "14","15", "16", "17", "18", "19", "20", "21","22", "23", "24", "25", "26", "27", "28","29", "30", "31", "32", "33", "34", "35","36", "37", "38", "39", "40", "41", "42", "43", "44", "45","46", "47", "48", "49", "50", "51", "52", "53", "54", "55","56", "57", "58", "59", "60", "61", "62", "63", "64", "65","66", "67", "68", "69", "70", "71", "72", "73", "Least\nDeprived\nTown"))+
#   labs(x = "Year", y = "Multiple Deprivation Rank")+
#   scale_x_discrete(labels = c("2010", "2015", "2019"))+
#   coord_cartesian(xlim = c(0.9, 3.2))+
#   geom_text( aes(label=(my_labs)),  color = "navyblue", fontface = "bold", nudge_x= -0.17, size=3)+
#   geom_text( aes(label=(my_labs2)),  color = "navyblue", fontface = "bold", nudge_y= -0.0 , size=3)+
#   geom_text( aes(label=(my_labs3)),  color = "navyblue", fontface = "bold", nudge_x= 0.27 , size=3)+
#   geom_line(size=0.4, color="grey", linetype="solid", arrow = arrow(ends = "last", length = unit(0.1 , "inches")))+
#   theme(axis.text.y = element_text(size = 10))
# ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "depplotsouthwest32.png", plot=effort3, width=6.4, height=5, dpi=600)
# #  geom_segment(aes(x = 0.4, xend = 0.4, y = 39.5, yend=38.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 34.6, yend=35.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 34.4, yend=33.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 30.5, yend=29.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 26.6, yend=27.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 25.5, yend=26.4), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 26.5, yend=27.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 21.5, yend=20.6), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 20.4, yend=19.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 14.5, yend=15.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 10.5, yend=9.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)+
# # geom_segment(aes(x = 0.4, xend = 0.4, y = 7.5, yend=8.5), arrow = arrow(length = unit(0.2, "cm"), type = "closed"), size=1)
# # 


growthwestmids <- ggplot(data=westmids, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_text(aes(x=0, y=100, label="British Town Average"))
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "westmidsfigure1.png", plot=growthwestmids, width=15, height=10, dpi=600)

# 
# growthscotlandz <- ggplot(data=westmids, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
#   geom_point()+
#   labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in Scottish Towns and Population Size")+
#   scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
#   theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
#   geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
#   stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

growthwestmidsservice <- ggplot(data=westmids, aes(x=services, y=publicchange, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Public Service Levels (Per Capita, 2018)", y = "Changes to Public Services (2011-18)")+
  scale_y_continuous(breaks=c(-2,2), labels = c("Service Decrease", "Service Increase"))+
  scale_x_continuous(breaks=c(-3,3), labels = c("Lower Service Provision", "Higher Service Provision"))+
  theme(axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black", size=2.7)+
  geom_text(aes(y=0, x=11, label="British Town Average"), size=3, angle=270)+
  geom_text(aes(y=6, x=0, label="British Town Average"), size=3)
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "westmidsfigure2.png", plot=growthwestmidsservice, width=13, height=10, dpi=600)

westmidsjobs <- ggplot(data=westmids, aes(x=jobdensity, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Job Density (2016)", y = "Population")+
  geom_vline(xintercept = 0.5, linetype="dashed", color="grey", size=0.8)+
  geom_vline(xintercept = 0.7, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 56.726, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black", size=3.5)+
  geom_text(aes(x=0.4, y=140000, label="Residential Towns"), size=4)+
  geom_text(aes(x=0.6, y=140000, label="Partially Residential Towns"), size=4)+
  geom_text(aes(x=0.8, y=140000, label="Working Towns"), size=4)
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "westmidsfigure3.png", plot=westmidsjobs, width=18, height=9.5, dpi=600)

wmdep <- ggplot(data=westmids, aes(x=severedepthreeorfour, y=housedep, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Severe Household Deprivation (%, 2011)", y = "Slight Household Deprivation (%, 2011)")+
  geom_vline(xintercept = 5.049, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 56.726, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black", size=3.4)+
  geom_text(aes(x=10, y=56.726, label="British Town Average"), size=4, angle=270)+
  geom_text(aes(x=5.049, y=71, label="British Town Average"), size=4)+
  geom_text(aes(x=4, y=70, label="1"), size=10)+
  geom_text(aes(x=8, y=70, label="2"), size=10)+
  geom_text(aes(x=8, y=46, label="3"), size=10)+
  geom_text(aes(x=4, y=46, label="4"), size=10)
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "westmidsfigure4.png", plot=wmdep, width=18, height=9.5, dpi=600)


###Comparing westmids to Britain###

data18 <- data18 %>% mutate(wm =ifelse(data18$Region=='West Midlands', "West Midlands Towns","Rest of Britain"))
data18$wm <- factor(data18$wm, levels = c("Rest of Britain", "West Midlands Towns"))

h <- lm(healths~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18)
b <- lm(buss~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18)
n <- lm(nurserys~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18)
s <- lm(schoolss~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18)
g <- lm(gpss~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18)
p <- lm(pboxs~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18)
output<-tab_model(h,b,n,s,g, p, show.ci=0.95,title="OLS Public Services" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Health Services", "Bus Stops", "Nurseries", "Schools","GPs",  "Post Boxes"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"West Midlands Dummy", "Job Density"), file = "table4wm")



# a <- plot_model(h,dot.size = 6,colors = "Dark2", terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey", axis.lim = c(-15,15), axis.title = "Relative Number of Health Providers", title = "Health Providers")
# b <- plot_model(b, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey", axis.lim = c(-40,40), axis.title = "Relative Number of Bus Stops", title = "Bus Stops")
# c <- plot_model(n,dot.size = 6,colors = "seagreen2", terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-4,4), axis.title = "Relative Number of Nursery Schools", title = "Nursery Schools")
# d <- plot_model(s, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-3,3), axis.title = "Relative Number of Schools", title = "Schools")
# e <- plot_model(g, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-1,1), axis.title = "Relative Number of Doctor's Surgeries", title = "Doctor's Surgeries")
# f <- plot_model(p, dot.size = 6, colors = "navyblue",terms = c("sw [South West Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-7,7), axis.title = "Relative Number of Post Boxes", title = "Post Boxes")
ff <- plot_model(n, dot.size = 6, colors = "navyblue",terms = c("wm [West Midlands Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-7,7), axis.title = "Relative Number of Nursery Schools", title = "Nursery Schools")
# ffff <- plot_model(h, dot.size = 6, colors = "navyblue",terms = c("sw [South West Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-7,7), axis.title = "Relative Number of Health Services", title = "Health-related Services")
# ffffff <- plot_model(b, dot.size = 6, colors = "navyblue",terms = c("sw [South West Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-30,30), axis.title = "Relative Number of Bus Stops", title = "Bus Stops")

# a <- a + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
# b <- b + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
# c <- c + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
# ffff <- ffff + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
ff <- ff + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
# f <- f + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
# ffffff <- ffffff + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")


ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "southwestfigure10.png", plot=ff, width=6, height=7.5, dpi=600)

# ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "southwestfigure11.png", plot=ffffff, width=6, height=3, dpi=600)

data18 <- data18 %>% mutate(mental=ifelse(data18$mentals>0, 1, 0))
data18 <- data18 %>% mutate(hospital=ifelse(data18$hospitals>0, 1, 0))
data18 <- data18 %>% mutate(station=ifelse(data18$stations>0, 1, 0))
data18 <- data18 %>% mutate(jobcent=ifelse(data18$jobcents>0, 1, 0))
data18 <- data18 %>% mutate(furthered=ifelse(data18$furthereds>0, 1, 0))
data18 <- data18 %>% mutate(police=ifelse(data18$polices>0, 1, 0))
data18 <- data18 %>% mutate(librari=ifelse(data18$librarys>1, 1, 0))



data18$mental <- factor( data18$mental)
data18$hospital <- factor( data18$hospital)
data18$furthered <- factor( data18$furthered)
data18$jobcent <- factor( data18$jobcent)
data18$station <- factor( data18$station)
data18$police <- factor( data18$police)
data18$country <- factor(data18$country)

a <- glm(mental~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18, family = "binomial")

b <- glm(hospital~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18, family = "binomial")

c <- glm(jobcent~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18, family = "binomial")
d <- glm(station~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18, family = "binomial")

e <- glm(furthered~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18, family = "binomial")

f <- glm(police~growth+distance+poplsoa+size+housedep+wm+jobdensity, data=data18, family = "binomial")
x<-tab_model(a,b,c,d,e, f, show.ci=0.95,title="Logistic Regression of Public Services availability" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Mental Health Practitioners", "Hospital", "Job Centre", "Train Station","Further Education College",  "Police Station"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"West Midlands Dummy", "Job Density"), file = "table5wm")

w<- plot_model(e,ci.lvl = 0 ,type = "pred", colors = c("navyblue", "blue"), terms = c("housedep", "wm"),  axis.title = c( "Households in Deprivation (%, 2011)", "Probability of Present FE College"), legend.title = "Region", title="Further Education Colleges")
y<- plot_model(c,ci.lvl = 0 ,type = "pred", colors = c("seagreen2", "seagreen4"), terms = c("housedep", "sw"), title = "Predicted Probability of a Town-based Job Centre", axis.title = c( "Households in Deprivation (%, 2011)", "Probability of Present Job Centre"), legend.title = "Region", title="Further Education Colleges")
z<- plot_model(f,ci.lvl = 0, type = "pred", colors = c("seagreen2", "seagreen4"), terms = c("workinggroup", "sw"), title = "Predicted Probability of a Town-based Police Station", axis.title = c( "Town Type", "Probability of Present Police Station"), legend.title = "Region")

w <- w+geom_line(size=2)
  # theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))
y <- y+geom_line(size=2)+coord_cartesian(ylim=c(0.25,1))
z <- z+geom_point(size=6)+coord_cartesian(ylim=c(0.25,1))

compare <- ggarrange(w,y, ncol=2, nrow=1)
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "westmidsfigure12.png", plot=w, width=7, height=5, dpi=600)


####covidwm####

coviddeaths <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\coviddeathsjune.csv")
names(coviddeaths)[names(coviddeaths)=="MSOA code"] <- "MSOA11CD"
msoabualookup <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\msoabualookup.csv")
msoapop18 <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\msoapop18.csv")


coviddeaths <-merge(coviddeaths, msoabualookup, by="MSOA11CD", all.x=TRUE)
coviddeaths <-merge(coviddeaths, msoapop18, by="MSOA11CD", all.x=TRUE)
coviddeaths <- aggregate(coviddeaths[-c(1,3,4,5,6,7,8,9,10,11)], coviddeaths["BUA11CD"], sum)

data18 <- merge(data18, coviddeaths, by="BUA11CD", all.x=TRUE)
data18$covidperten <- data18$coviddeathsmarjune/(data18$`All Ages`/10000)



westmids <- merge(westmids, coviddeaths, by="BUA11CD", all.x=TRUE)
westmids$covidperten <- westmids$coviddeathsmarjune/(westmids$`All Ages.y`/10000)


covidworkplot <- ggplot(data=westmids, aes(x=covidperten, y=growth, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Covid-related Deaths (Per 10,000 Capita, March-May 2020)", y = "Job Density ()")+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black", size=2.7)




###JSAchange###


jsachange <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\jsa1920.csv")
jsachange$ratechange <- jsachange$ratemay20-jsachange$ratemay19

data18 <- merge(data18, jsachange[-c(2,3,4,5)], by="bua11nm", all.x=TRUE)
westmids <- merge(westmids, jsachange[-c(2,3,4,5)], by="bua11nm", all.x=TRUE)

data18 <- data18 %>% mutate(Regionwm=ifelse(data18$Region=="West Midlands", "West Midlands",  "Other"))
data18$Regionwm <- factor(data18$Regionwm)

tourismplot <- ggplot(data=data18, aes(x=ratechange, y=housedep, group=BUA11CD, label=name))+
  geom_point(aes( colour=Regionwm, alpha=Regionwm))+
  #geom_text_repel(data=data18s %>% filter(Region=="South West" &tourism>15 ),
   #               aes(label=tourismsample$name), size=3, color="black")+
  labs(y = "Household Deprivation (%, 2011)", x = "Employment in food, drink, accomodation or tour guide industries\n(% of total employment, 2018)")+
  scale_alpha_discrete(range = c(0.5, 1))+
  scale_color_manual(values=c("#FF99CC", "#000066"))


ggplot(data=westmids, aes(x=ratechange, y=housedep, group=BUA11CD, label=name))+
  geom_point(aes( colour=Regionwm, alpha=Regionwm))+
  #geom_text_repel(data=data18s %>% filter(Region=="South West" &tourism>15 ),
  #               aes(label=tourismsample$name), size=3, color="black")+
  labs(y = "Household Deprivation (%, 2011)", x = "Employment in food, drink, accomodation or tour guide industries\n(% of total employment, 2018)")+
  scale_alpha_discrete(range = c(0.5, 1))+
  scale_color_manual(values=c("#FF99CC", "#000066"))

####claimantchane####
claimantchange <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\claimantcount1920.csv")
claimantchange$claimchange <- claimantchange$claimrate20-claimantchange$claimrate19
westmids <- merge(westmids, claimantchange[-c(2,3)], by="bua11nm", all.x=TRUE)



claimantchange$claimperchange <- (claimantchange$claimrate20/claimantchange$claimrate19)*100

scotclaim <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\scotclaim.csv")
intzonelookup <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\intlookup.csv")
intzonepop18 <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\interzonepop1664.csv")
intzonepop18 <- intzonepop18[which(intzonepop18$Year=="2018" & intzonepop18$Sex=="All"),]


scotclaim <-merge(scotclaim, intzonelookup, by="InterZone2011", all.x=TRUE)


intzonepop18 <-merge(intzonepop18, intzonelookup, by="InterZone2011", all.x=TRUE)
intzonepop18 <- aggregate(intzonepop18[-c(1,2,3,4,6)], intzonepop18["BUA11CD"], sum)
names(intzonepop18)[names(intzonepop18)=="Age1664"] <- "totpopint"

scotclaim <-merge(scotclaim, intzonepop18, by="BUA11CD", all.x=TRUE)
intzonepop18 <- read_csv("C:\\Users\\bengo\\OneDrive\\Documents\\interzonepop1664.csv")
intzonepop18 <- intzonepop18[which(intzonepop18$Year=="2018" & intzonepop18$Sex=="All"),]
scotclaim <-merge(scotclaim, intzonepop18, by="InterZone2011", all.x=TRUE)
scotclaim$claimratechange <- (scotclaim$)*(scotclaim$Age1664/scotclaim$totpopint)-

tourismscot <- aggregate(tourismscot[-c(1,2,3,4,5,6,7)], tourismscot["BUA11CD"], sum)



data18 <- merge(data18, claimantchange[-c(2,3)], by="bua11nm", all.x=TRUE)
data18 <- data18 %>% mutate(Regions=ifelse(data18$Region=="West Midlands", "North and Midlands", ifelse(data18$Region=="East Midlands", "North and Midlands", ifelse(data18$Region=="North West", "North and Midlands", ifelse(data18$Region=="North East", "North and Midlands", ifelse(data18$Region=="North East", "North and Midlands", ifelse(data18$Region=="Yorkshire and the Humber", "North and Midlands", ifelse(data18$Region=="South East", "South", ifelse(data18$Region=="South West", "South", ifelse(data18$Region=="London", "South",  "Wales"))))))))))
data18$Regions <- factor(data18$Regions)

dataeng <- data18[which(data18$country=="England"), ]
datascot <- data18[which(data18$country=="Scotland"), ]
datawales <- data18[which(data18$country=="Wales"), ]

Engsmaple <- dataeng[which(dataeng$claimchange>5),]
wasmaple <- datawales[which(datawales$claimchange>4),]
scotsmaple <- datascot[which(datascot$claimchange>5),]

engclaimantchangeplot <- ggplot(data=dataeng, aes(x=claimchange, y=housedep, label=name))+
  geom_point(colour = "#FF99CC")+
  geom_smooth(method= lm , color="#990099", se=FALSE)+
  geom_text_repel(data=dataeng %>% filter(claimchange>5 ),
                 aes(label=Engsmaple$name), size=3, color="black")+
  labs(y = "Household Deprivation (%, 2011)", x = "Change in Claimant rate between May 2019 and May 2020\n(% Points, Rate as a Proportion of Residents aged 16-64)", title = "Claimant Rate Change in English Towns", caption = "source: Census, 2011; ONS, Claimant Count (Experimental Statistics)")
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "engclaimantratechange.png", plot=engclaimantchangeplot, width=7, height=6, dpi=600)

walclaimantchangeplot <- ggplot(data=datawales, aes(x=claimchange, y=housedep, label=name))+
  geom_point(colour = "#FF99CC")+
  geom_smooth(method= lm , color="#990099", se=FALSE)+
  geom_text_repel(data=datawales %>% filter(claimchange>4 ),
                  aes(label=wasmaple$name), size=3, color="black")+
  labs(y = "Household Deprivation (%, 2011)", x = "Change in Claimant rate between May 2019 and May 2020\n(% Points, Rate as a Proportion of Residents aged 16-64)", title = "Claimant Rate Change in Welsh Towns", caption = "source: Census, 2011; ONS, Claimant Count (Experimental Statistics)")
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "walclaimantratechange.png", plot=walclaimantchangeplot, width=7, height=6, dpi=600)



claimantchangeplot <- ggplot(data=data18, aes(x=claimchange, y=housedep, label=name))+
  geom_point()+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black", size=3.4)+
  #geom_text_repel(data=data18s %>% filter(Region=="South West" &tourism>15 ),
  #               aes(label=tourismsample$name), size=3, color="black")+
  labs(y = "Household Deprivation (%, 2011)", x = "Change in Claimaint rate between May 2019 and May 2020 (% Points)")
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "claimantratechange.png", plot=claimantchangeplot, width=8, height=6, dpi=600)

why <- data18[complete.cases(data18),] 

tourismplot <- ggplot(data=why, aes(x=claimchange, y=housedep, group=BUA11CD))+
  geom_point()+
  geom_smooth(method=lm, color="red", fill="#69b3a2", se=TRUE)+
  labs(y = "Household Deprivation (%, 2011)", x = "Change in claimant count Between May 2019 and May 2020 (% point)")


tourismplot <- ggplot(data=why, aes(x=claimperchange, y=housedep, group=BUA11CD))+
  geom_point()+
  geom_smooth(method=lm, color="red", fill="#69b3a2", se=TRUE)+
  labs(y = "Household Deprivation (%, 2011)", x = "Change in claimant count Between May 2019 and May 2020 (% point)")




claimantchangeplot <- ggplot(data=westmids, aes(x=claimchange, y=housedep, label=name))+
  geom_point()+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black", size=3.4)+
  #geom_text_repel(data=data18s %>% filter(Region=="South West" &tourism>15 ),
  #               aes(label=tourismsample$name), size=3, color="black")+
  labs(y = "Household Deprivation (%, 2011)", x = "Change in Claimaint rate between May 2019 and May 2020 (% Points)")
ggsave(path = "C:\\Users\\bengo\\OneDrive\\Documents\\", filename = "wmidzclaimantchangeplot.png", plot=claimantchangeplot, width=9, height=6, dpi=600)


####wmcamap####

cas <- shapefile("caboundaries")
cas <- cas[which(cas$objectid==6),]
nuts1 <- shapefile("nuts")
nuts1 <- nuts1[which(nuts1$objectid==5),]
cas <- spTransform(cas, CRS("+init=epsg:4326"))
nuts1 <- spTransform(nuts1, CRS("+init=epsg:4326"))
names(cas)[names(cas)=="cauth19cd"] <- "nuts118cd"
names(cas)[names(cas)=="cauth19nm"] <- "nuts118nm"

cas2 <- rbind(cas, nuts1)
cas2$name <- c("Combined Authority", "West Midlands Region")


mymap <- leaflet(data=cas2) %>%
  addTiles() %>%  
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
                               opacity = 1.0, fillOpacity = 0.5,
                               fillColor = c("red", "yellow"), label = cas2$name, labelOptions = labelOptions(noHide = T, direction = c("right", "left"),style = list("color" = "blue", "font-weight" = "900"))) %>% 
  setView(lng=-2.2, lat=52.5275, zoom=9) 
mymap <- addMarkers(mymap, icon = icons, data=uktownswm, label =uktownswm$name, labelOptions = labelOptions(noHide = T, textOnly = T, direction = 'auto', style = list("font-style" = "italic", "font-weight" = "600")))  
   
mapshot(mymap, file="figurecawm4.png", zoom=20, vwidth=750, vheight=855)  

