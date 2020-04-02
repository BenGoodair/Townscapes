###load the packages###
library(htmlwidgets)
library(htmltools)

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

###public service index###

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


###North East####

###north east map###
nuts1 <- readOGR("C:\\Users\\bjg55\\Documents\\Towns\\2019Ring", "nuts1")
nuts1 <- nuts1[which(nuts1$objectid!=12), ]

nuts1 <- spTransform(nuts1, CRS("+init=epsg:27700"))
nuts1s <- as.data.frame(nuts1)
nuts1s <- nuts1s %>% mutate(NE = ifelse(nuts1s$objectid==1, 1, 0))
nuts1ss <- merge(nuts1, nuts1s, by="objectid")
# nuts1ss$NE <- factor(nuts1ss$NE)

map <- spplot(nuts1ss, "NE", col.regions = c("white", "seagreen2"), colorkey=FALSE)
# ggsave(map,"trypleaseggsave.png", width=7.64, height=5.5, dpi=100)
# dev.print(png, file="scripttry", width=764, height=550)
# print(spplot(nuts1ss, "NE", col.regions = c("white", "seagreen2"), colorkey=FALSE))
# png(file="scripttry", width=764, height=550, res=1200)
# dev.off()

northeast <- data18[which(data18$Region=='North East'), ]
northeastpanel <- masterpanel[which(masterpanel$Region=='North East'), ]
northeastlist <- masterlist[which(masterlist$Region=='North East'), ]

###school travel time map####

#  schooltime <- read_csv("school public transport time 2016.csv")
#  schooltime <- schooltime[which(schooltime$Region=="North East"), ]
#  lsoa <- readOGR("C:\\Users\\bjg55\\Documents\\Towns\\2019Ring", "lsoa")
#  lsoa <- spTransform(lsoa, CRS("+init=epsg:27700"))
#  lsoaschool <- merge(lsoa, schooltime, by="lsoa11cd")
#  lsoaschool <- sp.na.omit(lsoaschool)
#  cutseq <- seq(0,120, by=12) 
#  netowns <- merge(uktowns, region, by="BUA11CD")
#  netowns <- spTransform(netowns, CRS("+init=epsg:4326"))
# lsoaschool <- spTransform(lsoaschool, CRS("+init=epsg:4326"))

netowns <- netowns[which(netowns$Region=="North East"),]
netowns <- merge(netowns, popbua, by="BUA11CD")
netowns <- netowns[which(netowns$poplsoa>10000&netowns$poplsoa<175000),]

 # townlayer <- list( "sp.polygons",netowns, fill="gray")
 # schooltimemap <- spplot(lsoaschool, "SSPTt", alpha=0.5, sp.layout=list(townlayer), main=list(label="Travel Time to nearest Secondary School by Public Tranpsort (Minutes)", cex=1), col="transparent")
 # 
 # trying <- fortify(lsoaschool, region="lsoa11cd")
 # tryings <- merge(trying, lsoaschool@data, by.x="id", by.y="lsoa11cd")
 # 
 # tryingz <- fortify(netowns, region="BUA11CD")
 # tryingzz <- merge(tryingz, netowns@data, by.x="id", by.y="BUA11CD")
 # tryingsss <- tryings[which(tryings$SSPTt<30),]
 # 
 # plz <- ggplot()+
 #   geom_polygon(aes(long, lat, group=group, fill=tryingsss$SSPTt),tryingsss)+
 #   coord_quickmap()+
 #   geom_polygon(aes(long, lat, group=group),tryingzz, fill=NA, color="yellow", size=0.7)+
 #   labs(fill="Travel Time (Minutes)")+
 #   coord_cartesian(ylim=c(55.03,55.104), xlim=c(-1.62, -1.503))+
 #   geom_label(aes(x=-1.59, y=55.04, label="Cramlington", fontface="bold"))+
 #   geom_label(aes(x=-1.525, y=55.04, label="Seaton Dalevale",  fontface="bold"))+
 #   theme_nothing(legend=TRUE)
 #   
 # ggsave("Figure3.png", plot=plz, width=6.74, height=5.5, dpi=600)
 # 
   
 
 




growthnortheast <- ggplot(data=northeast, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_text(aes(x=0, y=70, label="British Town Average"))+
  geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
    labs(x = "Improvement Index", y = "Distance from Nearest City (km)")
  

#ggsave("Figure1.png", plot=growthnortheast, width=12, height=7.96, dpi=600)
growthnortheast <- ggplot(data=northeast, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_text(aes(x=0, y=70, label="British Town Average"))+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black")+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)")


ggsave("Figure1navy.png", plot=growthnortheast, width=12, height=7.96, dpi=600)

# growthnortheastz <- ggplot(data=northeast, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
#   geom_point()+
#   labs(x = "Improvement Index", y = "Population Size (2016)", title = "Improvement or Decline in North East Towns and Population Size")+
#   scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
#   theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
#   geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")

# growthnortheastservice <- ggplot(data=northeast, aes(x=services, y=growth, group=BUA11CD, label=name))+
#   geom_point()+
#   labs(x = "Public Service Index", y = "Improvement Index")+
#   scale_y_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
#   scale_x_continuous(breaks=c(-3,3), labels = c("Lower Service Provision", "Higher Service Provision"))+
#   theme(axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
#   theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
#   geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
#   geom_hline(yintercept = 0, linetype="dashed", color="grey", size=0.8)+
#   geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
#   geom_text(aes(y=0, x=10, label="British Town Average"), size=3, angle=270)+
#   geom_text(aes(y=6.7, x=0, label="British Town Average"), size=3)
# ggsave("Figure2.png", plot=growthnortheastservice, width=12, height=7.96, dpi=600)

growthnortheastservicechange <- ggplot(data=northeast, aes(x=services, y=publicchange, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Public Service Levels (Per Capita, 2018)", y = "Changes to Public Services (2011-18)")+
  scale_y_continuous(breaks=c(-2,2), labels = c("Service Decrease", "Service Increase"))+
  scale_x_continuous(breaks=c(-3,3), labels = c("Lower Service Provision", "Higher Service Provision"))+
  theme(axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
  geom_text(aes(y=0, x=11, label="British Town Average"), size=3, angle=270)+
  geom_text(aes(y=4, x=0, label="British Town Average"), size=3)
#ggsave("Figure2.png", plot=growthnortheastservicechange, width=10, height=7.96, dpi=600)
growthnortheastservicechange <- ggplot(data=northeast, aes(x=services, y=publicchange, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Public Service Levels (Per Capita, 2018)", y = "Changes to Public Services (2011-18)")+
  scale_y_continuous(breaks=c(-2,2), labels = c("Service Decrease", "Service Increase"))+
  scale_x_continuous(breaks=c(-3,3), labels = c("Lower Service Provision", "Higher Service Provision"))+
  theme(axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black")+
  geom_text(aes(y=0, x=11, label="British Town Average"), size=3, angle=270)+
  geom_text(aes(y=4, x=0, label="British Town Average"), size=3)
ggsave("Figure2navy.png", plot=growthnortheastservicechange, width=10, height=7.96, dpi=600)

nejobs <- ggplot(data=northeast, aes(x=jobdensity, y=housedep, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Job Density (2016)", y = "Household Deprivation (%, 2011)")+
  geom_vline(xintercept = 0.5, linetype="dashed", color="grey", size=0.8)+
  geom_vline(xintercept = 0.7, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 56.726, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
  geom_text(aes(x=0.4, y=70, label="Residential Towns"), size=4)+
  geom_text(aes(x=0.6, y=70, label="Partially Residential Towns"), size=4)+
  geom_text(aes(x=1.05, y=56.726, label="British Town Average"), size=4, angle=270)+
  geom_text(aes(x=0.8, y=70, label="Working Towns"), size=4)
#ggsave("Figure4.png", plot=nejobs, width=11, height=7.96, dpi=600)
nejobs <- ggplot(data=northeast, aes(x=jobdensity, y=housedep, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Job Density (2016)", y = "Household Deprivation (%, 2011)")+
  geom_vline(xintercept = 0.5, linetype="dashed", color="grey", size=0.8)+
  geom_vline(xintercept = 0.7, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 56.726, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black")+
  geom_text(aes(x=0.4, y=70, label="Residential Towns"), size=4)+
  geom_text(aes(x=0.6, y=70, label="Partially Residential Towns"), size=4)+
  geom_text(aes(x=1.05, y=56.726, label="British Town Average"), size=4, angle=270)+
  geom_text(aes(x=0.8, y=70, label="Working Towns"), size=4)
ggsave("Figure4navy.png", plot=nejobs, width=11, height=7.96, dpi=600)



northeast$busperpop <- northeast$buss/(northeast$poplsoa/1000)
northeast$healthsper <- northeast$healths/(northeast$poplsoa/1000)
northeast$pboxsper <- northeast$pboxs/(northeast$poplsoa/1000)
northeast$nurseryper <- northeast$nurserys/(northeast$poplsoa/1000)
northeast$schoolper <- northeast$schoolss/(northeast$poplsoa/1000)
northeast$wifiper <- northeast$wifis/(northeast$poplsoa/1000)
northeast$hallper <- northeast$halls/(northeast$poplsoa/1000)


northeast$train <- factor(northeast$stations)
northeast$declinegroup <- factor(northeast$declinegroup, levels = c("Declining", "Stagnant", "Improving"))
northeast$distancegroup <- factor(northeast$distancegroup, levels = c("Neighbouring", "Middle", "Isolated"))

busnortheast <- ggplot(data=northeast, aes(x=busperpop, y=growth, label=name))+
  geom_point()+
  labs(x = "Bus Stops (per 1000 capita)", y = "Population (2016)", title = "Bus Stops in North East Towns")+
  geom_label_repel()

northeast$id <- seq(1, nrow(northeast))
label_ne = northeast
numberofbar = nrow(label_ne)
angle = 90-360*(label_ne$id-0.5)/numberofbar
label_ne$hjust <- ifelse(angle< -90, 1, 0)
label_ne$angle <- ifelse(angle< -90, angle+180,angle)


t<- ggplot(northeast, aes(x=as.factor(id), y=busperpop))+
  geom_bar(stat="identity", fill=("seagreen2"))+
  ylim(-4,45)+
  theme_minimal()+
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), plot.margin = unit(rep(-10,5), "cm"))+
  coord_polar(start = 0)+
  geom_text(aes(x=0, y=24, label= "Bus Stops per 1000 Population in North East Towns"))+
  geom_text(data=label_ne, aes(x=id, y=busperpop+1, label=paste("",name," = ", round(busperpop, digits = 2),""), hjust=hjust), color="black", fontface="bold", size=3, angle= label_ne$angle, inherit.aes = FALSE)
  

busne <- ggplot(data=northeast, aes(x = northeast$distancegroup, y = busperpop))+
  geom_bar(aes(distancegroup, busperpop, fill = northeast$distancegroup), 
           position = "dodge", stat = "summary", fun.y = "mean")+
  labs(x = "Distance Group", y = "Bus Stops (per 1000 capita)", title = "Distance Groups")+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0,7))+
  scale_fill_brewer(palette = "Dark2")


busnes <- ggplot(data=northeast, aes(x = northeast$declinegroup, y = busperpop))+
  geom_bar(aes(declinegroup, busperpop, fill = northeast$declinegroup), 
           position = "dodge", stat = "summary", fun.y = "mean")+
  labs(x = "Improvement Group", y = "Bus Stops (per 1000 capita)", title = "Improvement Group")+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0,7))+
  scale_fill_brewer(palette = "Dark2")

busplots <- ggarrange(busne, busnes)

#busnes <- ggplot(data=northeast, aes(x = northeast$declinegroup, y = busperpop))+
 # geom_bar(aes(declinegroup, busperpop, fill = northeast$declinegroup), 
  #         position = "dodge", stat = "summary", fun.y = "mean")+
  #labs(x = "Improvement Group", y = "Bus Stops (per 1000 capita)", title = "Improvement Group")+
  #theme(legend.position = "none")+
  #coord_cartesian(ylim = c(0,7))+
  #facet_grid(.~distancegroup)



nurserytime <- ggplot(data=northeastpanel, aes(x=date, y=nurserys, group=BUA11CD))+
  geom_point()+
  geom_line()+
  labs(x = "Year", y = "Number of Nursery Schools", title = "North East Nursery School Provision")





data18 <- data18 %>% mutate(ne =ifelse(data18$Region=='North East', "North East","Rest of Britain"))
data18$ne <- factor(data18$ne, levels = c("Rest of Britain", "North East"))

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(ne) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = ne, y = busperpop, color = ne, fill = ne)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
 # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North East", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )
 
healthaverage <- ggplot(data18, aes(x = ne, y = healthsper, color = ne, fill = ne)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North East", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = ne, y = schoolper, color = ne, fill = ne)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North East", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = ne, y = busperpop, color = ne, fill = ne)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North East", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = ne, y = pboxsper, color = ne, fill = ne)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North East", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = ne, y = gpsper, color = ne, fill = ne)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North East", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)


data18 <- data18 %>% mutate(ne =ifelse(data18$Region=='North East', "North East Towns","Rest of Britain"))
data18$ne <- factor(data18$ne, levels = c("Rest of Britain", "North East Towns"))



comparehealthne <- lm(healths~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18)
comparebusne <- lm(buss~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18)
comparenurseryne <- lm(nurserys~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18)
compareschoolne <- lm(schoolss~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18)
comparegpne <- lm(gpss~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18)
comparepostboxne <- lm(pboxs~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18)
tab_model(comparehealthne,comparebusne,comparenurseryne,compareschoolne,comparegpne, comparepostboxne, show.ci=0.95,title="OLS Public Services" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Health Services", "Bus Stops", "Nurseries", "Schools","GPs",  "Post Boxes"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"North East Dummy", "Job Density"), file = "table4")


a <- plot_model(comparehealthne,dot.size = 6,colors = "Dark2", terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey", axis.lim = c(-15,15), axis.title = "Relative Number of Health Providers", title = "Health Providers")
b <- plot_model(comparebusne, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey", axis.lim = c(-40,40), axis.title = "Relative Number of Bus Stops", title = "Bus Stops")
c <- plot_model(comparenurseryne,dot.size = 6,colors = "seagreen2", terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-4,4), axis.title = "Relative Number of Nursery Schools", title = "")
d <- plot_model(compareschoolne, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-3,3), axis.title = "Relative Number of Schools", title = "Schools")
e <- plot_model(comparegpne, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-1,1), axis.title = "Relative Number of Doctor's Surgeries", title = "Doctor's Surgeries")
f <- plot_model(comparepostboxne, dot.size = 6, colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-8,8), axis.title = "Relative Number of Post Boxes", title = "Post Boxes")

a <- a + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
b <- b + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
c <- c + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
d <- d + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
e <- e + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
f <- f + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
ggsave("Figure10.png", plot=c, width=7.64, height=3.5, dpi=600)
c <- plot_model(comparenurseryne,dot.size = 6,colors = "navyblue", terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-4,4), axis.title = "Relative Number of Nursery Schools", title = "")
c <- c + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
ggsave("Figure10navy.png", plot=c, width=7.64, height=3.5, dpi=600)


averagecompared <- ggarrange(a,b,c,d,e,f,   ncol=3, nrow=2)


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
data18$librari <- factor( data18$librari)
data18$country <- factor(data18$country)

a <- glm(mental~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18, family = "binomial")

b <- glm(hospital~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18, family = "binomial")

c <- glm(jobcent~growth+distance+poplsoa+size+housedep+ne+workinggroup, data=data18, family = "binomial")

d <- glm(station~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18, family = "binomial")

e <- glm(furthered~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18, family = "binomial")

f <- glm(police~growth+distance+poplsoa+size+housedep+ne+workinggroup, data=data18, family = "binomial")

g <- glm(librari~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18, family = "binomial")

f <- glm(police~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18, family = "binomial")

tab_model(a,b,c,d,e, f, show.ci=0.95,title="Logistic Regression of Public Services availability" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Mental Health Practitioners", "Hospital", "Job Centre", "Train Station","Further Education College",  "Police Station"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"North East Dummy", "Job Density"), file = "table5")


# a <- glm(mental~growth+distance+poplsoa+size+housedep+workinggroup, data=data18, family = "binomial")
# 
# b <- glm(hospital~growth+distance+poplsoa+size+housedep+workinggroup, data=data18, family = "binomial")
# 
# c <- glm(jobcent~growth+distance+poplsoa+size+housedep+workinggroup, data=data18, family = "binomial")
# 
# d <- glm(station~growth+distance+poplsoa+size+housedep+workinggroup, data=data18, family = "binomial")
# 
# e <- glm(furthered~growth+distance+poplsoa+size+housedep+workinggroup, data=data18, family = "binomial")
# 
# f <- glm(police~growth+distance+poplsoa+size+housedep+workinggroup, data=data18, family = "binomial")
# 
# g <- glm(librari~growth+distance+poplsoa+size+housedep+workinggroup, data=data18, family = "binomial")

w<- plot_model(d, type = "pred", colors = "seagreen2", terms = c("growth [-5:5]"), title = "Predicted Probability of a Town-based Train Station", axis.title = c( "Economic and Social Improvement (2001-2011)", "Probability of Present Train Station"), legend.title = "Region")

w+theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=.96, label="Average Town Improvement"))+
  scale_x_continuous(breaks=c(-4,4), labels = c("Declining Towns", "Improving Towns"))

w<- plot_model(a, type = "pred", colors = "seagreen2", terms = c("housedep [35:80]"), title = "Predicted Probability of a Town-based Mental Health Provider", axis.title = c( "Households in Deprivation (%, 2011)", "Probability of Present Mental Health Service"), legend.title = "Region")
x<- plot_model(b, ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("growth [-5:5]", "ne"), title = "Predicted Probability of Hospital", axis.title = c( "Improvement Score", "Probability"), legend.title = "Region")
y<- plot_model(c,ci.lvl = 0,  colors = c("seagreen2", "seagreen3"),type = "pred", terms = c("growth [-5:5]", "ne"),title = "Predicted Probability of Job Centre", axis.title = c( "Improvement Score", "Probability"), legend.title = "Region")
z<- plot_model(d,ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("growth [-5:5]", "ne"),title = "Predicted Probability of Train Station", axis.title = c( "Improvement Score", "Probability"), legend.title = "Region")
t<- plot_model(e,ci.lvl = 0,  colors = c("seagreen2", "seagreen3"), type = "pred", terms = c("growth [-5:5]", "ne"),title = "Predicted Probability of Further Education", axis.title = c( "Improvement Score", "Probability"), legend.title = "Region")
u<- plot_model(f,ci.lvl = 0,dot.size=6  ,colors = c("seagreen2", "seagreen4"),  type = "pred", terms = c("workinggroup", "ne"),title = "", axis.title = c( "Town Type", "Probability"), legend.title = "Region")
v<- plot_model(g, type = "pred", terms = c("growth [-5:5]", "ne"))
x <- x+geom_line(size=2)+coord_cartesian(ylim=c(0.45,1))
y <- y+geom_line(size=2)+coord_cartesian(ylim=c(0.45,1))
z <- z+geom_line(size=2)+coord_cartesian(ylim=c(0.45,1))
t <- t+geom_line(size=2)+coord_cartesian(ylim=c(0.45,1))
v <- v+geom_line(size=2)+coord_cartesian(ylim=c(0.45,1))
ggsave("Figure12.png", plot=u, width=7.64, height=5, dpi=600)
u<- plot_model(f,ci.lvl = 0,dot.size=6  ,colors = c("navyblue", "blue"),  type = "pred", terms = c("workinggroup", "ne"),title = "", axis.title = c( "Town Type", "Probability"), legend.title = "Region")
ggsave("Figure12navy.png", plot=u, width=7.64, height=5, dpi=600)

averagecompared <-ggarrange(x,y,z,t, ncol=2, nrow=2)
v<- plot_model(g, type = "pred",ci.lvl = 0,  colors = "Dark2", terms = c("jobdensity [0:100]", "ne"),title = "Predicted Probability of Library", axis.title = c( "Job Density (2016)", "Probability"), legend.title = "Region")
v <- v+geom_line(size=2)+coord_cartesian(ylim=c(0.2,1))




x<- plot_model(b,ci.lvl = 0, colors = "Dark2",type = "pred", terms = c("distance [0:80]", "ne"), title = "Predicted Probability of Hospital", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
y<- plot_model(c,ci.lvl = 0, colors = "Dark2",type = "pred", terms = c("distance [0:80]", "ne"),title = "Predicted Probability of Job Centre", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
z<- plot_model(d, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("distance [0:80]", "ne"),title = "Predicted Probability of Train Station", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
t<- plot_model(e,ci.lvl = 0, colors = "Dark2", type = "pred", terms = c("distance [0:80]", "ne"),title = "Predicted Probability of Further Education", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
#u<- plot_model(f, type = "pred", terms = c("growth [-5:5]", "ne"))
#v<- plot_model(g, type = "pred", terms = c("growth [-5:5]", "ne"))

x <- x+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))
y <- y+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))
z <- z+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))
t <- t+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))

averagecompared <-ggarrange(x,y,z,t, ncol=2, nrow=2)

x<- plot_model(b,ci.lvl = 0, colors = "Dark2",type = "pred", terms = c("housedep [35:80]", "ne"), title = "Predicted Probability of Hospital", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
y<- plot_model(c,ci.lvl = 0, colors = c("seagreen2", "seagreen4"),type = "pred", terms = c("housedep [35:80]", "ne"),title = "", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
z<- plot_model(d, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("housedep [35:80]", "ne"),title = "Predicted Probability of Train Station", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
t<- plot_model(e,ci.lvl = 0, colors = "Dark2", type = "pred", terms = c("housedep [35:80]", "ne"),title = "Predicted Probability of Further Education", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
#u<- plot_model(f, type = "pred", terms = c("growth [-5:5]", "ne"))
#v<- plot_model(g, type = "pred", terms = c("growth [-5:5]", "ne"))
u<- plot_model(a,ci.lvl = 0, colors = "Dark2",type = "pred", terms = c("housedep [35:80]", "ne"), title = "Predicted Probability of Mental Health Service", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
u <- u+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))

x <- x+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))
y <- y+geom_line(size=2)+coord_cartesian(ylim=c(0.25,1))
z <- z+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))
t <- t+geom_line(size=2)+coord_cartesian(ylim=c(0.36,1))
ggsave("Figure11.png", plot=y, width=7.64, height=5.5, dpi=600)
y<- plot_model(c,ci.lvl = 0, colors = c("navyblue", "blue"),type = "pred", terms = c("housedep [35:80]", "ne"),title = "", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")

y <- y+geom_line(size=2)+coord_cartesian(ylim=c(0.25,1))
ggsave("Figure11navy.png", plot=y, width=7.64, height=5.5, dpi=600)

averagecompared <-ggarrange(x,y,z,t, ncol=2, nrow=2)

Health <- lm(healths~growth+distance+poplsoa+size+ne, data=data18)
Bus <- lm(buss~growth+distance+poplsoa+size+ne, data=data18)
Nurseries <- lm(nurserys~growth+distance+poplsoa+size+ne, data=data18)
Schools <- lm(schoolss~growth+distance+poplsoa+size+ne, data=data18)
GPs <- lm(gpss~growth+distance+poplsoa+size+ne, data=data18)
Postboxes <- lm(pboxs~growth+distance+poplsoa+size+ne, data=data18)

Healthy <- lm(healths~growth+distance+poplsoa+size, data=data18)
Busy <- lm(buss~growth+distance+poplsoa+size, data=data18)
Postboxesy <- lm(pboxs~growth+distance+poplsoa+size, data=data18)

disty <-  plot_coefs(Healthy, Busy, Postboxesy, colors = "Dark2", model.names = c("Health Providers", "Bus Stops", "Post Boxes") ,legend.title = "Service", coefs = c("North East" = "ne", "Distance from City" = "distance")) 
ney <- plot_coefs(Health, Bus, Postboxes,colors = "Dark2", model.names = c("Health Providers", "Bus Stops", "Post Boxes") ,legend.title = "Service", omit.coefs = c("poplsoa","size","distance", "(Intercept)", "growth"), coefs = c("North East" = "neNorth East Towns"))
little <- plot_coefs(Nurseries, Schools, GPs, scale = TRUE, model.names = c("Nursery Schools", "Schools", "GPs") ,legend.title = "Service", coefs = c("Improvement" = "growth", "Distance from City" = "distance") ,omit.coefs = c("poplsoa","size","neNorth East Towns", "(Intercept)"))

ney <- ney + labs(x="Expected Change in Service Numbers", title = "North East Compared to other Towns")  + geom_text(aes(x=0, y=1.5, label="Towns in Rest of Britain"), colour= "black")
disty <- disty + labs(x="Expected Change in Service Numbers (per 1km increase in distance from nearest city)", title = "Towns' Location and Service Provision")


try <- lm(services~growth+distance+poplsoa+size, data=data18)

northeast <- northeast %>% mutate(mental=ifelse(northeast$mentals>0, 1, 0))
northeast <- northeast %>% mutate(hospital=ifelse(northeast$hospitals>0, 1, 0))
northeast <- northeast %>% mutate(station=ifelse(northeast$stations>0, 1, 0))
northeast <- northeast %>% mutate(jobcent=ifelse(northeast$jobcents>0, 1, 0))
northeast <- northeast %>% mutate(furthered=ifelse(northeast$furthereds>0, 1, 0))
northeast <- northeast %>% mutate(police=ifelse(northeast$polices>0, 1, 0))
northeast <- northeast %>% mutate(librari=ifelse(northeast$librarys>1, 1, 0))



northeast$mental <- factor( northeast$mental)
northeast$hospital <- factor( northeast$hospital)
northeast$furthered <- factor( northeast$furthered)
northeast$jobcent <- factor( northeast$jobcent)
northeast$station <- factor( northeast$station)
northeast$police <- factor( northeast$police)
northeast$librari <- factor( northeast$librari)
northeast$country <- factor(northeast$country)

##trains##
d <- glm(station~growth+distance+poplsoa+size+housedep+workinggroup, data=northeast, family = "binomial")

z<- plot_model(d, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("housedep [35:80]"),title = "Predicted Probability of Train Stations", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
zz<- plot_model(d, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("distance [0:80]"),title = "Predicted Probability of Train Staions", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
zzz<- plot_model(d,ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("growth [-5:5]"),title = "Predicted Probability of a Train Station", axis.title = c( "Improvement Score", "Probability"), legend.title = "Region")
zzzz<- plot_model(d,ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("workinggroup"),title = "Predicted Probability of a Train Station", axis.title = c( "Town Type", "Probability"), legend.title = "Region")

z <- z+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zz <- zz+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zzz <- zzz+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zzzz <- zzzz+coord_cartesian(ylim=c(0.46,1))

averagecompared <-ggarrange(zzz,z, ncol=2, nrow=2)



comparehealthne <- lm(healths~growth+distance+poplsoa+size+housedep+workinggroup, data=northeast)
comparebusne <- lm(buss~growth+distance+poplsoa+size+housedep+workinggroup, data=northeast)
comparenurseryne <- lm(nurserys~growth+distance+poplsoa+size+housedep+workinggroup, data=northeast)
compareschoolne <- lm(schoolss~growth+distance+poplsoa+size+housedep+workinggroup, data=northeast)
comparegpne <- lm(gpss~growth+distance+poplsoa+size+housedep+workinggroup, data=northeast)
comparepostboxne <- lm(pboxs~growth+distance+poplsoa+size+housedep+workinggroup, data=northeast)
# 
# comparehealthne <- lm(healths~growth+distance+poplsoa+size+housedep+jobdensity, data=data18)
# comparebusne <- lm(buss~growth+distance+poplsoa+size+housedep+jobdensity, data=data18)
# comparenurseryne <- lm(nurserys~growth+distance+poplsoa+size+housedep+jobdensity, data=data18)
# compareschoolne <- lm(schoolss~growth+distance+poplsoa+size+housedep+jobdensity, data=data18)
# comparegpne <- lm(gpss~growth+distance+poplsoa+size+housedep+jobdensity, data=data18)
# comparepostboxne <- lm(pboxs~growth+distance+poplsoa+size+housedep+jobdensity, data=data18)


try <- effect_plot(comparepostboxne,cat.geom="line" , fill="seagreen2", plot.points = FALSE,  pred=growth, interval = TRUE,  x.label = "Improvement Index", y.label = "Predicted Post Boxes")+
  geom_line(colour="seagreen2", linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    title = "Estimated Number of all Post Boxes per Town",
    x = "Improvement Index",
    y = "Predicted number of Post Boxes")+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=50, label="Average Town Improvement"))+
  scale_x_continuous(breaks=c(-2.5,2.5), labels = c("Declining Towns", "Improving Towns"))
  
try <- effect_plot(comparehealthne,cat.geom="line" , pred=growth, interval = TRUE)+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    x = "Improvement Index",
    y = "Predicted number of Health Services")+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=85, label="Average Town Improvement"))+
  scale_x_continuous(breaks=c(-2.5,2.5), labels = c("Declining Towns", "Improving Towns"))
ggsave("Figure9.png", plot=try, width=7.64, height=5.5, dpi=600)
try <- effect_plot(comparehealthne,cat.geom="line" , pred=growth, interval = TRUE)+
  geom_line(colour="navyblue",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("navyblue")+
  scale_color_manual("navyblue")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "blue", colour="blue", size=1)+
  labs(
    x = "Improvement Index",
    y = "Predicted number of Health Services")+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=85, label="Average Town Improvement"))+
  scale_x_continuous(breaks=c(-2.5,2.5), labels = c("Declining Towns", "Improving Towns"))
ggsave("Figure9navy.png", plot=try, width=7.64, height=5.5, dpi=600)

try <- effect_plot(comparehealthne,main.title = NULL ,plot.points = FALSE ,fill="Black",colors = "Black", interval = TRUE,pred=housedep,  point.size = 2 )
try + labs( colour = "black")+
  theme_minimal() +
  
  labs(
  title = "Estimated Number of all Health-related Services per Town",
  x = "Households in Deprivation (%, 2011)",
  y = "Predicted Health providers"
)
  
try <- effect_plot(comparehealthne,cat.geom = "bar" ,main.title = NULL ,plot.points = FALSE ,fill="seagreen2",colors = "seagreen2", interval = FALSE,pred=workinggroup )
try + geom_col(fill="seagreen2")+
  labs(  title = "Estimated Number of all Health-related Services per Town",
         x = "Town Type",
         y = "Predicted number of Health Services")+
  theme_minimal() 
  
  labs(
    title = "Estimated Number of all Health-related Services per Town",
    x = "Households in Deprivation (%, 2011)",
    y = "Predicted Health providers"
  )

  




##trains##
d <- glm(station~growth+distance+poplsoa+size+housedep+ne+workinggroup, data=data18, family = "binomial")

z<- plot_model(d, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("housedep [35:80]", "ne"),title = "Predicted Probability of Train Stations", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
zz<- plot_model(d, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("distance [0:80]", "ne"),title = "Predicted Probability of Train Staions", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
zzz<- plot_model(d,ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("growth [-5:5]", "ne"),title = "Predicted Probability of a Train Station", axis.title = c( "Improvement Score", "Probability"), legend.title = "Region")
zzzz<- plot_model(d,ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("workinggroup", "ne"),title = "Predicted Probability of a Train Station", axis.title = c( "Town Type", "Probability"), legend.title = "Region")

z <- z+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zz <- zz+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zzz <- zzz+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zzzz <- zzzz+coord_cartesian(ylim=c(0.46,1))

averagecompared <-ggarrange(zzz,z,zz,zzzz, ncol=2, nrow=2)

b <- glm(hospital~growth+distance+poplsoa+size+housedep+ne+workinggroup, data=data18, family = "binomial")

z<- plot_model(b, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("housedep [35:80]", "ne"),title = "Predicted Probability of a Hospital", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
zz<- plot_model(b, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("distance [0:80]", "ne"),title = "Predicted Probability of a Hospital", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
zzz<- plot_model(b,ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("growth [-5:5]", "ne"),title = "Predicted Probability of a Hospital", axis.title = c( "Improvement Score", "Probability"), legend.title = "Region")
zzzz<- plot_model(b,ci.lvl = 0,  colors = "Dark2",type = "pred", terms = c("workinggroup", "ne"),title = "Predicted Probability of a Hospital", axis.title = c( "Town Type", "Probability"), legend.title = "Region")

z <- z+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zz <- zz+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zzz <- zzz+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zzzz <- zzzz+coord_cartesian(ylim=c(0.46,1))

averagecompared <-ggarrange(zzzz,zz,z,zzz, ncol=2, nrow=2)
a <- glm(mental~growth+distance+poplsoa+size+housedep+ne+jobdensity, data=data18, family = "binomial")
z<- plot_model(a, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("housedep [35:80]", "ne"),title = "Predicted Probability of a Meantal Health Provider", axis.title = c( "Households in Deprivation (%, 2011)", "Probability"), legend.title = "Region")
zz<- plot_model(a, ci.lvl = 0,colors = "Dark2",type = "pred", terms = c("distance [0:80]", "ne"),title = "Predicted Probability of a Meantal Health Provider", axis.title = c( "Distance from City (Km)", "Probability"), legend.title = "Region")
z <- z+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))
zz <- zz+geom_line(size=2)+coord_cartesian(ylim=c(0.46,1))

averagecompared <-ggarrange(zz,z, ncol=2, nrow=1)


##North West##

northwest <- data18[which(data18$Region=='North West'), ]
northwestpanel <- masterpanel[which(masterpanel$Region=='North West'), ]
northwestlist <- masterlist[which(masterlist$Region=='North West'), ]


growthnorthwest <- ggplot(data=northwest, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in North West Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=74, label="British Town Average"))

growthnorthwestz <- ggplot(data=northwest, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in North West Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))


gpnew <- ggplot(data=northwest, aes(x=growth, y=gpsper, group=BUA11CD))+
  geom_point()+
  labs(x = "Distance from a City (Km)", y = "GPs (per 100 population)", title = "GPs per capita and Distance from a City in North West Towns")+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(nw =ifelse(data18$Region=='North West', 1,0))
data18$nw <- factor(data18$nw)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(nw) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = nw, y = busperpop, color = nw, fill = nw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North West", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = nw, y = healthsper, color = nw, fill = nw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North West", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = nw, y = schoolper, color = nw, fill = nw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North West", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = nw, y = busperpop, color = nw, fill = nw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North West", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = nw, y = pboxsper, color = nw, fill = nw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North West", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = nw, y = gpsper, color = nw, fill = nw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("North West", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)

##Scotland###

scotland <- data18[which(data18$Region=='Scotland'), ]
scotlandpanel <- masterpanel[which(masterpanel$Region=='Scotland'), ]
scotlandlist <- masterlist[which(masterlist$Region=='Scotland'), ]


growthscotland <- ggplot(data=scotland, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in Scottish Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=74, label="British Town Average"))

growthscotlandz <- ggplot(data=scotland, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in Scottish Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

masterpanel <- masterpanel[complete.cases(masterpanel$Region), ]
masterpanel$busper <- (masterpanel$buss/masterpanel$poplsoa)*1000
mastpanels <- masterpanel[which(masterpanel$Region!='London'), ]

scotbustime <- ggplot(data=mastpanels, aes(x=date, y=busper, group=BUA11CD))+
   stat_smooth(aes(group = 1)) + 
  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3) + 
  labs(x = "Date", y= "Bus Stops (per 1000 population)", title="Bus Stops in British Regions' Towns")+ 
  facet_grid(. ~ Region)

scotlibtime <- ggplot(data=mastpanels, aes(x=date, y=librarys, group=BUA11CD))+
  stat_smooth(aes(group = 1)) + 
  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 17, size = 3) + 
  labs(x = "Date", y= "Bus Stops (per 1000 population)", title="Bus Stops in British Regions' Towns")+ 
  facet_grid(. ~ Region)



data18 <- data18 %>% mutate(sc =ifelse(data18$Region=='Scotland', 1,0))
data18$sc <- factor(data18$sc)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(sc) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = sc, y = busperpop, color = sc, fill = sc)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Scotland", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = sc, y = healthsper, color = sc, fill = sc)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Scotland", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = sc, y = schoolper, color = sc, fill = sc)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Scotland", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = sc, y = busperpop, color = sc, fill = sc)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Scotland", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = sc, y = pboxsper, color = sc, fill = sc)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Scotland", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = sc, y = gpsper, color = sc, fill = sc)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Scotland", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)


###wales###

wales <- data18[which(data18$Region=='Wales'), ]
walespanel <- masterpanel[which(masterpanel$Region=='Wales'), ]
waleslist <- masterlist[which(masterlist$Region=='Wales'), ]

growthwales <- ggplot(data=wales, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in Welsh Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=60, label="British Town Average"))

growthwalesz <- ggplot(data=wales, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in Welsh Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(wa =ifelse(data18$Region=='Wales', 1,0))
data18$wa <- factor(data18$wa)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(wa) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = wa, y = busperpop, color = wa, fill = wa)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Wales", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = wa, y = healthsper, color = wa, fill = wa)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Wales", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = wa, y = schoolper, color = wa, fill = wa)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Wales", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = wa, y = busperpop, color = wa, fill = wa)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Wales", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = wa, y = pboxsper, color = wa, fill = wa)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Wales", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = wa, y = gpsper, color = wa, fill = wa)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Wales", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)

wales <- wales%>% mutate(hosp =ifelse(wales$hospitals==0,1,0)) 
wales <- wales%>% mutate(lib =ifelse(wales$librarys==0,1,0)) 
wales <- wales%>% mutate(fir =ifelse(wales$fires==0,1,0)) 
wales <- wales%>% mutate(pol =ifelse(wales$polices==0,1,0)) 
wales <- wales%>% mutate(stat =ifelse(wales$stations==0,1,0)) 
wales <- wales%>% mutate(fur =ifelse(wales$furthereds==0,1,0)) 

wales$missing <- wales$hosp+wales$fir+wales$fur+wales$pol+wales$lib+wales$stat


###South west##
southwest <- data18[which(data18$Region=='South West'), ]
southwestpanel <- masterpanel[which(masterpanel$Region=='South West'), ]
southwestlist <- masterlist[which(masterlist$Region=='South West'), ]

growthsouthwest <- ggplot(data=southwest, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in South West Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=85, label="British Town Average"))

growthsouthwestz <- ggplot(data=southwest, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in South West Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(sw =ifelse(data18$Region=='South West', 1,0))
data18$sw <- factor(data18$sw)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(sw) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = sw, y = busperpop, color = sw, fill = sw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South West", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = sw, y = healthsper, color = sw, fill = sw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South West", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = sw, y = schoolper, color = sw, fill = sw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South West", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = sw, y = busperpop, color = sw, fill = sw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South West", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = sw, y = pboxsper, color = sw, fill = sw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South West", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = sw, y = gpsper, color = sw, fill = sw)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South West", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)


###South East##
southeast <- data18[which(data18$Region=='South East'), ]
southeastpanel <- masterpanel[which(masterpanel$Region=='South East'), ]
southeastlist <- masterlist[which(masterlist$Region=='South East'), ]



growthsoutheast <- ggplot(data=southeast, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in South East Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=60, label="British Town Average"))

growthsoutheastz <- ggplot(data=southeast, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in South East Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(se =ifelse(data18$Region=='South East', 1,0))
data18$se <- factor(data18$se)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(se) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = se, y = busperpop, color = se, fill = se)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South East", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = se, y = healthsper, color = se, fill = se)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South East", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = se, y = schoolper, color = se, fill = se)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South East", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = se, y = busperpop, color = se, fill = se)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South East", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = se, y = pboxsper, color = se, fill = se)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South East", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = se, y = gpsper, color = se, fill = se)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("South East", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)


###East of England###

easteng <- data18[which(data18$Region=='East of England'), ]
eastengpanel <- masterpanel[which(masterpanel$Region=='East of England'), ]
eastenglist <- masterlist[which(masterlist$Region=='East of England'), ]

growtheasteng <- ggplot(data=easteng, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in East of England Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=60, label="British Town Average"))

growtheastengz <- ggplot(data=easteng, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in East of England Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(ee =ifelse(data18$Region=='East of England', 1,0))
data18$ee <- factor(data18$ee)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(ee) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = ee, y = busperpop, color = ee, fill = ee)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East of England", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = ee, y = healthsper, color = ee, fill = ee)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East of England", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = ee, y = schoolper, color = ee, fill = ee)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East of England", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = ee, y = busperpop, color = ee, fill = ee)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East of England", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = ee, y = pboxsper, color = ee, fill = ee)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East of England", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = ee, y = gpsper, color = ee, fill = ee)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East of England", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)

###Yorkshire and the Humber###
yorkshire <- data18[which(data18$Region=='Yorkshire and The Humber'), ]
yorkshirepanel <- masterpanel[which(masterpanel$Region=='Yorkshire and The Humber'), ]
yorkshirelist <- masterlist[which(masterlist$Region=='Yorkshire and The Humber'), ]

growthyorkshires <- ggplot(data=yorkshire, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in Yorkshire and the Humber Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=60, label="British Town Average"))

growthyorkshirez <- ggplot(data=yorkshire, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in Yorkshire and the Humber Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(yo =ifelse(data18$Region=='Yorkshire and The Humber', 1,0))
data18$yo <- factor(data18$yo)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(yo) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = yo, y = busperpop, color = yo, fill = yo)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Yorkshire & \n The Humber", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = yo, y = healthsper, color = yo, fill = yo)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Yorkshire & \n The Humber", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = yo, y = schoolper, color = yo, fill = yo)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Yorkshire & \n The Humber", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = yo, y = busperpop, color = yo, fill = yo)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Yorkshire & \n The Humber", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = yo, y = pboxsper, color = yo, fill = yo)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Yorkshire & \n The Humber", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = yo, y = gpsper, color = yo, fill = yo)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("Yorkshire & \n The Humber", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)

###East Midlands###
eastmidlands <- data18[which(data18$Region=='East Midlands'), ]
eastmidlandspanel <- masterpanel[which(masterpanel$Region=='East Midlands'), ]
eastmidlandslist <- masterlist[which(masterlist$Region=='East Midlands'), ]

growtheastmidlands <- ggplot(data=eastmidlands, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in East Midlands Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=60, label="British Town Average"))

growtheastmidlandsz <- ggplot(data=eastmidlands, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in East Midlands Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  scale_y_continuous(breaks = c(0,50000, 100000), labels = c("0", "50000", "100000"))+
  coord_cartesian(ylim=c(0,140000))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(em =ifelse(data18$Region=='East Midlands', 1,0))
data18$em <- factor(data18$em)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(em) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = em, y = busperpop, color = em, fill = em)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East Midlands", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = em, y = healthsper, color = em, fill = em)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East Midlands", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = em, y = schoolper, color = em, fill = em)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East Midlands", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = em, y = busperpop, color = em, fill = em)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East Midlands", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = em, y = pboxsper, color = em, fill = em)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East Midlands", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = em, y = gpsper, color = em, fill = em)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("East Midlands", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)

###West Midlands###

westmidlands <- data18[which(data18$Region=='West Midlands'), ]
westmidlandspanel <- masterpanel[which(masterpanel$Region=='West Midlands'), ]
westmidlandslist <- masterlist[which(masterlist$Region=='West Midlands'), ]
growthwestmidlands <- ggplot(data=westmidlands, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)", title = "Improvement or Decline in West Midlands Towns and Distance from City")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=60, label="British Town Average"))

growthwestmidlandsz <- ggplot(data=westmidlands, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in West Midlands Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel()+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

data18 <- data18 %>% mutate(wm =ifelse(data18$Region=='West Midlands', 1,0))
data18$wm <- factor(data18$wm)

data18 <- data18[complete.cases(data18$bua11nm), ]


gd <- data18 %>% 
  group_by(wm) %>% 
  summarise(
    busperpop = mean(busperpop),
    healthsper = mean(healthsper),
    schoolper = mean(schoolper),
    pboxsper = mean(pboxsper),
    nurseryper = mean(nurseryper),
    schoolper = mean(schoolper),
    gpsper = mean(gpsper),
    hallper = mean(hallper)
    
  )

busaverage <- ggplot(data18, aes(x = wm, y = busperpop, color = wm, fill = wm)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  # geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("West Midlands", "Rest of Britain"))+
  labs(
    title = "Bus Stops",
    x = "Region",
    y = "Service per 1000 Population"
  )

healthaverage <- ggplot(data18, aes(x = wm, y = healthsper, color = wm, fill = wm)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("West Midlands", "Rest of Britain"))+
  labs(
    title = "Health Providers",
    x = "Region",
    y = "Service per 1000 Population"
  )

Schoolaverage <- ggplot(data18, aes(x = wm, y = schoolper, color = wm, fill = wm)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("West Midlands", "Rest of Britain"))+
  labs(
    title = "Schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

nurseryaverage <- ggplot(data18, aes(x = wm, y = busperpop, color = wm, fill = wm)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("West Midlands", "Rest of Britain"))+
  labs(
    title = "Nursery schools",
    x = "Region",
    y = "Service per 1000 Population"
  )

postboxesper <- ggplot(data18, aes(x = wm, y = pboxsper, color = wm, fill = wm)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("West Midlands", "Rest of Britain"))+
  labs(
    title = "Post Boxes",
    x = "Region",
    y = "Service per 1000 Population"
  )

gpaverage <- ggplot(data18, aes(x = wm, y = gpsper, color = wm, fill = wm)) +
  geom_bar(data = gd, stat = "identity", alpha = .3) +
  geom_point() +
  #geom_text_repel(aes(label=ifelse(ne==1, paste0(name), "")))+
  guides(color = "none", fill = "none") +
  theme_bw() +
  scale_x_discrete(breaks=c(1,0), labels = c("West Midlands", "Rest of Britain"))+
  labs(
    title = "GPs",
    x = "Region",
    y = "Service per 1000 Population"
  )


averagecompared <- ggarrange(Schoolaverage, nurseryaverage, healthaverage, postboxesper, gpaverage, busaverage, ncol=3, nrow=2)






###overall trends####


h <- lm(healths~growth+distance+poplsoa+size+housedep+workinggroup, data=data18)
b <- lm(buss~growth+distance+poplsoa+size+housedep+workinggroup, data=data18)
n <- lm(nurserys~growth+distance+poplsoa+size+housedep+workinggroup, data=data18)
s <- lm(schoolss~growth+distance+poplsoa+size+housedep+workinggroup, data=data18)
g <- lm(gpss~growth+distance+poplsoa+size+housedep+workinggroup, data=data18)
p <- lm(pboxs~growth+distance+poplsoa+size+housedep+workinggroup, data=data18)

hh <- effect_plot(p,cat.geom="line" , pred=growth, interval = TRUE)+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    title = "Estimated Number of Post Boxes per Town",
    x = "Improvement Index",
    y = "Predicted number of GPS")+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_text(aes(x=0, y=50, label="Average Town Improvement"))+
  scale_x_continuous(breaks=c(-2.5,2.5), labels = c("Declining Towns", "Improving Towns"))

pp <- effect_plot(p,cat.geom="line" , pred=housedep, interval = TRUE)+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    title = "Estimated Number of Post Boxes per Town",
    x = "Household Deprivation (%, 2016)",
    y = "Predicted number of Postboxes")





heal <- ggplot(data = data18, aes(x = distance, y = healthsper, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Health providers (per 1000 population)", title="Health Providers and Distance from City")
bu <- ggplot(data = data18, aes(x = distance, y = busperpop, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Distance from nearest City (km)", y= "Bus Stops (per 1000 population)", title="Bus Stops and Distance from City")
wif <- ggplot(data = data18, aes(x = distance, y = wifiper, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Distance from nearest City (km)", y= "Wifi Hotspots (per 1000 population)", title="Wifi Hotspots and Distance from City")
nurs <- ggplot(data = data18, aes(x = distance, y = nurseryper, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Nurseries (per 1000 population)", title="Nurseries and Distance from City")
po <- ggplot(data = data18, aes(x = distance, y = pboxsper, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Post Boxes (per 1000 population)", title="Post Boxes and Distance from City")
scho <- ggplot(data = data18, aes(x = distance, y = schoolper, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Schools (per 1000 population)", title="Schools and Distance from City")
hal <- ggplot(data = data18, aes(x = distance, y = hallper, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Distance from nearest City (km)", y= "Community Halls (per 1000 population)", title="Community Halls and Distance from City")
ggarrange(heal, bu, wif, nurs, scho, po, hal, ncol=2, nrow=4)

heal <- ggplot(data = data18, aes(x = growth, y = healthsper, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Health providers (per 1000 population)", title="Health Providers and Economic Growth") + scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
bu <- ggplot(data = data18, aes(x = growth, y = busperpop, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Growth Index", y= "Bus Stops (per 1000 population)", title="Bus Stops and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
wif <- ggplot(data = data18, aes(x = growth, y = wifiper, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+  labs(x = "Growth Index", y= "Wifi Hotspots (per 1000 population)", title="Wifi Hotspots and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
nurs <- ggplot(data = data18, aes(x = growth, y = nurseryper, group = BUA11CD))+geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Nurseries (per 1000 population)", title="Nurseries and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
po <- ggplot(data = data18, aes(x = growth, y = pboxsper, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Post Boxes (per 1000 population)", title="Post Boxes and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
scho <- ggplot(data = data18, aes(x = growth, y = schoolper, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Schools (per 1000 population)", title="Schools and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
hal <- ggplot(data = data18, aes(x = growth, y = hallper, group = BUA11CD))+ geom_point()+ stat_smooth(method = 'lm', se=FALSE,aes(group = 1))+ labs(x = "Growth Index", y= "Community Halls (per 1000 population)", title="Community Halls and Economic Growth")+ scale_x_continuous(breaks=c(-3,3),labels=c("Declining Towns", "Growing Towns") ) +theme(axis.line.x=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
ggarrange(heal, bu, wif, nurs, scho, po, hal, ncol=2, nrow=4)


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
data18$librari <- factor( data18$librari)
data18$country <- factor(data18$country)

h <- glm(mental~growth+distance+poplsoa+size+country, data=data18, family = "binomial")
h2 <- margins(h)

f <- glm(hospital~growth+distance+poplsoa+size+country, data=data18, family = "binomial")
f2 <- margins(f)

d <- glm(jobcent~growth+distance+poplsoa+size+country, data=data18, family = "binomial")
d2 <- margins(d)

s <- glm(station~growth+distance+poplsoa+size+country, data=data18, family = "binomial")
s2 <- margins(s)

a <- glm(furthered~growth+distance+poplsoa+size+country, data=data18, family = "binomial")
a2 <- margins(a)

q <- glm(police~growth+distance+poplsoa+size+country, data=data18, family = "binomial")
q2 <- margins(q)

p <- glm(librari~growth+distance+poplsoa+size+country, data=data18, family = "binomial")
p2 <- margins(q)

fine <- stargazer(q, a, s, d, f, h, p, style= "apsr", title="Results", align=TRUE, type="text")

plot_model(s, type = "pred", terms = "growth [-5:5]")
plot_model(f, type = "pred", terms = "growth [-5:5]")
plot_model(d, type = "pred", terms = "growth [-5:5]")




h2 <- summary(h2)
f2 <- summary(f2)
d2 <- summary(d2)
s2 <- summary(s2)
a2 <- summary(a2)
q2 <- summary(q2)
p2 <- summary(p2)

w <- ggplot(data = h2) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  coord_cartesian(ylim=c(-0.05,0.11)) +
  ggtitle("Mental Health")

x <- ggplot(data = f2) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  coord_cartesian(ylim=c(-0.05,0.11)) +
  ggtitle("Hospital")

y <- ggplot(data = d2) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  coord_cartesian(ylim=c(-0.05,0.11)) +
  ggtitle("Jobcentre")

z <- ggplot(data = s2) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  coord_cartesian(ylim=c(-0.05,0.11)) +
  ggtitle("Train Station")

t <- ggplot(data = a2) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  coord_cartesian(ylim=c(-0.05,0.11)) +
  ggtitle("Further Education")

r <- ggplot(data = q2) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  coord_cartesian(ylim=c(-0.05,0.11)) +
  ggtitle("GP")
m <- ggplot(data = p2) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  coord_cartesian(ylim=c(-0.05,0.11)) +
  ggtitle("GP")


ggarrange(w, x, y, z, t, r, m, ncol=4, nrow=2)



try <- ggplot(data = masterpanel, aes(x = date, y = healthsper, group = BUA11CD))
healthplot <- try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                          geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Health per person", title="Health Providers per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
try <- ggplot(data = masterpanel, aes(x = date, y = schoolper, group = BUA11CD))
schoolplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                         geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Schools per person", title="Schools per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = nurseryper, group = BUA11CD))
nurseryplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                          geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Nurseries per person", title="Children's Nurseries per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = pboxsper, group = BUA11CD))
pboxplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Post Boxes per person", title="Post Box per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = busperpop, group = BUA11CD))
busplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                      geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Bus stops per person", title="Bus stops per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = hallper, group = BUA11CD))
hallplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Community Centres/ Halls per person", title="Community Centres/ Halls per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = wifiper, group = BUA11CD))
wifiplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Wifi Hotspots per person", title="Wifi Hotspots per Person")+ facet_grid(. ~ declinegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

ggarrange(wifiplot, hallplot, busplot, pboxplot, nurseryplot, schoolplot, healthplot, ncol=4, nrow=2)


try <- ggplot(data = masterpanel, aes(x = date, y = healthsper, group = BUA11CD))
healthplot <- try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                          geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Health Providers per person", title="Health Providers per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))
try <- ggplot(data = masterpanel, aes(x = date, y = schoolper, group = BUA11CD))
schoolplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                         geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Schools per person", title="Schools per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = nurseryper, group = BUA11CD))
nurseryplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                          geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Nurseries per person", title="Nurseries per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = pboxsper, group = BUA11CD))
pboxplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Post Boxes per person", title="Post Box per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = busperpop, group = BUA11CD))
busplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                      geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Bus stops per person", title="Bus stops per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = hallper, group = BUA11CD))
hallplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Community Centres/ Halls per person", title="Community Centres/ Halls per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

try <- ggplot(data = masterpanel, aes(x = date, y = wifiper, group = BUA11CD))
wifiplot <-try + stat_smooth(method = 'loess', span = 1,aes(group = 1)) + stat_summary(aes(group = 1),
                                                                                       geom = "point", fun.y = mean, shape = 17, size = 3) + labs(x = "Date", y= "Number of Wifi Hotspots per person", title="Wifi Hotspots per Person")+ facet_grid(. ~ distancegroup)+ geom_hline(yintercept=0)+ scale_y_continuous(breaks=c(-0.5,0.0,0.5),labels=c(expression("Lower than Average          -0.5"*sigma), expression("Average            0"*sigma),expression("Higher than average       0.5"*sigma) ) )+ coord_cartesian(ylim=c(-0.5,0.5))+theme(axis.line.y=element_line(arrow = arrow(length = unit(3, 'mm'), ends="both")))

ggarrange(wifiplot, hallplot, busplot, pboxplot, nurseryplot, schoolplot, healthplot, ncol=4, nrow=2)








####pubs and post offices####



pub196 <- read_delim("pubs and bars.csv", 
                           "|", escape_double = FALSE, trim_ws = TRUE)
pub196 <- pub196[complete.cases(pub196$`Feature Easting`), ]
names(pub196)[names(pub196)=="PointX Classification Code"] <- "code"
pub196$code <- as.numeric(pub196$code)


pub196 <- pub196[which(pub196$code==1020034), ]

pub196 <- subset(pub196, "Feature Easting" != "" | "Feature Northing" != "")
pub196$poi_ID <- 1:nrow(pub196)
coords <- cbind(Easting = as.numeric(as.character(pub196$`Feature Easting`)),
                Northing = as.numeric(as.character(pub196$`Feature Northing`)))


pub196 <- SpatialPointsDataFrame(coords, data = data.frame(pub196$Name ,
                                                                 pub196$poi_ID), proj4string = CRS("+init=epsg:27700"))


pub196 <- over(uktowns, geometry(pub196), returnList = TRUE)
pub196 <- sapply(pub196, length)
pub196 <- spCbind(uktowns, pub196)


pub196 <- merge(pub196, popbua, by="BUA11CD")

pub196 <- pub196[which(pub196$poplsoa>10000& pub196$poplsoa<175000), ]

pub196 <- merge(pub196, region, by="BUA11CD")

pub196 <- pub196[which(pub196$Region!="Scotland"), ]




post196 <- read_delim("postoffices.csv", 
                     "|", escape_double = FALSE, trim_ws = TRUE)
post196 <- post196[complete.cases(post196$`Feature Easting`), ]
names(post196)[names(post196)=="PointX Classification Code"] <- "code"
post196$code <- as.numeric(post196$code)


post196 <- post196[which(post196$code==9480763), ]

post196 <- subset(post196, "Feature Easting" != "" | "Feature Northing" != "")
post196$poi_ID <- 1:nrow(post196)
coords <- cbind(Easting = as.numeric(as.character(post196$`Feature Easting`)),
                Northing = as.numeric(as.character(post196$`Feature Northing`)))


post196 <- SpatialPointsDataFrame(coords, data = data.frame(post196$Name ,
                                                           post196$poi_ID), proj4string = CRS("+init=epsg:27700"))


post196 <- over(uktowns, geometry(post196), returnList = TRUE)
post196 <- sapply(post196, length)
post196 <- spCbind(uktowns, post196)


post196 <- merge(post196, popbua, by="BUA11CD")

post196 <- post196[which(post196$poplsoa>10000& post196$poplsoa<175000), ]

post196 <- merge(post196, region, by="BUA11CD")

post196 <- post196[which(post196$Region!="Scotland"), ]


postandpub <- merge(post196, pub196, by="BUA11CD")



postandpubs <- merge(postandpub, uktowns, by="BUA11CD")



data18 <- data18 %>% mutate(decline19 = ifelse(data18$growth< -1.07 , "Declining", ifelse(data18$growth> 0.64, "Improving", "Stagnant")))
data18 <- data18 %>% mutate(deprived19 = ifelse(data18$housedep< 53.1 , "Affluent", ifelse(data18$housedep> 60, "Deprived", "Middling")))


postandpubz <- merge(postandpubs, data18, by="BUA11CD")

postandpubzz <- as.data.frame(postandpubz)

bins <- c(0,2,5,Inf)
palbins <- c(0,10,24,Inf)
palpost <- colorBin("YlOrRd", domain = postandpubz$post196, bins=bins)
palpub <- colorBin("YlOrRd", domain = postandpubz$pub196, bins=palbins)

postandpubz <- spTransform(postandpubz, CRS("+init=epsg:4326"))
postandpubz <- postandpubz[which(postandpubz$Region.x!="Wales"),]

  
# mymap <- leaflet(data=postandpubz) %>%
#   addTiles() %>%  
#   setView(lng=-1.3823, lat=53.0975, zoom=5)
# addPolygons(mymap, data = postandpubz,  fillColor = ~palpost(post196), weight= 2, opacity = 1, color = "black", dashArray = "3", fillOpacity = 0.8, highlight = highlightOptions(weight=5, color = "#666", dashArray = "", fillOpacity = "0.9", bringToFront = TRUE), group = "Post Office")%>%
#   addPolygons(mymap,data = postandpubz, fillColor = ~palpub(pub196), weight= 2, opacity = 1, color = "black", dashArray = "3", fillOpacity = 0.8, highlight = highlightOptions(weight=5, color = "#666", dashArray = "", fillOpacity = "0.9", bringToFront = TRUE), group = "Pubs, Inns and Bars")%>%
#   addLayersControl(
#     baseGroups = c("Post Offices", "Pubs, Inns and Bars"), data=postandpubz,
#     overlayGroups = "Legends", position= "topright",
#     options = layersControlOptions(collapsed = FALSE)
#   ) 

labels <- sprintf(
  "<strong>%s</strong><br/>%g Post Offices ",
  postandpubz$name, postandpubz$post196
) %>% lapply(htmltools::HTML)

labelspub <- sprintf(
  "<strong>%s</strong><br/>%g Pubs, Inns or Bars ",
  postandpubz$name, postandpubz$pub196
) %>% lapply(htmltools::HTML)


# map2 <- leaflet(postandpubz) %>%
#   setView(lng=-1.3823, lat=53.0975, zoom=8)%>%
#   addTiles()%>%
#   addPolygons(
#     fillColor = ~palpost(post196),
#     weight = 1,
#     opacity = 1,
#     color = "white",
#     dashArray = "3",
#     fillOpacity = 0.7,
#     highlight = highlightOptions(
#       weight = 5,
#       color = "#666",
#       dashArray = "",
#       fillOpacity = 0.7,
#       bringToFront = TRUE),
#     label = labels,
#     labelOptions = labelOptions(
#       style = list("font-weight" = "normal", padding = "3px 8px"),
#       textsize = "15px",
#       direction = "auto"), group = "Post Offices") %>%
#   addPolygons(
#     fillColor = ~palpub(pub196),
#     weight = 1,
#     opacity = 1,
#     color = "white",
#     dashArray = "3",
#     fillOpacity = 0.7,
#     highlight = highlightOptions(
#       weight = 5,
#       color = "#666",
#       dashArray = "",
#       fillOpacity = 0.7,
#       bringToFront = TRUE),
#     label = labelspub,
#     labelOptions = labelOptions(
#       style = list("font-weight" = "normal", padding = "3px 8px"),
#       textsize = "15px",
#       direction = "auto"), group = "Pubs, Inns and Bars") %>%
#       addLayersControl(
#         baseGroups = c("Post Offices", "Pubs, Inns and Bars"), data=postandpubz,
#         overlayGroups = "Legends", position= "topright",
#         options = layersControlOptions(collapsed = FALSE)
#       ) 




getColor <- function(postandpubz) {
  sapply(postandpubz$pub196, function(pub196) {
    if(pub196 >= 20) {
      "lightblue"
    } else if(pub196 >= 10) {
      "blue"
    } else {
      "darkblue"
    } })
}


getColorz <- function(postandpubz) {
  sapply(postandpubz$post196, function(post196) {
    if(post196 >= 5) {
      "lightblue"
    } else if(post196 >= 2) {
      "blue"
    } else {
      "darkblue"
    } })
}



iconz <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorz(postandpubz))

iconp <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColor(postandpubz))

postandpubz <- SpatialPointsDataFrame(gCentroid(postandpubz, byid=TRUE), postandpubz@data, match.ID=FALSE)

# postandpubz <- postandpubz[which(postandpubz$Region.x!="Wales"),]



mymap <- leaflet(data=postandpubz) %>%
  addTiles() %>%  
  setView(lng=-1.3823, lat=53.0975, zoom=7) 
  mymap <- addLegend(mymap, group = "Legends and Data Source",
          position = "bottomright",
          colors = c("lightblue", "blue", "darkblue", ""),
          labels = c("More than Five", "More than Two", "Fewer than Two", "Data is provided by Ordnance Survey as of June 2019: © Crown copyright and database rights 2019 Ordnance Survey (100025252). This material includes data licensed from PointX© Database Right/Copyright 2019."), opacity = 1,
          title = "Number of Post Offices") 
  mymap <-   addLegend(mymap, group = "Legends and Data Source",
            position = "bottomright",
            colors = c("lightblue", "blue", "darkblue", ""),
            labels = c("More than Twenty", "More than Ten", "Fewer than Ten", "Data is provided by Ordnance Survey as of June 2019: © Crown copyright and database rights 2019 Ordnance Survey (100025252). This material includes data licensed from PointX© Database Right/Copyright 2019."), opacity = 1,
            title = "Number of Pubs, Inns and Bars") 

mymap <-  addAwesomeMarkers(mymap, data=postandpubz,icon=iconz,
                  popup = labels, label= labels, group = "Post Offices")
  mymap <- addAwesomeMarkers(mymap, data=postandpubz, icon=iconp,
                     popup = labelspub, label= labelspub, group = "Pubs, Inns and Bars") 

mymap <- addLayersControl(mymap,
    baseGroups = c("Post Offices", "Pubs, Inns and Bars"),
   position= "topright", overlayGroups = "Legends and Data Source",
    options = layersControlOptions(collapsed = FALSE))
# 
# tag.map.data <- tags$style(HTML(".leaflet-control.map-title { transform: translate(-50%, 20%); position: fixed !important; left: 50%; text-align: center; paddling-left:2px; padding-right: 2px; background: rgba(255,255,255,0.75); font-weight: bold; fontsize: 0.5px;}"))
# 
# # rr <- tags$div(HTML('<a href="https://cran.r-project.org/"> <img border="0" alt="Data is provided by Ordnance Survey as of June 2019: © Crown copyright and database rights 2019 Ordnance Survey (100025252). This material includes data licensed from PointX© Database Right/Copyright 2019." src="/PathToImageR.jpeg" width="300", height="100"> </a>'))
# 
# title <- tags$div(tag.map.data, HTML("Data is provided by Ordnance Survey as of June 2019: © Crown copyright and database rights 2019 Ordnance Survey (100025252). This material includes data licensed from PointX© Database Right/Copyright 2019.", size = 0.5))
# title2 <- tags$p(tags$style("p {color: black; font-size: 10px; padding-left:2px; padding-right: 2px}"), tags$b("Data is provided by Ordnance Survey as of June 2019: © Crown copyright and database rights 2019 Ordnance Survey (100025252). <br> This material includes data licensed from PointX© Database Right/Copyright 2019."))
# 
# mymaps <- addControl(mymap, title2, position = "bottomleft", layerId = "Data", className = "title2")
# 
# mymapz <- addLayersControl(mymaps,
#                           baseGroups = c("Post Offices", "Pubs, Inns and Bars"),
#                           position= "topright", overlayGroups = c("Legends", "Data"),
#                           options = layersControlOptions(collapsed = FALSE))

    
mymap %>% hideGroup("Legends and Data Source")

postandpubzz$deprived19 <- factor(postandpubzz$deprived19, levels = c("Deprived", "Middling", "Affluent"))
postandpubzz$decline19 <- factor(postandpubzz$decline19, levels = c("Declining", "Stagnant", "Improving"))
postandpubzz <- postandpubzz[which(postandpubzz$Region.x!="Wales"),]
postandpubzz <- postandpubzz[complete.cases(postandpubzz$deprived19),]


x <- ggplot(postandpubzz)+
geom_bar( aes(decline19, post196, fill=as.factor(decline19)),
       stat= "summary", fun.y="mean") +
  theme_tufte()+
  xlab("Decline or Improvement")+
  ylab("Average number of Post Offices")+
  ggtitle("Average Number of Post Offices in Declining or Improving Towns")+
  theme(legend.position = "none")
ggsave(plot=x,"graphelection4.png", width=6, height=5.5, dpi=600)


y <- ggplot(postandpubzz)+
  geom_bar( aes(deprived19, post196, fill=as.factor(deprived19)),
            stat= "summary", fun.y="mean") +
  theme_tufte()+
  xlab("Decline or Improvement")+
  ylab("Average number of Post Offices")+
  ggtitle("Average Number of Post Offices in Deprived or Aflluent Towns")+
  theme(legend.position = "none")
ggsave(plot=y,"graphelection3.png", width=6, height=5.5, dpi=600)

z <- ggplot(postandpubzz)+
  geom_bar( aes(decline19, pub196, fill=as.factor(decline19)),
            stat= "summary", fun.y="mean") +
  theme_tufte()+
  xlab("Decline or Improvement")+
  ylab("Average number of Pubs")+
  ggtitle("Average Number of Pubs in Declining or Improving Towns")+
  theme(legend.position = "none")
ggsave(plot=z,"graphelection2.png", width=6, height=5.5, dpi=600)


j <- ggplot(postandpubzz)+
  geom_bar( aes(deprived19, pub196, fill=as.factor(deprived19)),
            stat= "summary", fun.y="mean") +
  theme_tufte()+
  xlab("Deprived or Affluent")+
  ylab("Average number of Pubs")+
  ggtitle("Average Number of Pubs in Deprived or Aflluent Towns")+
  theme(legend.position = "none")
  
ggsave(plot=j,"graphelection1.png", width=6, height=5.5, dpi=600)


postandpubzz$postper <- postandpubzz$post196/(postandpubzz$poplsoa.x/10000)
postandpubzz$pubper <- postandpubzz$pub196/(postandpubzz$poplsoa.x/10000)

z <- ggplot(postandpubzz)+
  geom_bar( aes(deprived19, pubper, fill=as.factor(decline19)),
            stat= "summary", fun.y="mean") +
  theme_tufte()+
  xlab("Decline or Improvement")+
  ylab("Average number of Pubs")+
  ggtitle("Average Number of Pubs in Declining or Improving Towns")+
  theme(legend.position = "none")

