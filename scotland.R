###load the packages###

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


##Scotland####


nuts1 <- readOGR("C:\\Users\\bjg55\\Documents\\Towns\\2019Ring", "nuts1")
nuts1 <- nuts1[which(nuts1$objectid!=12), ]

nuts1 <- spTransform(nuts1, CRS("+init=epsg:27700"))
nuts1s <- as.data.frame(nuts1)
nuts1s <- nuts1s %>% mutate(SC = ifelse(nuts1s$objectid==11, 1, 0))
nuts1ss <- merge(nuts1, nuts1s, by="objectid")
# nuts1ss$NE <- factor(nuts1ss$NE)

map <- spplot(nuts1ss, "SC", col.regions = c("white", "seagreen2"), colorkey=FALSE)


scotland <- data18[which(data18$Region=='Scotland'), ]
scotlandpanel <- masterpanel[which(masterpanel$Region=='Scotland'), ]
scotlandlist <- masterlist[which(masterlist$Region=='Scotland'), ]


growthscotland <- ggplot(data=scotland, aes(x=growth, y=distance, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Distance from Nearest City (km)")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black")+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_text(aes(x=0, y=134, label="British Town Average"))
#ggsave("Figure1scotversiontwo.png", plot=growthscotland, width=18.6, height=15, dpi=600)

growthscotlandz <- ggplot(data=scotland, aes(x=growth, y=poplsoa, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Improvement Index", y = "Population Size", title = "Improvement or Decline in Scottish Towns and Population Size")+
  scale_x_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
  stat_smooth(method = 'lm', se=FALSE,aes(group = 1))

growthscotlandservice <- ggplot(data=scotland, aes(x=services, y=growth, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Public Service Index", y = "Improvement Index", title = "Improvement and Level of Public Services")+
  scale_y_continuous(breaks=c(-2,2), labels = c("Declining Towns", "Improving Towns"))+
  scale_x_continuous(breaks=c(-3,3), labels = c("Lower Service Provision", "Higher Service Provision"))+
  theme(axis.line.y = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "seagreen2", color = "white", fontface = "bold", segment.color = "Black")+
  geom_text(aes(y=0, x=11, label="British Town Average"), size=3, angle=270)+
  geom_text(aes(y=6.7, x=0, label="British Town Average"), size=3)

growthnortheastservicechange <- ggplot(data=scotland, aes(x=services, y=publicchange, group=BUA11CD, label=name))+
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
  geom_text(aes(y=8, x=0, label="British Town Average"), size=3)
#ggsave("Figure2scot.png", plot=growthnortheastservicechange, width=14.2, height=12.2, dpi=600)

scotjobs <- ggplot(data=scotland, aes(x=jobdensity, y=housedep, group=BUA11CD, label=name))+
  geom_point()+
  labs(x = "Job Density (2016)", y = "Household Deprivation (%, 2011)")+
  geom_vline(xintercept = 0.5, linetype="dashed", color="grey", size=0.8)+
  geom_vline(xintercept = 0.7, linetype="dashed", color="grey", size=0.8)+
  geom_hline(yintercept = 56.726, linetype="dashed", color="grey", size=0.8)+
  geom_label_repel(fill = "navyblue", color = "white", fontface = "bold", segment.color = "Black")+
  geom_text(aes(x=0.4, y=74, label="Residential Towns"), size=4)+
  geom_text(aes(x=0.6, y=74, label="Partially Residential Towns"), size=4)+
  geom_text(aes(x=1.25, y=56.726, label="British Town Average"), size=4, angle=270)+
  geom_text(aes(x=0.8, y=74, label="Working Towns"), size=4)
#ggsave("Figure3scot.png", plot=scotjobs, width=15, height=14, dpi=600)



data18 <- data18 %>% mutate(sc =ifelse(data18$Region=='Scotland', 1,0))
data18$sc <- factor(data18$sc)

data18 <- data18[complete.cases(data18$bua11nm), ]

###Scottish relationships###
h <- lm(healths~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland)
b <- lm(buss~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland)
n <- lm(nurserys~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland)
s <- lm(schoolss~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland)
g <- lm(gpss~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland)
p <- lm(pboxs~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland)

hh <- effect_plot(g,cat.geom="line" , pred=growth, interval = TRUE)+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    title = "Estimated Number of Doctor's Surgeries per Town",
    x = "Improvement Index",
    y = "Predicted number of GPS")+
  theme(axis.line.x = element_line(arrow = arrow(length = unit(3, 'mm'), ends='both')))+
  geom_vline(xintercept = 0, linetype="dashed", color="black", size=0.8)+
  geom_text(aes(x=0, y=7, label="Average Town Improvement"))+
  scale_x_continuous(breaks=c(-2.5,2.5), labels = c("Declining Towns", "Improving Towns"))

pp <- effect_plot(p,cat.geom="line" , pred=housedep, interval = TRUE)+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    title = "Estimated Number of all Post Boxes per Town",
    x = "Household Deprivation (%, 2016)",
    y = "Predicted number of Postboxes")

gg <- effect_plot(g,cat.geom="line" , pred=housedep, interval = TRUE)+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    title = "Estimated Number of all Doctor's Surgeries per Town",
    x = "Household Deprivation (%, 2016)",
    y = "Predicted number of GPs")

nn <- effect_plot(n,cat.geom="line" , pred=housedep, interval = TRUE)+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)+
  labs(
    title = "Estimated Number of all Nursery Schools per Town",
    x = "Household Deprivation (%, 2016)",
    y = "Predicted number of Nursery Schools")

yep <- ggarrange(gg,pp,nn, ncol=1, nrow=3)

scotland <- scotland %>% mutate(mental=ifelse(scotland$mentals>0, 1, 0))
scotland <- scotland %>% mutate(hospital=ifelse(scotland$hospitals>0, 1, 0))
scotland <- scotland %>% mutate(station=ifelse(scotland$stations>0, 1, 0))
scotland <- scotland %>% mutate(jobcent=ifelse(scotland$jobcents>0, 1, 0))
scotland <- scotland %>% mutate(furthered=ifelse(scotland$furthereds>0, 1, 0))
scotland <- scotland %>% mutate(police=ifelse(scotland$polices>0, 1, 0))



scotland$mental <- factor( scotland$mental)
scotland$hospital <- factor( scotland$hospital)
scotland$furthered <- factor( scotland$furthered)
scotland$jobcent <- factor( scotland$jobcent)
scotland$station <- factor( scotland$station)
scotland$police <- factor( scotland$police)
scotland$country <- factor(scotland$country)



a <- glm(mental~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland, family = "binomial")

b <- glm(hospital~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland, family = "binomial")

c <- glm(jobcent~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland, family = "binomial")

d <- glm(station~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland, family = "binomial")

e <- glm(furthered~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland, family = "binomial")

f <- glm(police~growth+distance+poplsoa+size+housedep+workinggroup, data=scotland, family = "binomial")

w<- plot_model(d ,type = "pred",  terms = c("distance"), title = "Predicted Probability of a Town-based Train Station", axis.title = c( "Distance from Nearest City (km)", "Probability of Present Train Station"), legend.title = "Region")+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)

y<- plot_model(b,type = "pred", colors = "seagreen2", terms = c("distance"), title = "Predicted Probability of a Town-based Hospital", axis.title = c( "Distance from Nearest City (km)", "Probability of Present Hospital"), legend.title = "Region")+
  geom_line(colour="seagreen2",linetype = "solid", size=2)+
  theme_minimal()+
  scale_fill_manual("seagreen2")+
  scale_color_manual("seagreen2")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill = "seagreen1", colour="seagreen3", size=1)

###Comparing Scotland to Britain###

data18 <- data18 %>% mutate(sc =ifelse(data18$Region=='Scotland', "Scottish Towns","Rest of Britain"))
data18$sc <- factor(data18$sc, levels = c("Rest of Britain", "Scottish Towns"))

h <- lm(healths~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18)
b <- lm(buss~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18)
n <- lm(nurserys~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18)
s <- lm(schoolss~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18)
g <- lm(gpss~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18)
p <- lm(pboxs~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18)

output<-tab_model(h,b,n,s,g, p, show.ci=0.95,title="OLS Public Services" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Health Services", "Bus Stops", "Nurseries", "Schools","GPs",  "Post Boxes"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"Scotland Dummy", "Job Density"))
tab_model(h,b,n,s,g, p, show.ci=0.95,title="OLS Public Services" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Health Services", "Bus Stops", "Nurseries", "Schools","GPs",  "Post Boxes"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"Scotland Dummy", "Job Density"), file = "table4scot")


# a <- plot_model(h,dot.size = 6,colors = "Dark2", terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey", axis.lim = c(-15,15), axis.title = "Relative Number of Health Providers", title = "Health Providers")
# b <- plot_model(b, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey", axis.lim = c(-40,40), axis.title = "Relative Number of Bus Stops", title = "Bus Stops")
# c <- plot_model(n,dot.size = 6,colors = "seagreen2", terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-4,4), axis.title = "Relative Number of Nursery Schools", title = "Nursery Schools")
# d <- plot_model(s, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-3,3), axis.title = "Relative Number of Schools", title = "Schools")
# e <- plot_model(g, dot.size = 6,colors = "Dark2",terms = c("ne [North East Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE,vline.color = "Grey",axis.lim = c(-1,1), axis.title = "Relative Number of Doctor's Surgeries", title = "Doctor's Surgeries")
f <- plot_model(p, dot.size = 6, colors = "navyblue",terms = c("sc [Scottish Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-8,8), axis.title = "Relative Number of Post Boxes", title = "Post Boxes")
ff <- plot_model(n, dot.size = 6, colors = "navyblue",terms = c("sc [Scottish Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-8,8), axis.title = "Relative Number of Nursery Schools", title = "Nursery Schools")
ffff <- plot_model(h, dot.size = 6, colors = "navyblue",terms = c("sc [Scottish Towns]"), ci.lvl = 0, show.values = TRUE, show.p = TRUE, vline.color = "Grey",axis.lim = c(-8,8), axis.title = "Relative Number of Health Services", title = "Health-related Services")

# a <- a + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
# b <- b + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
# c <- c + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
ffff <- ffff + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
ff <- ff + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")
f <- f + geom_text(aes(x=1.5, y=0, label="Towns in Rest of Britain"), colour= "black")


averagecompared <- ggarrange(ff,ffff,f,   ncol=1, nrow=3)
#ggsave("Figure9scot.png", plot=averagecompared, width=5.5, height=7.96, dpi=600)

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

a <- glm(mental~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18, family = "binomial")

b <- glm(hospital~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18, family = "binomial")

c <- glm(jobcent~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18, family = "binomial")

d <- glm(station~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18, family = "binomial")

e <- glm(furthered~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18, family = "binomial")

f <- glm(police~growth+distance+poplsoa+size+housedep+sc+jobdensity, data=data18, family = "binomial")
x<-tab_model(a,b,c,d,e, f, show.ci=0.95,title="Logistic Regression of Public Services availability" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Mental Health Practitioners", "Hospital", "Job Centre", "Train Station","Further Education College",  "Police Station"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"Scotland Dummy", "Job Density"))
tab_model(a,b,c,d,e, f, show.ci=0.95,title="Logistic Regression of Public Services availability" ,p.style = c("numeric"), auto.label=FALSE, dv.labels = c("Mental Health Practitioners", "Hospital", "Job Centre", "Train Station","Further Education College",  "Police Station"), pred.labels = c("Intercept",  "Improvement Index", "Distance from City", "Population", "Land Area","Household Deprivation" ,"Scotland Dummy", "Job Density"), file = "table5scot")

w<- plot_model(a,ci.lvl = 0 ,type = "pred", colors = c("navyblue", "blue"), terms = c("housedep", "sc"),title = "", axis.title = c( "Households in Deprivation (%, 2011)", "Probability of Present Mental Health Service"), legend.title = "Region")
y<- plot_model(c,ci.lvl = 0 ,type = "pred", colors = c("navyblue", "blue"), terms = c("housedep", "sc"), title = "Predicted Probability of a Town-based Job Centre", axis.title = c( "Households in Deprivation (%, 2011)", "Probability of Present Job Centre"), legend.title = "Region")
z<- plot_model(f,ci.lvl = 0, type = "pred", colors = c("navyblue", "blue"), terms = c("jobdensity", "sc"), title = "Predicted Probability of a Town-based Police Station", axis.title = c( "Town Type", "Probability of Present Police Station"), legend.title = "Region")

w <- w+geom_line(size=2)+coord_cartesian(ylim=c(0.25,1))
y <- y+geom_line(size=2)+coord_cartesian(ylim=c(0.25,1))
z <- z+geom_point(size=6)+coord_cartesian(ylim=c(0.25,1))
#ggsave("Figure10scot.png", plot=w, width=7.64, height=5.5, dpi=600)

compare <- ggarrange(w,y, ncol=2, nrow=1)

####flufff####

scotlandlist$allservicenumbers11 <- scotlandlist$bus11+scotlandlist$GPs11+scotlandlist$hospitals11+scotlandlist$furthered11+scotlandlist$station11+scotlandlist$schools11+scotlandlist$nursery11+scotlandlist$library11+scotlandlist$jobcent11+scotlandlist$police11+scotlandlist$mental11+scotlandlist$health11+scotlandlist$halls11+scotlandlist$fire11
scotlandlist$allservicenumbers18 <- scotlandlist$bus18+scotlandlist$GPs18+scotlandlist$hospitals18+scotlandlist$furthered18+scotlandlist$station18+scotlandlist$school18+scotlandlist$nursery18+scotlandlist$library18+scotlandlist$jobcent18+scotlandlist$police18+scotlandlist$mental18+scotlandlist$health18+scotlandlist$hall18+scotlandlist$fire18
sum(scotlandlist$fire18)
sum(scotlandlist$fire11)



