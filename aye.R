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
proj4string(buas) <- crs(world)
collar <- (brewer.pal( 10, "RdBu"))
citylayer <- list( "sp.polygons",citybua, fill= "gray77", col="transparent")
theme_set(theme_minimal())

###UK geography###

scotbua <- readOGR("C:\\Users\\bjg55\\Documents\\towns\\2019Ring", "Settlements2010")
scotbua2 <- readOGR("C:\\Users\\bjg55\\Documents\\towns\\2019Ring", "localities")
proj4string(scotbua) <- CRS("+init=epsg:27700")
proj4string(scotbua2) <- CRS("+init=epsg:27700")

scotbua <- spTransform(scotbua, CRS("+init=epsg:27700"))
scotbua2 <- spTransform(scotbua2, CRS("+init=epsg:27700"))

scotpop <- read_csv("scotpop.csv")
scotpop11 <- read_csv("scot11pops.csv")

names(scotbua)[names(scotbua)=="S10Code"] <- "BUA11CD"
names(scotbua)[names(scotbua)=="Shape_Area"] <- "st_areasha"
names(scotbua)[names(scotbua)=="Shape_Le_1"] <- "st_lengths"
names(scotbua)[names(scotbua)=="S10Name"] <- "bua11nm"
names(scotbua2)[names(scotbua2)=="code"] <- "BUA11CD"
names(scotbua2)[names(scotbua2)=="Shape_Area"] <- "st_areasha"
names(scotbua2)[names(scotbua2)=="Shape_Le_1"] <- "st_lengths"
names(scotbua2)[names(scotbua2)=="name"] <- "bua11nm"

scotbua <- merge(scotbua, scotpop11, by="bua11nm")
scotbua2 <- merge(scotbua2, scotpop, by="BUA11CD", all.x=TRUE)

scotbua <- scotbua[which(scotbua$pop11>10000 & scotbua$pop11<175000), ]
scotbua2 <- scotbua2[which(scotbua2$poplsoa>10000 & scotbua2$poplsoa<175000), ]
scotbua2 <- spTransform(scotbua2, CRS("+init=epsg:4326"))
scotbua <- spTransform(scotbua, CRS("+init=epsg:4326"))
plz <- ggplot()+
     geom_polygon(aes(long, lat, group=group),scotbua, fill=NA, color="black", size=0.07)+
     coord_quickmap()+
     geom_polygon(aes(long, lat, group=group),scotbua2, fill=NA, color="yellow", size=0.07)+
     theme_nothing(legend=TRUE)





