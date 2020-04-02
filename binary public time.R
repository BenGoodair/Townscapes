###load the packages###

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
citybua <- merge(buas, popbua, by="BUA11CD")
citybua <- citybua[which(citybua$poplsoa>175000), ]
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

### download the data###
POI_2013 <- read_csv("D:/poi/POI_2013.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2013 <- POI_2013 [complete.cases(POI_2013$feature_easting), ]
POI_2013 <- POI_2013 [complete.cases(POI_2013$feature_northing), ]
names(POI_2013)[names(POI_2013)=="pointx_classification_code"] <- "code"



mental13 <- POI_2013[which(POI_2013$code==5280372), ]
mental13 <- subset(mental13, "feature_easting" != "" | "feature_northing" != "")
mental13$poi_ID <- 1:nrow(mental13)
coords <- cbind(Easting = as.numeric(as.character(mental13$feature_easting)),
                Northing = as.numeric(as.character(mental13$feature_northing)))
mental13sp <- SpatialPointsDataFrame(coords, data = data.frame(mental13$name,
                                                            mental13$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental13 <- over(uktowns, geometry(mental13sp), returnList = TRUE)
mental13 <- sapply(mental13, length)
mental13 <- spCbind(uktowns, mental13)
mental13 <- merge(mental13, popbua, by="BUA11CD")
mental13 <- mental13[which(mental13$poplsoa<175000&mental13$poplsoa>10000), ]
mental13cent <- SpatialPointsDataFrame(gCentroid(mental13, byid=TRUE), mental13@data, match.ID=FALSE)
mental13 <- mental13[-c(2,3,4,5,7)]


station13 <- POI_2013[which(POI_2013$code==10570738), ]
station13 <- subset(station13, "feature_easting" != "" | "feature_northing" != "")
station13$poi_ID <- 1:nrow(station13)
coords <- cbind(Easting = as.numeric(as.character(station13$feature_easting)),
                Northing = as.numeric(as.character(station13$feature_northing)))
station13sp <- SpatialPointsDataFrame(coords, data = data.frame(station13$name,
                                                            station13$poi_ID), proj4string = CRS("+init=epsg:27700"))


station13 <- over(uktowns, geometry(station13sp), returnList = TRUE)
station13 <- sapply(station13, length)
station13 <- spCbind(uktowns, station13)
station13 <- merge(station13, popbua, by="BUA11CD")
station13 <- station13[which(station13$poplsoa<175000&station13$poplsoa>10000), ]
station13cent <- SpatialPointsDataFrame(gCentroid(station13, byid=TRUE), station13@data, match.ID=FALSE)
station13 <- station13[-c(2,3,4,5,7)]

jobcent13 <- POI_2013[which(POI_2013$code==6330418), ]
jobcent13 <- subset(jobcent13, "feature_easting" != "" | "feature_northing" != "")
jobcent13$poi_ID <- 1:nrow(jobcent13)
coords <- cbind(Easting = as.numeric(as.character(jobcent13$feature_easting)),
                Northing = as.numeric(as.character(jobcent13$feature_northing)))
jobcent13sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent13$name,
                                                            jobcent13$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent13 <- over(uktowns, geometry(jobcent13sp), returnList = TRUE)
jobcent13 <- sapply(jobcent13, length)
jobcent13 <- spCbind(uktowns, jobcent13)
jobcent13 <- merge(jobcent13, popbua, by="BUA11CD")
jobcent13 <- jobcent13[which(jobcent13$poplsoa<175000&jobcent13$poplsoa>10000), ]
jobcent13cent <- SpatialPointsDataFrame(gCentroid(jobcent13, byid=TRUE), jobcent13@data, match.ID=FALSE)
jobcent13 <- jobcent13[-c(2,3,4,5,7)]

furthered13 <- POI_2013[which(POI_2013$code==5310376), ]
furthered13 <- subset(furthered13, "feature_easting" != "" | "feature_northing" != "")
furthered13$poi_ID <- 1:nrow(furthered13)
coords <- cbind(Easting = as.numeric(as.character(furthered13$feature_easting)),
                Northing = as.numeric(as.character(furthered13$feature_northing)))
furthered13sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered13$name,
                                                            furthered13$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered13 <- over(uktowns, geometry(furthered13sp), returnList = TRUE)
furthered13 <- sapply(furthered13, length)
furthered13 <- spCbind(uktowns, furthered13)
furthered13 <- merge(furthered13, popbua, by="BUA11CD")
furthered13 <- furthered13[which(furthered13$poplsoa<175000&furthered13$poplsoa>10000), ]
furthered13cent <- SpatialPointsDataFrame(gCentroid(furthered13, byid=TRUE), furthered13@data, match.ID=FALSE)
furthered13 <- furthered13[-c(2,3,4,5,7)]

GPs13 <- POI_2013[which(POI_2013$code==5280369), ]
GPs13 <- subset(GPs13, "feature_easting" != "" | "feature_northing" != "")
GPs13$poi_ID <- 1:nrow(GPs13)
coords <- cbind(Easting = as.numeric(as.character(GPs13$feature_easting)),
                Northing = as.numeric(as.character(GPs13$feature_northing)))
GPs13sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs13$name,
                                                            GPs13$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs13 <- over(uktowns, geometry(GPs13sp), returnList = TRUE)
GPs13 <- sapply(GPs13, length)
GPs13 <- spCbind(uktowns, GPs13)
GPs13 <- merge(GPs13, popbua, by="BUA11CD")
GPs13 <- GPs13[which(GPs13$poplsoa<175000&GPs13$poplsoa>10000), ]
GPs13cent <- SpatialPointsDataFrame(gCentroid(GPs13, byid=TRUE), GPs13@data, match.ID=FALSE)
GPs13 <- GPs13[-c(2,3,4,5,7)]

hospitals13 <- POI_2013[which(POI_2013$code==5280371), ]
hospitals13 <- subset(hospitals13, "feature_easting" != "" | "feature_northing" != "")
hospitals13$poi_ID <- 1:nrow(hospitals13)
coords <- cbind(Easting = as.numeric(as.character(hospitals13$feature_easting)),
                Northing = as.numeric(as.character(hospitals13$feature_northing)))
hospitals13sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals13$name,
                                                            hospitals13$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals13 <- over(uktowns, geometry(hospitals13sp), returnList = TRUE)
hospitals13 <- sapply(hospitals13, length)
hospitals13 <- spCbind(uktowns, hospitals13)
hospitals13 <- merge(hospitals13, popbua, by="BUA11CD")
hospitals13 <- hospitals13[which(hospitals13$poplsoa<175000&hospitals13$poplsoa>10000), ]
hospitals13cent <- SpatialPointsDataFrame(gCentroid(hospitals13, byid=TRUE), hospitals13@data, match.ID=FALSE)
hospitals13 <- hospitals13[-c(2,3,4,5,7)]

hmrc13 <- POI_2013[which(POI_2013$code==5280371), ]
hmrc13 <- subset(hmrc13, "feature_easting" != "" | "feature_northing" != "")
hmrc13$poi_ID <- 1:nrow(hmrc13)
coords <- cbind(Easting = as.numeric(as.character(hmrc13$feature_easting)),
                Northing = as.numeric(as.character(hmrc13$feature_northing)))
hmrc13sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc13$name,
                                                                  hmrc13$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc13 <- over(uktowns, geometry(hmrc13sp), returnList = TRUE)
hmrc13 <- sapply(hmrc13, length)
hmrc13 <- spCbind(uktowns, hmrc13)
hmrc13 <- merge(hmrc13, popbua, by="BUA11CD")
hmrc13 <- hmrc13[which(hmrc13$poplsoa<175000&hmrc13$poplsoa>10000), ]
hmrc13cent <- SpatialPointsDataFrame(gCentroid(hmrc13, byid=TRUE), hmrc13@data, match.ID=FALSE)
hmrc13 <- hmrc13[-c(2,3,4,5,7)]



hospital <-  as.data.frame(hospitals13)
furthered <-  as.data.frame(furthered13)
hmrc <-  as.data.frame(hmrc13)
gps <-  as.data.frame(GPs13)
mental <-  as.data.frame(mental13)
station <-  as.data.frame(station13)
jobcent <-  as.data.frame(jobcent13)


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2012 <- read_csv("D:/poi/POI_2012.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2012 <- POI_2012 [complete.cases(POI_2012$feature_easting), ]
POI_2012 <- POI_2012 [complete.cases(POI_2012$feature_northing), ]
names(POI_2012)[names(POI_2012)=="pointx_classification_code"] <- "code"



mental12 <- POI_2012[which(POI_2012$code==5280372), ]
mental12 <- subset(mental12, "feature_easting" != "" | "feature_northing" != "")
mental12$poi_ID <- 1:nrow(mental12)
coords <- cbind(Easting = as.numeric(as.character(mental12$feature_easting)),
                Northing = as.numeric(as.character(mental12$feature_northing)))
mental12sp <- SpatialPointsDataFrame(coords, data = data.frame(mental12$name,
                                                               mental12$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental12 <- over(uktowns, geometry(mental12sp), returnList = TRUE)
mental12 <- sapply(mental12, length)
mental12 <- spCbind(uktowns, mental12)
mental12 <- merge(mental12, popbua, by="BUA11CD")
mental12 <- mental12[which(mental12$poplsoa<175000&mental12$poplsoa>10000), ]
mental12cent <- SpatialPointsDataFrame(gCentroid(mental12, byid=TRUE), mental12@data, match.ID=FALSE)
mental12 <- mental12[-c(2,3,4,5,7)]


station12 <- POI_2012[which(POI_2012$code==10570738), ]
station12 <- subset(station12, "feature_easting" != "" | "feature_northing" != "")
station12$poi_ID <- 1:nrow(station12)
coords <- cbind(Easting = as.numeric(as.character(station12$feature_easting)),
                Northing = as.numeric(as.character(station12$feature_northing)))
station12sp <- SpatialPointsDataFrame(coords, data = data.frame(station12$name,
                                                                station12$poi_ID), proj4string = CRS("+init=epsg:27700"))


station12 <- over(uktowns, geometry(station12sp), returnList = TRUE)
station12 <- sapply(station12, length)
station12 <- spCbind(uktowns, station12)
station12 <- merge(station12, popbua, by="BUA11CD")
station12 <- station12[which(station12$poplsoa<175000&station12$poplsoa>10000), ]
station12cent <- SpatialPointsDataFrame(gCentroid(station12, byid=TRUE), station12@data, match.ID=FALSE)
station12 <- station12[-c(2,3,4,5,7)]

jobcent12 <- POI_2012[which(POI_2012$code==6330418), ]
jobcent12 <- subset(jobcent12, "feature_easting" != "" | "feature_northing" != "")
jobcent12$poi_ID <- 1:nrow(jobcent12)
coords <- cbind(Easting = as.numeric(as.character(jobcent12$feature_easting)),
                Northing = as.numeric(as.character(jobcent12$feature_northing)))
jobcent12sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent12$name,
                                                                jobcent12$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent12 <- over(uktowns, geometry(jobcent12sp), returnList = TRUE)
jobcent12 <- sapply(jobcent12, length)
jobcent12 <- spCbind(uktowns, jobcent12)
jobcent12 <- merge(jobcent12, popbua, by="BUA11CD")
jobcent12 <- jobcent12[which(jobcent12$poplsoa<175000&jobcent12$poplsoa>10000), ]
jobcent12cent <- SpatialPointsDataFrame(gCentroid(jobcent12, byid=TRUE), jobcent12@data, match.ID=FALSE)
jobcent12 <- jobcent12[-c(2,3,4,5,7)]

furthered12 <- POI_2012[which(POI_2012$code==5310376), ]
furthered12 <- subset(furthered12, "feature_easting" != "" | "feature_northing" != "")
furthered12$poi_ID <- 1:nrow(furthered12)
coords <- cbind(Easting = as.numeric(as.character(furthered12$feature_easting)),
                Northing = as.numeric(as.character(furthered12$feature_northing)))
furthered12sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered12$name,
                                                                  furthered12$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered12 <- over(uktowns, geometry(furthered12sp), returnList = TRUE)
furthered12 <- sapply(furthered12, length)
furthered12 <- spCbind(uktowns, furthered12)
furthered12 <- merge(furthered12, popbua, by="BUA11CD")
furthered12 <- furthered12[which(furthered12$poplsoa<175000&furthered12$poplsoa>10000), ]
furthered12cent <- SpatialPointsDataFrame(gCentroid(furthered12, byid=TRUE), furthered12@data, match.ID=FALSE)
furthered12 <- furthered12[-c(2,3,4,5,7)]

GPs12 <- POI_2012[which(POI_2012$code==5280369), ]
GPs12 <- subset(GPs12, "feature_easting" != "" | "feature_northing" != "")
GPs12$poi_ID <- 1:nrow(GPs12)
coords <- cbind(Easting = as.numeric(as.character(GPs12$feature_easting)),
                Northing = as.numeric(as.character(GPs12$feature_northing)))
GPs12sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs12$name,
                                                            GPs12$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs12 <- over(uktowns, geometry(GPs12sp), returnList = TRUE)
GPs12 <- sapply(GPs12, length)
GPs12 <- spCbind(uktowns, GPs12)
GPs12 <- merge(GPs12, popbua, by="BUA11CD")
GPs12 <- GPs12[which(GPs12$poplsoa<175000&GPs12$poplsoa>10000), ]
GPs12cent <- SpatialPointsDataFrame(gCentroid(GPs12, byid=TRUE), GPs12@data, match.ID=FALSE)
GPs12 <- GPs12[-c(2,3,4,5,7)]

hospitals12 <- POI_2012[which(POI_2012$code==5280371), ]
hospitals12 <- subset(hospitals12, "feature_easting" != "" | "feature_northing" != "")
hospitals12$poi_ID <- 1:nrow(hospitals12)
coords <- cbind(Easting = as.numeric(as.character(hospitals12$feature_easting)),
                Northing = as.numeric(as.character(hospitals12$feature_northing)))
hospitals12sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals12$name,
                                                                  hospitals12$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals12 <- over(uktowns, geometry(hospitals12sp), returnList = TRUE)
hospitals12 <- sapply(hospitals12, length)
hospitals12 <- spCbind(uktowns, hospitals12)
hospitals12 <- merge(hospitals12, popbua, by="BUA11CD")
hospitals12 <- hospitals12[which(hospitals12$poplsoa<175000&hospitals12$poplsoa>10000), ]
hospitals12cent <- SpatialPointsDataFrame(gCentroid(hospitals12, byid=TRUE), hospitals12@data, match.ID=FALSE)
hospitals12 <- hospitals12[-c(2,3,4,5,7)]

hmrc12 <- POI_2012[which(POI_2012$code==5280371), ]
hmrc12 <- subset(hmrc12, "feature_easting" != "" | "feature_northing" != "")
hmrc12$poi_ID <- 1:nrow(hmrc12)
coords <- cbind(Easting = as.numeric(as.character(hmrc12$feature_easting)),
                Northing = as.numeric(as.character(hmrc12$feature_northing)))
hmrc12sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc12$name,
                                                             hmrc12$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc12 <- over(uktowns, geometry(hmrc12sp), returnList = TRUE)
hmrc12 <- sapply(hmrc12, length)
hmrc12 <- spCbind(uktowns, hmrc12)
hmrc12 <- merge(hmrc12, popbua, by="BUA11CD")
hmrc12 <- hmrc12[which(hmrc12$poplsoa<175000&hmrc12$poplsoa>10000), ]
hmrc12cent <- SpatialPointsDataFrame(gCentroid(hmrc12, byid=TRUE), hmrc12@data, match.ID=FALSE)
hmrc12 <- hmrc12[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals12, by="BUA11CD")
furthered <-  merge(furthered, furthered12, by="BUA11CD")
hmrc <-  merge(hmrc, hmrc12, by="BUA11CD")
gps <-  merge(gps, GPs12, by="BUA11CD")
mental <-  merge(mental, mental12, "BUA11CD")
station <-  merge(station, station12, by="BUA11CD")
jobcent <-  merge(jobcent, jobcent12, by="BUA11CD")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])





POI_2008 <- read_csv("D:/poi/POI_2008.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2008 <- POI_2008 [complete.cases(POI_2008$easting), ]
POI_2008 <- POI_2008 [complete.cases(POI_2008$northing), ]
names(POI_2008)[names(POI_2008)=="pointx_classification_code"] <- "code"



mental08 <- POI_2008[which(POI_2008$code==5280372), ]
mental08 <- subset(mental08, "easting" != "" | "northing" != "")
mental08$poi_ID <- 1:nrow(mental08)
coords <- cbind(Easting = as.numeric(as.character(mental08$easting)),
                Northing = as.numeric(as.character(mental08$northing)))
mental08sp <- SpatialPointsDataFrame(coords, data = data.frame(mental08$name,
                                                               mental08$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental08 <- over(uktowns, geometry(mental08sp), returnList = TRUE)
mental08 <- sapply(mental08, length)
mental08 <- spCbind(uktowns, mental08)
mental08 <- merge(mental08, popbua, by="BUA11CD")
mental08 <- mental08[which(mental08$poplsoa<175000&mental08$poplsoa>10000), ]
mental08cent <- SpatialPointsDataFrame(gCentroid(mental08, byid=TRUE), mental08@data, match.ID=FALSE)
mental08 <- mental08[-c(2,3,4,5,7)]


station08 <- POI_2008[which(POI_2008$code==10570738), ]
station08 <- subset(station08, "easting" != "" | "northing" != "")
station08$poi_ID <- 1:nrow(station08)
coords <- cbind(Easting = as.numeric(as.character(station08$easting)),
                Northing = as.numeric(as.character(station08$northing)))
station08sp <- SpatialPointsDataFrame(coords, data = data.frame(station08$name,
                                                                station08$poi_ID), proj4string = CRS("+init=epsg:27700"))


station08 <- over(uktowns, geometry(station08sp), returnList = TRUE)
station08 <- sapply(station08, length)
station08 <- spCbind(uktowns, station08)
station08 <- merge(station08, popbua, by="BUA11CD")
station08 <- station08[which(station08$poplsoa<175000&station08$poplsoa>10000), ]
station08cent <- SpatialPointsDataFrame(gCentroid(station08, byid=TRUE), station08@data, match.ID=FALSE)
station08 <- station08[-c(2,3,4,5,7)]

jobcent08 <- POI_2008[which(POI_2008$code==6330418), ]
jobcent08 <- subset(jobcent08, "easting" != "" | "northing" != "")
jobcent08$poi_ID <- 1:nrow(jobcent08)
coords <- cbind(Easting = as.numeric(as.character(jobcent08$easting)),
                Northing = as.numeric(as.character(jobcent08$northing)))
jobcent08sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent08$name,
                                                                jobcent08$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent08 <- over(uktowns, geometry(jobcent08sp), returnList = TRUE)
jobcent08 <- sapply(jobcent08, length)
jobcent08 <- spCbind(uktowns, jobcent08)
jobcent08 <- merge(jobcent08, popbua, by="BUA11CD")
jobcent08 <- jobcent08[which(jobcent08$poplsoa<175000&jobcent08$poplsoa>10000), ]
jobcent08cent <- SpatialPointsDataFrame(gCentroid(jobcent08, byid=TRUE), jobcent08@data, match.ID=FALSE)
jobcent08 <- jobcent08[-c(2,3,4,5,7)]

furthered08 <- POI_2008[which(POI_2008$code==5310376), ]
furthered08 <- subset(furthered08, "easting" != "" | "northing" != "")
furthered08$poi_ID <- 1:nrow(furthered08)
coords <- cbind(Easting = as.numeric(as.character(furthered08$easting)),
                Northing = as.numeric(as.character(furthered08$northing)))
furthered08sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered08$name,
                                                                  furthered08$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered08 <- over(uktowns, geometry(furthered08sp), returnList = TRUE)
furthered08 <- sapply(furthered08, length)
furthered08 <- spCbind(uktowns, furthered08)
furthered08 <- merge(furthered08, popbua, by="BUA11CD")
furthered08 <- furthered08[which(furthered08$poplsoa<175000&furthered08$poplsoa>10000), ]
furthered08cent <- SpatialPointsDataFrame(gCentroid(furthered08, byid=TRUE), furthered08@data, match.ID=FALSE)
furthered08 <- furthered08[-c(2,3,4,5,7)]

GPs08 <- POI_2008[which(POI_2008$code==5280369), ]
GPs08 <- subset(GPs08, "easting" != "" | "northing" != "")
GPs08$poi_ID <- 1:nrow(GPs08)
coords <- cbind(Easting = as.numeric(as.character(GPs08$easting)),
                Northing = as.numeric(as.character(GPs08$northing)))
GPs08sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs08$name,
                                                            GPs08$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs08 <- over(uktowns, geometry(GPs08sp), returnList = TRUE)
GPs08 <- sapply(GPs08, length)
GPs08 <- spCbind(uktowns, GPs08)
GPs08 <- merge(GPs08, popbua, by="BUA11CD")
GPs08 <- GPs08[which(GPs08$poplsoa<175000&GPs08$poplsoa>10000), ]
GPs08cent <- SpatialPointsDataFrame(gCentroid(GPs08, byid=TRUE), GPs08@data, match.ID=FALSE)
GPs08 <- GPs08[-c(2,3,4,5,7)]

hospitals08 <- POI_2008[which(POI_2008$code==5280371), ]
hospitals08 <- subset(hospitals08, "easting" != "" | "northing" != "")
hospitals08$poi_ID <- 1:nrow(hospitals08)
coords <- cbind(Easting = as.numeric(as.character(hospitals08$easting)),
                Northing = as.numeric(as.character(hospitals08$northing)))
hospitals08sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals08$name,
                                                                  hospitals08$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals08 <- over(uktowns, geometry(hospitals08sp), returnList = TRUE)
hospitals08 <- sapply(hospitals08, length)
hospitals08 <- spCbind(uktowns, hospitals08)
hospitals08 <- merge(hospitals08, popbua, by="BUA11CD")
hospitals08 <- hospitals08[which(hospitals08$poplsoa<175000&hospitals08$poplsoa>10000), ]
hospitals08cent <- SpatialPointsDataFrame(gCentroid(hospitals08, byid=TRUE), hospitals08@data, match.ID=FALSE)
hospitals08 <- hospitals08[-c(2,3,4,5,7)]

hmrc08 <- POI_2008[which(POI_2008$code==5280371), ]
hmrc08 <- subset(hmrc08, "easting" != "" | "northing" != "")
hmrc08$poi_ID <- 1:nrow(hmrc08)
coords <- cbind(Easting = as.numeric(as.character(hmrc08$easting)),
                Northing = as.numeric(as.character(hmrc08$northing)))
hmrc08sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc08$name,
                                                             hmrc08$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc08 <- over(uktowns, geometry(hmrc08sp), returnList = TRUE)
hmrc08 <- sapply(hmrc08, length)
hmrc08 <- spCbind(uktowns, hmrc08)
hmrc08 <- merge(hmrc08, popbua, by="BUA11CD")
hmrc08 <- hmrc08[which(hmrc08$poplsoa<175000&hmrc08$poplsoa>10000), ]
hmrc08cent <- SpatialPointsDataFrame(gCentroid(hmrc08, byid=TRUE), hmrc08@data, match.ID=FALSE)
hmrc08 <- hmrc08[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals08, by="BUA11CD")
furthered <-  merge(furthered, furthered08, by="BUA11CD")
hmrc <-  merge(hmrc, hmrc08, by="BUA11CD")
gps <-  merge(gps, GPs08, by="BUA11CD")
mental <-  merge(mental, mental08, "BUA11CD")
station <-  merge(station, station08, by="BUA11CD")
jobcent <-  merge(jobcent, jobcent08, by="BUA11CD")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])




POI_2009 <- read_csv("D:/poi/POI_2009.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2009 <- POI_2009 [complete.cases(POI_2009$easting), ]
POI_2009 <- POI_2009 [complete.cases(POI_2009$northing), ]
names(POI_2009)[names(POI_2009)=="pointx_classification_code"] <- "code"



mental09 <- POI_2009[which(POI_2009$code==5280372), ]
mental09 <- subset(mental09, "easting" != "" | "northing" != "")
mental09$poi_ID <- 1:nrow(mental09)
coords <- cbind(Easting = as.numeric(as.character(mental09$easting)),
                Northing = as.numeric(as.character(mental09$northing)))
mental09sp <- SpatialPointsDataFrame(coords, data = data.frame(mental09$name,
                                                               mental09$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental09 <- over(uktowns, geometry(mental09sp), returnList = TRUE)
mental09 <- sapply(mental09, length)
mental09 <- spCbind(uktowns, mental09)
mental09 <- merge(mental09, popbua, by="BUA11CD")
mental09 <- mental09[which(mental09$poplsoa<175000&mental09$poplsoa>10000), ]
mental09cent <- SpatialPointsDataFrame(gCentroid(mental09, byid=TRUE), mental09@data, match.ID=FALSE)
mental09 <- mental09[-c(2,3,4,5,7)]


station09 <- POI_2009[which(POI_2009$code==10570738), ]
station09 <- subset(station09, "easting" != "" | "northing" != "")
station09$poi_ID <- 1:nrow(station09)
coords <- cbind(Easting = as.numeric(as.character(station09$easting)),
                Northing = as.numeric(as.character(station09$northing)))
station09sp <- SpatialPointsDataFrame(coords, data = data.frame(station09$name,
                                                                station09$poi_ID), proj4string = CRS("+init=epsg:27700"))


station09 <- over(uktowns, geometry(station09sp), returnList = TRUE)
station09 <- sapply(station09, length)
station09 <- spCbind(uktowns, station09)
station09 <- merge(station09, popbua, by="BUA11CD")
station09 <- station09[which(station09$poplsoa<175000&station09$poplsoa>10000), ]
station09cent <- SpatialPointsDataFrame(gCentroid(station09, byid=TRUE), station09@data, match.ID=FALSE)
station09 <- station09[-c(2,3,4,5,7)]

jobcent09 <- POI_2009[which(POI_2009$code==6330418), ]
jobcent09 <- subset(jobcent09, "easting" != "" | "northing" != "")
jobcent09$poi_ID <- 1:nrow(jobcent09)
coords <- cbind(Easting = as.numeric(as.character(jobcent09$easting)),
                Northing = as.numeric(as.character(jobcent09$northing)))
jobcent09sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent09$name,
                                                                jobcent09$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent09 <- over(uktowns, geometry(jobcent09sp), returnList = TRUE)
jobcent09 <- sapply(jobcent09, length)
jobcent09 <- spCbind(uktowns, jobcent09)
jobcent09 <- merge(jobcent09, popbua, by="BUA11CD")
jobcent09 <- jobcent09[which(jobcent09$poplsoa<175000&jobcent09$poplsoa>10000), ]
jobcent09cent <- SpatialPointsDataFrame(gCentroid(jobcent09, byid=TRUE), jobcent09@data, match.ID=FALSE)
jobcent09 <- jobcent09[-c(2,3,4,5,7)]

furthered09 <- POI_2009[which(POI_2009$code==5310376), ]
furthered09 <- subset(furthered09, "easting" != "" | "northing" != "")
furthered09$poi_ID <- 1:nrow(furthered09)
coords <- cbind(Easting = as.numeric(as.character(furthered09$easting)),
                Northing = as.numeric(as.character(furthered09$northing)))
furthered09sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered09$name,
                                                                  furthered09$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered09 <- over(uktowns, geometry(furthered09sp), returnList = TRUE)
furthered09 <- sapply(furthered09, length)
furthered09 <- spCbind(uktowns, furthered09)
furthered09 <- merge(furthered09, popbua, by="BUA11CD")
furthered09 <- furthered09[which(furthered09$poplsoa<175000&furthered09$poplsoa>10000), ]
furthered09cent <- SpatialPointsDataFrame(gCentroid(furthered09, byid=TRUE), furthered09@data, match.ID=FALSE)
furthered09 <- furthered09[-c(2,3,4,5,7)]

GPs09 <- POI_2009[which(POI_2009$code==5280369), ]
GPs09 <- subset(GPs09, "easting" != "" | "northing" != "")
GPs09$poi_ID <- 1:nrow(GPs09)
coords <- cbind(Easting = as.numeric(as.character(GPs09$easting)),
                Northing = as.numeric(as.character(GPs09$northing)))
GPs09sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs09$name,
                                                            GPs09$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs09 <- over(uktowns, geometry(GPs09sp), returnList = TRUE)
GPs09 <- sapply(GPs09, length)
GPs09 <- spCbind(uktowns, GPs09)
GPs09 <- merge(GPs09, popbua, by="BUA11CD")
GPs09 <- GPs09[which(GPs09$poplsoa<175000&GPs09$poplsoa>10000), ]
GPs09cent <- SpatialPointsDataFrame(gCentroid(GPs09, byid=TRUE), GPs09@data, match.ID=FALSE)
GPs09 <- GPs09[-c(2,3,4,5,7)]

hospitals09 <- POI_2009[which(POI_2009$code==5280371), ]
hospitals09 <- subset(hospitals09, "easting" != "" | "northing" != "")
hospitals09$poi_ID <- 1:nrow(hospitals09)
coords <- cbind(Easting = as.numeric(as.character(hospitals09$easting)),
                Northing = as.numeric(as.character(hospitals09$northing)))
hospitals09sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals09$name,
                                                                  hospitals09$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals09 <- over(uktowns, geometry(hospitals09sp), returnList = TRUE)
hospitals09 <- sapply(hospitals09, length)
hospitals09 <- spCbind(uktowns, hospitals09)
hospitals09 <- merge(hospitals09, popbua, by="BUA11CD")
hospitals09 <- hospitals09[which(hospitals09$poplsoa<175000&hospitals09$poplsoa>10000), ]
hospitals09cent <- SpatialPointsDataFrame(gCentroid(hospitals09, byid=TRUE), hospitals09@data, match.ID=FALSE)
hospitals09 <- hospitals09[-c(2,3,4,5,7)]

hmrc09 <- POI_2009[which(POI_2009$code==5280371), ]
hmrc09 <- subset(hmrc09, "easting" != "" | "northing" != "")
hmrc09$poi_ID <- 1:nrow(hmrc09)
coords <- cbind(Easting = as.numeric(as.character(hmrc09$easting)),
                Northing = as.numeric(as.character(hmrc09$northing)))
hmrc09sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc09$name,
                                                             hmrc09$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc09 <- over(uktowns, geometry(hmrc09sp), returnList = TRUE)
hmrc09 <- sapply(hmrc09, length)
hmrc09 <- spCbind(uktowns, hmrc09)
hmrc09 <- merge(hmrc09, popbua, by="BUA11CD")
hmrc09 <- hmrc09[which(hmrc09$poplsoa<175000&hmrc09$poplsoa>10000), ]
hmrc09cent <- SpatialPointsDataFrame(gCentroid(hmrc09, byid=TRUE), hmrc09@data, match.ID=FALSE)
hmrc09 <- hmrc09[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals09, by="BUA11CD")
furthered <-  merge(furthered, furthered09, by="BUA11CD")
hmrc <-  merge(hmrc, hmrc09, by="BUA11CD")
gps <-  merge(gps, GPs09, by="BUA11CD")
mental <-  merge(mental, mental09, "BUA11CD")
station <-  merge(station, station09, by="BUA11CD")
jobcent <-  merge(jobcent, jobcent09, by="BUA11CD")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])





POI_2010 <- read_csv("D:/poi/POI_2010.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2010 <- POI_2010 [complete.cases(POI_2010$easting), ]
POI_2010 <- POI_2010 [complete.cases(POI_2010$northing), ]
names(POI_2010)[names(POI_2010)=="pointx_classification_code"] <- "code"



mental10 <- POI_2010[which(POI_2010$code==5280372), ]
mental10 <- subset(mental10, "easting" != "" | "northing" != "")
mental10$poi_ID <- 1:nrow(mental10)
coords <- cbind(Easting = as.numeric(as.character(mental10$easting)),
                Northing = as.numeric(as.character(mental10$northing)))
mental10sp <- SpatialPointsDataFrame(coords, data = data.frame(mental10$name,
                                                               mental10$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental10 <- over(uktowns, geometry(mental10sp), returnList = TRUE)
mental10 <- sapply(mental10, length)
mental10 <- spCbind(uktowns, mental10)
mental10 <- merge(mental10, popbua, by="BUA11CD")
mental10 <- mental10[which(mental10$poplsoa<175000&mental10$poplsoa>10000), ]
mental10cent <- SpatialPointsDataFrame(gCentroid(mental10, byid=TRUE), mental10@data, match.ID=FALSE)
mental10 <- mental10[-c(2,3,4,5,7)]


station10 <- POI_2010[which(POI_2010$code==10570738), ]
station10 <- subset(station10, "easting" != "" | "northing" != "")
station10$poi_ID <- 1:nrow(station10)
coords <- cbind(Easting = as.numeric(as.character(station10$easting)),
                Northing = as.numeric(as.character(station10$northing)))
station10sp <- SpatialPointsDataFrame(coords, data = data.frame(station10$name,
                                                                station10$poi_ID), proj4string = CRS("+init=epsg:27700"))


station10 <- over(uktowns, geometry(station10sp), returnList = TRUE)
station10 <- sapply(station10, length)
station10 <- spCbind(uktowns, station10)
station10 <- merge(station10, popbua, by="BUA11CD")
station10 <- station10[which(station10$poplsoa<175000&station10$poplsoa>10000), ]
station10cent <- SpatialPointsDataFrame(gCentroid(station10, byid=TRUE), station10@data, match.ID=FALSE)
station10 <- station10[-c(2,3,4,5,7)]

jobcent10 <- POI_2010[which(POI_2010$code==6330418), ]
jobcent10 <- subset(jobcent10, "easting" != "" | "northing" != "")
jobcent10$poi_ID <- 1:nrow(jobcent10)
coords <- cbind(Easting = as.numeric(as.character(jobcent10$easting)),
                Northing = as.numeric(as.character(jobcent10$northing)))
jobcent10sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent10$name,
                                                                jobcent10$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent10 <- over(uktowns, geometry(jobcent10sp), returnList = TRUE)
jobcent10 <- sapply(jobcent10, length)
jobcent10 <- spCbind(uktowns, jobcent10)
jobcent10 <- merge(jobcent10, popbua, by="BUA11CD")
jobcent10 <- jobcent10[which(jobcent10$poplsoa<175000&jobcent10$poplsoa>10000), ]
jobcent10cent <- SpatialPointsDataFrame(gCentroid(jobcent10, byid=TRUE), jobcent10@data, match.ID=FALSE)
jobcent10 <- jobcent10[-c(2,3,4,5,7)]

furthered10 <- POI_2010[which(POI_2010$code==5310376), ]
furthered10 <- subset(furthered10, "easting" != "" | "northing" != "")
furthered10$poi_ID <- 1:nrow(furthered10)
coords <- cbind(Easting = as.numeric(as.character(furthered10$easting)),
                Northing = as.numeric(as.character(furthered10$northing)))
furthered10sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered10$name,
                                                                  furthered10$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered10 <- over(uktowns, geometry(furthered10sp), returnList = TRUE)
furthered10 <- sapply(furthered10, length)
furthered10 <- spCbind(uktowns, furthered10)
furthered10 <- merge(furthered10, popbua, by="BUA11CD")
furthered10 <- furthered10[which(furthered10$poplsoa<175000&furthered10$poplsoa>10000), ]
furthered10cent <- SpatialPointsDataFrame(gCentroid(furthered10, byid=TRUE), furthered10@data, match.ID=FALSE)
furthered10 <- furthered10[-c(2,3,4,5,7)]

GPs10 <- POI_2010[which(POI_2010$code==5280369), ]
GPs10 <- subset(GPs10, "easting" != "" | "northing" != "")
GPs10$poi_ID <- 1:nrow(GPs10)
coords <- cbind(Easting = as.numeric(as.character(GPs10$easting)),
                Northing = as.numeric(as.character(GPs10$northing)))
GPs10sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs10$name,
                                                            GPs10$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs10 <- over(uktowns, geometry(GPs10sp), returnList = TRUE)
GPs10 <- sapply(GPs10, length)
GPs10 <- spCbind(uktowns, GPs10)
GPs10 <- merge(GPs10, popbua, by="BUA11CD")
GPs10 <- GPs10[which(GPs10$poplsoa<175000&GPs10$poplsoa>10000), ]
GPs10cent <- SpatialPointsDataFrame(gCentroid(GPs10, byid=TRUE), GPs10@data, match.ID=FALSE)
GPs10 <- GPs10[-c(2,3,4,5,7)]

hospitals10 <- POI_2010[which(POI_2010$code==5280371), ]
hospitals10 <- subset(hospitals10, "easting" != "" | "northing" != "")
hospitals10$poi_ID <- 1:nrow(hospitals10)
coords <- cbind(Easting = as.numeric(as.character(hospitals10$easting)),
                Northing = as.numeric(as.character(hospitals10$northing)))
hospitals10sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals10$name,
                                                                  hospitals10$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals10 <- over(uktowns, geometry(hospitals10sp), returnList = TRUE)
hospitals10 <- sapply(hospitals10, length)
hospitals10 <- spCbind(uktowns, hospitals10)
hospitals10 <- merge(hospitals10, popbua, by="BUA11CD")
hospitals10 <- hospitals10[which(hospitals10$poplsoa<175000&hospitals10$poplsoa>10000), ]
hospitals10cent <- SpatialPointsDataFrame(gCentroid(hospitals10, byid=TRUE), hospitals10@data, match.ID=FALSE)
hospitals10 <- hospitals10[-c(2,3,4,5,7)]

hmrc10 <- POI_2010[which(POI_2010$code==5280371), ]
hmrc10 <- subset(hmrc10, "easting" != "" | "northing" != "")
hmrc10$poi_ID <- 1:nrow(hmrc10)
coords <- cbind(Easting = as.numeric(as.character(hmrc10$easting)),
                Northing = as.numeric(as.character(hmrc10$northing)))
hmrc10sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc10$name,
                                                             hmrc10$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc10 <- over(uktowns, geometry(hmrc10sp), returnList = TRUE)
hmrc10 <- sapply(hmrc10, length)
hmrc10 <- spCbind(uktowns, hmrc10)
hmrc10 <- merge(hmrc10, popbua, by="BUA11CD")
hmrc10 <- hmrc10[which(hmrc10$poplsoa<175000&hmrc10$poplsoa>10000), ]
hmrc10cent <- SpatialPointsDataFrame(gCentroid(hmrc10, byid=TRUE), hmrc10@data, match.ID=FALSE)
hmrc10 <- hmrc10[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals10, by="BUA11CD")
furthered <-  merge(furthered, furthered10, by="BUA11CD")
hmrc <-  merge(hmrc, hmrc10, by="BUA11CD")
gps <-  merge(gps, GPs10, by="BUA11CD")
mental <-  merge(mental, mental10, "BUA11CD")
station <-  merge(station, station10, by="BUA11CD")
jobcent <-  merge(jobcent, jobcent10, by="BUA11CD")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2011 <- read_csv("D:/poi/POI_2011.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2011 <- POI_2011 [complete.cases(POI_2011$feature_easting), ]
POI_2011 <- POI_2011 [complete.cases(POI_2011$feature_northing), ]
names(POI_2011)[names(POI_2011)=="pointx_classification_code"] <- "code"



mental11 <- POI_2011[which(POI_2011$code==5280372), ]
mental11 <- subset(mental11, "feature_easting" != "" | "feature_northing" != "")
mental11$poi_ID <- 1:nrow(mental11)
coords <- cbind(Easting = as.numeric(as.character(mental11$feature_easting)),
                Northing = as.numeric(as.character(mental11$feature_northing)))
mental11sp <- SpatialPointsDataFrame(coords, data = data.frame(mental11$name,
                                                               mental11$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental11 <- over(uktowns, geometry(mental11sp), returnList = TRUE)
mental11 <- sapply(mental11, length)
mental11 <- spCbind(uktowns, mental11)
mental11 <- merge(mental11, popbua, by="BUA11CD")
mental11 <- mental11[which(mental11$poplsoa<175000&mental11$poplsoa>10000), ]
mental11cent <- SpatialPointsDataFrame(gCentroid(mental11, byid=TRUE), mental11@data, match.ID=FALSE)
mental11 <- mental11[-c(2,3,4,5,7)]


station11 <- POI_2011[which(POI_2011$code==10570738), ]
station11 <- subset(station11, "feature_easting" != "" | "feature_northing" != "")
station11$poi_ID <- 1:nrow(station11)
coords <- cbind(Easting = as.numeric(as.character(station11$feature_easting)),
                Northing = as.numeric(as.character(station11$feature_northing)))
station11sp <- SpatialPointsDataFrame(coords, data = data.frame(station11$name,
                                                                station11$poi_ID), proj4string = CRS("+init=epsg:27700"))


station11 <- over(uktowns, geometry(station11sp), returnList = TRUE)
station11 <- sapply(station11, length)
station11 <- spCbind(uktowns, station11)
station11 <- merge(station11, popbua, by="BUA11CD")
station11 <- station11[which(station11$poplsoa<175000&station11$poplsoa>10000), ]
station11cent <- SpatialPointsDataFrame(gCentroid(station11, byid=TRUE), station11@data, match.ID=FALSE)
station11 <- station11[-c(2,3,4,5,7)]

jobcent11 <- POI_2011[which(POI_2011$code==6330418), ]
jobcent11 <- subset(jobcent11, "feature_easting" != "" | "feature_northing" != "")
jobcent11$poi_ID <- 1:nrow(jobcent11)
coords <- cbind(Easting = as.numeric(as.character(jobcent11$feature_easting)),
                Northing = as.numeric(as.character(jobcent11$feature_northing)))
jobcent11sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent11$name,
                                                                jobcent11$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent11 <- over(uktowns, geometry(jobcent11sp), returnList = TRUE)
jobcent11 <- sapply(jobcent11, length)
jobcent11 <- spCbind(uktowns, jobcent11)
jobcent11 <- merge(jobcent11, popbua, by="BUA11CD")
jobcent11 <- jobcent11[which(jobcent11$poplsoa<175000&jobcent11$poplsoa>10000), ]
jobcent11cent <- SpatialPointsDataFrame(gCentroid(jobcent11, byid=TRUE), jobcent11@data, match.ID=FALSE)
jobcent11 <- jobcent11[-c(2,3,4,5,7)]

furthered11 <- POI_2011[which(POI_2011$code==5310376), ]
furthered11 <- subset(furthered11, "feature_easting" != "" | "feature_northing" != "")
furthered11$poi_ID <- 1:nrow(furthered11)
coords <- cbind(Easting = as.numeric(as.character(furthered11$feature_easting)),
                Northing = as.numeric(as.character(furthered11$feature_northing)))
furthered11sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered11$name,
                                                                  furthered11$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered11 <- over(uktowns, geometry(furthered11sp), returnList = TRUE)
furthered11 <- sapply(furthered11, length)
furthered11 <- spCbind(uktowns, furthered11)
furthered11 <- merge(furthered11, popbua, by="BUA11CD")
furthered11 <- furthered11[which(furthered11$poplsoa<175000&furthered11$poplsoa>10000), ]
furthered11cent <- SpatialPointsDataFrame(gCentroid(furthered11, byid=TRUE), furthered11@data, match.ID=FALSE)
furthered11 <- furthered11[-c(2,3,4,5,7)]

GPs11 <- POI_2011[which(POI_2011$code==5280369), ]
GPs11 <- subset(GPs11, "feature_easting" != "" | "feature_northing" != "")
GPs11$poi_ID <- 1:nrow(GPs11)
coords <- cbind(Easting = as.numeric(as.character(GPs11$feature_easting)),
                Northing = as.numeric(as.character(GPs11$feature_northing)))
GPs11sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs11$name,
                                                            GPs11$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs11 <- over(uktowns, geometry(GPs11sp), returnList = TRUE)
GPs11 <- sapply(GPs11, length)
GPs11 <- spCbind(uktowns, GPs11)
GPs11 <- merge(GPs11, popbua, by="BUA11CD")
GPs11 <- GPs11[which(GPs11$poplsoa<175000&GPs11$poplsoa>10000), ]
GPs11cent <- SpatialPointsDataFrame(gCentroid(GPs11, byid=TRUE), GPs11@data, match.ID=FALSE)
GPs11 <- GPs11[-c(2,3,4,5,7)]

hospitals11 <- POI_2011[which(POI_2011$code==5280371), ]
hospitals11 <- subset(hospitals11, "feature_easting" != "" | "feature_northing" != "")
hospitals11$poi_ID <- 1:nrow(hospitals11)
coords <- cbind(Easting = as.numeric(as.character(hospitals11$feature_easting)),
                Northing = as.numeric(as.character(hospitals11$feature_northing)))
hospitals11sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals11$name,
                                                                  hospitals11$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals11 <- over(uktowns, geometry(hospitals11sp), returnList = TRUE)
hospitals11 <- sapply(hospitals11, length)
hospitals11 <- spCbind(uktowns, hospitals11)
hospitals11 <- merge(hospitals11, popbua, by="BUA11CD")
hospitals11 <- hospitals11[which(hospitals11$poplsoa<175000&hospitals11$poplsoa>10000), ]
hospitals11cent <- SpatialPointsDataFrame(gCentroid(hospitals11, byid=TRUE), hospitals11@data, match.ID=FALSE)
hospitals11 <- hospitals11[-c(2,3,4,5,7)]

hmrc11 <- POI_2011[which(POI_2011$code==5280371), ]
hmrc11 <- subset(hmrc11, "feature_easting" != "" | "feature_northing" != "")
hmrc11$poi_ID <- 1:nrow(hmrc11)
coords <- cbind(Easting = as.numeric(as.character(hmrc11$feature_easting)),
                Northing = as.numeric(as.character(hmrc11$feature_northing)))
hmrc11sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc11$name,
                                                             hmrc11$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc11 <- over(uktowns, geometry(hmrc11sp), returnList = TRUE)
hmrc11 <- sapply(hmrc11, length)
hmrc11 <- spCbind(uktowns, hmrc11)
hmrc11 <- merge(hmrc11, popbua, by="BUA11CD")
hmrc11 <- hmrc11[which(hmrc11$poplsoa<175000&hmrc11$poplsoa>10000), ]
hmrc11cent <- SpatialPointsDataFrame(gCentroid(hmrc11, byid=TRUE), hmrc11@data, match.ID=FALSE)
hmrc11 <- hmrc11[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals11, by="BUA11CD")
furthered <-  merge(furthered, furthered11, by="BUA11CD")
hmrc <-  merge(hmrc, hmrc11, by="BUA11CD")
gps <-  merge(gps, GPs11, by="BUA11CD")
mental <-  merge(mental, mental11, "BUA11CD")
station <-  merge(station, station11, by="BUA11CD")
jobcent <-  merge(jobcent, jobcent11, by="BUA11CD")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


###2014 onwards###


health_2014 <- read_csv("health_2014.csv", 
                     col_types = cols('PointX Classification Code' = col_number()))
health_2014 <- health_2014 [complete.cases(health_2014$`Feature Easting`), ]
health_2014 <- health_2014 [complete.cases(health_2014$`Feature Northing`), ]
names(health_2014)[names(health_2014)=="PointX Classification Code"] <- "code"

schools_2014 <- read_csv("school_2014.csv", 
                        col_types = cols('PointX Classification Code' = col_number()))
schools_2014 <- schools_2014 [complete.cases(schools_2014$`Feature Easting`), ]
schools_2014 <- schools_2014 [complete.cases(schools_2014$`Feature Northing`), ]
names(schools_2014)[names(schools_2014)=="PointX Classification Code"] <- "code"

gov_2014 <- read_csv("gov_2014.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Easting`), ]
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Northing`), ]
names(gov_2014)[names(gov_2014)=="PointX Classification Code"] <- "code"

station_2014 <- read_csv("station_2014.csv", 
                     col_types = cols('PointX Classification Code' = col_number()))
station_2014 <- station_2014 [complete.cases(station_2014$`Feature Easting`), ]
station_2014 <- station_2014 [complete.cases(station_2014$`Feature Northing`), ]
names(station_2014)[names(station_2014)=="PointX Classification Code"] <- "code"


mental14 <- health_2014[which(health_2014$code==05280372), ]
mental14 <- subset(mental14, "Feature Easting" != "" | "Feature Northing" != "")
mental14$poi_ID <- 1:nrow(mental14)
coords <- cbind(Easting = as.numeric(as.character(mental14$`Feature Easting`)),
                Northing = as.numeric(as.character(mental14$`Feature Northing`)))
mental14sp <- SpatialPointsDataFrame(coords, data = data.frame(mental14$Name,
                                                               mental14$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental14 <- over(uktowns, geometry(mental14sp), returnList = TRUE)
mental14 <- sapply(mental14, length)
mental14 <- spCbind(uktowns, mental14)
mental14 <- merge(mental14, popbua, by="BUA14CD")
mental14 <- mental14[which(mental14$poplsoa<175000&mental14$poplsoa>10000), ]
mental14cent <- SpatialPointsDataFrame(gCentroid(mental14, byid=TRUE), mental14@data, match.ID=FALSE)
mental14 <- mental14[-c(2,3,4,5,7)]


station14 <- station_2014[which(station_2014$code==010570738), ]
station14 <- subset(station14, "Feature Easting" != "" | "Feature Northing" != "")
station14$poi_ID <- 1:nrow(station14)
coords <- cbind(Easting = as.numeric(as.character(station14$`Feature Easting`)),
                Northing = as.numeric(as.character(station14$`Feature Northing`)))
station14sp <- SpatialPointsDataFrame(coords, data = data.frame(station14$Name,
                                                                station14$poi_ID), proj4string = CRS("+init=epsg:27700"))


station14 <- over(uktowns, geometry(station14sp), returnList = TRUE)
station14 <- sapply(station14, length)
station14 <- spCbind(uktowns, station14)
station14 <- merge(station14, popbua, by="BUA14CD")
station14 <- station14[which(station14$poplsoa<175000&station14$poplsoa>10000), ]
station14cent <- SpatialPointsDataFrame(gCentroid(station14, byid=TRUE), station14@data, match.ID=FALSE)
station14 <- station14[-c(2,3,4,5,7)]

jobcent14 <- gov_2014[which(gov_2014$code==06330418), ]
jobcent14 <- subset(jobcent14, "Feature Easting" != "" | "Feature Northing" != "")
jobcent14$poi_ID <- 1:nrow(jobcent14)
coords <- cbind(Easting = as.numeric(as.character(jobcent14$`Feature Easting`)),
                Northing = as.numeric(as.character(jobcent14$`Feature Northing`)))
jobcent14sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent14$Name,
                                                                jobcent14$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent14 <- over(uktowns, geometry(jobcent14sp), returnList = TRUE)
jobcent14 <- sapply(jobcent14, length)
jobcent14 <- spCbind(uktowns, jobcent14)
jobcent14 <- merge(jobcent14, popbua, by="BUA14CD")
jobcent14 <- jobcent14[which(jobcent14$poplsoa<175000&jobcent14$poplsoa>10000), ]
jobcent14cent <- SpatialPointsDataFrame(gCentroid(jobcent14, byid=TRUE), jobcent14@data, match.ID=FALSE)
jobcent14 <- jobcent14[-c(2,3,4,5,7)]

furthered14 <- schools_2014[which(schools_2014$code==05310376), ]
furthered14 <- subset(furthered14, "Feature Easting" != "" | "Feature Northing" != "")
furthered14$poi_ID <- 1:nrow(furthered14)
coords <- cbind(Easting = as.numeric(as.character(furthered14$`Feature Easting`)),
                Northing = as.numeric(as.character(furthered14$`Feature Northing`)))
furthered14sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered14$Name,
                                                                  furthered14$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered14 <- over(uktowns, geometry(furthered14sp), returnList = TRUE)
furthered14 <- sapply(furthered14, length)
furthered14 <- spCbind(uktowns, furthered14)
furthered14 <- merge(furthered14, popbua, by="BUA14CD")
furthered14 <- furthered14[which(furthered14$poplsoa<175000&furthered14$poplsoa>10000), ]
furthered14cent <- SpatialPointsDataFrame(gCentroid(furthered14, byid=TRUE), furthered14@data, match.ID=FALSE)
furthered14 <- furthered14[-c(2,3,4,5,7)]

GPs14 <- health_2014[which(health_2014$code==05280369), ]
GPs14 <- subset(GPs14, "Feature Easting" != "" | "Feature Northing" != "")
GPs14$poi_ID <- 1:nrow(GPs14)
coords <- cbind(Easting = as.numeric(as.character(GPs14$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs14$`Feature Northing`)))
GPs14sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs14$Name,
                                                            GPs14$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs14 <- over(uktowns, geometry(GPs14sp), returnList = TRUE)
GPs14 <- sapply(GPs14, length)
GPs14 <- spCbind(uktowns, GPs14)
GPs14 <- merge(GPs14, popbua, by="BUA14CD")
GPs14 <- GPs14[which(GPs14$poplsoa<175000&GPs14$poplsoa>10000), ]
GPs14cent <- SpatialPointsDataFrame(gCentroid(GPs14, byid=TRUE), GPs14@data, match.ID=FALSE)
GPs14 <- GPs14[-c(2,3,4,5,7)]

hospitals14 <- health_2014[which(health_2014$code==05280371), ]
hospitals14 <- subset(hospitals14, "Feature Easting" != "" | "Feature Northing" != "")
hospitals14$poi_ID <- 1:nrow(hospitals14)
coords <- cbind(Easting = as.numeric(as.character(hospitals14$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals14$`Feature Northing`)))
hospitals14sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals14$Name,
                                                                  hospitals14$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals14 <- over(uktowns, geometry(hospitals14sp), returnList = TRUE)
hospitals14 <- sapply(hospitals14, length)
hospitals14 <- spCbind(uktowns, hospitals14)
hospitals14 <- merge(hospitals14, popbua, by="BUA14CD")
hospitals14 <- hospitals14[which(hospitals14$poplsoa<175000&hospitals14$poplsoa>10000), ]
hospitals14cent <- SpatialPointsDataFrame(gCentroid(hospitals14, byid=TRUE), hospitals14@data, match.ID=FALSE)
hospitals14 <- hospitals14[-c(2,3,4,5,7)]

hmrc14 <- gov_2014[which(gov_2014$code==05280371), ]
hmrc14 <- subset(hmrc14, "Feature Easting" != "" | "Feature Northing" != "")
hmrc14$poi_ID <- 1:nrow(hmrc14)
coords <- cbind(Easting = as.numeric(as.character(hmrc14$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc14$`Feature Northing`)))
hmrc14sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc14$Name,
                                                             hmrc14$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc14 <- over(uktowns, geometry(hmrc14sp), returnList = TRUE)
hmrc14 <- sapply(hmrc14, length)
hmrc14 <- spCbind(uktowns, hmrc14)
hmrc14 <- merge(hmrc14, popbua, by="BUA14CD")
hmrc14 <- hmrc14[which(hmrc14$poplsoa<175000&hmrc14$poplsoa>10000), ]
hmrc14cent <- SpatialPointsDataFrame(gCentroid(hmrc14, byid=TRUE), hmrc14@data, match.ID=FALSE)
hmrc14 <- hmrc14[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals14, by="BUA14CD")
furthered <-  merge(furthered, furthered14, by="BUA14CD")
hmrc <-  merge(hmrc, hmrc14, by="BUA14CD")
gps <-  merge(gps, GPs14, by="BUA14CD")
mental <-  merge(mental, mental14, "BUA14CD")
station <-  merge(station, station14, by="BUA14CD")
jobcent <-  merge(jobcent, jobcent14, by="BUA14CD")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


health_2015 <- read_csv("health_2015.csv", 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2015 <- health_2015 [complete.cases(health_2015$`Feature Easting`), ]
health_2015 <- health_2015 [complete.cases(health_2015$`Feature Northing`), ]
names(health_2015)[names(health_2015)=="PointX Classification Code"] <- "code"

schools_2015 <- read_csv("school_2015.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2015 <- schools_2015 [complete.cases(schools_2015$`Feature Easting`), ]
schools_2015 <- schools_2015 [complete.cases(schools_2015$`Feature Northing`), ]
names(schools_2015)[names(schools_2015)=="PointX Classification Code"] <- "code"

gov_2015 <- read_csv("gov_2015.csv", 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Easting`), ]
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Northing`), ]
names(gov_2015)[names(gov_2015)=="PointX Classification Code"] <- "code"

station_2015 <- read_csv("station_2015.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
station_2015 <- station_2015 [complete.cases(station_2015$`Feature Easting`), ]
station_2015 <- station_2015 [complete.cases(station_2015$`Feature Northing`), ]
names(station_2015)[names(station_2015)=="PointX Classification Code"] <- "code"


mental15 <- health_2015[which(health_2015$code==05280372), ]
mental15 <- subset(mental15, "Feature Easting" != "" | "Feature Northing" != "")
mental15$poi_ID <- 1:nrow(mental15)
coords <- cbind(Easting = as.numeric(as.character(mental15$`Feature Easting`)),
                Northing = as.numeric(as.character(mental15$`Feature Northing`)))
mental15sp <- SpatialPointsDataFrame(coords, data = data.frame(mental15$Name,
                                                               mental15$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental15 <- over(uktowns, geometry(mental15sp), returnList = TRUE)
mental15 <- sapply(mental15, length)
mental15 <- spCbind(uktowns, mental15)
mental15 <- merge(mental15, popbua, by="BUA15CD")
mental15 <- mental15[which(mental15$poplsoa<175000&mental15$poplsoa>10000), ]
mental15cent <- SpatialPointsDataFrame(gCentroid(mental15, byid=TRUE), mental15@data, match.ID=FALSE)
mental15 <- mental15[-c(2,3,4,5,7)]


station15 <- station_2015[which(station_2015$code==010570738), ]
station15 <- subset(station15, "Feature Easting" != "" | "Feature Northing" != "")
station15$poi_ID <- 1:nrow(station15)
coords <- cbind(Easting = as.numeric(as.character(station15$`Feature Easting`)),
                Northing = as.numeric(as.character(station15$`Feature Northing`)))
station15sp <- SpatialPointsDataFrame(coords, data = data.frame(station15$Name,
                                                                station15$poi_ID), proj4string = CRS("+init=epsg:27700"))


station15 <- over(uktowns, geometry(station15sp), returnList = TRUE)
station15 <- sapply(station15, length)
station15 <- spCbind(uktowns, station15)
station15 <- merge(station15, popbua, by="BUA15CD")
station15 <- station15[which(station15$poplsoa<175000&station15$poplsoa>10000), ]
station15cent <- SpatialPointsDataFrame(gCentroid(station15, byid=TRUE), station15@data, match.ID=FALSE)
station15 <- station15[-c(2,3,4,5,7)]

jobcent15 <- gov_2015[which(gov_2015$code==06330418), ]
jobcent15 <- subset(jobcent15, "Feature Easting" != "" | "Feature Northing" != "")
jobcent15$poi_ID <- 1:nrow(jobcent15)
coords <- cbind(Easting = as.numeric(as.character(jobcent15$`Feature Easting`)),
                Northing = as.numeric(as.character(jobcent15$`Feature Northing`)))
jobcent15sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent15$Name,
                                                                jobcent15$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent15 <- over(uktowns, geometry(jobcent15sp), returnList = TRUE)
jobcent15 <- sapply(jobcent15, length)
jobcent15 <- spCbind(uktowns, jobcent15)
jobcent15 <- merge(jobcent15, popbua, by="BUA15CD")
jobcent15 <- jobcent15[which(jobcent15$poplsoa<175000&jobcent15$poplsoa>10000), ]
jobcent15cent <- SpatialPointsDataFrame(gCentroid(jobcent15, byid=TRUE), jobcent15@data, match.ID=FALSE)
jobcent15 <- jobcent15[-c(2,3,4,5,7)]

furthered15 <- schools_2015[which(schools_2015$code==05310376), ]
furthered15 <- subset(furthered15, "Feature Easting" != "" | "Feature Northing" != "")
furthered15$poi_ID <- 1:nrow(furthered15)
coords <- cbind(Easting = as.numeric(as.character(furthered15$`Feature Easting`)),
                Northing = as.numeric(as.character(furthered15$`Feature Northing`)))
furthered15sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered15$Name,
                                                                  furthered15$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered15 <- over(uktowns, geometry(furthered15sp), returnList = TRUE)
furthered15 <- sapply(furthered15, length)
furthered15 <- spCbind(uktowns, furthered15)
furthered15 <- merge(furthered15, popbua, by="BUA15CD")
furthered15 <- furthered15[which(furthered15$poplsoa<175000&furthered15$poplsoa>10000), ]
furthered15cent <- SpatialPointsDataFrame(gCentroid(furthered15, byid=TRUE), furthered15@data, match.ID=FALSE)
furthered15 <- furthered15[-c(2,3,4,5,7)]

GPs15 <- health_2015[which(health_2015$code==05280369), ]
GPs15 <- subset(GPs15, "Feature Easting" != "" | "Feature Northing" != "")
GPs15$poi_ID <- 1:nrow(GPs15)
coords <- cbind(Easting = as.numeric(as.character(GPs15$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs15$`Feature Northing`)))
GPs15sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs15$Name,
                                                            GPs15$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs15 <- over(uktowns, geometry(GPs15sp), returnList = TRUE)
GPs15 <- sapply(GPs15, length)
GPs15 <- spCbind(uktowns, GPs15)
GPs15 <- merge(GPs15, popbua, by="BUA15CD")
GPs15 <- GPs15[which(GPs15$poplsoa<175000&GPs15$poplsoa>10000), ]
GPs15cent <- SpatialPointsDataFrame(gCentroid(GPs15, byid=TRUE), GPs15@data, match.ID=FALSE)
GPs15 <- GPs15[-c(2,3,4,5,7)]

hospitals15 <- health_2015[which(health_2015$code==05280371), ]
hospitals15 <- subset(hospitals15, "Feature Easting" != "" | "Feature Northing" != "")
hospitals15$poi_ID <- 1:nrow(hospitals15)
coords <- cbind(Easting = as.numeric(as.character(hospitals15$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals15$`Feature Northing`)))
hospitals15sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals15$Name,
                                                                  hospitals15$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals15 <- over(uktowns, geometry(hospitals15sp), returnList = TRUE)
hospitals15 <- sapply(hospitals15, length)
hospitals15 <- spCbind(uktowns, hospitals15)
hospitals15 <- merge(hospitals15, popbua, by="BUA15CD")
hospitals15 <- hospitals15[which(hospitals15$poplsoa<175000&hospitals15$poplsoa>10000), ]
hospitals15cent <- SpatialPointsDataFrame(gCentroid(hospitals15, byid=TRUE), hospitals15@data, match.ID=FALSE)
hospitals15 <- hospitals15[-c(2,3,4,5,7)]

hmrc15 <- gov_2015[which(gov_2015$code==05280371), ]
hmrc15 <- subset(hmrc15, "Feature Easting" != "" | "Feature Northing" != "")
hmrc15$poi_ID <- 1:nrow(hmrc15)
coords <- cbind(Easting = as.numeric(as.character(hmrc15$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc15$`Feature Northing`)))
hmrc15sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc15$Name,
                                                             hmrc15$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc15 <- over(uktowns, geometry(hmrc15sp), returnList = TRUE)
hmrc15 <- sapply(hmrc15, length)
hmrc15 <- spCbind(uktowns, hmrc15)
hmrc15 <- merge(hmrc15, popbua, by="BUA15CD")
hmrc15 <- hmrc15[which(hmrc15$poplsoa<175000&hmrc15$poplsoa>10000), ]
hmrc15cent <- SpatialPointsDataFrame(gCentroid(hmrc15, byid=TRUE), hmrc15@data, match.ID=FALSE)
hmrc15 <- hmrc15[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals15, by="BUA15CD")
furthered <-  merge(furthered, furthered15, by="BUA15CD")
hmrc <-  merge(hmrc, hmrc15, by="BUA15CD")
gps <-  merge(gps, GPs15, by="BUA15CD")
mental <-  merge(mental, mental15, "BUA15CD")
station <-  merge(station, station15, by="BUA15CD")
jobcent <-  merge(jobcent, jobcent15, by="BUA15CD")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


health_2016 <- read_csv("health_2016.csv", 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2016 <- health_2016 [complete.cases(health_2016$`Feature Easting`), ]
health_2016 <- health_2016 [complete.cases(health_2016$`Feature Northing`), ]
names(health_2016)[names(health_2016)=="PointX Classification Code"] <- "code"

schools_2016 <- read_csv("school_2016.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2016 <- schools_2016 [complete.cases(schools_2016$`Feature Easting`), ]
schools_2016 <- schools_2016 [complete.cases(schools_2016$`Feature Northing`), ]
names(schools_2016)[names(schools_2016)=="PointX Classification Code"] <- "code"

gov_2016 <- read_csv("gov_2016.csv", 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Easting`), ]
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Northing`), ]
names(gov_2016)[names(gov_2016)=="PointX Classification Code"] <- "code"

station_2016 <- read_csv("station_2016.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
station_2016 <- station_2016 [complete.cases(station_2016$`Feature Easting`), ]
station_2016 <- station_2016 [complete.cases(station_2016$`Feature Northing`), ]
names(station_2016)[names(station_2016)=="PointX Classification Code"] <- "code"


mental16 <- health_2016[which(health_2016$code==05280372), ]
mental16 <- subset(mental16, "Feature Easting" != "" | "Feature Northing" != "")
mental16$poi_ID <- 1:nrow(mental16)
coords <- cbind(Easting = as.numeric(as.character(mental16$`Feature Easting`)),
                Northing = as.numeric(as.character(mental16$`Feature Northing`)))
mental16sp <- SpatialPointsDataFrame(coords, data = data.frame(mental16$Name,
                                                               mental16$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental16 <- over(uktowns, geometry(mental16sp), returnList = TRUE)
mental16 <- sapply(mental16, length)
mental16 <- spCbind(uktowns, mental16)
mental16 <- merge(mental16, popbua, by="BUA16CD")
mental16 <- mental16[which(mental16$poplsoa<175000&mental16$poplsoa>10000), ]
mental16cent <- SpatialPointsDataFrame(gCentroid(mental16, byid=TRUE), mental16@data, match.ID=FALSE)
mental16 <- mental16[-c(2,3,4,5,7)]


station16 <- station_2016[which(station_2016$code==010570738), ]
station16 <- subset(station16, "Feature Easting" != "" | "Feature Northing" != "")
station16$poi_ID <- 1:nrow(station16)
coords <- cbind(Easting = as.numeric(as.character(station16$`Feature Easting`)),
                Northing = as.numeric(as.character(station16$`Feature Northing`)))
station16sp <- SpatialPointsDataFrame(coords, data = data.frame(station16$Name,
                                                                station16$poi_ID), proj4string = CRS("+init=epsg:27700"))


station16 <- over(uktowns, geometry(station16sp), returnList = TRUE)
station16 <- sapply(station16, length)
station16 <- spCbind(uktowns, station16)
station16 <- merge(station16, popbua, by="BUA16CD")
station16 <- station16[which(station16$poplsoa<175000&station16$poplsoa>10000), ]
station16cent <- SpatialPointsDataFrame(gCentroid(station16, byid=TRUE), station16@data, match.ID=FALSE)
station16 <- station16[-c(2,3,4,5,7)]

jobcent16 <- gov_2016[which(gov_2016$code==06330418), ]
jobcent16 <- subset(jobcent16, "Feature Easting" != "" | "Feature Northing" != "")
jobcent16$poi_ID <- 1:nrow(jobcent16)
coords <- cbind(Easting = as.numeric(as.character(jobcent16$`Feature Easting`)),
                Northing = as.numeric(as.character(jobcent16$`Feature Northing`)))
jobcent16sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent16$Name,
                                                                jobcent16$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent16 <- over(uktowns, geometry(jobcent16sp), returnList = TRUE)
jobcent16 <- sapply(jobcent16, length)
jobcent16 <- spCbind(uktowns, jobcent16)
jobcent16 <- merge(jobcent16, popbua, by="BUA16CD")
jobcent16 <- jobcent16[which(jobcent16$poplsoa<175000&jobcent16$poplsoa>10000), ]
jobcent16cent <- SpatialPointsDataFrame(gCentroid(jobcent16, byid=TRUE), jobcent16@data, match.ID=FALSE)
jobcent16 <- jobcent16[-c(2,3,4,5,7)]

furthered16 <- schools_2016[which(schools_2016$code==05310376), ]
furthered16 <- subset(furthered16, "Feature Easting" != "" | "Feature Northing" != "")
furthered16$poi_ID <- 1:nrow(furthered16)
coords <- cbind(Easting = as.numeric(as.character(furthered16$`Feature Easting`)),
                Northing = as.numeric(as.character(furthered16$`Feature Northing`)))
furthered16sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered16$Name,
                                                                  furthered16$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered16 <- over(uktowns, geometry(furthered16sp), returnList = TRUE)
furthered16 <- sapply(furthered16, length)
furthered16 <- spCbind(uktowns, furthered16)
furthered16 <- merge(furthered16, popbua, by="BUA16CD")
furthered16 <- furthered16[which(furthered16$poplsoa<175000&furthered16$poplsoa>10000), ]
furthered16cent <- SpatialPointsDataFrame(gCentroid(furthered16, byid=TRUE), furthered16@data, match.ID=FALSE)
furthered16 <- furthered16[-c(2,3,4,5,7)]

GPs16 <- health_2016[which(health_2016$code==05280369), ]
GPs16 <- subset(GPs16, "Feature Easting" != "" | "Feature Northing" != "")
GPs16$poi_ID <- 1:nrow(GPs16)
coords <- cbind(Easting = as.numeric(as.character(GPs16$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs16$`Feature Northing`)))
GPs16sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs16$Name,
                                                            GPs16$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs16 <- over(uktowns, geometry(GPs16sp), returnList = TRUE)
GPs16 <- sapply(GPs16, length)
GPs16 <- spCbind(uktowns, GPs16)
GPs16 <- merge(GPs16, popbua, by="BUA16CD")
GPs16 <- GPs16[which(GPs16$poplsoa<175000&GPs16$poplsoa>10000), ]
GPs16cent <- SpatialPointsDataFrame(gCentroid(GPs16, byid=TRUE), GPs16@data, match.ID=FALSE)
GPs16 <- GPs16[-c(2,3,4,5,7)]

hospitals16 <- health_2016[which(health_2016$code==05280371), ]
hospitals16 <- subset(hospitals16, "Feature Easting" != "" | "Feature Northing" != "")
hospitals16$poi_ID <- 1:nrow(hospitals16)
coords <- cbind(Easting = as.numeric(as.character(hospitals16$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals16$`Feature Northing`)))
hospitals16sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals16$Name,
                                                                  hospitals16$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals16 <- over(uktowns, geometry(hospitals16sp), returnList = TRUE)
hospitals16 <- sapply(hospitals16, length)
hospitals16 <- spCbind(uktowns, hospitals16)
hospitals16 <- merge(hospitals16, popbua, by="BUA16CD")
hospitals16 <- hospitals16[which(hospitals16$poplsoa<175000&hospitals16$poplsoa>10000), ]
hospitals16cent <- SpatialPointsDataFrame(gCentroid(hospitals16, byid=TRUE), hospitals16@data, match.ID=FALSE)
hospitals16 <- hospitals16[-c(2,3,4,5,7)]

hmrc16 <- gov_2016[which(gov_2016$code==05280371), ]
hmrc16 <- subset(hmrc16, "Feature Easting" != "" | "Feature Northing" != "")
hmrc16$poi_ID <- 1:nrow(hmrc16)
coords <- cbind(Easting = as.numeric(as.character(hmrc16$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc16$`Feature Northing`)))
hmrc16sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc16$Name,
                                                             hmrc16$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc16 <- over(uktowns, geometry(hmrc16sp), returnList = TRUE)
hmrc16 <- sapply(hmrc16, length)
hmrc16 <- spCbind(uktowns, hmrc16)
hmrc16 <- merge(hmrc16, popbua, by="BUA16CD")
hmrc16 <- hmrc16[which(hmrc16$poplsoa<175000&hmrc16$poplsoa>10000), ]
hmrc16cent <- SpatialPointsDataFrame(gCentroid(hmrc16, byid=TRUE), hmrc16@data, match.ID=FALSE)
hmrc16 <- hmrc16[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals16, by="BUA16CD")
furthered <-  merge(furthered, furthered16, by="BUA16CD")
hmrc <-  merge(hmrc, hmrc16, by="BUA16CD")
gps <-  merge(gps, GPs16, by="BUA16CD")
mental <-  merge(mental, mental16, "BUA16CD")
station <-  merge(station, station16, by="BUA16CD")
jobcent <-  merge(jobcent, jobcent16, by="BUA16CD")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


health_2017 <- read_csv("health_2017.csv", 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2017 <- health_2017 [complete.cases(health_2017$`Feature Easting`), ]
health_2017 <- health_2017 [complete.cases(health_2017$`Feature Northing`), ]
names(health_2017)[names(health_2017)=="PointX Classification Code"] <- "code"

schools_2017 <- read_csv("school_2017.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2017 <- schools_2017 [complete.cases(schools_2017$`Feature Easting`), ]
schools_2017 <- schools_2017 [complete.cases(schools_2017$`Feature Northing`), ]
names(schools_2017)[names(schools_2017)=="PointX Classification Code"] <- "code"

gov_2017 <- read_csv("gov_2017.csv", 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Easting`), ]
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Northing`), ]
names(gov_2017)[names(gov_2017)=="PointX Classification Code"] <- "code"

station_2017 <- read_csv("station_2017.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
station_2017 <- station_2017 [complete.cases(station_2017$`Feature Easting`), ]
station_2017 <- station_2017 [complete.cases(station_2017$`Feature Northing`), ]
names(station_2017)[names(station_2017)=="PointX Classification Code"] <- "code"


mental17 <- health_2017[which(health_2017$code==05280372), ]
mental17 <- subset(mental17, "Feature Easting" != "" | "Feature Northing" != "")
mental17$poi_ID <- 1:nrow(mental17)
coords <- cbind(Easting = as.numeric(as.character(mental17$`Feature Easting`)),
                Northing = as.numeric(as.character(mental17$`Feature Northing`)))
mental17sp <- SpatialPointsDataFrame(coords, data = data.frame(mental17$Name,
                                                               mental17$poi_ID), proj4string = CRS("+init=epsg:27700"))


mental17 <- over(uktowns, geometry(mental17sp), returnList = TRUE)
mental17 <- sapply(mental17, length)
mental17 <- spCbind(uktowns, mental17)
mental17 <- merge(mental17, popbua, by="BUA17CD")
mental17 <- mental17[which(mental17$poplsoa<175000&mental17$poplsoa>10000), ]
mental17cent <- SpatialPointsDataFrame(gCentroid(mental17, byid=TRUE), mental17@data, match.ID=FALSE)
mental17 <- mental17[-c(2,3,4,5,7)]


station17 <- station_2017[which(station_2017$code==010570738), ]
station17 <- subset(station17, "Feature Easting" != "" | "Feature Northing" != "")
station17$poi_ID <- 1:nrow(station17)
coords <- cbind(Easting = as.numeric(as.character(station17$`Feature Easting`)),
                Northing = as.numeric(as.character(station17$`Feature Northing`)))
station17sp <- SpatialPointsDataFrame(coords, data = data.frame(station17$Name,
                                                                station17$poi_ID), proj4string = CRS("+init=epsg:27700"))


station17 <- over(uktowns, geometry(station17sp), returnList = TRUE)
station17 <- sapply(station17, length)
station17 <- spCbind(uktowns, station17)
station17 <- merge(station17, popbua, by="BUA17CD")
station17 <- station17[which(station17$poplsoa<175000&station17$poplsoa>10000), ]
station17cent <- SpatialPointsDataFrame(gCentroid(station17, byid=TRUE), station17@data, match.ID=FALSE)
station17 <- station17[-c(2,3,4,5,7)]

jobcent17 <- gov_2017[which(gov_2017$code==06330418), ]
jobcent17 <- subset(jobcent17, "Feature Easting" != "" | "Feature Northing" != "")
jobcent17$poi_ID <- 1:nrow(jobcent17)
coords <- cbind(Easting = as.numeric(as.character(jobcent17$`Feature Easting`)),
                Northing = as.numeric(as.character(jobcent17$`Feature Northing`)))
jobcent17sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent17$Name,
                                                                jobcent17$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent17 <- over(uktowns, geometry(jobcent17sp), returnList = TRUE)
jobcent17 <- sapply(jobcent17, length)
jobcent17 <- spCbind(uktowns, jobcent17)
jobcent17 <- merge(jobcent17, popbua, by="BUA17CD")
jobcent17 <- jobcent17[which(jobcent17$poplsoa<175000&jobcent17$poplsoa>10000), ]
jobcent17cent <- SpatialPointsDataFrame(gCentroid(jobcent17, byid=TRUE), jobcent17@data, match.ID=FALSE)
jobcent17 <- jobcent17[-c(2,3,4,5,7)]

furthered17 <- schools_2017[which(schools_2017$code==05310376), ]
furthered17 <- subset(furthered17, "Feature Easting" != "" | "Feature Northing" != "")
furthered17$poi_ID <- 1:nrow(furthered17)
coords <- cbind(Easting = as.numeric(as.character(furthered17$`Feature Easting`)),
                Northing = as.numeric(as.character(furthered17$`Feature Northing`)))
furthered17sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered17$Name,
                                                                  furthered17$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered17 <- over(uktowns, geometry(furthered17sp), returnList = TRUE)
furthered17 <- sapply(furthered17, length)
furthered17 <- spCbind(uktowns, furthered17)
furthered17 <- merge(furthered17, popbua, by="BUA17CD")
furthered17 <- furthered17[which(furthered17$poplsoa<175000&furthered17$poplsoa>10000), ]
furthered17cent <- SpatialPointsDataFrame(gCentroid(furthered17, byid=TRUE), furthered17@data, match.ID=FALSE)
furthered17 <- furthered17[-c(2,3,4,5,7)]

GPs17 <- health_2017[which(health_2017$code==05280369), ]
GPs17 <- subset(GPs17, "Feature Easting" != "" | "Feature Northing" != "")
GPs17$poi_ID <- 1:nrow(GPs17)
coords <- cbind(Easting = as.numeric(as.character(GPs17$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs17$`Feature Northing`)))
GPs17sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs17$Name,
                                                            GPs17$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs17 <- over(uktowns, geometry(GPs17sp), returnList = TRUE)
GPs17 <- sapply(GPs17, length)
GPs17 <- spCbind(uktowns, GPs17)
GPs17 <- merge(GPs17, popbua, by="BUA17CD")
GPs17 <- GPs17[which(GPs17$poplsoa<175000&GPs17$poplsoa>10000), ]
GPs17cent <- SpatialPointsDataFrame(gCentroid(GPs17, byid=TRUE), GPs17@data, match.ID=FALSE)
GPs17 <- GPs17[-c(2,3,4,5,7)]

hospitals17 <- health_2017[which(health_2017$code==05280371), ]
hospitals17 <- subset(hospitals17, "Feature Easting" != "" | "Feature Northing" != "")
hospitals17$poi_ID <- 1:nrow(hospitals17)
coords <- cbind(Easting = as.numeric(as.character(hospitals17$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals17$`Feature Northing`)))
hospitals17sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals17$Name,
                                                                  hospitals17$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals17 <- over(uktowns, geometry(hospitals17sp), returnList = TRUE)
hospitals17 <- sapply(hospitals17, length)
hospitals17 <- spCbind(uktowns, hospitals17)
hospitals17 <- merge(hospitals17, popbua, by="BUA17CD")
hospitals17 <- hospitals17[which(hospitals17$poplsoa<175000&hospitals17$poplsoa>10000), ]
hospitals17cent <- SpatialPointsDataFrame(gCentroid(hospitals17, byid=TRUE), hospitals17@data, match.ID=FALSE)
hospitals17 <- hospitals17[-c(2,3,4,5,7)]

hmrc17 <- gov_2017[which(gov_2017$code==05280371), ]
hmrc17 <- subset(hmrc17, "Feature Easting" != "" | "Feature Northing" != "")
hmrc17$poi_ID <- 1:nrow(hmrc17)
coords <- cbind(Easting = as.numeric(as.character(hmrc17$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc17$`Feature Northing`)))
hmrc17sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc17$Name,
                                                             hmrc17$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc17 <- over(uktowns, geometry(hmrc17sp), returnList = TRUE)
hmrc17 <- sapply(hmrc17, length)
hmrc17 <- spCbind(uktowns, hmrc17)
hmrc17 <- merge(hmrc17, popbua, by="BUA17CD")
hmrc17 <- hmrc17[which(hmrc17$poplsoa<175000&hmrc17$poplsoa>10000), ]
hmrc17cent <- SpatialPointsDataFrame(gCentroid(hmrc17, byid=TRUE), hmrc17@data, match.ID=FALSE)
hmrc17 <- hmrc17[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals17, by="BUA17CD")
furthered <-  merge(furthered, furthered17, by="BUA17CD")
hmrc <-  merge(hmrc, hmrc17, by="BUA17CD")
gps <-  merge(gps, GPs17, by="BUA17CD")
mental <-  merge(mental, mental17, "BUA17CD")
station <-  merge(station, station17, by="BUA17CD")
jobcent <-  merge(jobcent, jobcent17, by="BUA17CD")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


health_2018 <- read_csv("health_2018.csv", 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Easting`), ]
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Northing`), ]
names(health_2018)[names(health_2018)=="PointX Classification Code"] <- "code"

schools_2018 <- read_csv("school_2018.csv", 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2018 <- schools_2018 [complete.cases(schools_2018$`Feature Easting`), ]
schools_2018 <- schools_2018 [complete.cases(schools_2018$`Feature Northing`), ]
names(schools_2018)[names(schools_2018)=="PointX Classification Code"] <- "code"

gov_2018 <- read_csv("gov_2018.csv", 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"

station_2018 <- read_csv("station_2018.csv", 
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


mental18 <- over(uktowns, geometry(mental18sp), returnList = TRUE)
mental18 <- sapply(mental18, length)
mental18 <- spCbind(uktowns, mental18)
mental18 <- merge(mental18, popbua, by="BUA18CD")
mental18 <- mental18[which(mental18$poplsoa<175000&mental18$poplsoa>10000), ]
mental18cent <- SpatialPointsDataFrame(gCentroid(mental18, byid=TRUE), mental18@data, match.ID=FALSE)
mental18 <- mental18[-c(2,3,4,5,7)]


station18 <- station_2018[which(station_2018$code==010570738), ]
station18 <- subset(station18, "Feature Easting" != "" | "Feature Northing" != "")
station18$poi_ID <- 1:nrow(station18)
coords <- cbind(Easting = as.numeric(as.character(station18$`Feature Easting`)),
                Northing = as.numeric(as.character(station18$`Feature Northing`)))
station18sp <- SpatialPointsDataFrame(coords, data = data.frame(station18$Name,
                                                                station18$poi_ID), proj4string = CRS("+init=epsg:27700"))


station18 <- over(uktowns, geometry(station18sp), returnList = TRUE)
station18 <- sapply(station18, length)
station18 <- spCbind(uktowns, station18)
station18 <- merge(station18, popbua, by="BUA18CD")
station18 <- station18[which(station18$poplsoa<175000&station18$poplsoa>10000), ]
station18cent <- SpatialPointsDataFrame(gCentroid(station18, byid=TRUE), station18@data, match.ID=FALSE)
station18 <- station18[-c(2,3,4,5,7)]

jobcent18 <- gov_2018[which(gov_2018$code==06330418), ]
jobcent18 <- subset(jobcent18, "Feature Easting" != "" | "Feature Northing" != "")
jobcent18$poi_ID <- 1:nrow(jobcent18)
coords <- cbind(Easting = as.numeric(as.character(jobcent18$`Feature Easting`)),
                Northing = as.numeric(as.character(jobcent18$`Feature Northing`)))
jobcent18sp <- SpatialPointsDataFrame(coords, data = data.frame(jobcent18$Name,
                                                                jobcent18$poi_ID), proj4string = CRS("+init=epsg:27700"))


jobcent18 <- over(uktowns, geometry(jobcent18sp), returnList = TRUE)
jobcent18 <- sapply(jobcent18, length)
jobcent18 <- spCbind(uktowns, jobcent18)
jobcent18 <- merge(jobcent18, popbua, by="BUA18CD")
jobcent18 <- jobcent18[which(jobcent18$poplsoa<175000&jobcent18$poplsoa>10000), ]
jobcent18cent <- SpatialPointsDataFrame(gCentroid(jobcent18, byid=TRUE), jobcent18@data, match.ID=FALSE)
jobcent18 <- jobcent18[-c(2,3,4,5,7)]

furthered18 <- schools_2018[which(schools_2018$code==05310376), ]
furthered18 <- subset(furthered18, "Feature Easting" != "" | "Feature Northing" != "")
furthered18$poi_ID <- 1:nrow(furthered18)
coords <- cbind(Easting = as.numeric(as.character(furthered18$`Feature Easting`)),
                Northing = as.numeric(as.character(furthered18$`Feature Northing`)))
furthered18sp <- SpatialPointsDataFrame(coords, data = data.frame(furthered18$Name,
                                                                  furthered18$poi_ID), proj4string = CRS("+init=epsg:27700"))


furthered18 <- over(uktowns, geometry(furthered18sp), returnList = TRUE)
furthered18 <- sapply(furthered18, length)
furthered18 <- spCbind(uktowns, furthered18)
furthered18 <- merge(furthered18, popbua, by="BUA18CD")
furthered18 <- furthered18[which(furthered18$poplsoa<175000&furthered18$poplsoa>10000), ]
furthered18cent <- SpatialPointsDataFrame(gCentroid(furthered18, byid=TRUE), furthered18@data, match.ID=FALSE)
furthered18 <- furthered18[-c(2,3,4,5,7)]

GPs18 <- health_2018[which(health_2018$code==05280369), ]
GPs18 <- subset(GPs18, "Feature Easting" != "" | "Feature Northing" != "")
GPs18$poi_ID <- 1:nrow(GPs18)
coords <- cbind(Easting = as.numeric(as.character(GPs18$`Feature Easting`)),
                Northing = as.numeric(as.character(GPs18$`Feature Northing`)))
GPs18sp <- SpatialPointsDataFrame(coords, data = data.frame(GPs18$Name,
                                                            GPs18$poi_ID), proj4string = CRS("+init=epsg:27700"))


GPs18 <- over(uktowns, geometry(GPs18sp), returnList = TRUE)
GPs18 <- sapply(GPs18, length)
GPs18 <- spCbind(uktowns, GPs18)
GPs18 <- merge(GPs18, popbua, by="BUA18CD")
GPs18 <- GPs18[which(GPs18$poplsoa<175000&GPs18$poplsoa>10000), ]
GPs18cent <- SpatialPointsDataFrame(gCentroid(GPs18, byid=TRUE), GPs18@data, match.ID=FALSE)
GPs18 <- GPs18[-c(2,3,4,5,7)]

hospitals18 <- health_2018[which(health_2018$code==05280371), ]
hospitals18 <- subset(hospitals18, "Feature Easting" != "" | "Feature Northing" != "")
hospitals18$poi_ID <- 1:nrow(hospitals18)
coords <- cbind(Easting = as.numeric(as.character(hospitals18$`Feature Easting`)),
                Northing = as.numeric(as.character(hospitals18$`Feature Northing`)))
hospitals18sp <- SpatialPointsDataFrame(coords, data = data.frame(hospitals18$Name,
                                                                  hospitals18$poi_ID), proj4string = CRS("+init=epsg:27700"))


hospitals18 <- over(uktowns, geometry(hospitals18sp), returnList = TRUE)
hospitals18 <- sapply(hospitals18, length)
hospitals18 <- spCbind(uktowns, hospitals18)
hospitals18 <- merge(hospitals18, popbua, by="BUA18CD")
hospitals18 <- hospitals18[which(hospitals18$poplsoa<175000&hospitals18$poplsoa>10000), ]
hospitals18cent <- SpatialPointsDataFrame(gCentroid(hospitals18, byid=TRUE), hospitals18@data, match.ID=FALSE)
hospitals18 <- hospitals18[-c(2,3,4,5,7)]

hmrc18 <- gov_2018[which(gov_2018$code==05280371), ]
hmrc18 <- subset(hmrc18, "Feature Easting" != "" | "Feature Northing" != "")
hmrc18$poi_ID <- 1:nrow(hmrc18)
coords <- cbind(Easting = as.numeric(as.character(hmrc18$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc18$`Feature Northing`)))
hmrc18sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc18$Name,
                                                             hmrc18$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc18 <- over(uktowns, geometry(hmrc18sp), returnList = TRUE)
hmrc18 <- sapply(hmrc18, length)
hmrc18 <- spCbind(uktowns, hmrc18)
hmrc18 <- merge(hmrc18, popbua, by="BUA18CD")
hmrc18 <- hmrc18[which(hmrc18$poplsoa<175000&hmrc18$poplsoa>10000), ]
hmrc18cent <- SpatialPointsDataFrame(gCentroid(hmrc18, byid=TRUE), hmrc18@data, match.ID=FALSE)
hmrc18 <- hmrc18[-c(2,3,4,5,7)]



hospital <-  merge(hospital, hospitals18, by="BUA18CD")
furthered <-  merge(furthered, furthered18, by="BUA18CD")
hmrc <-  merge(hmrc, hmrc18, by="BUA18CD")
gps <-  merge(gps, GPs18, by="BUA18CD")
mental <-  merge(mental, mental18, "BUA18CD")
station <-  merge(station, station18, by="BUA18CD")
jobcent <-  merge(jobcent, jobcent18, by="BUA18CD")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])





hospital$"2018-06-01" <- hospital$hospital18/(hospital$poplsoa/1000)
hospital$"2017-06-01" <- hospital$hospital17/(hospital$poplsoa/1000)
hospital$"2016-06-01" <- hospital$hospital16/(hospital$poplsoa/1000)
hospital$"2015-06-01" <- hospital$hospital15/(hospital$poplsoa/1000)
hospital$"2014-09-01" <- hospital$hospital14/(hospital$poplsoa/1000)
hospital$"2013-06-01" <- hospital$hospital13/(hospital$poplsoa/1000)
hospital$"2012-06-01" <- hospital$hospital12/(hospital$poplsoa/1000)
hospital$"2011-06-01" <- hospital$hospital11/(hospital$poplsoa/1000)
hospital$"2010-09-01" <- hospital$hospital10/(hospital$poplsoa/1000)
hospital$"2009-06-01" <- hospital$hospital09/(hospital$poplsoa/1000)
hospital$"2008-03-01" <- hospital$hospital08/(hospital$poplsoa/1000)

jobcent$"2018-06-01" <- jobcent$jobcent18/(jobcent$poplsoa/1000)
jobcent$"2017-06-01" <- jobcent$jobcent17/(jobcent$poplsoa/1000)
jobcent$"2016-06-01" <- jobcent$jobcent16/(jobcent$poplsoa/1000)
jobcent$"2015-06-01" <- jobcent$jobcent15/(jobcent$poplsoa/1000)
jobcent$"2014-09-01" <- jobcent$jobcent14/(jobcent$poplsoa/1000)
jobcent$"2013-06-01" <- jobcent$jobcent13/(jobcent$poplsoa/1000)
jobcent$"2012-06-01" <- jobcent$jobcent12/(jobcent$poplsoa/1000)
jobcent$"2011-06-01" <- jobcent$jobcent11/(jobcent$poplsoa/1000)
jobcent$"2010-09-01" <- jobcent$jobcent10/(jobcent$poplsoa/1000)
jobcent$"2009-06-01" <- jobcent$jobcent09/(jobcent$poplsoa/1000)
jobcent$"2008-03-01" <- jobcent$jobcent08/(jobcent$poplsoa/1000)

hmrc$"2018-06-01" <- hmrc$hmrc18/(hmrc$poplsoa/1000)
hmrc$"2017-06-01" <- hmrc$hmrc17/(hmrc$poplsoa/1000)
hmrc$"2016-06-01" <- hmrc$hmrc16/(hmrc$poplsoa/1000)
hmrc$"2015-06-01" <- hmrc$hmrc15/(hmrc$poplsoa/1000)
hmrc$"2014-09-01" <- hmrc$hmrc14/(hmrc$poplsoa/1000)
hmrc$"2013-06-01" <- hmrc$hmrc13/(hmrc$poplsoa/1000)
hmrc$"2012-06-01" <- hmrc$hmrc12/(hmrc$poplsoa/1000)
hmrc$"2011-06-01" <- hmrc$hmrc11/(hmrc$poplsoa/1000)
hmrc$"2010-09-01" <- hmrc$hmrc10/(hmrc$poplsoa/1000)
hmrc$"2009-06-01" <- hmrc$hmrc09/(hmrc$poplsoa/1000)
hmrc$"2008-03-01" <- hmrc$hmrc08/(hmrc$poplsoa/1000)

mental$"2018-06-01" <- mental$mental18/(mental$poplsoa/1000)
mental$"2017-06-01" <- mental$mental17/(mental$poplsoa/1000)
mental$"2016-06-01" <- mental$mental16/(mental$poplsoa/1000)
mental$"2015-06-01" <- mental$mental15/(mental$poplsoa/1000)
mental$"2014-09-01" <- mental$mental14/(mental$poplsoa/1000)
mental$"2013-06-01" <- mental$mental13/(mental$poplsoa/1000)
mental$"2012-06-01" <- mental$mental12/(mental$poplsoa/1000)
mental$"2011-06-01" <- mental$mental11/(mental$poplsoa/1000)
mental$"2010-09-01" <- mental$mental10/(mental$poplsoa/1000)
mental$"2009-06-01" <- mental$mental09/(mental$poplsoa/1000)
mental$"2008-03-01" <- mental$mental08/(mental$poplsoa/1000)

gps$"2018-06-01" <- gps$gps18/(gps$poplsoa/1000)
gps$"2017-06-01" <- gps$gps17/(gps$poplsoa/1000)
gps$"2016-06-01" <- gps$gps16/(gps$poplsoa/1000)
gps$"2015-06-01" <- gps$gps15/(gps$poplsoa/1000)
gps$"2014-09-01" <- gps$gps14/(gps$poplsoa/1000)
gps$"2013-06-01" <- gps$gps13/(gps$poplsoa/1000)
gps$"2012-06-01" <- gps$gps12/(gps$poplsoa/1000)
gps$"2011-06-01" <- gps$gps11/(gps$poplsoa/1000)
gps$"2010-09-01" <- gps$gps10/(gps$poplsoa/1000)
gps$"2009-06-01" <- gps$gps09/(gps$poplsoa/1000)
gps$"2008-03-01" <- gps$gps08/(gps$poplsoa/1000)

station$"2018-06-01" <- station$station18/(station$poplsoa/1000)
station$"2017-06-01" <- station$station17/(station$poplsoa/1000)
station$"2016-06-01" <- station$station16/(station$poplsoa/1000)
station$"2015-06-01" <- station$station15/(station$poplsoa/1000)
station$"2014-09-01" <- station$station14/(station$poplsoa/1000)
station$"2013-06-01" <- station$station13/(station$poplsoa/1000)
station$"2012-06-01" <- station$station12/(station$poplsoa/1000)
station$"2011-06-01" <- station$station11/(station$poplsoa/1000)
station$"2010-09-01" <- station$station10/(station$poplsoa/1000)
station$"2009-06-01" <- station$station09/(station$poplsoa/1000)
station$"2008-03-01" <- station$station08/(station$poplsoa/1000)

furthered$"2018-06-01" <- furthered$furthered18/(furthered$poplsoa/1000)
furthered$"2017-06-01" <- furthered$furthered17/(furthered$poplsoa/1000)
furthered$"2016-06-01" <- furthered$furthered16/(furthered$poplsoa/1000)
furthered$"2015-06-01" <- furthered$furthered15/(furthered$poplsoa/1000)
furthered$"2014-09-01" <- furthered$furthered14/(furthered$poplsoa/1000)
furthered$"2013-06-01" <- furthered$furthered13/(furthered$poplsoa/1000)
furthered$"2012-06-01" <- furthered$furthered12/(furthered$poplsoa/1000)
furthered$"2011-06-01" <- furthered$furthered11/(furthered$poplsoa/1000)
furthered$"2010-09-01" <- furthered$furthered10/(furthered$poplsoa/1000)
furthered$"2009-06-01" <- furthered$furthered09/(furthered$poplsoa/1000)
furthered$"2008-03-01" <- furthered$furthered08/(furthered$poplsoa/1000)





POI_2010 <- read_csv("D:/poi/POI_2010.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2010 <- POI_2010 [complete.cases(POI_2010$easting), ]
POI_2010 <- POI_2010 [complete.cases(POI_2010$northing), ]
names(POI_2010)[names(POI_2010)=="pointx_classification_code"] <- "code"



police10 <- POI_2010[which(POI_2010$code==6330422), ]
police10 <- subset(police10, "easting" != "" | "northing" != "")
police10$poi_ID <- 1:nrow(police10)
coords <- cbind(Easting = as.numeric(as.character(police10$easting)),
                Northing = as.numeric(as.character(police10$northing)))
police10sp <- SpatialPointsDataFrame(coords, data = data.frame(police10$name,
                                                               police10$poi_ID), proj4string = CRS("+init=epsg:27700"))


police10 <- over(uktowns, geometry(police10sp), returnList = TRUE)
police10 <- sapply(police10, length)
police10 <- spCbind(uktowns, police10)
police10 <- merge(police10, popbua, by="BUA11CD")
police10 <- police10[which(police10$poplsoa<175000&police10$poplsoa>10000), ]
police10cent <- SpatialPointsDataFrame(gCentroid(police10, byid=TRUE), police10@data, match.ID=FALSE)
police10 <- police10[-c(2,3,4,5,7)]

write.csv(police10, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police10.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


POI_2011 <- read_csv("D:/poi/POI_2011.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2011 <- as.data.frame(POI_2011)
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_easting'), ]
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_northing'), ]
names(POI_2011)[names(POI_2011)=="pointx_classification_code"] <- "code"


police11 <- POI_2011[which(POI_2011$code==6330422), ]
police11 <- subset(police11, "`feature_easting`" != "" | "'feature_northing'" != "")
police11$poi_ID <- 1:nrow(police11)
coords <- cbind(Easting = as.numeric(as.character(police11$'feature_easting')),
                Northing = as.numeric(as.character(police11$'feature_northing')))
police11sp <- SpatialPointsDataFrame(coords, data = data.frame(police11$name,
                                                               police11$poi_ID), proj4string = CRS("+init=epsg:27700"))


police11 <- over(uktowns, geometry(police11sp), returnList = TRUE)
police11 <- sapply(police11, length)
police11 <- spCbind(uktowns, police11)
police11 <- merge(police11, popbua, by="BUA11CD")
police11 <- police11[which(police11$poplsoa<175000&police11$poplsoa>10000), ]
police11cent <- SpatialPointsDataFrame(gCentroid(police11, byid=TRUE), police11@data, match.ID=FALSE)
police11 <- police11[-c(2,3,4,5,7)]

write.csv(police11, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police11.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2009 <- read_csv("D:/poi/POI_2009.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2009 <- POI_2009 [complete.cases(POI_2009$easting), ]
POI_2009 <- POI_2009 [complete.cases(POI_2009$northing), ]
names(POI_2009)[names(POI_2009)=="pointx_classification_code"] <- "code"



police09 <- POI_2009[which(POI_2009$code==6330422), ]
police09 <- subset(police09, "easting" != "" | "northing" != "")
police09$poi_ID <- 1:nrow(police09)
coords <- cbind(Easting = as.numeric(as.character(police09$easting)),
                Northing = as.numeric(as.character(police09$northing)))
police09sp <- SpatialPointsDataFrame(coords, data = data.frame(police09$name,
                                                               police09$poi_ID), proj4string = CRS("+init=epsg:27700"))


police09 <- over(uktowns, geometry(police09sp), returnList = TRUE)
police09 <- sapply(police09, length)
police09 <- spCbind(uktowns, police09)
police09 <- merge(police09, popbua, by="BUA11CD")
police09 <- police09[which(police09$poplsoa<175000&police09$poplsoa>10000), ]
police09cent <- SpatialPointsDataFrame(gCentroid(police09, byid=TRUE), police09@data, match.ID=FALSE)
police09 <- police09[-c(2,3,4,5,7)]

write.csv(police09, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police09.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2008 <- read_csv("D:/poi/POI_2008.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2008 <- POI_2008 [complete.cases(POI_2008$easting), ]
POI_2008 <- POI_2008 [complete.cases(POI_2008$northing), ]
names(POI_2008)[names(POI_2008)=="pointx_classification_code"] <- "code"



police08 <- POI_2008[which(POI_2008$code==6330422), ]
police08 <- subset(police08, "easting" != "" | "northing" != "")
police08$poi_ID <- 1:nrow(police08)
coords <- cbind(Easting = as.numeric(as.character(police08$easting)),
                Northing = as.numeric(as.character(police08$northing)))
police08sp <- SpatialPointsDataFrame(coords, data = data.frame(police08$name,
                                                               police08$poi_ID), proj4string = CRS("+init=epsg:27700"))


police08 <- over(uktowns, geometry(police08sp), returnList = TRUE)
police08 <- sapply(police08, length)
police08 <- spCbind(uktowns, police08)
police08 <- merge(police08, popbua, by="BUA11CD")
police08 <- police08[which(police08$poplsoa<175000&police08$poplsoa>10000), ]
police08cent <- SpatialPointsDataFrame(gCentroid(police08, byid=TRUE), police08@data, match.ID=FALSE)
police08 <- police08[-c(2,3,4,5,7)]

write.csv(police08, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police08.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2013 <- read_csv("D:/poi/POI_2013.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Easting`), ]
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Northing`), ]
names(POI_2013)[names(POI_2013)=="pointx_classification_code"] <- "code"



police13 <- POI_2013[which(POI_2013$code==6330422), ]
police13 <- subset(police13, "`feature_easting`" != "" | "'feature_northing'" != "")
police13$poi_ID <- 1:nrow(police13)
coords <- cbind(Easting = as.numeric(as.character(police13$'feature_easting')),
                Northing = as.numeric(as.character(police13$'feature_northing')))
police13sp <- SpatialPointsDataFrame(coords, data = data.frame(police13$name,
                                                               police13$poi_ID), proj4string = CRS("+init=epsg:27700"))


police13 <- over(uktowns, geometry(police13sp), returnList = TRUE)
police13 <- sapply(police13, length)
police13 <- spCbind(uktowns, police13)
police13 <- merge(police13, popbua, by="BUA11CD")
police13 <- police13[which(police13$poplsoa<175000&police13$poplsoa>10000), ]
police13cent <- SpatialPointsDataFrame(gCentroid(police13, byid=TRUE), police13@data, match.ID=FALSE)
police13 <- police13[-c(2,3,4,5,7)]

write.csv(police13, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police13.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2012 <- read_csv("D:/poi/POI_2012.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2012 <- POI_2012 [complete.cases(POI_2012$`feature_easting`), ]
POI_2012 <- POI_2012 [complete.cases(POI_2012$`Feature Northing`), ]
names(POI_2012)[names(POI_2012)=="pointx_classification_code"] <- "code"



police12 <- POI_2012[which(POI_2012$code==6330422), ]
police12 <- subset(police12, "`feature_easting`" != "" | "'feature_northing'" != "")
police12$poi_ID <- 1:nrow(police12)
coords <- cbind(Easting = as.numeric(as.character(police12$`feature_easting`)),
                Northing = as.numeric(as.character(police12$'feature_northing')))
police12sp <- SpatialPointsDataFrame(coords, data = data.frame(police12$name,
                                                               police12$poi_ID), proj4string = CRS("+init=epsg:27700"))


police12 <- over(uktowns, geometry(police12sp), returnList = TRUE)
police12 <- sapply(police12, length)
police12 <- spCbind(uktowns, police12)
police12 <- merge(police12, popbua, by="BUA11CD")
police12 <- police12[which(police12$poplsoa<175000&police12$poplsoa>10000), ]
police12cent <- SpatialPointsDataFrame(gCentroid(police12, byid=TRUE), police12@data, match.ID=FALSE)
police12 <- police12[-c(2,3,4,5,7)]


write.csv(police12, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police12.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


gov_2014 <- read_delim("gov_2014.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Easting`), ]
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Northing`), ]
names(gov_2014)[names(gov_2014)=="PointX Classification Code"] <- "code"

gov_2015 <- read_delim("gov_2015.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Easting`), ]
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Northing`), ]
names(gov_2015)[names(gov_2015)=="PointX Classification Code"] <- "code"

gov_2016 <- read_delim("gov_2016.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Easting`), ]
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Northing`), ]
names(gov_2016)[names(gov_2016)=="PointX Classification Code"] <- "code"

gov_2017 <- read_delim("gov_2017.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Easting`), ]
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Northing`), ]
names(gov_2017)[names(gov_2017)=="PointX Classification Code"] <- "code"

gov_2018 <- read_delim("gov_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"


police15 <- gov_2015[which(gov_2015$code==06330422), ]
police15 <- subset(police15, "Feature Easting" != "" | "Feature Northing" != "")
police15$poi_ID <- 1:nrow(police15)
coords <- cbind(Easting = as.numeric(as.character(police15$`Feature Easting`)),
                Northing = as.numeric(as.character(police15$`Feature Northing`)))
police15sp <- SpatialPointsDataFrame(coords, data = data.frame(police15$Name,
                                                             police15$poi_ID), proj4string = CRS("+init=epsg:27700"))


police15 <- over(uktowns, geometry(police15sp), returnList = TRUE)
police15 <- sapply(police15, length)
police15 <- spCbind(uktowns, police15)
police15 <- merge(police15, popbua, by="BUA11CD")
police15 <- police15[which(police15$poplsoa<175000&police15$poplsoa>10000), ]
police15cent <- SpatialPointsDataFrame(gCentroid(police15, byid=TRUE), police15@data, match.ID=FALSE)
police15 <- police15[-c(2,3,4,5,7)]


police <-  as.data.frame(police15)

write.csv(police, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police15.csv")


police14 <- gov_2014[which(gov_2014$code==06330422), ]
police14 <- subset(police14, "Feature Easting" != "" | "Feature Northing" != "")
police14$poi_ID <- 1:nrow(police14)
coords <- cbind(Easting = as.numeric(as.character(police14$`Feature Easting`)),
                Northing = as.numeric(as.character(police14$`Feature Northing`)))
police14sp <- SpatialPointsDataFrame(coords, data = data.frame(police14$Name,
                                                               police14$poi_ID), proj4string = CRS("+init=epsg:27700"))


police14 <- over(uktowns, geometry(police14sp), returnList = TRUE)
police14 <- sapply(police14, length)
police14 <- spCbind(uktowns, police14)
police14 <- merge(police14, popbua, by="BUA11CD")
police14 <- police14[which(police14$poplsoa<175000&police14$poplsoa>10000), ]
police14cent <- SpatialPointsDataFrame(gCentroid(police14, byid=TRUE), police14@data, match.ID=FALSE)
police14 <- police14[-c(2,3,4,5,7)]


police <-  as.data.frame(police14)

write.csv(police, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police14.csv")


police16 <- gov_2016[which(gov_2016$code==06330422), ]
police16 <- subset(police16, "Feature Easting" != "" | "Feature Northing" != "")
police16$poi_ID <- 1:nrow(police16)
coords <- cbind(Easting = as.numeric(as.character(police16$`Feature Easting`)),
                Northing = as.numeric(as.character(police16$`Feature Northing`)))
police16sp <- SpatialPointsDataFrame(coords, data = data.frame(police16$Name,
                                                               police16$poi_ID), proj4string = CRS("+init=epsg:27700"))


police16 <- over(uktowns, geometry(police16sp), returnList = TRUE)
police16 <- sapply(police16, length)
police16 <- spCbind(uktowns, police16)
police16 <- merge(police16, popbua, by="BUA11CD")
police16 <- police16[which(police16$poplsoa<175000&police16$poplsoa>10000), ]
police16cent <- SpatialPointsDataFrame(gCentroid(police16, byid=TRUE), police16@data, match.ID=FALSE)
police16 <- police16[-c(2,3,4,5,7)]


police <-  as.data.frame(police16)

write.csv(police, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police16.csv")


police17 <- gov_2017[which(gov_2017$code==06330422), ]
police17 <- subset(police17, "Feature Easting" != "" | "Feature Northing" != "")
police17$poi_ID <- 1:nrow(police17)
coords <- cbind(Easting = as.numeric(as.character(police17$`Feature Easting`)),
                Northing = as.numeric(as.character(police17$`Feature Northing`)))
police17sp <- SpatialPointsDataFrame(coords, data = data.frame(police17$Name,
                                                               police17$poi_ID), proj4string = CRS("+init=epsg:27700"))


police17 <- over(uktowns, geometry(police17sp), returnList = TRUE)
police17 <- sapply(police17, length)
police17 <- spCbind(uktowns, police17)
police17 <- merge(police17, popbua, by="BUA11CD")
police17 <- police17[which(police17$poplsoa<175000&police17$poplsoa>10000), ]
police17cent <- SpatialPointsDataFrame(gCentroid(police17, byid=TRUE), police17@data, match.ID=FALSE)
police17 <- police17[-c(2,3,4,5,7)]


police <-  as.data.frame(police17)

write.csv(police, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police17.csv")


police18 <- gov_2018[which(gov_2018$code==06330422), ]
police18 <- subset(police18, "Feature Easting" != "" | "Feature Northing" != "")
police18$poi_ID <- 1:nrow(police18)
coords <- cbind(Easting = as.numeric(as.character(police18$`Feature Easting`)),
                Northing = as.numeric(as.character(police18$`Feature Northing`)))
police18sp <- SpatialPointsDataFrame(coords, data = data.frame(police18$Name,
                                                               police18$poi_ID), proj4string = CRS("+init=epsg:27700"))


police18 <- over(uktowns, geometry(police18sp), returnList = TRUE)
police18 <- sapply(police18, length)
police18 <- spCbind(uktowns, police18)
police18 <- merge(police18, popbua, by="BUA11CD")
police18 <- police18[which(police18$poplsoa<175000&police18$poplsoa>10000), ]
police18cent <- SpatialPointsDataFrame(gCentroid(police18, byid=TRUE), police18@data, match.ID=FALSE)
police18 <- police18[-c(2,3,4,5,7)]


police <-  as.data.frame(police18)

write.csv(police, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\police18.csv")







POI_2010 <- read_csv("D:/poi/POI_2010.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2010 <- POI_2010 [complete.cases(POI_2010$easting), ]
POI_2010 <- POI_2010 [complete.cases(POI_2010$northing), ]
names(POI_2010)[names(POI_2010)=="pointx_classification_code"] <- "code"



fire10 <- POI_2010[which(POI_2010$code==6330414), ]
fire10 <- subset(fire10, "easting" != "" | "northing" != "")
fire10$poi_ID <- 1:nrow(fire10)
coords <- cbind(Easting = as.numeric(as.character(fire10$easting)),
                Northing = as.numeric(as.character(fire10$northing)))
fire10sp <- SpatialPointsDataFrame(coords, data = data.frame(fire10$name,
                                                               fire10$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire10 <- over(uktowns, geometry(fire10sp), returnList = TRUE)
fire10 <- sapply(fire10, length)
fire10 <- spCbind(uktowns, fire10)
fire10 <- merge(fire10, popbua, by="BUA11CD")
fire10 <- fire10[which(fire10$poplsoa<175000&fire10$poplsoa>10000), ]
fire10cent <- SpatialPointsDataFrame(gCentroid(fire10, byid=TRUE), fire10@data, match.ID=FALSE)
fire10 <- fire10[-c(2,3,4,5,7)]

write.csv(fire10, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire10.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


POI_2011 <- read_csv("D:/poi/POI_2011.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2011 <- as.data.frame(POI_2011)
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_easting'), ]
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_northing'), ]
names(POI_2011)[names(POI_2011)=="pointx_classification_code"] <- "code"


fire11 <- POI_2011[which(POI_2011$code==6330414), ]
fire11 <- subset(fire11, "`feature_easting`" != "" | "'feature_northing'" != "")
fire11$poi_ID <- 1:nrow(fire11)
coords <- cbind(Easting = as.numeric(as.character(fire11$'feature_easting')),
                Northing = as.numeric(as.character(fire11$'feature_northing')))
fire11sp <- SpatialPointsDataFrame(coords, data = data.frame(fire11$name,
                                                               fire11$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire11 <- over(uktowns, geometry(fire11sp), returnList = TRUE)
fire11 <- sapply(fire11, length)
fire11 <- spCbind(uktowns, fire11)
fire11 <- merge(fire11, popbua, by="BUA11CD")
fire11 <- fire11[which(fire11$poplsoa<175000&fire11$poplsoa>10000), ]
fire11cent <- SpatialPointsDataFrame(gCentroid(fire11, byid=TRUE), fire11@data, match.ID=FALSE)
fire11 <- fire11[-c(2,3,4,5,7)]

write.csv(fire11, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire11.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2009 <- read_csv("D:/poi/POI_2009.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2009 <- POI_2009 [complete.cases(POI_2009$easting), ]
POI_2009 <- POI_2009 [complete.cases(POI_2009$northing), ]
names(POI_2009)[names(POI_2009)=="pointx_classification_code"] <- "code"



fire09 <- POI_2009[which(POI_2009$code==6330414), ]
fire09 <- subset(fire09, "easting" != "" | "northing" != "")
fire09$poi_ID <- 1:nrow(fire09)
coords <- cbind(Easting = as.numeric(as.character(fire09$easting)),
                Northing = as.numeric(as.character(fire09$northing)))
fire09sp <- SpatialPointsDataFrame(coords, data = data.frame(fire09$name,
                                                               fire09$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire09 <- over(uktowns, geometry(fire09sp), returnList = TRUE)
fire09 <- sapply(fire09, length)
fire09 <- spCbind(uktowns, fire09)
fire09 <- merge(fire09, popbua, by="BUA11CD")
fire09 <- fire09[which(fire09$poplsoa<175000&fire09$poplsoa>10000), ]
fire09cent <- SpatialPointsDataFrame(gCentroid(fire09, byid=TRUE), fire09@data, match.ID=FALSE)
fire09 <- fire09[-c(2,3,4,5,7)]

write.csv(fire09, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire09.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2008 <- read_csv("D:/poi/POI_2008.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2008 <- POI_2008 [complete.cases(POI_2008$easting), ]
POI_2008 <- POI_2008 [complete.cases(POI_2008$northing), ]
names(POI_2008)[names(POI_2008)=="pointx_classification_code"] <- "code"



fire08 <- POI_2008[which(POI_2008$code==6330414), ]
fire08 <- subset(fire08, "easting" != "" | "northing" != "")
fire08$poi_ID <- 1:nrow(fire08)
coords <- cbind(Easting = as.numeric(as.character(fire08$easting)),
                Northing = as.numeric(as.character(fire08$northing)))
fire08sp <- SpatialPointsDataFrame(coords, data = data.frame(fire08$name,
                                                               fire08$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire08 <- over(uktowns, geometry(fire08sp), returnList = TRUE)
fire08 <- sapply(fire08, length)
fire08 <- spCbind(uktowns, fire08)
fire08 <- merge(fire08, popbua, by="BUA11CD")
fire08 <- fire08[which(fire08$poplsoa<175000&fire08$poplsoa>10000), ]
fire08cent <- SpatialPointsDataFrame(gCentroid(fire08, byid=TRUE), fire08@data, match.ID=FALSE)
fire08 <- fire08[-c(2,3,4,5,7)]

write.csv(fire08, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire08.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2013 <- read_csv("D:/poi/POI_2013.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Easting`), ]
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Northing`), ]
names(POI_2013)[names(POI_2013)=="pointx_classification_code"] <- "code"



fire13 <- POI_2013[which(POI_2013$code==6330414), ]
fire13 <- subset(fire13, "`feature_easting`" != "" | "'feature_northing'" != "")
fire13$poi_ID <- 1:nrow(fire13)
coords <- cbind(Easting = as.numeric(as.character(fire13$'feature_easting')),
                Northing = as.numeric(as.character(fire13$'feature_northing')))
fire13sp <- SpatialPointsDataFrame(coords, data = data.frame(fire13$name,
                                                               fire13$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire13 <- over(uktowns, geometry(fire13sp), returnList = TRUE)
fire13 <- sapply(fire13, length)
fire13 <- spCbind(uktowns, fire13)
fire13 <- merge(fire13, popbua, by="BUA11CD")
fire13 <- fire13[which(fire13$poplsoa<175000&fire13$poplsoa>10000), ]
fire13cent <- SpatialPointsDataFrame(gCentroid(fire13, byid=TRUE), fire13@data, match.ID=FALSE)
fire13 <- fire13[-c(2,3,4,5,7)]

write.csv(fire13, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire13.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2012 <- read_csv("D:/poi/POI_2012.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2012 <- POI_2012 [complete.cases(POI_2012$`feature_easting`), ]
POI_2012 <- POI_2012 [complete.cases(POI_2012$`Feature Northing`), ]
names(POI_2012)[names(POI_2012)=="pointx_classification_code"] <- "code"



fire12 <- POI_2012[which(POI_2012$code==6330414), ]
fire12 <- subset(fire12, "`feature_easting`" != "" | "'feature_northing'" != "")
fire12$poi_ID <- 1:nrow(fire12)
coords <- cbind(Easting = as.numeric(as.character(fire12$`feature_easting`)),
                Northing = as.numeric(as.character(fire12$'feature_northing')))
fire12sp <- SpatialPointsDataFrame(coords, data = data.frame(fire12$name,
                                                               fire12$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire12 <- over(uktowns, geometry(fire12sp), returnList = TRUE)
fire12 <- sapply(fire12, length)
fire12 <- spCbind(uktowns, fire12)
fire12 <- merge(fire12, popbua, by="BUA11CD")
fire12 <- fire12[which(fire12$poplsoa<175000&fire12$poplsoa>10000), ]
fire12cent <- SpatialPointsDataFrame(gCentroid(fire12, byid=TRUE), fire12@data, match.ID=FALSE)
fire12 <- fire12[-c(2,3,4,5,7)]


write.csv(fire12, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire12.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


gov_2014 <- read_delim("gov_2014.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Easting`), ]
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Northing`), ]
names(gov_2014)[names(gov_2014)=="PointX Classification Code"] <- "code"

gov_2015 <- read_delim("gov_2015.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Easting`), ]
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Northing`), ]
names(gov_2015)[names(gov_2015)=="PointX Classification Code"] <- "code"

gov_2016 <- read_delim("gov_2016.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Easting`), ]
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Northing`), ]
names(gov_2016)[names(gov_2016)=="PointX Classification Code"] <- "code"

gov_2017 <- read_delim("gov_2017.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Easting`), ]
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Northing`), ]
names(gov_2017)[names(gov_2017)=="PointX Classification Code"] <- "code"

gov_2018 <- read_delim("gov_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Easting`), ]
gov_2018 <- gov_2018 [complete.cases(gov_2018$`Feature Northing`), ]
names(gov_2018)[names(gov_2018)=="PointX Classification Code"] <- "code"


fire15 <- gov_2015[which(gov_2015$code==06330414), ]
fire15 <- subset(fire15, "Feature Easting" != "" | "Feature Northing" != "")
fire15$poi_ID <- 1:nrow(fire15)
coords <- cbind(Easting = as.numeric(as.character(fire15$`Feature Easting`)),
                Northing = as.numeric(as.character(fire15$`Feature Northing`)))
fire15sp <- SpatialPointsDataFrame(coords, data = data.frame(fire15$Name,
                                                               fire15$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire15 <- over(uktowns, geometry(fire15sp), returnList = TRUE)
fire15 <- sapply(fire15, length)
fire15 <- spCbind(uktowns, fire15)
fire15 <- merge(fire15, popbua, by="BUA11CD")
fire15 <- fire15[which(fire15$poplsoa<175000&fire15$poplsoa>10000), ]
fire15cent <- SpatialPointsDataFrame(gCentroid(fire15, byid=TRUE), fire15@data, match.ID=FALSE)
fire15 <- fire15[-c(2,3,4,5,7)]


fire <-  as.data.frame(fire15)

write.csv(fire, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire15.csv")


fire14 <- gov_2014[which(gov_2014$code==06330414), ]
fire14 <- subset(fire14, "Feature Easting" != "" | "Feature Northing" != "")
fire14$poi_ID <- 1:nrow(fire14)
coords <- cbind(Easting = as.numeric(as.character(fire14$`Feature Easting`)),
                Northing = as.numeric(as.character(fire14$`Feature Northing`)))
fire14sp <- SpatialPointsDataFrame(coords, data = data.frame(fire14$Name,
                                                               fire14$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire14 <- over(uktowns, geometry(fire14sp), returnList = TRUE)
fire14 <- sapply(fire14, length)
fire14 <- spCbind(uktowns, fire14)
fire14 <- merge(fire14, popbua, by="BUA11CD")
fire14 <- fire14[which(fire14$poplsoa<175000&fire14$poplsoa>10000), ]
fire14cent <- SpatialPointsDataFrame(gCentroid(fire14, byid=TRUE), fire14@data, match.ID=FALSE)
fire14 <- fire14[-c(2,3,4,5,7)]


fire <-  as.data.frame(fire14)

write.csv(fire, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire14.csv")


fire16 <- gov_2016[which(gov_2016$code==06330414), ]
fire16 <- subset(fire16, "Feature Easting" != "" | "Feature Northing" != "")
fire16$poi_ID <- 1:nrow(fire16)
coords <- cbind(Easting = as.numeric(as.character(fire16$`Feature Easting`)),
                Northing = as.numeric(as.character(fire16$`Feature Northing`)))
fire16sp <- SpatialPointsDataFrame(coords, data = data.frame(fire16$Name,
                                                               fire16$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire16 <- over(uktowns, geometry(fire16sp), returnList = TRUE)
fire16 <- sapply(fire16, length)
fire16 <- spCbind(uktowns, fire16)
fire16 <- merge(fire16, popbua, by="BUA11CD")
fire16 <- fire16[which(fire16$poplsoa<175000&fire16$poplsoa>10000), ]
fire16cent <- SpatialPointsDataFrame(gCentroid(fire16, byid=TRUE), fire16@data, match.ID=FALSE)
fire16 <- fire16[-c(2,3,4,5,7)]


fire <-  as.data.frame(fire16)

write.csv(fire, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire16.csv")


fire17 <- gov_2017[which(gov_2017$code==06330414), ]
fire17 <- subset(fire17, "Feature Easting" != "" | "Feature Northing" != "")
fire17$poi_ID <- 1:nrow(fire17)
coords <- cbind(Easting = as.numeric(as.character(fire17$`Feature Easting`)),
                Northing = as.numeric(as.character(fire17$`Feature Northing`)))
fire17sp <- SpatialPointsDataFrame(coords, data = data.frame(fire17$Name,
                                                               fire17$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire17 <- over(uktowns, geometry(fire17sp), returnList = TRUE)
fire17 <- sapply(fire17, length)
fire17 <- spCbind(uktowns, fire17)
fire17 <- merge(fire17, popbua, by="BUA11CD")
fire17 <- fire17[which(fire17$poplsoa<175000&fire17$poplsoa>10000), ]
fire17cent <- SpatialPointsDataFrame(gCentroid(fire17, byid=TRUE), fire17@data, match.ID=FALSE)
fire17 <- fire17[-c(2,3,4,5,7)]


fire <-  as.data.frame(fire17)

write.csv(fire, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire17.csv")


fire18 <- gov_2018[which(gov_2018$code==06330414), ]
fire18 <- subset(fire18, "Feature Easting" != "" | "Feature Northing" != "")
fire18$poi_ID <- 1:nrow(fire18)
coords <- cbind(Easting = as.numeric(as.character(fire18$`Feature Easting`)),
                Northing = as.numeric(as.character(fire18$`Feature Northing`)))
fire18sp <- SpatialPointsDataFrame(coords, data = data.frame(fire18$Name,
                                                               fire18$poi_ID), proj4string = CRS("+init=epsg:27700"))


fire18 <- over(uktowns, geometry(fire18sp), returnList = TRUE)
fire18 <- sapply(fire18, length)
fire18 <- spCbind(uktowns, fire18)
fire18 <- merge(fire18, popbua, by="BUA11CD")
fire18 <- fire18[which(fire18$poplsoa<175000&fire18$poplsoa>10000), ]
fire18cent <- SpatialPointsDataFrame(gCentroid(fire18, byid=TRUE), fire18@data, match.ID=FALSE)
fire18 <- fire18[-c(2,3,4,5,7)]


fire <-  as.data.frame(fire18)

write.csv(fire, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\fire18.csv")














POI_2010 <- read_csv("D:/poi/POI_2010.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2010 <- POI_2010 [complete.cases(POI_2010$easting), ]
POI_2010 <- POI_2010 [complete.cases(POI_2010$northing), ]
names(POI_2010)[names(POI_2010)=="pointx_classification_code"] <- "code"



library10 <- POI_2010[which(POI_2010$code==6340458), ]
library10 <- subset(library10, "easting" != "" | "northing" != "")
library10$poi_ID <- 1:nrow(library10)
coords <- cbind(Easting = as.numeric(as.character(library10$easting)),
                Northing = as.numeric(as.character(library10$northing)))
library10sp <- SpatialPointsDataFrame(coords, data = data.frame(library10$name,
                                                             library10$poi_ID), proj4string = CRS("+init=epsg:27700"))


library10 <- over(uktowns, geometry(library10sp), returnList = TRUE)
library10 <- sapply(library10, length)
library10 <- spCbind(uktowns, library10)
library10 <- merge(library10, popbua, by="BUA11CD")
library10 <- library10[which(library10$poplsoa<175000&library10$poplsoa>10000), ]
library10cent <- SpatialPointsDataFrame(gCentroid(library10, byid=TRUE), library10@data, match.ID=FALSE)
library10 <- library10[-c(2,3,4,5,7)]

write.csv(library10, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library10.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


POI_2011 <- read_csv("D:/poi/POI_2011.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2011 <- as.data.frame(POI_2011)
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_easting'), ]
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_northing'), ]
names(POI_2011)[names(POI_2011)=="pointx_classification_code"] <- "code"


library11 <- POI_2011[which(POI_2011$code==6340458), ]
library11 <- subset(library11, "`feature_easting`" != "" | "'feature_northing'" != "")
library11$poi_ID <- 1:nrow(library11)
coords <- cbind(Easting = as.numeric(as.character(library11$'feature_easting')),
                Northing = as.numeric(as.character(library11$'feature_northing')))
library11sp <- SpatialPointsDataFrame(coords, data = data.frame(library11$name,
                                                             library11$poi_ID), proj4string = CRS("+init=epsg:27700"))


library11 <- over(uktowns, geometry(library11sp), returnList = TRUE)
library11 <- sapply(library11, length)
library11 <- spCbind(uktowns, library11)
library11 <- merge(library11, popbua, by="BUA11CD")
library11 <- library11[which(library11$poplsoa<175000&library11$poplsoa>10000), ]
library11cent <- SpatialPointsDataFrame(gCentroid(library11, byid=TRUE), library11@data, match.ID=FALSE)
library11 <- library11[-c(2,3,4,5,7)]

write.csv(library11, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library11.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2009 <- read_csv("D:/poi/POI_2009.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2009 <- POI_2009 [complete.cases(POI_2009$easting), ]
POI_2009 <- POI_2009 [complete.cases(POI_2009$northing), ]
names(POI_2009)[names(POI_2009)=="pointx_classification_code"] <- "code"



library09 <- POI_2009[which(POI_2009$code==6340458), ]
library09 <- subset(library09, "easting" != "" | "northing" != "")
library09$poi_ID <- 1:nrow(library09)
coords <- cbind(Easting = as.numeric(as.character(library09$easting)),
                Northing = as.numeric(as.character(library09$northing)))
library09sp <- SpatialPointsDataFrame(coords, data = data.frame(library09$name,
                                                             library09$poi_ID), proj4string = CRS("+init=epsg:27700"))


library09 <- over(uktowns, geometry(library09sp), returnList = TRUE)
library09 <- sapply(library09, length)
library09 <- spCbind(uktowns, library09)
library09 <- merge(library09, popbua, by="BUA11CD")
library09 <- library09[which(library09$poplsoa<175000&library09$poplsoa>10000), ]
library09cent <- SpatialPointsDataFrame(gCentroid(library09, byid=TRUE), library09@data, match.ID=FALSE)
library09 <- library09[-c(2,3,4,5,7)]

write.csv(library09, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library09.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2008 <- read_csv("D:/poi/POI_2008.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2008 <- POI_2008 [complete.cases(POI_2008$easting), ]
POI_2008 <- POI_2008 [complete.cases(POI_2008$northing), ]
names(POI_2008)[names(POI_2008)=="pointx_classification_code"] <- "code"



library08 <- POI_2008[which(POI_2008$code==6340458), ]
library08 <- subset(library08, "easting" != "" | "northing" != "")
library08$poi_ID <- 1:nrow(library08)
coords <- cbind(Easting = as.numeric(as.character(library08$easting)),
                Northing = as.numeric(as.character(library08$northing)))
library08sp <- SpatialPointsDataFrame(coords, data = data.frame(library08$name,
                                                             library08$poi_ID), proj4string = CRS("+init=epsg:27700"))


library08 <- over(uktowns, geometry(library08sp), returnList = TRUE)
library08 <- sapply(library08, length)
library08 <- spCbind(uktowns, library08)
library08 <- merge(library08, popbua, by="BUA11CD")
library08 <- library08[which(library08$poplsoa<175000&library08$poplsoa>10000), ]
library08cent <- SpatialPointsDataFrame(gCentroid(library08, byid=TRUE), library08@data, match.ID=FALSE)
library08 <- library08[-c(2,3,4,5,7)]

write.csv(library08, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library08.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2013 <- read_csv("D:/poi/POI_2013.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Easting`), ]
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Northing`), ]
names(POI_2013)[names(POI_2013)=="pointx_classification_code"] <- "code"



library13 <- POI_2013[which(POI_2013$code==6340458), ]
library13 <- subset(library13, "`feature_easting`" != "" | "'feature_northing'" != "")
library13$poi_ID <- 1:nrow(library13)
coords <- cbind(Easting = as.numeric(as.character(library13$'feature_easting')),
                Northing = as.numeric(as.character(library13$'feature_northing')))
library13sp <- SpatialPointsDataFrame(coords, data = data.frame(library13$name,
                                                             library13$poi_ID), proj4string = CRS("+init=epsg:27700"))


library13 <- over(uktowns, geometry(library13sp), returnList = TRUE)
library13 <- sapply(library13, length)
library13 <- spCbind(uktowns, library13)
library13 <- merge(library13, popbua, by="BUA11CD")
library13 <- library13[which(library13$poplsoa<175000&library13$poplsoa>10000), ]
library13cent <- SpatialPointsDataFrame(gCentroid(library13, byid=TRUE), library13@data, match.ID=FALSE)
library13 <- library13[-c(2,3,4,5,7)]

write.csv(library13, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library13.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2012 <- read_csv("D:/poi/POI_2012.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2012 <- POI_2012 [complete.cases(POI_2012$`feature_easting`), ]
POI_2012 <- POI_2012 [complete.cases(POI_2012$`Feature Northing`), ]
names(POI_2012)[names(POI_2012)=="pointx_classification_code"] <- "code"



library12 <- POI_2012[which(POI_2012$code==6340458), ]
library12 <- subset(library12, "`feature_easting`" != "" | "'feature_northing'" != "")
library12$poi_ID <- 1:nrow(library12)
coords <- cbind(Easting = as.numeric(as.character(library12$`feature_easting`)),
                Northing = as.numeric(as.character(library12$'feature_northing')))
library12sp <- SpatialPointsDataFrame(coords, data = data.frame(library12$name,
                                                             library12$poi_ID), proj4string = CRS("+init=epsg:27700"))


library12 <- over(uktowns, geometry(library12sp), returnList = TRUE)
library12 <- sapply(library12, length)
library12 <- spCbind(uktowns, library12)
library12 <- merge(library12, popbua, by="BUA11CD")
library12 <- library12[which(library12$poplsoa<175000&library12$poplsoa>10000), ]
library12cent <- SpatialPointsDataFrame(gCentroid(library12, byid=TRUE), library12@data, match.ID=FALSE)
library12 <- library12[-c(2,3,4,5,7)]


write.csv(library12, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library12.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


infra2014 <- read_delim("infra14.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
infra2014 <- infra2014 [complete.cases(infra2014$`Feature Easting`), ]
infra2014 <- infra2014 [complete.cases(infra2014$`Feature Northing`), ]
names(infra2014)[names(infra2014)=="PointX Classification Code"] <- "code"

infra2015 <- read_delim("infra15.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
infra2015 <- infra2015 [complete.cases(infra2015$`Feature Easting`), ]
infra2015 <- infra2015 [complete.cases(infra2015$`Feature Northing`), ]
names(infra2015)[names(infra2015)=="PointX Classification Code"] <- "code"

infra2016 <- read_delim("infra16.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
infra2016 <- infra2016 [complete.cases(infra2016$`Feature Easting`), ]
infra2016 <- infra2016 [complete.cases(infra2016$`Feature Northing`), ]
names(infra2016)[names(infra2016)=="PointX Classification Code"] <- "code"

infra2017 <- read_delim("infra17.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
infra2017 <- infra2017 [complete.cases(infra2017$`Feature Easting`), ]
infra2017 <- infra2017 [complete.cases(infra2017$`Feature Northing`), ]
names(infra2017)[names(infra2017)=="PointX Classification Code"] <- "code"

infra2018 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                       col_types = cols('PointX Classification Code' = col_number()))
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Easting`), ]
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Northing`), ]
names(infra2018)[names(infra2018)=="PointX Classification Code"] <- "code"


library15 <- infra2015[which(infra2015$code==06340458), ]
library15 <- subset(library15, "Feature Easting" != "" | "Feature Northing" != "")
library15$poi_ID <- 1:nrow(library15)
coords <- cbind(Easting = as.numeric(as.character(library15$`Feature Easting`)),
                Northing = as.numeric(as.character(library15$`Feature Northing`)))
library15sp <- SpatialPointsDataFrame(coords, data = data.frame(library15$Name,
                                                             library15$poi_ID), proj4string = CRS("+init=epsg:27700"))


library15 <- over(uktowns, geometry(library15sp), returnList = TRUE)
library15 <- sapply(library15, length)
library15 <- spCbind(uktowns, library15)
library15 <- merge(library15, popbua, by="BUA11CD")
library15 <- library15[which(library15$poplsoa<175000&library15$poplsoa>10000), ]
library15cent <- SpatialPointsDataFrame(gCentroid(library15, byid=TRUE), library15@data, match.ID=FALSE)
library15 <- library15[-c(2,3,4,5,7)]


library <-  as.data.frame(library15)

write.csv(library, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library15.csv")


library14 <- infra2014[which(infra2014$code==06340458), ]
library14 <- subset(library14, "Feature Easting" != "" | "Feature Northing" != "")
library14$poi_ID <- 1:nrow(library14)
coords <- cbind(Easting = as.numeric(as.character(library14$`Feature Easting`)),
                Northing = as.numeric(as.character(library14$`Feature Northing`)))
library14sp <- SpatialPointsDataFrame(coords, data = data.frame(library14$Name,
                                                             library14$poi_ID), proj4string = CRS("+init=epsg:27700"))


library14 <- over(uktowns, geometry(library14sp), returnList = TRUE)
library14 <- sapply(library14, length)
library14 <- spCbind(uktowns, library14)
library14 <- merge(library14, popbua, by="BUA11CD")
library14 <- library14[which(library14$poplsoa<175000&library14$poplsoa>10000), ]
library14cent <- SpatialPointsDataFrame(gCentroid(library14, byid=TRUE), library14@data, match.ID=FALSE)
library14 <- library14[-c(2,3,4,5,7)]


library <-  as.data.frame(library14)

write.csv(library, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library14.csv")


library16 <- infra2016[which(infra2016$code==06340458), ]
library16 <- subset(library16, "Feature Easting" != "" | "Feature Northing" != "")
library16$poi_ID <- 1:nrow(library16)
coords <- cbind(Easting = as.numeric(as.character(library16$`Feature Easting`)),
                Northing = as.numeric(as.character(library16$`Feature Northing`)))
library16sp <- SpatialPointsDataFrame(coords, data = data.frame(library16$Name,
                                                             library16$poi_ID), proj4string = CRS("+init=epsg:27700"))


library16 <- over(uktowns, geometry(library16sp), returnList = TRUE)
library16 <- sapply(library16, length)
library16 <- spCbind(uktowns, library16)
library16 <- merge(library16, popbua, by="BUA11CD")
library16 <- library16[which(library16$poplsoa<175000&library16$poplsoa>10000), ]
library16cent <- SpatialPointsDataFrame(gCentroid(library16, byid=TRUE), library16@data, match.ID=FALSE)
library16 <- library16[-c(2,3,4,5,7)]


library <-  as.data.frame(library16)

write.csv(library, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library16.csv")


library17 <- infra2017[which(infra2017$code==06340458), ]
library17 <- subset(library17, "Feature Easting" != "" | "Feature Northing" != "")
library17$poi_ID <- 1:nrow(library17)
coords <- cbind(Easting = as.numeric(as.character(library17$`Feature Easting`)),
                Northing = as.numeric(as.character(library17$`Feature Northing`)))
library17sp <- SpatialPointsDataFrame(coords, data = data.frame(library17$Name,
                                                             library17$poi_ID), proj4string = CRS("+init=epsg:27700"))


library17 <- over(uktowns, geometry(library17sp), returnList = TRUE)
library17 <- sapply(library17, length)
library17 <- spCbind(uktowns, library17)
library17 <- merge(library17, popbua, by="BUA11CD")
library17 <- library17[which(library17$poplsoa<175000&library17$poplsoa>10000), ]
library17cent <- SpatialPointsDataFrame(gCentroid(library17, byid=TRUE), library17@data, match.ID=FALSE)
library17 <- library17[-c(2,3,4,5,7)]


library <-  as.data.frame(library17)

write.csv(library, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library17.csv")


library18 <- infra2018[which(infra2018$code==06340458), ]
library18 <- subset(library18, "Feature Easting" != "" | "Feature Northing" != "")
library18$poi_ID <- 1:nrow(library18)
coords <- cbind(Easting = as.numeric(as.character(library18$`Feature Easting`)),
                Northing = as.numeric(as.character(library18$`Feature Northing`)))
library18sp <- SpatialPointsDataFrame(coords, data = data.frame(library18$Name,
                                                             library18$poi_ID), proj4string = CRS("+init=epsg:27700"))


library18 <- over(uktowns, geometry(library18sp), returnList = TRUE)
library18 <- sapply(library18, length)
library18 <- spCbind(uktowns, library18)
library18 <- merge(library18, popbua, by="BUA11CD")
library18 <- library18[which(library18$poplsoa<175000&library18$poplsoa>10000), ]
library18cent <- SpatialPointsDataFrame(gCentroid(library18, byid=TRUE), library18@data, match.ID=FALSE)
library18 <- library18[-c(2,3,4,5,7)]


library <-  as.data.frame(library18)

write.csv(library, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\library18.csv")

















POI_2010 <- read_csv("D:/poi/POI_2010.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2010 <- POI_2010 [complete.cases(POI_2010$easting), ]
POI_2010 <- POI_2010 [complete.cases(POI_2010$northing), ]
names(POI_2010)[names(POI_2010)=="pointx_classification_code"] <- "code"



toilet10 <- POI_2010[which(POI_2010$code==6340461), ]
toilet10 <- subset(toilet10, "easting" != "" | "northing" != "")
toilet10$poi_ID <- 1:nrow(toilet10)
coords <- cbind(Easting = as.numeric(as.character(toilet10$easting)),
                Northing = as.numeric(as.character(toilet10$northing)))
toilet10sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet10$name,
                                                                toilet10$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet10 <- over(uktowns, geometry(toilet10sp), returnList = TRUE)
toilet10 <- sapply(toilet10, length)
toilet10 <- spCbind(uktowns, toilet10)
toilet10 <- merge(toilet10, popbua, by="BUA11CD")
toilet10 <- toilet10[which(toilet10$poplsoa<175000&toilet10$poplsoa>10000), ]
toilet10cent <- SpatialPointsDataFrame(gCentroid(toilet10, byid=TRUE), toilet10@data, match.ID=FALSE)
toilet10 <- toilet10[-c(2,3,4,5,7)]

write.csv(toilet10, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet10.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


POI_2011 <- read_csv("D:/poi/POI_2011.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2011 <- as.data.frame(POI_2011)
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_easting'), ]
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_northing'), ]
names(POI_2011)[names(POI_2011)=="pointx_classification_code"] <- "code"


toilet11 <- POI_2011[which(POI_2011$code==6340461), ]
toilet11 <- subset(toilet11, "`feature_easting`" != "" | "'feature_northing'" != "")
toilet11$poi_ID <- 1:nrow(toilet11)
coords <- cbind(Easting = as.numeric(as.character(toilet11$'feature_easting')),
                Northing = as.numeric(as.character(toilet11$'feature_northing')))
toilet11sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet11$name,
                                                                toilet11$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet11 <- over(uktowns, geometry(toilet11sp), returnList = TRUE)
toilet11 <- sapply(toilet11, length)
toilet11 <- spCbind(uktowns, toilet11)
toilet11 <- merge(toilet11, popbua, by="BUA11CD")
toilet11 <- toilet11[which(toilet11$poplsoa<175000&toilet11$poplsoa>10000), ]
toilet11cent <- SpatialPointsDataFrame(gCentroid(toilet11, byid=TRUE), toilet11@data, match.ID=FALSE)
toilet11 <- toilet11[-c(2,3,4,5,7)]

write.csv(toilet11, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet11.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2009 <- read_csv("D:/poi/POI_2009.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2009 <- POI_2009 [complete.cases(POI_2009$easting), ]
POI_2009 <- POI_2009 [complete.cases(POI_2009$northing), ]
names(POI_2009)[names(POI_2009)=="pointx_classification_code"] <- "code"



toilet09 <- POI_2009[which(POI_2009$code==6340461), ]
toilet09 <- subset(toilet09, "easting" != "" | "northing" != "")
toilet09$poi_ID <- 1:nrow(toilet09)
coords <- cbind(Easting = as.numeric(as.character(toilet09$easting)),
                Northing = as.numeric(as.character(toilet09$northing)))
toilet09sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet09$name,
                                                                toilet09$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet09 <- over(uktowns, geometry(toilet09sp), returnList = TRUE)
toilet09 <- sapply(toilet09, length)
toilet09 <- spCbind(uktowns, toilet09)
toilet09 <- merge(toilet09, popbua, by="BUA11CD")
toilet09 <- toilet09[which(toilet09$poplsoa<175000&toilet09$poplsoa>10000), ]
toilet09cent <- SpatialPointsDataFrame(gCentroid(toilet09, byid=TRUE), toilet09@data, match.ID=FALSE)
toilet09 <- toilet09[-c(2,3,4,5,7)]

write.csv(toilet09, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet09.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2008 <- read_csv("D:/poi/POI_2008.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2008 <- POI_2008 [complete.cases(POI_2008$easting), ]
POI_2008 <- POI_2008 [complete.cases(POI_2008$northing), ]
names(POI_2008)[names(POI_2008)=="pointx_classification_code"] <- "code"



toilet08 <- POI_2008[which(POI_2008$code==6340461), ]
toilet08 <- subset(toilet08, "easting" != "" | "northing" != "")
toilet08$poi_ID <- 1:nrow(toilet08)
coords <- cbind(Easting = as.numeric(as.character(toilet08$easting)),
                Northing = as.numeric(as.character(toilet08$northing)))
toilet08sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet08$name,
                                                                toilet08$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet08 <- over(uktowns, geometry(toilet08sp), returnList = TRUE)
toilet08 <- sapply(toilet08, length)
toilet08 <- spCbind(uktowns, toilet08)
toilet08 <- merge(toilet08, popbua, by="BUA11CD")
toilet08 <- toilet08[which(toilet08$poplsoa<175000&toilet08$poplsoa>10000), ]
toilet08cent <- SpatialPointsDataFrame(gCentroid(toilet08, byid=TRUE), toilet08@data, match.ID=FALSE)
toilet08 <- toilet08[-c(2,3,4,5,7)]

write.csv(toilet08, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet08.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2013 <- read_csv("D:/poi/POI_2013.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Easting`), ]
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Northing`), ]
names(POI_2013)[names(POI_2013)=="pointx_classification_code"] <- "code"



toilet13 <- POI_2013[which(POI_2013$code==6340461), ]
toilet13 <- subset(toilet13, "`feature_easting`" != "" | "'feature_northing'" != "")
toilet13$poi_ID <- 1:nrow(toilet13)
coords <- cbind(Easting = as.numeric(as.character(toilet13$'feature_easting')),
                Northing = as.numeric(as.character(toilet13$'feature_northing')))
toilet13sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet13$name,
                                                                toilet13$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet13 <- over(uktowns, geometry(toilet13sp), returnList = TRUE)
toilet13 <- sapply(toilet13, length)
toilet13 <- spCbind(uktowns, toilet13)
toilet13 <- merge(toilet13, popbua, by="BUA11CD")
toilet13 <- toilet13[which(toilet13$poplsoa<175000&toilet13$poplsoa>10000), ]
toilet13cent <- SpatialPointsDataFrame(gCentroid(toilet13, byid=TRUE), toilet13@data, match.ID=FALSE)
toilet13 <- toilet13[-c(2,3,4,5,7)]

write.csv(toilet13, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet13.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2012 <- read_csv("D:/poi/POI_2012.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2012 <- POI_2012 [complete.cases(POI_2012$`feature_easting`), ]
POI_2012 <- POI_2012 [complete.cases(POI_2012$`Feature Northing`), ]
names(POI_2012)[names(POI_2012)=="pointx_classification_code"] <- "code"



toilet12 <- POI_2012[which(POI_2012$code==6340461), ]
toilet12 <- subset(toilet12, "`feature_easting`" != "" | "'feature_northing'" != "")
toilet12$poi_ID <- 1:nrow(toilet12)
coords <- cbind(Easting = as.numeric(as.character(toilet12$`feature_easting`)),
                Northing = as.numeric(as.character(toilet12$'feature_northing')))
toilet12sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet12$name,
                                                                toilet12$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet12 <- over(uktowns, geometry(toilet12sp), returnList = TRUE)
toilet12 <- sapply(toilet12, length)
toilet12 <- spCbind(uktowns, toilet12)
toilet12 <- merge(toilet12, popbua, by="BUA11CD")
toilet12 <- toilet12[which(toilet12$poplsoa<175000&toilet12$poplsoa>10000), ]
toilet12cent <- SpatialPointsDataFrame(gCentroid(toilet12, byid=TRUE), toilet12@data, match.ID=FALSE)
toilet12 <- toilet12[-c(2,3,4,5,7)]


write.csv(toilet12, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet12.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


infra2014 <- read_delim("infra14.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
infra2014 <- infra2014 [complete.cases(infra2014$`Feature Easting`), ]
infra2014 <- infra2014 [complete.cases(infra2014$`Feature Northing`), ]
names(infra2014)[names(infra2014)=="PointX Classification Code"] <- "code"

infra2015 <- read_delim("infra15.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
infra2015 <- infra2015 [complete.cases(infra2015$`Feature Easting`), ]
infra2015 <- infra2015 [complete.cases(infra2015$`Feature Northing`), ]
names(infra2015)[names(infra2015)=="PointX Classification Code"] <- "code"

infra2016 <- read_delim("infra16.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
infra2016 <- infra2016 [complete.cases(infra2016$`Feature Easting`), ]
infra2016 <- infra2016 [complete.cases(infra2016$`Feature Northing`), ]
names(infra2016)[names(infra2016)=="PointX Classification Code"] <- "code"

infra2017 <- read_delim("infra17.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
infra2017 <- infra2017 [complete.cases(infra2017$`Feature Easting`), ]
infra2017 <- infra2017 [complete.cases(infra2017$`Feature Northing`), ]
names(infra2017)[names(infra2017)=="PointX Classification Code"] <- "code"

infra2018 <- read_delim("infra18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Easting`), ]
infra2018 <- infra2018 [complete.cases(infra2018$`Feature Northing`), ]
names(infra2018)[names(infra2018)=="PointX Classification Code"] <- "code"


toilet15 <- infra2015[which(infra2015$code==06340461), ]
toilet15 <- subset(toilet15, "Feature Easting" != "" | "Feature Northing" != "")
toilet15$poi_ID <- 1:nrow(toilet15)
coords <- cbind(Easting = as.numeric(as.character(toilet15$`Feature Easting`)),
                Northing = as.numeric(as.character(toilet15$`Feature Northing`)))
toilet15sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet15$Name,
                                                                toilet15$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet15 <- over(uktowns, geometry(toilet15sp), returnList = TRUE)
toilet15 <- sapply(toilet15, length)
toilet15 <- spCbind(uktowns, toilet15)
toilet15 <- merge(toilet15, popbua, by="BUA11CD")
toilet15 <- toilet15[which(toilet15$poplsoa<175000&toilet15$poplsoa>10000), ]
toilet15cent <- SpatialPointsDataFrame(gCentroid(toilet15, byid=TRUE), toilet15@data, match.ID=FALSE)
toilet15 <- toilet15[-c(2,3,4,5,7)]


toilet <-  as.data.frame(toilet15)

write.csv(toilet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet15.csv")


toilet14 <- infra2014[which(infra2014$code==06340461), ]
toilet14 <- subset(toilet14, "Feature Easting" != "" | "Feature Northing" != "")
toilet14$poi_ID <- 1:nrow(toilet14)
coords <- cbind(Easting = as.numeric(as.character(toilet14$`Feature Easting`)),
                Northing = as.numeric(as.character(toilet14$`Feature Northing`)))
toilet14sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet14$Name,
                                                                toilet14$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet14 <- over(uktowns, geometry(toilet14sp), returnList = TRUE)
toilet14 <- sapply(toilet14, length)
toilet14 <- spCbind(uktowns, toilet14)
toilet14 <- merge(toilet14, popbua, by="BUA11CD")
toilet14 <- toilet14[which(toilet14$poplsoa<175000&toilet14$poplsoa>10000), ]
toilet14cent <- SpatialPointsDataFrame(gCentroid(toilet14, byid=TRUE), toilet14@data, match.ID=FALSE)
toilet14 <- toilet14[-c(2,3,4,5,7)]


toilet <-  as.data.frame(toilet14)

write.csv(toilet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet14.csv")


toilet16 <- infra2016[which(infra2016$code==06340461), ]
toilet16 <- subset(toilet16, "Feature Easting" != "" | "Feature Northing" != "")
toilet16$poi_ID <- 1:nrow(toilet16)
coords <- cbind(Easting = as.numeric(as.character(toilet16$`Feature Easting`)),
                Northing = as.numeric(as.character(toilet16$`Feature Northing`)))
toilet16sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet16$Name,
                                                                toilet16$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet16 <- over(uktowns, geometry(toilet16sp), returnList = TRUE)
toilet16 <- sapply(toilet16, length)
toilet16 <- spCbind(uktowns, toilet16)
toilet16 <- merge(toilet16, popbua, by="BUA11CD")
toilet16 <- toilet16[which(toilet16$poplsoa<175000&toilet16$poplsoa>10000), ]
toilet16cent <- SpatialPointsDataFrame(gCentroid(toilet16, byid=TRUE), toilet16@data, match.ID=FALSE)
toilet16 <- toilet16[-c(2,3,4,5,7)]


toilet <-  as.data.frame(toilet16)

write.csv(toilet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet16.csv")


toilet17 <- infra2017[which(infra2017$code==06340461), ]
toilet17 <- subset(toilet17, "Feature Easting" != "" | "Feature Northing" != "")
toilet17$poi_ID <- 1:nrow(toilet17)
coords <- cbind(Easting = as.numeric(as.character(toilet17$`Feature Easting`)),
                Northing = as.numeric(as.character(toilet17$`Feature Northing`)))
toilet17sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet17$Name,
                                                                toilet17$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet17 <- over(uktowns, geometry(toilet17sp), returnList = TRUE)
toilet17 <- sapply(toilet17, length)
toilet17 <- spCbind(uktowns, toilet17)
toilet17 <- merge(toilet17, popbua, by="BUA11CD")
toilet17 <- toilet17[which(toilet17$poplsoa<175000&toilet17$poplsoa>10000), ]
toilet17cent <- SpatialPointsDataFrame(gCentroid(toilet17, byid=TRUE), toilet17@data, match.ID=FALSE)
toilet17 <- toilet17[-c(2,3,4,5,7)]


toilet <-  as.data.frame(toilet17)

write.csv(toilet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet17.csv")


toilet18 <- infra2018[which(infra2018$code==06340461), ]
toilet18 <- subset(toilet18, "Feature Easting" != "" | "Feature Northing" != "")
toilet18$poi_ID <- 1:nrow(toilet18)
coords <- cbind(Easting = as.numeric(as.character(toilet18$`Feature Easting`)),
                Northing = as.numeric(as.character(toilet18$`Feature Northing`)))
toilet18sp <- SpatialPointsDataFrame(coords, data = data.frame(toilet18$Name,
                                                                toilet18$poi_ID), proj4string = CRS("+init=epsg:27700"))


toilet18 <- over(uktowns, geometry(toilet18sp), returnList = TRUE)
toilet18 <- sapply(toilet18, length)
toilet18 <- spCbind(uktowns, toilet18)
toilet18 <- merge(toilet18, popbua, by="BUA11CD")
toilet18 <- toilet18[which(toilet18$poplsoa<175000&toilet18$poplsoa>10000), ]
toilet18cent <- SpatialPointsDataFrame(gCentroid(toilet18, byid=TRUE), toilet18@data, match.ID=FALSE)
toilet18 <- toilet18[-c(2,3,4,5,7)]


toilet <-  as.data.frame(toilet18)

write.csv(toilet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\toilet18.csv")






















POI_2010 <- read_csv("D:/poi/POI_2010.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2010 <- POI_2010 [complete.cases(POI_2010$easting), ]
POI_2010 <- POI_2010 [complete.cases(POI_2010$northing), ]
names(POI_2010)[names(POI_2010)=="pointx_classification_code"] <- "code"



gp10 <- POI_2010[which(POI_2010$code==5280369), ]
gp10 <- subset(gp10, "easting" != "" | "northing" != "")
gp10$poi_ID <- 1:nrow(gp10)
coords <- cbind(Easting = as.numeric(as.character(gp10$easting)),
                Northing = as.numeric(as.character(gp10$northing)))
gp10sp <- SpatialPointsDataFrame(coords, data = data.frame(gp10$name,
                                                               gp10$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp10 <- over(uktowns, geometry(gp10sp), returnList = TRUE)
gp10 <- sapply(gp10, length)
gp10 <- spCbind(uktowns, gp10)
gp10 <- merge(gp10, popbua, by="BUA11CD")
gp10 <- gp10[which(gp10$poplsoa<175000&gp10$poplsoa>10000), ]
gp10cent <- SpatialPointsDataFrame(gCentroid(gp10, byid=TRUE), gp10@data, match.ID=FALSE)
gp10 <- gp10[-c(2,3,4,5,7)]

write.csv(gp10, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps10.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


POI_2011 <- read_csv("D:/poi/POI_2011.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2011 <- as.data.frame(POI_2011)
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_easting'), ]
POI_2011 <- POI_2011[complete.cases(POI_2011$'feature_northing'), ]
names(POI_2011)[names(POI_2011)=="pointx_classification_code"] <- "code"


gp11 <- POI_2011[which(POI_2011$code==5280369), ]
gp11 <- subset(gp11, "`feature_easting`" != "" | "'feature_northing'" != "")
gp11$poi_ID <- 1:nrow(gp11)
coords <- cbind(Easting = as.numeric(as.character(gp11$'feature_easting')),
                Northing = as.numeric(as.character(gp11$'feature_northing')))
gp11sp <- SpatialPointsDataFrame(coords, data = data.frame(gp11$name,
                                                               gp11$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp11 <- over(uktowns, geometry(gp11sp), returnList = TRUE)
gp11 <- sapply(gp11, length)
gp11 <- spCbind(uktowns, gp11)
gp11 <- merge(gp11, popbua, by="BUA11CD")
gp11 <- gp11[which(gp11$poplsoa<175000&gp11$poplsoa>10000), ]
gp11cent <- SpatialPointsDataFrame(gCentroid(gp11, byid=TRUE), gp11@data, match.ID=FALSE)
gp11 <- gp11[-c(2,3,4,5,7)]

write.csv(gp11, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps11.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2009 <- read_csv("D:/poi/POI_2009.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2009 <- POI_2009 [complete.cases(POI_2009$easting), ]
POI_2009 <- POI_2009 [complete.cases(POI_2009$northing), ]
names(POI_2009)[names(POI_2009)=="pointx_classification_code"] <- "code"



gp09 <- POI_2009[which(POI_2009$code==5280369), ]
gp09 <- subset(gp09, "easting" != "" | "northing" != "")
gp09$poi_ID <- 1:nrow(gp09)
coords <- cbind(Easting = as.numeric(as.character(gp09$easting)),
                Northing = as.numeric(as.character(gp09$northing)))
gp09sp <- SpatialPointsDataFrame(coords, data = data.frame(gp09$name,
                                                               gp09$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp09 <- over(uktowns, geometry(gp09sp), returnList = TRUE)
gp09 <- sapply(gp09, length)
gp09 <- spCbind(uktowns, gp09)
gp09 <- merge(gp09, popbua, by="BUA11CD")
gp09 <- gp09[which(gp09$poplsoa<175000&gp09$poplsoa>10000), ]
gp09cent <- SpatialPointsDataFrame(gCentroid(gp09, byid=TRUE), gp09@data, match.ID=FALSE)
gp09 <- gp09[-c(2,3,4,5,7)]

write.csv(gp09, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps09.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2008 <- read_csv("D:/poi/POI_2008.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2008 <- POI_2008 [complete.cases(POI_2008$easting), ]
POI_2008 <- POI_2008 [complete.cases(POI_2008$northing), ]
names(POI_2008)[names(POI_2008)=="pointx_classification_code"] <- "code"



gp08 <- POI_2008[which(POI_2008$code==5280369), ]
gp08 <- subset(gp08, "easting" != "" | "northing" != "")
gp08$poi_ID <- 1:nrow(gp08)
coords <- cbind(Easting = as.numeric(as.character(gp08$easting)),
                Northing = as.numeric(as.character(gp08$northing)))
gp08sp <- SpatialPointsDataFrame(coords, data = data.frame(gp08$name,
                                                               gp08$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp08 <- over(uktowns, geometry(gp08sp), returnList = TRUE)
gp08 <- sapply(gp08, length)
gp08 <- spCbind(uktowns, gp08)
gp08 <- merge(gp08, popbua, by="BUA11CD")
gp08 <- gp08[which(gp08$poplsoa<175000&gp08$poplsoa>10000), ]
gp08cent <- SpatialPointsDataFrame(gCentroid(gp08, byid=TRUE), gp08@data, match.ID=FALSE)
gp08 <- gp08[-c(2,3,4,5,7)]

write.csv(gp08, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps08.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2013 <- read_csv("D:/poi/POI_2013.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Easting`), ]
POI_2013 <- POI_2013 [complete.cases(POI_2013$`Feature Northing`), ]
names(POI_2013)[names(POI_2013)=="pointx_classification_code"] <- "code"



gp13 <- POI_2013[which(POI_2013$code==5280369), ]
gp13 <- subset(gp13, "`feature_easting`" != "" | "'feature_northing'" != "")
gp13$poi_ID <- 1:nrow(gp13)
coords <- cbind(Easting = as.numeric(as.character(gp13$'feature_easting')),
                Northing = as.numeric(as.character(gp13$'feature_northing')))
gp13sp <- SpatialPointsDataFrame(coords, data = data.frame(gp13$name,
                                                               gp13$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp13 <- over(uktowns, geometry(gp13sp), returnList = TRUE)
gp13 <- sapply(gp13, length)
gp13 <- spCbind(uktowns, gp13)
gp13 <- merge(gp13, popbua, by="BUA11CD")
gp13 <- gp13[which(gp13$poplsoa<175000&gp13$poplsoa>10000), ]
gp13cent <- SpatialPointsDataFrame(gCentroid(gp13, byid=TRUE), gp13@data, match.ID=FALSE)
gp13 <- gp13[-c(2,3,4,5,7)]

write.csv(gp13, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps13.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])



POI_2012 <- read_csv("D:/poi/POI_2012.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2012 <- POI_2012 [complete.cases(POI_2012$`feature_easting`), ]
POI_2012 <- POI_2012 [complete.cases(POI_2012$`Feature Northing`), ]
names(POI_2012)[names(POI_2012)=="pointx_classification_code"] <- "code"



gp12 <- POI_2012[which(POI_2012$code==5280369), ]
gp12 <- subset(gp12, "`feature_easting`" != "" | "'feature_northing'" != "")
gp12$poi_ID <- 1:nrow(gp12)
coords <- cbind(Easting = as.numeric(as.character(gp12$`feature_easting`)),
                Northing = as.numeric(as.character(gp12$'feature_northing')))
gp12sp <- SpatialPointsDataFrame(coords, data = data.frame(gp12$name,
                                                               gp12$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp12 <- over(uktowns, geometry(gp12sp), returnList = TRUE)
gp12 <- sapply(gp12, length)
gp12 <- spCbind(uktowns, gp12)
gp12 <- merge(gp12, popbua, by="BUA11CD")
gp12 <- gp12[which(gp12$poplsoa<175000&gp12$poplsoa>10000), ]
gp12cent <- SpatialPointsDataFrame(gCentroid(gp12, byid=TRUE), gp12@data, match.ID=FALSE)
gp12 <- gp12[-c(2,3,4,5,7)]


write.csv(gp12, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps12.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" )])


health_2014 <- read_delim("health_2014.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2014 <- health_2014 [complete.cases(health_2014$`Feature Easting`), ]
health_2014 <- health_2014 [complete.cases(health_2014$`Feature Northing`), ]
names(health_2014)[names(health_2014)=="PointX Classification Code"] <- "code"

health_2015 <- read_delim("health_2015.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2015 <- health_2015 [complete.cases(health_2015$`Feature Easting`), ]
health_2015 <- health_2015 [complete.cases(health_2015$`Feature Northing`), ]
names(health_2015)[names(health_2015)=="PointX Classification Code"] <- "code"

health_2016 <- read_delim("health_2016.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2016 <- health_2016 [complete.cases(health_2016$`Feature Easting`), ]
health_2016 <- health_2016 [complete.cases(health_2016$`Feature Northing`), ]
names(health_2016)[names(health_2016)=="PointX Classification Code"] <- "code"

health_2017 <- read_delim("health_2017.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2017 <- health_2017 [complete.cases(health_2017$`Feature Easting`), ]
health_2017 <- health_2017 [complete.cases(health_2017$`Feature Northing`), ]
names(health_2017)[names(health_2017)=="PointX Classification Code"] <- "code"

health_2018 <- read_delim("health_2018.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Easting`), ]
health_2018 <- health_2018 [complete.cases(health_2018$`Feature Northing`), ]
names(health_2018)[names(health_2018)=="PointX Classification Code"] <- "code"


gp15 <- health_2015[which(health_2015$code==05280369), ]
gp15 <- subset(gp15, "Feature Easting" != "" | "Feature Northing" != "")
gp15$poi_ID <- 1:nrow(gp15)
coords <- cbind(Easting = as.numeric(as.character(gp15$`Feature Easting`)),
                Northing = as.numeric(as.character(gp15$`Feature Northing`)))
gp15sp <- SpatialPointsDataFrame(coords, data = data.frame(gp15$Name,
                                                               gp15$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp15 <- over(uktowns, geometry(gp15sp), returnList = TRUE)
gp15 <- sapply(gp15, length)
gp15 <- spCbind(uktowns, gp15)
gp15 <- merge(gp15, popbua, by="BUA11CD")
gp15 <- gp15[which(gp15$poplsoa<175000&gp15$poplsoa>10000), ]
gp15cent <- SpatialPointsDataFrame(gCentroid(gp15, byid=TRUE), gp15@data, match.ID=FALSE)
gp15 <- gp15[-c(2,3,4,5,7)]


gp <-  as.data.frame(gp15)

write.csv(gp, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps15.csv")


gp14 <- health_2014[which(health_2014$code==05280369), ]
gp14 <- subset(gp14, "Feature Easting" != "" | "Feature Northing" != "")
gp14$poi_ID <- 1:nrow(gp14)
coords <- cbind(Easting = as.numeric(as.character(gp14$`Feature Easting`)),
                Northing = as.numeric(as.character(gp14$`Feature Northing`)))
gp14sp <- SpatialPointsDataFrame(coords, data = data.frame(gp14$Name,
                                                               gp14$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp14 <- over(uktowns, geometry(gp14sp), returnList = TRUE)
gp14 <- sapply(gp14, length)
gp14 <- spCbind(uktowns, gp14)
gp14 <- merge(gp14, popbua, by="BUA11CD")
gp14 <- gp14[which(gp14$poplsoa<175000&gp14$poplsoa>10000), ]
gp14cent <- SpatialPointsDataFrame(gCentroid(gp14, byid=TRUE), gp14@data, match.ID=FALSE)
gp14 <- gp14[-c(2,3,4,5,7)]


gp <-  as.data.frame(gp14)

write.csv(gp, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps14.csv")


gp16 <- health_2016[which(health_2016$code==05280369), ]
gp16 <- subset(gp16, "Feature Easting" != "" | "Feature Northing" != "")
gp16$poi_ID <- 1:nrow(gp16)
coords <- cbind(Easting = as.numeric(as.character(gp16$`Feature Easting`)),
                Northing = as.numeric(as.character(gp16$`Feature Northing`)))
gp16sp <- SpatialPointsDataFrame(coords, data = data.frame(gp16$Name,
                                                               gp16$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp16 <- over(uktowns, geometry(gp16sp), returnList = TRUE)
gp16 <- sapply(gp16, length)
gp16 <- spCbind(uktowns, gp16)
gp16 <- merge(gp16, popbua, by="BUA11CD")
gp16 <- gp16[which(gp16$poplsoa<175000&gp16$poplsoa>10000), ]
gp16cent <- SpatialPointsDataFrame(gCentroid(gp16, byid=TRUE), gp16@data, match.ID=FALSE)
gp16 <- gp16[-c(2,3,4,5,7)]


gp <-  as.data.frame(gp16)

write.csv(gp, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps16.csv")


gp17 <- health_2017[which(health_2017$code==05280369), ]
gp17 <- subset(gp17, "Feature Easting" != "" | "Feature Northing" != "")
gp17$poi_ID <- 1:nrow(gp17)
coords <- cbind(Easting = as.numeric(as.character(gp17$`Feature Easting`)),
                Northing = as.numeric(as.character(gp17$`Feature Northing`)))
gp17sp <- SpatialPointsDataFrame(coords, data = data.frame(gp17$Name,
                                                               gp17$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp17 <- over(uktowns, geometry(gp17sp), returnList = TRUE)
gp17 <- sapply(gp17, length)
gp17 <- spCbind(uktowns, gp17)
gp17 <- merge(gp17, popbua, by="BUA11CD")
gp17 <- gp17[which(gp17$poplsoa<175000&gp17$poplsoa>10000), ]
gp17cent <- SpatialPointsDataFrame(gCentroid(gp17, byid=TRUE), gp17@data, match.ID=FALSE)
gp17 <- gp17[-c(2,3,4,5,7)]


gp <-  as.data.frame(gp17)

write.csv(gp, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps17.csv")


gp18 <- health_2018[which(health_2018$code==05280369), ]
gp18 <- subset(gp18, "Feature Easting" != "" | "Feature Northing" != "")
gp18$poi_ID <- 1:nrow(gp18)
coords <- cbind(Easting = as.numeric(as.character(gp18$`Feature Easting`)),
                Northing = as.numeric(as.character(gp18$`Feature Northing`)))
gp18sp <- SpatialPointsDataFrame(coords, data = data.frame(gp18$Name,
                                                               gp18$poi_ID), proj4string = CRS("+init=epsg:27700"))


gp18 <- over(uktowns, geometry(gp18sp), returnList = TRUE)
gp18 <- sapply(gp18, length)
gp18 <- spCbind(uktowns, gp18)
gp18 <- merge(gp18, popbua, by="BUA11CD")
gp18 <- gp18[which(gp18$poplsoa<175000&gp18$poplsoa>10000), ]
gp18cent <- SpatialPointsDataFrame(gCentroid(gp18, byid=TRUE), gp18@data, match.ID=FALSE)
gp18 <- gp18[-c(2,3,4,5,7)]


gp <-  as.data.frame(gp18)

write.csv(gp, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps18.csv")


