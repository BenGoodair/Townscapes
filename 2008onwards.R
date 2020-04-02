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


bet13 <- POI_2013[which(POI_2013$code==4220279), ]
bet13 <- subset(bet13, "feature_easting" != "" | "feature_northing" != "")
bet13$poi_ID <- 1:nrow(bet13)
coords <- cbind(Easting = as.numeric(as.character(bet13$feature_easting)),
                Northing = as.numeric(as.character(bet13$feature_northing)))
bet13sp <- SpatialPointsDataFrame(coords, data = data.frame(bet13$name,
                                                            bet13$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet13 <- over(uktowns, geometry(bet13sp), returnList = TRUE)
bet13 <- sapply(bet13, length)
bet13 <- spCbind(uktowns, bet13)
bet13 <- merge(bet13, popbua, by="BUA11CD")
bet13 <- bet13[which(bet13$poplsoa<175000&bet13$poplsoa>10000), ]
bet13cent <- SpatialPointsDataFrame(gCentroid(bet13, byid=TRUE), bet13@data, match.ID=FALSE)
bet13 <- bet13[-c(2,3,4,5,7)]

cine13 <- POI_2013[which(POI_2013$code==4250308 | POI_2013$code==4250315), ]
cine13 <- subset(cine13, "feature_easting" != "" | "feature_northing" != "")
cine13$poi_ID <- 1:nrow(cine13)
coords <- cbind(Easting = as.numeric(as.character(cine13$feature_easting)),
                Northing = as.numeric(as.character(cine13$feature_northing)))
cine13sp <- SpatialPointsDataFrame(coords, data = data.frame(cine13$name,
                                                             cine13$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine13 <- over(uktowns, geometry(cine13sp), returnList = TRUE)
cine13 <- sapply(cine13, length)
cine13 <- spCbind(uktowns, cine13)
cine13 <- merge(cine13, popbua, by="BUA11CD")
cine13 <- cine13[which(cine13$poplsoa<175000&cine13$poplsoa>10000), ]
cine13cent <- SpatialPointsDataFrame(gCentroid(cine13, byid=TRUE), cine13@data, match.ID=FALSE)
cine13 <- cine13[-c(2,3,4,5,7)]

night13 <- POI_2013[which(POI_2013$code==4250312 | POI_2013$code==4250311), ]
night13 <- subset(night13, "feature_easting" != "" | "feature_northing" != "")
night13$poi_ID <- 1:nrow(night13)
coords <- cbind(Easting = as.numeric(as.character(night13$feature_easting)),
                Northing = as.numeric(as.character(night13$feature_northing)))
night13sp <- SpatialPointsDataFrame(coords, data = data.frame(night13$name,
                                                              night13$poi_ID), proj4string = CRS("+init=epsg:27700"))


night13 <- over(uktowns, geometry(night13sp), returnList = TRUE)
night13 <- sapply(night13, length)
night13 <- spCbind(uktowns, night13)
night13 <- merge(night13, popbua, by="BUA11CD")
night13 <- night13[which(night13$poplsoa<175000&night13$poplsoa>10000), ]
night13cent <- SpatialPointsDataFrame(gCentroid(night13, byid=TRUE), night13@data, match.ID=FALSE)
night13 <- night13[-c(2,3,4,5,7)]


sport13 <- POI_2013[which(POI_2013$code>=4240000 & POI_2013$code<4250000), ]
sport13 <- subset(sport13, "feature_easting" != "" | "feature_northing" != "")
sport13$poi_ID <- 1:nrow(sport13)
coords <- cbind(Easting = as.numeric(as.character(sport13$feature_easting)),
                Northing = as.numeric(as.character(sport13$feature_northing)))
sport13sp <- SpatialPointsDataFrame(coords, data = data.frame(sport13$name,
                                                              sport13$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport13 <- over(uktowns, geometry(sport13sp), returnList = TRUE)
sport13 <- sapply(sport13, length)
sport13 <- spCbind(uktowns, sport13)
sport13 <- merge(sport13, popbua, by="BUA11CD")
sport13 <- sport13[which(sport13$poplsoa<175000&sport13$poplsoa>10000), ]
sport13cent <- SpatialPointsDataFrame(gCentroid(sport13, byid=TRUE), sport13@data, match.ID=FALSE)
sport13 <- sport13[-c(2,3,4,5,7)]


creative13 <- POI_2013[which(POI_2013$code==5320384 | POI_2013$code==5320394| POI_2013$code==5320395| POI_2013$code==5320396| POI_2013$code==5320389), ]
creative13 <- subset(creative13, "feature_easting" != "" | "feature_northing" != "")
creative13$poi_ID <- 1:nrow(creative13)
coords <- cbind(Easting = as.numeric(as.character(creative13$feature_easting)),
                Northing = as.numeric(as.character(creative13$feature_northing)))
creative13sp <- SpatialPointsDataFrame(coords, data = data.frame(creative13$name,
                                                                 creative13$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative13 <- over(uktowns, geometry(creative13sp), returnList = TRUE)
creative13 <- sapply(creative13, length)
creative13 <- spCbind(uktowns, creative13)
creative13 <- merge(creative13, popbua, by="BUA11CD")
creative13 <- creative13[which(creative13$poplsoa<175000&creative13$poplsoa>10000), ]
creative13cent <- SpatialPointsDataFrame(gCentroid(creative13, byid=TRUE), creative13@data, match.ID=FALSE)
creative13 <- creative13[-c(2,3,4,5,7)]





bet <-  as.data.frame(bet13)
creative <-  as.data.frame(creative13)
sport <-  as.data.frame(sport13)
night <-  as.data.frame(night13)
cine <-  as.data.frame(cine13)




##rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])

POI_2012 <- read_csv("D:/poi/POI_2012.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2012 <- POI_2012[complete.cases(POI_2012$feature_easting), ]
POI_2012 <- POI_2012[complete.cases(POI_2012$feature_northing), ]
names(POI_2012)[names(POI_2012)=="pointx_classification_code"] <- "code"


bet12 <- POI_2012[which(POI_2012$code==4220279), ]
bet12 <- subset(bet12, "feature_easting" != "" | "feature_northing" != "")
bet12$poi_ID <- 1:nrow(bet12)
coords <- cbind(Easting = as.numeric(as.character(bet12$feature_easting)),
                Northing = as.numeric(as.character(bet12$feature_northing)))
bet12sp <- SpatialPointsDataFrame(coords, data = data.frame(bet12$name,
                                                            bet12$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet12 <- over(uktowns, geometry(bet12sp), returnList = TRUE)
bet12 <- sapply(bet12, length)
bet12 <- spCbind(uktowns, bet12)
bet12 <- merge(bet12, popbua, by="BUA11CD")
bet12 <- bet12[which(bet12$poplsoa<175000&bet12$poplsoa>10000), ]
bet12cent <- SpatialPointsDataFrame(gCentroid(bet12, byid=TRUE), bet12@data, match.ID=FALSE)
bet12 <- bet12[-c(2,3,4,5,7)]

cine12 <- POI_2012[which(POI_2012$code==4250308 | POI_2012$code==4250315), ]
cine12 <- subset(cine12, "feature_easting" != "" | "feature_northing" != "")
cine12$poi_ID <- 1:nrow(cine12)
coords <- cbind(Easting = as.numeric(as.character(cine12$feature_easting)),
                Northing = as.numeric(as.character(cine12$feature_northing)))
cine12sp <- SpatialPointsDataFrame(coords, data = data.frame(cine12$name,
                                                             cine12$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine12 <- over(uktowns, geometry(cine12sp), returnList = TRUE)
cine12 <- sapply(cine12, length)
cine12 <- spCbind(uktowns, cine12)
cine12 <- merge(cine12, popbua, by="BUA11CD")
cine12 <- cine12[which(cine12$poplsoa<175000&cine12$poplsoa>10000), ]
cine12cent <- SpatialPointsDataFrame(gCentroid(cine12, byid=TRUE), cine12@data, match.ID=FALSE)
cine12 <- cine12[-c(2,3,4,5,7)]

night12 <- POI_2012[which(POI_2012$code==4250312 | POI_2012$code==4250311), ]
night12 <- subset(night12, "feature_easting" != "" | "feature_northing" != "")
night12$poi_ID <- 1:nrow(night12)
coords <- cbind(Easting = as.numeric(as.character(night12$feature_easting)),
                Northing = as.numeric(as.character(night12$feature_northing)))
night12sp <- SpatialPointsDataFrame(coords, data = data.frame(night12$name,
                                                              night12$poi_ID), proj4string = CRS("+init=epsg:27700"))


night12 <- over(uktowns, geometry(night12sp), returnList = TRUE)
night12 <- sapply(night12, length)
night12 <- spCbind(uktowns, night12)
night12 <- merge(night12, popbua, by="BUA11CD")
night12 <- night12[which(night12$poplsoa<175000&night12$poplsoa>10000), ]
night12cent <- SpatialPointsDataFrame(gCentroid(night12, byid=TRUE), night12@data, match.ID=FALSE)
night12 <- night12[-c(2,3,4,5,7)]


sport12 <- POI_2012[which(POI_2012$code>=4240000 & POI_2012$code< 4250000), ]
sport12 <- subset(sport12, "feature_easting" != "" | "feature_northing" != "")
sport12$poi_ID <- 1:nrow(sport12)
coords <- cbind(Easting = as.numeric(as.character(sport12$feature_easting)),
                Northing = as.numeric(as.character(sport12$feature_northing)))
sport12sp <- SpatialPointsDataFrame(coords, data = data.frame(sport12$name,
                                                              sport12$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport12 <- over(uktowns, geometry(sport12sp), returnList = TRUE)
sport12 <- sapply(sport12, length)
sport12 <- spCbind(uktowns, sport12)
sport12 <- merge(sport12, popbua, by="BUA11CD")
sport12 <- sport12[which(sport12$poplsoa<175000&sport12$poplsoa>10000), ]
sport12cent <- SpatialPointsDataFrame(gCentroid(sport12, byid=TRUE), sport12@data, match.ID=FALSE)
sport12 <- sport12[-c(2,3,4,5,7)]


creative12 <- POI_2012[which(POI_2012$code==5320384 | POI_2012$code==5320394| POI_2012$code==5320395| POI_2012$code==5320396| POI_2012$code==5320389), ]
creative12 <- subset(creative12, "feature_easting" != "" | "feature_northing" != "")
creative12$poi_ID <- 1:nrow(creative12)
coords <- cbind(Easting = as.numeric(as.character(creative12$feature_easting)),
                Northing = as.numeric(as.character(creative12$feature_northing)))
creative12sp <- SpatialPointsDataFrame(coords, data = data.frame(creative12$name,
                                                                 creative12$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative12 <- over(uktowns, geometry(creative12sp), returnList = TRUE)
creative12 <- sapply(creative12, length)
creative12 <- spCbind(uktowns, creative12)
creative12 <- merge(creative12, popbua, by="BUA11CD")
creative12 <- creative12[which(creative12$poplsoa<175000&creative12$poplsoa>10000), ]
creative12cent <- SpatialPointsDataFrame(gCentroid(creative12, byid=TRUE), creative12@data, match.ID=FALSE)
creative12 <- creative12[-c(2,3,4,5,7)]


bet <- merge(bet, bet12, by="BUA11CD")
creative <- merge(creative, creative12, by="BUA11CD")
sport <- merge(sport, sport12, by="BUA11CD")
night <- merge(night, night12, by="BUA11CD")
cine <- merge(cine, cine12, by="BUA11CD")


##rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])

POI_2011 <- read_csv("D:/poi/POI_2011.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2011 <- POI_2011[complete.cases(POI_2011$feature_easting), ]
POI_2011 <- POI_2011[complete.cases(POI_2011$feature_northing), ]
names(POI_2011)[names(POI_2011)=="pointx_classification_code"] <- "code"


bet11 <- POI_2011[which(POI_2011$code==4220279), ]
bet11 <- subset(bet11, "feature_easting" != "" | "feature_northing" != "")
bet11$poi_ID <- 1:nrow(bet11)
coords <- cbind(Easting = as.numeric(as.character(bet11$feature_easting)),
                Northing = as.numeric(as.character(bet11$feature_northing)))
bet11sp <- SpatialPointsDataFrame(coords, data = data.frame(bet11$name,
                                                            bet11$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet11 <- over(uktowns, geometry(bet11sp), returnList = TRUE)
bet11 <- sapply(bet11, length)
bet11 <- spCbind(uktowns, bet11)
bet11 <- merge(bet11, popbua, by="BUA11CD")
bet11 <- bet11[which(bet11$poplsoa<175000&bet11$poplsoa>10000), ]
bet11cent <- SpatialPointsDataFrame(gCentroid(bet11, byid=TRUE), bet11@data, match.ID=FALSE)
bet11 <- bet11[-c(2,3,4,5,7)]

cine11 <- POI_2011[which(POI_2011$code==4250308 | POI_2011$code==4250315), ]
cine11 <- subset(cine11, "feature_easting" != "" | "feature_northing" != "")
cine11$poi_ID <- 1:nrow(cine11)
coords <- cbind(Easting = as.numeric(as.character(cine11$feature_easting)),
                Northing = as.numeric(as.character(cine11$feature_northing)))
cine11sp <- SpatialPointsDataFrame(coords, data = data.frame(cine11$name,
                                                             cine11$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine11 <- over(uktowns, geometry(cine11sp), returnList = TRUE)
cine11 <- sapply(cine11, length)
cine11 <- spCbind(uktowns, cine11)
cine11 <- merge(cine11, popbua, by="BUA11CD")
cine11 <- cine11[which(cine11$poplsoa<175000&cine11$poplsoa>10000), ]
cine11cent <- SpatialPointsDataFrame(gCentroid(cine11, byid=TRUE), cine11@data, match.ID=FALSE)
cine11 <- cine11[-c(2,3,4,5,7)]

night11 <- POI_2011[which(POI_2011$code==4250312 | POI_2011$code==4250311), ]
night11 <- subset(night11, "feature_easting" != "" | "feature_northing" != "")
night11$poi_ID <- 1:nrow(night11)
coords <- cbind(Easting = as.numeric(as.character(night11$feature_easting)),
                Northing = as.numeric(as.character(night11$feature_northing)))
night11sp <- SpatialPointsDataFrame(coords, data = data.frame(night11$name,
                                                              night11$poi_ID), proj4string = CRS("+init=epsg:27700"))


night11 <- over(uktowns, geometry(night11sp), returnList = TRUE)
night11 <- sapply(night11, length)
night11 <- spCbind(uktowns, night11)
night11 <- merge(night11, popbua, by="BUA11CD")
night11 <- night11[which(night11$poplsoa<175000&night11$poplsoa>10000), ]
night11cent <- SpatialPointsDataFrame(gCentroid(night11, byid=TRUE), night11@data, match.ID=FALSE)
night11 <- night11[-c(2,3,4,5,7)]


sport11 <- POI_2011[which(POI_2011$code>=4240000 & POI_2011$code< 4250000), ]
sport11 <- subset(sport11, "feature_easting" != "" | "feature_northing" != "")
sport11$poi_ID <- 1:nrow(sport11)
coords <- cbind(Easting = as.numeric(as.character(sport11$feature_easting)),
                Northing = as.numeric(as.character(sport11$feature_northing)))
sport11sp <- SpatialPointsDataFrame(coords, data = data.frame(sport11$name,
                                                              sport11$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport11 <- over(uktowns, geometry(sport11sp), returnList = TRUE)
sport11 <- sapply(sport11, length)
sport11 <- spCbind(uktowns, sport11)
sport11 <- merge(sport11, popbua, by="BUA11CD")
sport11 <- sport11[which(sport11$poplsoa<175000&sport11$poplsoa>10000), ]
sport11cent <- SpatialPointsDataFrame(gCentroid(sport11, byid=TRUE), sport11@data, match.ID=FALSE)
sport11 <- sport11[-c(2,3,4,5,7)]


creative11 <- POI_2011[which(POI_2011$code==5320384 | POI_2011$code==5320394| POI_2011$code==5320395| POI_2011$code==5320396| POI_2011$code==5320389), ]
creative11 <- subset(creative11, "feature_easting" != "" | "feature_northing" != "")
creative11$poi_ID <- 1:nrow(creative11)
coords <- cbind(Easting = as.numeric(as.character(creative11$feature_easting)),
                Northing = as.numeric(as.character(creative11$feature_northing)))
creative11sp <- SpatialPointsDataFrame(coords, data = data.frame(creative11$name,
                                                                 creative11$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative11 <- over(uktowns, geometry(creative11sp), returnList = TRUE)
creative11 <- sapply(creative11, length)
creative11 <- spCbind(uktowns, creative11)
creative11 <- merge(creative11, popbua, by="BUA11CD")
creative11 <- creative11[which(creative11$poplsoa<175000&creative11$poplsoa>10000), ]
creative11cent <- SpatialPointsDataFrame(gCentroid(creative11, byid=TRUE), creative11@data, match.ID=FALSE)
creative11 <- creative11[-c(2,3,4,5,7)]


bet <-  as.data.frame(bet11)
creative <-  as.data.frame(creative11)
sport <-  as.data.frame(sport11)
night <-  as.data.frame(night11)
cine <-  as.data.frame(cine11)


write.csv(creative, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\creative11.csv")
write.csv(night, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\night11.csv")
write.csv(bet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\bet11.csv")
write.csv(cine, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\cine11.csv")
write.csv(sport, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\sport11.csv")


bet <- merge(bet, bet11, by="BUA11CD")
creative <- merge(creative, creative11, by="BUA11CD")
sport <- merge(sport, sport11, by="BUA11CD")
night <- merge(night, night11, by="BUA11CD")
cine <- merge(cine, cine11, by="BUA11CD")



rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])

POI_2010 <- read_csv("D:/poi/POI_2010.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2010 <- POI_2010[complete.cases(POI_2010$easting), ]
POI_2010 <- POI_2010[complete.cases(POI_2010$northing), ]
names(POI_2010)[names(POI_2010)=="pointx_classification_code"] <- "code"


bet10 <- POI_2010[which(POI_2010$code==4220279), ]
bet10 <- subset(bet10, "easting" != "" | "northing" != "")
bet10$poi_ID <- 1:nrow(bet10)
coords <- cbind(Easting = as.numeric(as.character(bet10$easting)),
                Northing = as.numeric(as.character(bet10$northing)))
bet10sp <- SpatialPointsDataFrame(coords, data = data.frame(bet10$name,
                                                            bet10$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet10 <- over(uktowns, geometry(bet10sp), returnList = TRUE)
bet10 <- sapply(bet10, length)
bet10 <- spCbind(uktowns, bet10)
bet10 <- merge(bet10, popbua, by="BUA11CD")
bet10 <- bet10[which(bet10$poplsoa<175000&bet10$poplsoa>10000), ]
bet10cent <- SpatialPointsDataFrame(gCentroid(bet10, byid=TRUE), bet10@data, match.ID=FALSE)
bet10 <- bet10[-c(2,3,4,5,7)]

cine10 <- POI_2010[which(POI_2010$code==4250308 | POI_2010$code==4250315), ]
cine10 <- subset(cine10, "easting" != "" | "northing" != "")
cine10$poi_ID <- 1:nrow(cine10)
coords <- cbind(Easting = as.numeric(as.character(cine10$easting)),
                Northing = as.numeric(as.character(cine10$northing)))
cine10sp <- SpatialPointsDataFrame(coords, data = data.frame(cine10$name,
                                                             cine10$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine10 <- over(uktowns, geometry(cine10sp), returnList = TRUE)
cine10 <- sapply(cine10, length)
cine10 <- spCbind(uktowns, cine10)
cine10 <- merge(cine10, popbua, by="BUA11CD")
cine10 <- cine10[which(cine10$poplsoa<175000&cine10$poplsoa>10000), ]
cine10cent <- SpatialPointsDataFrame(gCentroid(cine10, byid=TRUE), cine10@data, match.ID=FALSE)
cine10 <- cine10[-c(2,3,4,5,7)]

night10 <- POI_2010[which(POI_2010$code==4250312 | POI_2010$code==4250311), ]
night10 <- subset(night10, "easting" != "" | "northing" != "")
night10$poi_ID <- 1:nrow(night10)
coords <- cbind(Easting = as.numeric(as.character(night10$easting)),
                Northing = as.numeric(as.character(night10$northing)))
night10sp <- SpatialPointsDataFrame(coords, data = data.frame(night10$name,
                                                              night10$poi_ID), proj4string = CRS("+init=epsg:27700"))


night10 <- over(uktowns, geometry(night10sp), returnList = TRUE)
night10 <- sapply(night10, length)
night10 <- spCbind(uktowns, night10)
night10 <- merge(night10, popbua, by="BUA11CD")
night10 <- night10[which(night10$poplsoa<175000&night10$poplsoa>10000), ]
night10cent <- SpatialPointsDataFrame(gCentroid(night10, byid=TRUE), night10@data, match.ID=FALSE)
night10 <- night10[-c(2,3,4,5,7)]


sport10 <- POI_2010[which(POI_2010$code>=4240000 & POI_2010$code<4250000), ]
sport10 <- subset(sport10, "easting" != "" | "northing" != "")
sport10$poi_ID <- 1:nrow(sport10)
coords <- cbind(Easting = as.numeric(as.character(sport10$easting)),
                Northing = as.numeric(as.character(sport10$northing)))
sport10sp <- SpatialPointsDataFrame(coords, data = data.frame(sport10$name,
                                                              sport10$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport10 <- over(uktowns, geometry(sport10sp), returnList = TRUE)
sport10 <- sapply(sport10, length)
sport10 <- spCbind(uktowns, sport10)
sport10 <- merge(sport10, popbua, by="BUA11CD")
sport10 <- sport10[which(sport10$poplsoa<175000&sport10$poplsoa>10000), ]
sport10cent <- SpatialPointsDataFrame(gCentroid(sport10, byid=TRUE), sport10@data, match.ID=FALSE)
sport10 <- sport10[-c(2,3,4,5,7)]


creative10 <- POI_2010[which(POI_2010$code==5320384 | POI_2010$code==5320394| POI_2010$code==5320395| POI_2010$code==5320396| POI_2010$code==5320389), ]
creative10 <- subset(creative10, "easting" != "" | "northing" != "")
creative10$poi_ID <- 1:nrow(creative10)
coords <- cbind(Easting = as.numeric(as.character(creative10$easting)),
                Northing = as.numeric(as.character(creative10$northing)))
creative10sp <- SpatialPointsDataFrame(coords, data = data.frame(creative10$name,
                                                                 creative10$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative10 <- over(uktowns, geometry(creative10sp), returnList = TRUE)
creative10 <- sapply(creative10, length)
creative10 <- spCbind(uktowns, creative10)
creative10 <- merge(creative10, popbua, by="BUA11CD")
creative10 <- creative10[which(creative10$poplsoa<175000&creative10$poplsoa>10000), ]
creative10cent <- SpatialPointsDataFrame(gCentroid(creative10, byid=TRUE), creative10@data, match.ID=FALSE)
creative10 <- creative10[-c(2,3,4,5,7)]


bet <-  as.data.frame(bet10)
creative <-  as.data.frame(creative10)
sport <-  as.data.frame(sport10)
night <-  as.data.frame(night10)
cine <-  as.data.frame(cine10)

write.csv(creative, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\creative10.csv")
write.csv(night, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\night10.csv")
write.csv(bet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\bet10.csv")
write.csv(cine, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\cine10.csv")
write.csv(sport, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\sport10.csv")


bet <- merge(bet, bet11, by="BUA11CD")
creative <- merge(creative, creative11, by="BUA11CD")
sport <- merge(sport, sport11, by="BUA11CD")
night <- merge(night, night11, by="BUA11CD")
cine <- merge(cine, cine11, by="BUA11CD")


bet <- merge(bet, bet10, by="BUA11CD")
creative <- merge(creative, creative10, by="BUA11CD")
sport <- merge(sport, sport10, by="BUA11CD")
night <- merge(night, night10, by="BUA11CD")
cine <- merge(cine, cine10, by="BUA11CD")




rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])


POI_2009 <- read_csv("D:/poi/POI_2009.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2009 <- POI_2009[complete.cases(POI_2009$easting), ]
POI_2009 <- POI_2009[complete.cases(POI_2009$northing), ]
names(POI_2009)[names(POI_2009)=="pointx_classification_code"] <- "code"

bet09 <- POI_2009[which(POI_2009$code==4220279), ]
bet09 <- subset(bet09, "easting" != "" | "northing" != "")
bet09$poi_ID <- 1:nrow(bet09)
coords <- cbind(Easting = as.numeric(as.character(bet09$easting)),
                Northing = as.numeric(as.character(bet09$northing)))
bet09sp <- SpatialPointsDataFrame(coords, data = data.frame(bet09$name,
                                                            bet09$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet09 <- over(uktowns, geometry(bet09sp), returnList = TRUE)
bet09 <- sapply(bet09, length)
bet09 <- spCbind(uktowns, bet09)
bet09 <- merge(bet09, popbua, by="BUA11CD")
bet09 <- bet09[which(bet09$poplsoa<175000&bet09$poplsoa>10000), ]
bet09cent <- SpatialPointsDataFrame(gCentroid(bet09, byid=TRUE), bet09@data, match.ID=FALSE)
bet09 <- bet09[-c(2,3,4,5,7)]

cine09 <- POI_2009[which(POI_2009$code==4250308 | POI_2009$code==4250315), ]
cine09 <- subset(cine09, "easting" != "" | "northing" != "")
cine09$poi_ID <- 1:nrow(cine09)
coords <- cbind(Easting = as.numeric(as.character(cine09$easting)),
                Northing = as.numeric(as.character(cine09$northing)))
cine09sp <- SpatialPointsDataFrame(coords, data = data.frame(cine09$name,
                                                             cine09$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine09 <- over(uktowns, geometry(cine09sp), returnList = TRUE)
cine09 <- sapply(cine09, length)
cine09 <- spCbind(uktowns, cine09)
cine09 <- merge(cine09, popbua, by="BUA11CD")
cine09 <- cine09[which(cine09$poplsoa<175000&cine09$poplsoa>10000), ]
cine09cent <- SpatialPointsDataFrame(gCentroid(cine09, byid=TRUE), cine09@data, match.ID=FALSE)
cine09 <- cine09[-c(2,3,4,5,7)]

night09 <- POI_2009[which(POI_2009$code==4250312 | POI_2009$code==4250311), ]
night09 <- subset(night09, "easting" != "" | "northing" != "")
night09$poi_ID <- 1:nrow(night09)
coords <- cbind(Easting = as.numeric(as.character(night09$easting)),
                Northing = as.numeric(as.character(night09$northing)))
night09sp <- SpatialPointsDataFrame(coords, data = data.frame(night09$name,
                                                              night09$poi_ID), proj4string = CRS("+init=epsg:27700"))


night09 <- over(uktowns, geometry(night09sp), returnList = TRUE)
night09 <- sapply(night09, length)
night09 <- spCbind(uktowns, night09)
night09 <- merge(night09, popbua, by="BUA11CD")
night09 <- night09[which(night09$poplsoa<175000&night09$poplsoa>10000), ]
night09cent <- SpatialPointsDataFrame(gCentroid(night09, byid=TRUE), night09@data, match.ID=FALSE)
night09 <- night09[-c(2,3,4,5,7)]


sport09 <- POI_2009[which(POI_2009$code>=4240000 & POI_2009$code< 4250000), ]
sport09 <- subset(sport09, "easting" != "" | "northing" != "")
sport09$poi_ID <- 1:nrow(sport09)
coords <- cbind(Easting = as.numeric(as.character(sport09$easting)),
                Northing = as.numeric(as.character(sport09$northing)))
sport09sp <- SpatialPointsDataFrame(coords, data = data.frame(sport09$name,
                                                              sport09$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport09 <- over(uktowns, geometry(sport09sp), returnList = TRUE)
sport09 <- sapply(sport09, length)
sport09 <- spCbind(uktowns, sport09)
sport09 <- merge(sport09, popbua, by="BUA11CD")
sport09 <- sport09[which(sport09$poplsoa<175000&sport09$poplsoa>10000), ]
sport09cent <- SpatialPointsDataFrame(gCentroid(sport09, byid=TRUE), sport09@data, match.ID=FALSE)
sport09 <- sport09[-c(2,3,4,5,7)]


creative09 <- POI_2009[which(POI_2009$code==5320384 | POI_2009$code==5320394| POI_2009$code==5320395| POI_2009$code==5320396| POI_2009$code==5320389), ]
creative09 <- subset(creative09, "easting" != "" | "northing" != "")
creative09$poi_ID <- 1:nrow(creative09)
coords <- cbind(Easting = as.numeric(as.character(creative09$easting)),
                Northing = as.numeric(as.character(creative09$northing)))
creative09sp <- SpatialPointsDataFrame(coords, data = data.frame(creative09$name,
                                                                 creative09$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative09 <- over(uktowns, geometry(creative09sp), returnList = TRUE)
creative09 <- sapply(creative09, length)
creative09 <- spCbind(uktowns, creative09)
creative09 <- merge(creative09, popbua, by="BUA11CD")
creative09 <- creative09[which(creative09$poplsoa<175000&creative09$poplsoa>10000), ]
creative09cent <- SpatialPointsDataFrame(gCentroid(creative09, byid=TRUE), creative09@data, match.ID=FALSE)
creative09 <- creative09[-c(2,3,4,5,7)]



bet <-  as.data.frame(bet09)
creative <-  as.data.frame(creative09)
sport <-  as.data.frame(sport09)
night <-  as.data.frame(night09)
cine <-  as.data.frame(cine09)

write.csv(creative, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\creative09.csv")
write.csv(night, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\night09.csv")
write.csv(bet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\bet09.csv")
write.csv(cine, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\cine09.csv")
write.csv(sport, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\sport09.csv")



bet <- merge(bet, bet09, by="BUA11CD")
creative <- merge(creative, creative09, by="BUA11CD")
sport <- merge(sport, sport09, by="BUA11CD")
night <- merge(night, night09, by="BUA11CD")
cine <- merge(cine, cine09, by="BUA11CD")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])


POI_2008 <- read_csv("D:/poi/POI_2008.csv", 
                     col_types = cols(pointx_classification_code = col_number()))
POI_2008 <- POI_2008[complete.cases(POI_2008$easting), ]
POI_2008 <- POI_2008[complete.cases(POI_2008$northing), ]
names(POI_2008)[names(POI_2008)=="pointx_classification_code"] <- "code"


bet08 <- POI_2008[which(POI_2008$code==4220279), ]
bet08 <- subset(bet08, "easting" != "" | "northing" != "")
bet08$poi_ID <- 1:nrow(bet08)
coords <- cbind(Easting = as.numeric(as.character(bet08$easting)),
                Northing = as.numeric(as.character(bet08$northing)))
bet08sp <- SpatialPointsDataFrame(coords, data = data.frame(bet08$name,
                                                            bet08$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet08 <- over(uktowns, geometry(bet08sp), returnList = TRUE)
bet08 <- sapply(bet08, length)
bet08 <- spCbind(uktowns, bet08)
bet08 <- merge(bet08, popbua, by="BUA11CD")
bet08 <- bet08[which(bet08$poplsoa<175000&bet08$poplsoa>10000), ]
bet08cent <- SpatialPointsDataFrame(gCentroid(bet08, byid=TRUE), bet08@data, match.ID=FALSE)
bet08 <- bet08[-c(2,3,4,5,7)]

cine08 <- POI_2008[which(POI_2008$code==4250308 | POI_2008$code==4250315), ]
cine08 <- subset(cine08, "easting" != "" | "northing" != "")
cine08$poi_ID <- 1:nrow(cine08)
coords <- cbind(Easting = as.numeric(as.character(cine08$easting)),
                Northing = as.numeric(as.character(cine08$northing)))
cine08sp <- SpatialPointsDataFrame(coords, data = data.frame(cine08$name,
                                                             cine08$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine08 <- over(uktowns, geometry(cine08sp), returnList = TRUE)
cine08 <- sapply(cine08, length)
cine08 <- spCbind(uktowns, cine08)
cine08 <- merge(cine08, popbua, by="BUA11CD")
cine08 <- cine08[which(cine08$poplsoa<175000&cine08$poplsoa>10000), ]
cine08cent <- SpatialPointsDataFrame(gCentroid(cine08, byid=TRUE), cine08@data, match.ID=FALSE)
cine08 <- cine08[-c(2,3,4,5,7)]

night08 <- POI_2008[which(POI_2008$code==4250312 | POI_2008$code==4250311), ]
night08 <- subset(night08, "easting" != "" | "northing" != "")
night08$poi_ID <- 1:nrow(night08)
coords <- cbind(Easting = as.numeric(as.character(night08$easting)),
                Northing = as.numeric(as.character(night08$northing)))
night08sp <- SpatialPointsDataFrame(coords, data = data.frame(night08$name,
                                                              night08$poi_ID), proj4string = CRS("+init=epsg:27700"))


night08 <- over(uktowns, geometry(night08sp), returnList = TRUE)
night08 <- sapply(night08, length)
night08 <- spCbind(uktowns, night08)
night08 <- merge(night08, popbua, by="BUA11CD")
night08 <- night08[which(night08$poplsoa<175000&night08$poplsoa>10000), ]
night08cent <- SpatialPointsDataFrame(gCentroid(night08, byid=TRUE), night08@data, match.ID=FALSE)
night08 <- night08[-c(2,3,4,5,7)]


sport08 <- POI_2008[which(POI_2008$code>=4240000 & POI_2008$code< 4250000), ]
sport08 <- subset(sport08, "easting" != "" | "northing" != "")
sport08$poi_ID <- 1:nrow(sport08)
coords <- cbind(Easting = as.numeric(as.character(sport08$easting)),
                Northing = as.numeric(as.character(sport08$northing)))
sport08sp <- SpatialPointsDataFrame(coords, data = data.frame(sport08$name,
                                                              sport08$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport08 <- over(uktowns, geometry(sport08sp), returnList = TRUE)
sport08 <- sapply(sport08, length)
sport08 <- spCbind(uktowns, sport08)
sport08 <- merge(sport08, popbua, by="BUA11CD")
sport08 <- sport08[which(sport08$poplsoa<175000&sport08$poplsoa>10000), ]
sport08cent <- SpatialPointsDataFrame(gCentroid(sport08, byid=TRUE), sport08@data, match.ID=FALSE)
sport08 <- sport08[-c(2,3,4,5,7)]


creative08 <- POI_2008[which(POI_2008$code==5320384 | POI_2008$code==5320394| POI_2008$code==5320395| POI_2008$code==5320396| POI_2008$code==5320389), ]
creative08 <- subset(creative08, "easting" != "" | "northing" != "")
creative08$poi_ID <- 1:nrow(creative08)
coords <- cbind(Easting = as.numeric(as.character(creative08$easting)),
                Northing = as.numeric(as.character(creative08$northing)))
creative08sp <- SpatialPointsDataFrame(coords, data = data.frame(creative08$name,
                                                                 creative08$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative08 <- over(uktowns, geometry(creative08sp), returnList = TRUE)
creative08 <- sapply(creative08, length)
creative08 <- spCbind(uktowns, creative08)
creative08 <- merge(creative08, popbua, by="BUA11CD")
creative08 <- creative08[which(creative08$poplsoa<175000&creative08$poplsoa>10000), ]
creative08cent <- SpatialPointsDataFrame(gCentroid(creative08, byid=TRUE), creative08@data, match.ID=FALSE)
creative08 <- creative08[-c(2,3,4,5,7)]

bet <-  as.data.frame(bet08)
creative <-  as.data.frame(creative08)
sport <-  as.data.frame(sport08)
night <-  as.data.frame(night08)
cine <-  as.data.frame(cine08)


write.csv(creative, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\creative08.csv")
write.csv(night, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\night08.csv")
write.csv(bet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\bet08.csv")
write.csv(cine, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\cine08.csv")
write.csv(sport, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\sport08.csv")




rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])


###14 onwards###



creative15 <- read_delim("educatjun15.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
creative15 <- creative15[complete.cases(creative15$`Feature Easting`), ]
creative15 <- creative15[complete.cases(creative15$`Feature Northing`), ]
names(creative15)[names(creative15)=="PointX Classification Code"] <- "code"
creative15 <- creative15[which(creative15$code==05320384 | creative15$code==05320394| creative15$code==05320395| creative15$code==05320396| creative15$code==05320389), ]
creative15 <- subset(creative15, "Feature Easting" != "" | "Feature Northing" != "")
creative15$poi_ID <- 1:nrow(creative15)
coords <- cbind(Easting = as.numeric(as.character(creative15$"Feature Easting")),
                Northing = as.numeric(as.character(creative15$"Feature Northing")))
creative15sp <- SpatialPointsDataFrame(coords, data = data.frame(creative15$Name,
                                                                   creative15$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative15 <- over(uktowns, geometry(creative15sp), returnList = TRUE)
creative15 <- sapply(creative15, length)
creative15 <- spCbind(uktowns, creative15)
creative15 <- merge(creative15, popbua, by="BUA11CD")
creative15 <- creative15[which(creative15$poplsoa<175000&creative15$poplsoa>10000), ]
creative15cent <- SpatialPointsDataFrame(gCentroid(creative15, byid=TRUE), creative15@data, match.ID=FALSE)
creative15 <- creative15[-c(2,3,4,5,7)]


creative18 <- read_delim("educatjun18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
creative18 <- creative18[complete.cases(creative18$`Feature Easting`), ]
creative18 <- creative18[complete.cases(creative18$`Feature Northing`), ]
names(creative18)[names(creative18)=="PointX Classification Code"] <- "code"
creative18 <- creative18[which(creative18$code==05320384 | creative18$code==05320394| creative18$code==05320395| creative18$code==05320396| creative18$code==05320389), ]
creative18 <- subset(creative18, "Feature Easting" != "" | "Feature Northing" != "")
creative18$poi_ID <- 1:nrow(creative18)
coords <- cbind(Easting = as.numeric(as.character(creative18$"Feature Easting")),
                Northing = as.numeric(as.character(creative18$"Feature Northing")))
creative18sp <- SpatialPointsDataFrame(coords, data = data.frame(creative18$Name,
                                                                 creative18$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative18 <- over(uktowns, geometry(creative18sp), returnList = TRUE)
creative18 <- sapply(creative18, length)
creative18 <- spCbind(uktowns, creative18)
creative18 <- merge(creative18, popbua, by="BUA11CD")
creative18 <- creative18[which(creative18$poplsoa<175000&creative18$poplsoa>10000), ]
creative18cent <- SpatialPointsDataFrame(gCentroid(creative18, byid=TRUE), creative18@data, match.ID=FALSE)
creative18 <- creative18[-c(2,3,4,5,7)]


creative17 <- read_delim("educatjun17.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
creative17 <- creative17[complete.cases(creative17$`Feature Easting`), ]
creative17 <- creative17[complete.cases(creative17$`Feature Northing`), ]
names(creative17)[names(creative17)=="PointX Classification Code"] <- "code"
creative17 <- creative17[which(creative17$code==05320384 | creative17$code==05320394| creative17$code==05320395| creative17$code==05320396| creative17$code==05320389), ]
creative17 <- subset(creative17, "Feature Easting" != "" | "Feature Northing" != "")
creative17$poi_ID <- 1:nrow(creative17)
coords <- cbind(Easting = as.numeric(as.character(creative17$"Feature Easting")),
                Northing = as.numeric(as.character(creative17$"Feature Northing")))
creative17sp <- SpatialPointsDataFrame(coords, data = data.frame(creative17$Name,
                                                                 creative17$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative17 <- over(uktowns, geometry(creative17sp), returnList = TRUE)
creative17 <- sapply(creative17, length)
creative17 <- spCbind(uktowns, creative17)
creative17 <- merge(creative17, popbua, by="BUA11CD")
creative17 <- creative17[which(creative17$poplsoa<175000&creative17$poplsoa>10000), ]
creative17cent <- SpatialPointsDataFrame(gCentroid(creative17, byid=TRUE), creative17@data, match.ID=FALSE)
creative17 <- creative17[-c(2,3,4,5,7)]


creative16 <- read_delim("educatjun16.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
creative16 <- creative16[complete.cases(creative16$`Feature Easting`), ]
creative16 <- creative16[complete.cases(creative16$`Feature Northing`), ]
names(creative16)[names(creative16)=="PointX Classification Code"] <- "code"
creative16 <- creative16[which(creative16$code==05320384 | creative16$code==05320394| creative16$code==05320395| creative16$code==05320396| creative16$code==05320389), ]
creative16 <- subset(creative16, "Feature Easting" != "" | "Feature Northing" != "")
creative16$poi_ID <- 1:nrow(creative16)
coords <- cbind(Easting = as.numeric(as.character(creative16$"Feature Easting")),
                Northing = as.numeric(as.character(creative16$"Feature Northing")))
creative16sp <- SpatialPointsDataFrame(coords, data = data.frame(creative16$Name,
                                                                 creative16$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative16 <- over(uktowns, geometry(creative16sp), returnList = TRUE)
creative16 <- sapply(creative16, length)
creative16 <- spCbind(uktowns, creative16)
creative16 <- merge(creative16, popbua, by="BUA11CD")
creative16 <- creative16[which(creative16$poplsoa<175000&creative16$poplsoa>10000), ]
creative16cent <- SpatialPointsDataFrame(gCentroid(creative16, byid=TRUE), creative16@data, match.ID=FALSE)
creative16 <- creative16[-c(2,3,4,5,7)]

creative14 <- read_delim("educatsept14.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
creative14 <- creative14[complete.cases(creative14$`Feature Easting`), ]
creative14 <- creative14[complete.cases(creative14$`Feature Northing`), ]
names(creative14)[names(creative14)=="PointX Classification Code"] <- "code"
creative14 <- creative14[which(creative14$code==05320384 | creative14$code==05320394| creative14$code==05320395| creative14$code==05320396| creative14$code==05320389), ]
creative14 <- subset(creative14, "Feature Easting" != "" | "Feature Northing" != "")
creative14$poi_ID <- 1:nrow(creative14)
coords <- cbind(Easting = as.numeric(as.character(creative14$"Feature Easting")),
                Northing = as.numeric(as.character(creative14$"Feature Northing")))
creative14sp <- SpatialPointsDataFrame(coords, data = data.frame(creative14$Name,
                                                                 creative14$poi_ID), proj4string = CRS("+init=epsg:27700"))


creative14 <- over(uktowns, geometry(creative14sp), returnList = TRUE)
creative14 <- sapply(creative14, length)
creative14 <- spCbind(uktowns, creative14)
creative14 <- merge(creative14, popbua, by="BUA11CD")
creative14 <- creative14[which(creative14$poplsoa<175000&creative14$poplsoa>10000), ]
creative14cent <- SpatialPointsDataFrame(gCentroid(creative14, byid=TRUE), creative14@data, match.ID=FALSE)
creative14 <- creative14[-c(2,3,4,5,7)]



creative <- merge(creative14, creative15, by="BUA11CD")
creative <- merge(creative, creative16, by="BUA11CD")
creative <- merge(creative, creative17, by="BUA11CD")
creative <- merge(creative, creative18, by="BUA11CD")



write.csv(creative, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\creative148.csv")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])


night15 <- read_delim("discojun15.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
night15 <- night15[complete.cases(night15$`Feature Easting`), ]
night15 <- night15[complete.cases(night15$`Feature Northing`), ]
names(night15)[names(night15)=="PointX Classification Code"] <- "code"
night15 <- night15[which(night15$code==04250312 | night15$code==04250311), ]
night15 <- subset(night15, "Feature Easting" != "" | "Feature Northing" != "")
night15$poi_ID <- 1:nrow(night15)
coords <- cbind(Easting = as.numeric(as.character(night15$"Feature Easting")),
                Northing = as.numeric(as.character(night15$"Feature Northing")))
night15sp <- SpatialPointsDataFrame(coords, data = data.frame(night15$Name,
                                                                 night15$poi_ID), proj4string = CRS("+init=epsg:27700"))


night15 <- over(uktowns, geometry(night15sp), returnList = TRUE)
night15 <- sapply(night15, length)
night15 <- spCbind(uktowns, night15)
night15 <- merge(night15, popbua, by="BUA11CD")
night15 <- night15[which(night15$poplsoa<175000&night15$poplsoa>10000), ]
night15cent <- SpatialPointsDataFrame(gCentroid(night15, byid=TRUE), night15@data, match.ID=FALSE)
night15 <- night15[-c(2,3,4,5,7)]


night18 <- read_delim("discojun18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
night18 <- night18[complete.cases(night18$`Feature Easting`), ]
night18 <- night18[complete.cases(night18$`Feature Northing`), ]
names(night18)[names(night18)=="PointX Classification Code"] <- "code"
night18 <- night18[which(night18$code==04250312 | night18$code==04250311), ]
night18 <- subset(night18, "Feature Easting" != "" | "Feature Northing" != "")
night18$poi_ID <- 1:nrow(night18)
coords <- cbind(Easting = as.numeric(as.character(night18$"Feature Easting")),
                Northing = as.numeric(as.character(night18$"Feature Northing")))
night18sp <- SpatialPointsDataFrame(coords, data = data.frame(night18$Name,
                                                                 night18$poi_ID), proj4string = CRS("+init=epsg:27700"))


night18 <- over(uktowns, geometry(night18sp), returnList = TRUE)
night18 <- sapply(night18, length)
night18 <- spCbind(uktowns, night18)
night18 <- merge(night18, popbua, by="BUA11CD")
night18 <- night18[which(night18$poplsoa<175000&night18$poplsoa>10000), ]
night18cent <- SpatialPointsDataFrame(gCentroid(night18, byid=TRUE), night18@data, match.ID=FALSE)
night18 <- night18[-c(2,3,4,5,7)]


night17 <- read_delim("discojun17.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
night17 <- night17[complete.cases(night17$`Feature Easting`), ]
night17 <- night17[complete.cases(night17$`Feature Northing`), ]
names(night17)[names(night17)=="PointX Classification Code"] <- "code"
night17 <- night17[which(night17$code==04250312 | night17$code==04250311), ]
night17 <- subset(night17, "Feature Easting" != "" | "Feature Northing" != "")
night17$poi_ID <- 1:nrow(night17)
coords <- cbind(Easting = as.numeric(as.character(night17$"Feature Easting")),
                Northing = as.numeric(as.character(night17$"Feature Northing")))
night17sp <- SpatialPointsDataFrame(coords, data = data.frame(night17$Name,
                                                                 night17$poi_ID), proj4string = CRS("+init=epsg:27700"))


night17 <- over(uktowns, geometry(night17sp), returnList = TRUE)
night17 <- sapply(night17, length)
night17 <- spCbind(uktowns, night17)
night17 <- merge(night17, popbua, by="BUA11CD")
night17 <- night17[which(night17$poplsoa<175000&night17$poplsoa>10000), ]
night17cent <- SpatialPointsDataFrame(gCentroid(night17, byid=TRUE), night17@data, match.ID=FALSE)
night17 <- night17[-c(2,3,4,5,7)]


night16 <- read_delim("discojun16.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
night16 <- night16[complete.cases(night16$`Feature Easting`), ]
night16 <- night16[complete.cases(night16$`Feature Northing`), ]
names(night16)[names(night16)=="PointX Classification Code"] <- "code"
night16 <- night16[which(night16$code==04250312 | night16$code==04250311), ]
night16 <- subset(night16, "Feature Easting" != "" | "Feature Northing" != "")
night16$poi_ID <- 1:nrow(night16)
coords <- cbind(Easting = as.numeric(as.character(night16$"Feature Easting")),
                Northing = as.numeric(as.character(night16$"Feature Northing")))
night16sp <- SpatialPointsDataFrame(coords, data = data.frame(night16$Name,
                                                                 night16$poi_ID), proj4string = CRS("+init=epsg:27700"))


night16 <- over(uktowns, geometry(night16sp), returnList = TRUE)
night16 <- sapply(night16, length)
night16 <- spCbind(uktowns, night16)
night16 <- merge(night16, popbua, by="BUA11CD")
night16 <- night16[which(night16$poplsoa<175000&night16$poplsoa>10000), ]
night16cent <- SpatialPointsDataFrame(gCentroid(night16, byid=TRUE), night16@data, match.ID=FALSE)
night16 <- night16[-c(2,3,4,5,7)]

night14 <- read_delim("discosept14.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
night14 <- night14[complete.cases(night14$`Feature Easting`), ]
night14 <- night14[complete.cases(night14$`Feature Northing`), ]
names(night14)[names(night14)=="PointX Classification Code"] <- "code"
night14 <- night14[which(night14$code==04250312 | night14$code==04250311), ]
night14 <- subset(night14, "Feature Easting" != "" | "Feature Northing" != "")
night14$poi_ID <- 1:nrow(night14)
coords <- cbind(Easting = as.numeric(as.character(night14$"Feature Easting")),
                Northing = as.numeric(as.character(night14$"Feature Northing")))
night14sp <- SpatialPointsDataFrame(coords, data = data.frame(night14$Name,
                                                                 night14$poi_ID), proj4string = CRS("+init=epsg:27700"))


night14 <- over(uktowns, geometry(night14sp), returnList = TRUE)
night14 <- sapply(night14, length)
night14 <- spCbind(uktowns, night14)
night14 <- merge(night14, popbua, by="BUA11CD")
night14 <- night14[which(night14$poplsoa<175000&night14$poplsoa>10000), ]
night14cent <- SpatialPointsDataFrame(gCentroid(night14, byid=TRUE), night14@data, match.ID=FALSE)
night14 <- night14[-c(2,3,4,5,7)]


night <- merge(night14, night15, by="BUA11CD")
night <- merge(night, night16, by="BUA11CD")
night <- merge(night, night17, by="BUA11CD")
night <- merge(night, night18, by="BUA11CD")





write.csv(night, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\night148.csv")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])

cine15 <- read_delim("cinejun15.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
cine15 <- cine15[complete.cases(cine15$`Feature Easting`), ]
cine15 <- cine15[complete.cases(cine15$`Feature Northing`), ]
names(cine15)[names(cine15)=="PointX Classification Code"] <- "code"
cine15 <- cine15[which(cine15$code==04250308 | cine15$code==04250315), ]
cine15 <- subset(cine15, "Feature Easting" != "" | "Feature Northing" != "")
cine15$poi_ID <- 1:nrow(cine15)
coords <- cbind(Easting = as.numeric(as.character(cine15$"Feature Easting")),
                Northing = as.numeric(as.character(cine15$"Feature Northing")))
cine15sp <- SpatialPointsDataFrame(coords, data = data.frame(cine15$Name,
                                                                 cine15$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine15 <- over(uktowns, geometry(cine15sp), returnList = TRUE)
cine15 <- sapply(cine15, length)
cine15 <- spCbind(uktowns, cine15)
cine15 <- merge(cine15, popbua, by="BUA11CD")
cine15 <- cine15[which(cine15$poplsoa<175000&cine15$poplsoa>10000), ]
cine15cent <- SpatialPointsDataFrame(gCentroid(cine15, byid=TRUE), cine15@data, match.ID=FALSE)
cine15 <- cine15[-c(2,3,4,5,7)]


cine18 <- read_delim("cinejun18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
cine18 <- cine18[complete.cases(cine18$`Feature Easting`), ]
cine18 <- cine18[complete.cases(cine18$`Feature Northing`), ]
names(cine18)[names(cine18)=="PointX Classification Code"] <- "code"
cine18 <- cine18[which(cine18$code==04250308 | cine18$code==04250315), ]
cine18 <- subset(cine18, "Feature Easting" != "" | "Feature Northing" != "")
cine18$poi_ID <- 1:nrow(cine18)
coords <- cbind(Easting = as.numeric(as.character(cine18$"Feature Easting")),
                Northing = as.numeric(as.character(cine18$"Feature Northing")))
cine18sp <- SpatialPointsDataFrame(coords, data = data.frame(cine18$Name,
                                                                 cine18$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine18 <- over(uktowns, geometry(cine18sp), returnList = TRUE)
cine18 <- sapply(cine18, length)
cine18 <- spCbind(uktowns, cine18)
cine18 <- merge(cine18, popbua, by="BUA11CD")
cine18 <- cine18[which(cine18$poplsoa<175000&cine18$poplsoa>10000), ]
cine18cent <- SpatialPointsDataFrame(gCentroid(cine18, byid=TRUE), cine18@data, match.ID=FALSE)
cine18 <- cine18[-c(2,3,4,5,7)]


cine17 <- read_delim("cinejun17.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
cine17 <- cine17[complete.cases(cine17$`Feature Easting`), ]
cine17 <- cine17[complete.cases(cine17$`Feature Northing`), ]
names(cine17)[names(cine17)=="PointX Classification Code"] <- "code"
cine17 <- cine17[which(cine17$code==04250308 | cine17$code==04250315), ]
cine17 <- subset(cine17, "Feature Easting" != "" | "Feature Northing" != "")
cine17$poi_ID <- 1:nrow(cine17)
coords <- cbind(Easting = as.numeric(as.character(cine17$"Feature Easting")),
                Northing = as.numeric(as.character(cine17$"Feature Northing")))
cine17sp <- SpatialPointsDataFrame(coords, data = data.frame(cine17$Name,
                                                                 cine17$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine17 <- over(uktowns, geometry(cine17sp), returnList = TRUE)
cine17 <- sapply(cine17, length)
cine17 <- spCbind(uktowns, cine17)
cine17 <- merge(cine17, popbua, by="BUA11CD")
cine17 <- cine17[which(cine17$poplsoa<175000&cine17$poplsoa>10000), ]
cine17cent <- SpatialPointsDataFrame(gCentroid(cine17, byid=TRUE), cine17@data, match.ID=FALSE)
cine17 <- cine17[-c(2,3,4,5,7)]


cine16 <- read_delim("cinejun16.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
cine16 <- cine16[complete.cases(cine16$`Feature Easting`), ]
cine16 <- cine16[complete.cases(cine16$`Feature Northing`), ]
names(cine16)[names(cine16)=="PointX Classification Code"] <- "code"
cine16 <- cine16[which(cine16$code==04250308 | cine16$code==04250315), ]
cine16 <- subset(cine16, "Feature Easting" != "" | "Feature Northing" != "")
cine16$poi_ID <- 1:nrow(cine16)
coords <- cbind(Easting = as.numeric(as.character(cine16$"Feature Easting")),
                Northing = as.numeric(as.character(cine16$"Feature Northing")))
cine16sp <- SpatialPointsDataFrame(coords, data = data.frame(cine16$Name,
                                                                 cine16$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine16 <- over(uktowns, geometry(cine16sp), returnList = TRUE)
cine16 <- sapply(cine16, length)
cine16 <- spCbind(uktowns, cine16)
cine16 <- merge(cine16, popbua, by="BUA11CD")
cine16 <- cine16[which(cine16$poplsoa<175000&cine16$poplsoa>10000), ]
cine16cent <- SpatialPointsDataFrame(gCentroid(cine16, byid=TRUE), cine16@data, match.ID=FALSE)
cine16 <- cine16[-c(2,3,4,5,7)]

cine14 <- read_delim("cinesept14.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
cine14 <- cine14[complete.cases(cine14$`Feature Easting`), ]
cine14 <- cine14[complete.cases(cine14$`Feature Northing`), ]
names(cine14)[names(cine14)=="PointX Classification Code"] <- "code"
cine14 <- cine14[which(cine14$code==04250308 | cine14$code==04250315), ]
cine14 <- subset(cine14, "Feature Easting" != "" | "Feature Northing" != "")
cine14$poi_ID <- 1:nrow(cine14)
coords <- cbind(Easting = as.numeric(as.character(cine14$"Feature Easting")),
                Northing = as.numeric(as.character(cine14$"Feature Northing")))
cine14sp <- SpatialPointsDataFrame(coords, data = data.frame(cine14$Name,
                                                                 cine14$poi_ID), proj4string = CRS("+init=epsg:27700"))


cine14 <- over(uktowns, geometry(cine14sp), returnList = TRUE)
cine14 <- sapply(cine14, length)
cine14 <- spCbind(uktowns, cine14)
cine14 <- merge(cine14, popbua, by="BUA11CD")
cine14 <- cine14[which(cine14$poplsoa<175000&cine14$poplsoa>10000), ]
cine14cent <- SpatialPointsDataFrame(gCentroid(cine14, byid=TRUE), cine14@data, match.ID=FALSE)
cine14 <- cine14[-c(2,3,4,5,7)]


cine <- merge(cine14, cine15, by="BUA11CD")
cine <- merge(cine, cine16, by="BUA11CD")
cine <- merge(cine, cine17, by="BUA11CD")
cine <- merge(cine, cine18, by="BUA11CD")




write.csv(cine, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\cine148.csv")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])

sport15 <- read_delim("sportjun15.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
sport15 <- sport15[complete.cases(sport15$`Feature Easting`), ]
sport15 <- sport15[complete.cases(sport15$`Feature Northing`), ]
names(sport15)[names(sport15)=="PointX Classification Code"] <- "code"
sport15 <- sport15[which(sport15$code>=04240000 & sport15$code< 04250000), ]
sport15 <- subset(sport15, "Feature Easting" != "" | "Feature Northing" != "")
sport15$poi_ID <- 1:nrow(sport15)
coords <- cbind(Easting = as.numeric(as.character(sport15$"Feature Easting")),
                Northing = as.numeric(as.character(sport15$"Feature Northing")))
sport15sp <- SpatialPointsDataFrame(coords, data = data.frame(sport15$Name,
                                                                 sport15$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport15 <- over(uktowns, geometry(sport15sp), returnList = TRUE)
sport15 <- sapply(sport15, length)
sport15 <- spCbind(uktowns, sport15)
sport15 <- merge(sport15, popbua, by="BUA11CD")
sport15 <- sport15[which(sport15$poplsoa<175000&sport15$poplsoa>10000), ]
sport15cent <- SpatialPointsDataFrame(gCentroid(sport15, byid=TRUE), sport15@data, match.ID=FALSE)
sport15 <- sport15[-c(2,3,4,5,7)]


sport18 <- read_delim("sportjun18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
sport18 <- sport18[complete.cases(sport18$`Feature Easting`), ]
sport18 <- sport18[complete.cases(sport18$`Feature Northing`), ]
names(sport18)[names(sport18)=="PointX Classification Code"] <- "code"
sport18 <- sport18[which(sport18$code>=04240000 & sport18$code< 04250000), ]
sport18 <- subset(sport18, "Feature Easting" != "" | "Feature Northing" != "")
sport18$poi_ID <- 1:nrow(sport18)
coords <- cbind(Easting = as.numeric(as.character(sport18$"Feature Easting")),
                Northing = as.numeric(as.character(sport18$"Feature Northing")))
sport18sp <- SpatialPointsDataFrame(coords, data = data.frame(sport18$Name,
                                                                 sport18$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport18 <- over(uktowns, geometry(sport18sp), returnList = TRUE)
sport18 <- sapply(sport18, length)
sport18 <- spCbind(uktowns, sport18)
sport18 <- merge(sport18, popbua, by="BUA11CD")
sport18 <- sport18[which(sport18$poplsoa<175000&sport18$poplsoa>10000), ]
sport18cent <- SpatialPointsDataFrame(gCentroid(sport18, byid=TRUE), sport18@data, match.ID=FALSE)
sport18 <- sport18[-c(2,3,4,5,7)]


sport17 <- read_delim("sportjun17.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
sport17 <- sport17[complete.cases(sport17$`Feature Easting`), ]
sport17 <- sport17[complete.cases(sport17$`Feature Northing`), ]
names(sport17)[names(sport17)=="PointX Classification Code"] <- "code"
sport17 <- sport17[which(sport17$code>=04240000 & sport17$code< 04250000), ]
sport17 <- subset(sport17, "Feature Easting" != "" | "Feature Northing" != "")
sport17$poi_ID <- 1:nrow(sport17)
coords <- cbind(Easting = as.numeric(as.character(sport17$"Feature Easting")),
                Northing = as.numeric(as.character(sport17$"Feature Northing")))
sport17sp <- SpatialPointsDataFrame(coords, data = data.frame(sport17$Name,
                                                                 sport17$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport17 <- over(uktowns, geometry(sport17sp), returnList = TRUE)
sport17 <- sapply(sport17, length)
sport17 <- spCbind(uktowns, sport17)
sport17 <- merge(sport17, popbua, by="BUA11CD")
sport17 <- sport17[which(sport17$poplsoa<175000&sport17$poplsoa>10000), ]
sport17cent <- SpatialPointsDataFrame(gCentroid(sport17, byid=TRUE), sport17@data, match.ID=FALSE)
sport17 <- sport17[-c(2,3,4,5,7)]


sport16 <- read_delim("sportjun16.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
sport16 <- sport16[complete.cases(sport16$`Feature Easting`), ]
sport16 <- sport16[complete.cases(sport16$`Feature Northing`), ]
names(sport16)[names(sport16)=="PointX Classification Code"] <- "code"
sport16 <- sport16[which(sport16$code>=04240000 & sport16$code< 04250000), ]
sport16 <- subset(sport16, "Feature Easting" != "" | "Feature Northing" != "")
sport16$poi_ID <- 1:nrow(sport16)
coords <- cbind(Easting = as.numeric(as.character(sport16$"Feature Easting")),
                Northing = as.numeric(as.character(sport16$"Feature Northing")))
sport16sp <- SpatialPointsDataFrame(coords, data = data.frame(sport16$Name,
                                                                 sport16$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport16 <- over(uktowns, geometry(sport16sp), returnList = TRUE)
sport16 <- sapply(sport16, length)
sport16 <- spCbind(uktowns, sport16)
sport16 <- merge(sport16, popbua, by="BUA11CD")
sport16 <- sport16[which(sport16$poplsoa<175000&sport16$poplsoa>10000), ]
sport16cent <- SpatialPointsDataFrame(gCentroid(sport16, byid=TRUE), sport16@data, match.ID=FALSE)
sport16 <- sport16[-c(2,3,4,5,7)]

sport14 <- read_delim("sportsept14.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
sport14 <- sport14[complete.cases(sport14$`Feature Easting`), ]
sport14 <- sport14[complete.cases(sport14$`Feature Northing`), ]
names(sport14)[names(sport14)=="PointX Classification Code"] <- "code"
sport14 <- sport14[which(sport14$code>=04240000 & sport14$code< 04250000), ]
sport14 <- subset(sport14, "Feature Easting" != "" | "Feature Northing" != "")
sport14$poi_ID <- 1:nrow(sport14)
coords <- cbind(Easting = as.numeric(as.character(sport14$"Feature Easting")),
                Northing = as.numeric(as.character(sport14$"Feature Northing")))
sport14sp <- SpatialPointsDataFrame(coords, data = data.frame(sport14$Name,
                                                                 sport14$poi_ID), proj4string = CRS("+init=epsg:27700"))


sport14 <- over(uktowns, geometry(sport14sp), returnList = TRUE)
sport14 <- sapply(sport14, length)
sport14 <- spCbind(uktowns, sport14)
sport14 <- merge(sport14, popbua, by="BUA11CD")
sport14 <- sport14[which(sport14$poplsoa<175000&sport14$poplsoa>10000), ]
sport14cent <- SpatialPointsDataFrame(gCentroid(sport14, byid=TRUE), sport14@data, match.ID=FALSE)
sport14 <- sport14[-c(2,3,4,5,7)]


sport <- merge(sport14, sport15, by="BUA11CD")
sport <- merge(sport, sport16, by="BUA11CD")
sport <- merge(sport, sport17, by="BUA11CD")
sport <- merge(sport, sport18, by="BUA11CD")




write.csv(sport, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\sport148.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])

bet15 <- read_delim("gambjun15.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
bet15 <- bet15[complete.cases(bet15$`Feature Easting`), ]
bet15 <- bet15[complete.cases(bet15$`Feature Northing`), ]
names(bet15)[names(bet15)=="PointX Classification Code"] <- "code"
bet15 <- bet15[which(bet15$code==04220279), ]
bet15 <- subset(bet15, "Feature Easting" != "" | "Feature Northing" != "")
bet15$poi_ID <- 1:nrow(bet15)
coords <- cbind(Easting = as.numeric(as.character(bet15$"Feature Easting")),
                Northing = as.numeric(as.character(bet15$"Feature Northing")))
bet15sp <- SpatialPointsDataFrame(coords, data = data.frame(bet15$Name,
                                                                 bet15$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet15 <- over(uktowns, geometry(bet15sp), returnList = TRUE)
bet15 <- sapply(bet15, length)
bet15 <- spCbind(uktowns, bet15)
bet15 <- merge(bet15, popbua, by="BUA11CD")
bet15 <- bet15[which(bet15$poplsoa<175000&bet15$poplsoa>10000), ]
bet15cent <- SpatialPointsDataFrame(gCentroid(bet15, byid=TRUE), bet15@data, match.ID=FALSE)
bet15 <- bet15[-c(2,3,4,5,7)]


bet18 <- read_delim("gambjun18.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
bet18 <- bet18[complete.cases(bet18$`Feature Easting`), ]
bet18 <- bet18[complete.cases(bet18$`Feature Northing`), ]
names(bet18)[names(bet18)=="PointX Classification Code"] <- "code"
bet18 <- bet18[which(bet18$code==04220279), ]
bet18 <- subset(bet18, "Feature Easting" != "" | "Feature Northing" != "")
bet18$poi_ID <- 1:nrow(bet18)
coords <- cbind(Easting = as.numeric(as.character(bet18$"Feature Easting")),
                Northing = as.numeric(as.character(bet18$"Feature Northing")))
bet18sp <- SpatialPointsDataFrame(coords, data = data.frame(bet18$Name,
                                                                 bet18$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet18 <- over(uktowns, geometry(bet18sp), returnList = TRUE)
bet18 <- sapply(bet18, length)
bet18 <- spCbind(uktowns, bet18)
bet18 <- merge(bet18, popbua, by="BUA11CD")
bet18 <- bet18[which(bet18$poplsoa<175000&bet18$poplsoa>10000), ]
bet18cent <- SpatialPointsDataFrame(gCentroid(bet18, byid=TRUE), bet18@data, match.ID=FALSE)
bet18 <- bet18[-c(2,3,4,5,7)]


bet17 <- read_delim("gambjun17.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
bet17 <- bet17[complete.cases(bet17$`Feature Easting`), ]
bet17 <- bet17[complete.cases(bet17$`Feature Northing`), ]
names(bet17)[names(bet17)=="PointX Classification Code"] <- "code"
bet17 <- bet17[which(bet17$code==04220279), ]
bet17 <- subset(bet17, "Feature Easting" != "" | "Feature Northing" != "")
bet17$poi_ID <- 1:nrow(bet17)
coords <- cbind(Easting = as.numeric(as.character(bet17$"Feature Easting")),
                Northing = as.numeric(as.character(bet17$"Feature Northing")))
bet17sp <- SpatialPointsDataFrame(coords, data = data.frame(bet17$Name,
                                                                 bet17$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet17 <- over(uktowns, geometry(bet17sp), returnList = TRUE)
bet17 <- sapply(bet17, length)
bet17 <- spCbind(uktowns, bet17)
bet17 <- merge(bet17, popbua, by="BUA11CD")
bet17 <- bet17[which(bet17$poplsoa<175000&bet17$poplsoa>10000), ]
bet17cent <- SpatialPointsDataFrame(gCentroid(bet17, byid=TRUE), bet17@data, match.ID=FALSE)
bet17 <- bet17[-c(2,3,4,5,7)]


bet16 <- read_delim("gambjun16.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
bet16 <- bet16[complete.cases(bet16$`Feature Easting`), ]
bet16 <- bet16[complete.cases(bet16$`Feature Northing`), ]
names(bet16)[names(bet16)=="PointX Classification Code"] <- "code"
bet16 <- bet16[which(bet16$code==04220279), ]
bet16 <- subset(bet16, "Feature Easting" != "" | "Feature Northing" != "")
bet16$poi_ID <- 1:nrow(bet16)
coords <- cbind(Easting = as.numeric(as.character(bet16$"Feature Easting")),
                Northing = as.numeric(as.character(bet16$"Feature Northing")))
bet16sp <- SpatialPointsDataFrame(coords, data = data.frame(bet16$Name,
                                                                 bet16$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet16 <- over(uktowns, geometry(bet16sp), returnList = TRUE)
bet16 <- sapply(bet16, length)
bet16 <- spCbind(uktowns, bet16)
bet16 <- merge(bet16, popbua, by="BUA11CD")
bet16 <- bet16[which(bet16$poplsoa<175000&bet16$poplsoa>10000), ]
bet16cent <- SpatialPointsDataFrame(gCentroid(bet16, byid=TRUE), bet16@data, match.ID=FALSE)
bet16 <- bet16[-c(2,3,4,5,7)]

bet14 <- read_delim("gambsept14.csv",                            "|", escape_double = FALSE, trim_ws = TRUE,                          col_types = cols('PointX Classification Code' = col_number()))
bet14 <- bet14[complete.cases(bet14$`Feature Easting`), ]
bet14 <- bet14[complete.cases(bet14$`Feature Northing`), ]
names(bet14)[names(bet14)=="PointX Classification Code"] <- "code"
bet14 <- bet14[which(bet14$code==04220279), ]
bet14 <- subset(bet14, "Feature Easting" != "" | "Feature Northing" != "")
bet14$poi_ID <- 1:nrow(bet14)
coords <- cbind(Easting = as.numeric(as.character(bet14$"Feature Easting")),
                Northing = as.numeric(as.character(bet14$"Feature Northing")))
bet14sp <- SpatialPointsDataFrame(coords, data = data.frame(bet14$Name,
                                                                 bet14$poi_ID), proj4string = CRS("+init=epsg:27700"))


bet14 <- over(uktowns, geometry(bet14sp), returnList = TRUE)
bet14 <- sapply(bet14, length)
bet14 <- spCbind(uktowns, bet14)
bet14 <- merge(bet14, popbua, by="BUA11CD")
bet14 <- bet14[which(bet14$poplsoa<175000&bet14$poplsoa>10000), ]
bet14cent <- SpatialPointsDataFrame(gCentroid(bet14, byid=TRUE), bet14@data, match.ID=FALSE)
bet14 <- bet14[-c(2,3,4,5,7)]


bet <- merge(bet14, bet15, by="BUA11CD")
bet <- merge(bet, bet16, by="BUA11CD")
bet <- merge(bet, bet17, by="BUA11CD")
bet <- merge(bet, bet18, by="BUA11CD")





write.csv(bet, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\bet148.csv")




rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine" )])


bet <- merge(bet, popbua, by="BUA11CD")
creative <- merge(creative, popbua, by="BUA11CD")
sport <- merge(sport, popbua, by="BUA11CD")
night <- merge(night, popbua, by="BUA11CD")
cine <- merge(cine, popbua, by="BUA11CD")



bet$"2018-06-01" <- bet$bet18/(bet$poplsoa/1000)
bet$"2017-06-01" <- bet$bet17/(bet$poplsoa/1000)
bet$"2016-06-01" <- bet$bet16/(bet$poplsoa/1000)
bet$"2015-06-01" <- bet$bet15/(bet$poplsoa/1000)
bet$"2014-09-01" <- bet$bet14/(bet$poplsoa/1000)
bet$"2013-06-01" <- bet$bet13/(bet$poplsoa/1000)
bet$"2012-06-01" <- bet$bet12/(bet$poplsoa/1000)
bet$"2011-06-01" <- bet$bet11/(bet$poplsoa/1000)
bet$"2010-09-01" <- bet$bet10/(bet$poplsoa/1000)
bet$"2009-06-01" <- bet$bet09/(bet$poplsoa/1000)
bet$"2008-03-01" <- bet$bet08/(bet$poplsoa/1000)


creative$"2018-06-01" <- creative$creative18/(creative$poplsoa/1000)
creative$"2017-06-01" <- creative$creative17/(creative$poplsoa/1000)
creative$"2016-06-01" <- creative$creative16/(creative$poplsoa/1000)
creative$"2015-06-01" <- creative$creative15/(creative$poplsoa/1000)
creative$"2014-09-01" <- creative$creative14/(creative$poplsoa/1000)
creative$"2013-06-01" <- creative$creative13/(creative$poplsoa/1000)
creative$"2012-06-01" <- creative$creative12/(creative$poplsoa/1000)
creative$"2011-06-01" <- creative$creative11/(creative$poplsoa/1000)
creative$"2010-09-01" <- creative$creative10/(creative$poplsoa/1000)
creative$"2009-06-01" <- creative$creative09/(creative$poplsoa/1000)
creative$"2008-03-01" <- creative$creative08/(creative$poplsoa/1000)

sport$"2018-06-01" <- sport$sport18/(sport$poplsoa/1000)
sport$"2017-06-01" <- sport$sport17/(sport$poplsoa/1000)
sport$"2016-06-01" <- sport$sport16/(sport$poplsoa/1000)
sport$"2015-06-01" <- sport$sport15/(sport$poplsoa/1000)
sport$"2014-09-01" <- sport$sport14/(sport$poplsoa/1000)
sport$"2013-06-01" <- sport$sport13/(sport$poplsoa/1000)
sport$"2012-06-01" <- sport$sport12/(sport$poplsoa/1000)
sport$"2011-06-01" <- sport$sport11/(sport$poplsoa/1000)
sport$"2010-09-01" <- sport$sport10/(sport$poplsoa/1000)
sport$"2009-06-01" <- sport$sport09/(sport$poplsoa/1000)
sport$"2008-03-01" <- sport$sport08/(sport$poplsoa/1000)

night$"2018-06-01" <- night$night18/(night$poplsoa/1000)
night$"2017-06-01" <- night$night17/(night$poplsoa/1000)
night$"2016-06-01" <- night$night16/(night$poplsoa/1000)
night$"2015-06-01" <- night$night15/(night$poplsoa/1000)
night$"2014-09-01" <- night$night14/(night$poplsoa/1000)
night$"2013-06-01" <- night$night13/(night$poplsoa/1000)
night$"2012-06-01" <- night$night12/(night$poplsoa/1000)
night$"2011-06-01" <- night$night11/(night$poplsoa/1000)
night$"2010-09-01" <- night$night10/(night$poplsoa/1000)
night$"2009-06-01" <- night$night09/(night$poplsoa/1000)
night$"2008-03-01" <- night$night08/(night$poplsoa/1000)

cine$"2018-06-01" <- cine$cine18/(cine$poplsoa/1000)
cine$"2017-06-01" <- cine$cine17/(cine$poplsoa/1000)
cine$"2016-06-01" <- cine$cine16/(cine$poplsoa/1000)
cine$"2015-06-01" <- cine$cine15/(cine$poplsoa/1000)
cine$"2014-09-01" <- cine$cine14/(cine$poplsoa/1000)
cine$"2013-06-01" <- cine$cine13/(cine$poplsoa/1000)
cine$"2012-06-01" <- cine$cine12/(cine$poplsoa/1000)
cine$"2011-06-01" <- cine$cine11/(cine$poplsoa/1000)
cine$"2010-09-01" <- cine$cine10/(cine$poplsoa/1000)
cine$"2009-06-01" <- cine$cine09/(cine$poplsoa/1000)
cine$"2008-03-01" <- cine$cine08/(cine$poplsoa/1000)

cines <- cine[-c(2,3,4,5,6,7,8,9,10,11,12)]
nights <- night[-c(2,3,4,5,6,7,8,9,10,11,12)]
bets <- bet[-c(2,3,4,5,6,7,8,9,10,11,12)]
sports <- sport[-c(2,3,4,5,6,7,8,9,10,11,12)]
creatives <- creative[-c(2,3,4,5,6,7,8,9,10,11,12)]

cines <- as.data.frame(cines)

cines <- reshape(cines, idvar = "BUA11CD", ids = cines$BUA11CD,
                 times = names(cines), timevar = "date",
                 varying = list(names(cines)),v.names="cines", new.row.names = 1:((dim(cines)[2])*(dim(cines)[1])),direction = "long")

cines$date <- as.Date(cines$date)

nights <- as.data.frame(nights)

nights <- reshape(nights, idvar = "BUA11CD", ids = nights$BUA11CD,
                 times = names(nights), timevar = "date",
                 varying = list(names(nights)),v.names="nights", new.row.names = 1:((dim(nights)[2])*(dim(nights)[1])),direction = "long")

nights$date <- as.Date(nights$date)

bets <- as.data.frame(bets)

bets <- reshape(bets, idvar = "BUA11CD", ids = bets$BUA11CD,
                 times = names(bets), timevar = "date",
                 varying = list(names(bets)),v.names="bets", new.row.names = 1:((dim(bets)[2])*(dim(bets)[1])),direction = "long")

bets$date <- as.Date(bets$date)

sports <- as.data.frame(sports)

sports <- reshape(sports, idvar = "BUA11CD", ids = sports$BUA11CD,
                 times = names(sports), timevar = "date",
                 varying = list(names(sports)),v.names="sports", new.row.names = 1:((dim(sports)[2])*(dim(sports)[1])),direction = "long")

sports$date <- as.Date(sports$date)

creatives <- as.data.frame(creatives)

creatives <- reshape(creatives, idvar = "BUA11CD", ids = creatives$BUA11CD,
                 times = names(creatives), timevar = "date",
                 varying = list(names(creatives)),v.names="creatives", new.row.names = 1:((dim(creatives)[2])*(dim(creatives)[1])),direction = "long")

creatives$date <- as.Date(creatives$date)


artsreal <- merge(bets, sports, by=c("BUA11CD", "date"))
artsreal <- merge(artsreal, nights, by=c("BUA11CD", "date"))
artsreal <- merge(artsreal, cines, by=c("BUA11CD", "date"))
artsreal <- merge(artsreal, creatives, by=c("BUA11CD", "date"))







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
station13 <- POI_2013[which(POI_2013$code==10570738), ]
write.csv(station13, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\whystationwhy.csv")


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

hmrc13 <- POI_2013[which(POI_2013$code==6330417), ]
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




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital13.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered13.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc13.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental13.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station13.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps13.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent13.csv")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" , "creatives", "nights", "artsreal", "sports", "cines", "bets" )])



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

hmrc12 <- POI_2012[which(POI_2012$code==6330417), ]
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


hospital <-  as.data.frame(hospitals12)
furthered <-  as.data.frame(furthered12)
hmrc <-  as.data.frame(hmrc12)
gps <-  as.data.frame(GPs12)
mental <-  as.data.frame(mental12)
station <-  as.data.frame(station12)
jobcent <-  as.data.frame(jobcent12)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital12.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered12.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc12.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental12.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station12.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps12.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent12.csv")




rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" , "creatives", "nights", "artsreal", "sports", "cines", "bets" )])





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

hmrc08 <- POI_2008[which(POI_2008$code==6330417), ]
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


hospital <-  as.data.frame(hospitals08)
furthered <-  as.data.frame(furthered08)
hmrc <-  as.data.frame(hmrc08)
gps <-  as.data.frame(GPs08)
mental <-  as.data.frame(mental08)
station <-  as.data.frame(station08)
jobcent <-  as.data.frame(jobcent08)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital08.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered08.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc08.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental08.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station08.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps08.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent08.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])




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

hmrc09 <- POI_2009[which(POI_2009$code==6330417), ]
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




hospital <-  as.data.frame(hospitals09)
furthered <-  as.data.frame(furthered09)
hmrc <-  as.data.frame(hmrc09)
gps <-  as.data.frame(GPs09)
mental <-  as.data.frame(mental09)
station <-  as.data.frame(station09)
jobcent <-  as.data.frame(jobcent09)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital09.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered09.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc09.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental09.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station09.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps09.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent09.csv")


rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])





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

hmrc10 <- POI_2010[which(POI_2010$code==6330417), ]
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




hospital <-  as.data.frame(hospitals10)
furthered <-  as.data.frame(furthered10)
hmrc <-  as.data.frame(hmrc10)
gps <-  as.data.frame(GPs10)
mental <-  as.data.frame(mental10)
station <-  as.data.frame(station10)
jobcent <-  as.data.frame(jobcent10)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital10.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered10.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc10.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental10.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station10.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps10.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent10.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])



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

hmrc11 <- POI_2011[which(POI_2011$code==6330417), ]
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




hospital <-  as.data.frame(hospitals11)
furthered <-  as.data.frame(furthered11)
hmrc <-  as.data.frame(hmrc11)
gps <-  as.data.frame(GPs11)
mental <-  as.data.frame(mental11)
station <-  as.data.frame(station11)
jobcent <-  as.data.frame(jobcent11)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital11.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered11.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc11.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental11.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station11.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps11.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent11.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])


###2014 onwards###


health_2014 <- read_delim("health_2014.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2014 <- health_2014 [complete.cases(health_2014$`Feature Easting`), ]
health_2014 <- health_2014 [complete.cases(health_2014$`Feature Northing`), ]
names(health_2014)[names(health_2014)=="PointX Classification Code"] <- "code"

schools_2014 <- read_delim("schools_2014.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2014 <- schools_2014 [complete.cases(schools_2014$`Feature Easting`), ]
schools_2014 <- schools_2014 [complete.cases(schools_2014$`Feature Northing`), ]
names(schools_2014)[names(schools_2014)=="PointX Classification Code"] <- "code"

gov_2014 <- read_delim("gov_2014.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Easting`), ]
gov_2014 <- gov_2014 [complete.cases(gov_2014$`Feature Northing`), ]
names(gov_2014)[names(gov_2014)=="PointX Classification Code"] <- "code"

station_2014 <- read_delim("station_2014.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
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
mental14 <- merge(mental14, popbua, by="BUA11CD")
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
station14 <- merge(station14, popbua, by="BUA11CD")
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
jobcent14 <- merge(jobcent14, popbua, by="BUA11CD")
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
furthered14 <- merge(furthered14, popbua, by="BUA11CD")
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
GPs14 <- merge(GPs14, popbua, by="BUA11CD")
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
hospitals14 <- merge(hospitals14, popbua, by="BUA11CD")
hospitals14 <- hospitals14[which(hospitals14$poplsoa<175000&hospitals14$poplsoa>10000), ]
hospitals14cent <- SpatialPointsDataFrame(gCentroid(hospitals14, byid=TRUE), hospitals14@data, match.ID=FALSE)
hospitals14 <- hospitals14[-c(2,3,4,5,7)]

hmrc14 <- gov_2014[which(gov_2014$code==06330417), ]
hmrc14 <- subset(hmrc14, "Feature Easting" != "" | "Feature Northing" != "")
hmrc14$poi_ID <- 1:nrow(hmrc14)
coords <- cbind(Easting = as.numeric(as.character(hmrc14$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc14$`Feature Northing`)))
hmrc14sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc14$Name,
                                                             hmrc14$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc14 <- over(uktowns, geometry(hmrc14sp), returnList = TRUE)
hmrc14 <- sapply(hmrc14, length)
hmrc14 <- spCbind(uktowns, hmrc14)
hmrc14 <- merge(hmrc14, popbua, by="BUA11CD")
hmrc14 <- hmrc14[which(hmrc14$poplsoa<175000&hmrc14$poplsoa>10000), ]
hmrc14cent <- SpatialPointsDataFrame(gCentroid(hmrc14, byid=TRUE), hmrc14@data, match.ID=FALSE)
hmrc14 <- hmrc14[-c(2,3,4,5,7)]




hospital <-  as.data.frame(hospitals14)
furthered <-  as.data.frame(furthered14)
hmrc <-  as.data.frame(hmrc14)
gps <-  as.data.frame(GPs14)
mental <-  as.data.frame(mental14)
station <-  as.data.frame(station14)
jobcent <-  as.data.frame(jobcent14)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital14.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered14.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc14.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental14.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station14.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps14.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent14.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])


health_2015 <- read_delim("health_2015.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2015 <- health_2015 [complete.cases(health_2015$`Feature Easting`), ]
health_2015 <- health_2015 [complete.cases(health_2015$`Feature Northing`), ]
names(health_2015)[names(health_2015)=="PointX Classification Code"] <- "code"

schools_2015 <- read_delim("schools_2015.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2015 <- schools_2015 [complete.cases(schools_2015$`Feature Easting`), ]
schools_2015 <- schools_2015 [complete.cases(schools_2015$`Feature Northing`), ]
names(schools_2015)[names(schools_2015)=="PointX Classification Code"] <- "code"

gov_2015 <- read_delim("gov_2015.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Easting`), ]
gov_2015 <- gov_2015 [complete.cases(gov_2015$`Feature Northing`), ]
names(gov_2015)[names(gov_2015)=="PointX Classification Code"] <- "code"

station_2015 <- read_delim("station_2015.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
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
mental15 <- merge(mental15, popbua, by="BUA11CD")
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
station15 <- merge(station15, popbua, by="BUA11CD")
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
jobcent15 <- merge(jobcent15, popbua, by="BUA11CD")
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
furthered15 <- merge(furthered15, popbua, by="BUA11CD")
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
GPs15 <- merge(GPs15, popbua, by="BUA11CD")
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
hospitals15 <- merge(hospitals15, popbua, by="BUA11CD")
hospitals15 <- hospitals15[which(hospitals15$poplsoa<175000&hospitals15$poplsoa>10000), ]
hospitals15cent <- SpatialPointsDataFrame(gCentroid(hospitals15, byid=TRUE), hospitals15@data, match.ID=FALSE)
hospitals15 <- hospitals15[-c(2,3,4,5,7)]

hmrc15 <- gov_2015[which(gov_2015$code==06330417), ]
hmrc15 <- subset(hmrc15, "Feature Easting" != "" | "Feature Northing" != "")
hmrc15$poi_ID <- 1:nrow(hmrc15)
coords <- cbind(Easting = as.numeric(as.character(hmrc15$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc15$`Feature Northing`)))
hmrc15sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc15$Name,
                                                             hmrc15$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc15 <- over(uktowns, geometry(hmrc15sp), returnList = TRUE)
hmrc15 <- sapply(hmrc15, length)
hmrc15 <- spCbind(uktowns, hmrc15)
hmrc15 <- merge(hmrc15, popbua, by="BUA11CD")
hmrc15 <- hmrc15[which(hmrc15$poplsoa<175000&hmrc15$poplsoa>10000), ]
hmrc15cent <- SpatialPointsDataFrame(gCentroid(hmrc15, byid=TRUE), hmrc15@data, match.ID=FALSE)
hmrc15 <- hmrc15[-c(2,3,4,5,7)]


hospital <-  as.data.frame(hospitals15)
furthered <-  as.data.frame(furthered15)
hmrc <-  as.data.frame(hmrc15)
gps <-  as.data.frame(GPs15)
mental <-  as.data.frame(mental15)
station <-  as.data.frame(station15)
jobcent <-  as.data.frame(jobcent15)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital15.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered15.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc15.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental15.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station15.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps15.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent15.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc" , "creatives", "nights", "artsreal", "sports", "cines", "bets" )])


health_2016 <- read_delim("health_2016.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2016 <- health_2016 [complete.cases(health_2016$`Feature Easting`), ]
health_2016 <- health_2016 [complete.cases(health_2016$`Feature Northing`), ]
names(health_2016)[names(health_2016)=="PointX Classification Code"] <- "code"

schools_2016 <- read_delim("schools_2016.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2016 <- schools_2016 [complete.cases(schools_2016$`Feature Easting`), ]
schools_2016 <- schools_2016 [complete.cases(schools_2016$`Feature Northing`), ]
names(schools_2016)[names(schools_2016)=="PointX Classification Code"] <- "code"

gov_2016 <- read_delim("gov_2016.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Easting`), ]
gov_2016 <- gov_2016 [complete.cases(gov_2016$`Feature Northing`), ]
names(gov_2016)[names(gov_2016)=="PointX Classification Code"] <- "code"

station_2016 <- read_delim("station_2016.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
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
mental16 <- merge(mental16, popbua, by="BUA11CD")
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
station16 <- merge(station16, popbua, by="BUA11CD")
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
jobcent16 <- merge(jobcent16, popbua, by="BUA11CD")
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
furthered16 <- merge(furthered16, popbua, by="BUA11CD")
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
GPs16 <- merge(GPs16, popbua, by="BUA11CD")
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
hospitals16 <- merge(hospitals16, popbua, by="BUA11CD")
hospitals16 <- hospitals16[which(hospitals16$poplsoa<175000&hospitals16$poplsoa>10000), ]
hospitals16cent <- SpatialPointsDataFrame(gCentroid(hospitals16, byid=TRUE), hospitals16@data, match.ID=FALSE)
hospitals16 <- hospitals16[-c(2,3,4,5,7)]

hmrc16 <- gov_2016[which(gov_2016$code==06330417), ]
hmrc16 <- subset(hmrc16, "Feature Easting" != "" | "Feature Northing" != "")
hmrc16$poi_ID <- 1:nrow(hmrc16)
coords <- cbind(Easting = as.numeric(as.character(hmrc16$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc16$`Feature Northing`)))
hmrc16sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc16$Name,
                                                             hmrc16$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc16 <- over(uktowns, geometry(hmrc16sp), returnList = TRUE)
hmrc16 <- sapply(hmrc16, length)
hmrc16 <- spCbind(uktowns, hmrc16)
hmrc16 <- merge(hmrc16, popbua, by="BUA11CD")
hmrc16 <- hmrc16[which(hmrc16$poplsoa<175000&hmrc16$poplsoa>10000), ]
hmrc16cent <- SpatialPointsDataFrame(gCentroid(hmrc16, byid=TRUE), hmrc16@data, match.ID=FALSE)
hmrc16 <- hmrc16[-c(2,3,4,5,7)]



hospital <-  as.data.frame(hospitals16)
furthered <-  as.data.frame(furthered16)
hmrc <-  as.data.frame(hmrc16)
gps <-  as.data.frame(GPs16)
mental <-  as.data.frame(mental16)
station <-  as.data.frame(station16)
jobcent <-  as.data.frame(jobcent16)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital16.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered16.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc16.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental16.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station16.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps16.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent16.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])


health_2017 <- read_delim("health_2017.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                        col_types = cols('PointX Classification Code' = col_number()))
health_2017 <- health_2017 [complete.cases(health_2017$`Feature Easting`), ]
health_2017 <- health_2017 [complete.cases(health_2017$`Feature Northing`), ]
names(health_2017)[names(health_2017)=="PointX Classification Code"] <- "code"

schools_2017 <- read_delim("schools_2017.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                         col_types = cols('PointX Classification Code' = col_number()))
schools_2017 <- schools_2017 [complete.cases(schools_2017$`Feature Easting`), ]
schools_2017 <- schools_2017 [complete.cases(schools_2017$`Feature Northing`), ]
names(schools_2017)[names(schools_2017)=="PointX Classification Code"] <- "code"

gov_2017 <- read_delim("gov_2017.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
                     col_types = cols('PointX Classification Code' = col_number()))
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Easting`), ]
gov_2017 <- gov_2017 [complete.cases(gov_2017$`Feature Northing`), ]
names(gov_2017)[names(gov_2017)=="PointX Classification Code"] <- "code"

station_2017 <- read_delim("station_2017.csv",                            "|", escape_double = FALSE, trim_ws = TRUE, 
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
mental17 <- merge(mental17, popbua, by="BUA11CD")
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
station17 <- merge(station17, popbua, by="BUA11CD")
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
jobcent17 <- merge(jobcent17, popbua, by="BUA11CD")
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
furthered17 <- merge(furthered17, popbua, by="BUA11CD")
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
GPs17 <- merge(GPs17, popbua, by="BUA11CD")
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
hospitals17 <- merge(hospitals17, popbua, by="BUA11CD")
hospitals17 <- hospitals17[which(hospitals17$poplsoa<175000&hospitals17$poplsoa>10000), ]
hospitals17cent <- SpatialPointsDataFrame(gCentroid(hospitals17, byid=TRUE), hospitals17@data, match.ID=FALSE)
hospitals17 <- hospitals17[-c(2,3,4,5,7)]

hmrc17 <- gov_2017[which(gov_2017$code==06330417), ]
hmrc17 <- subset(hmrc17, "Feature Easting" != "" | "Feature Northing" != "")
hmrc17$poi_ID <- 1:nrow(hmrc17)
coords <- cbind(Easting = as.numeric(as.character(hmrc17$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc17$`Feature Northing`)))
hmrc17sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc17$Name,
                                                             hmrc17$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc17 <- over(uktowns, geometry(hmrc17sp), returnList = TRUE)
hmrc17 <- sapply(hmrc17, length)
hmrc17 <- spCbind(uktowns, hmrc17)
hmrc17 <- merge(hmrc17, popbua, by="BUA11CD")
hmrc17 <- hmrc17[which(hmrc17$poplsoa<175000&hmrc17$poplsoa>10000), ]
hmrc17cent <- SpatialPointsDataFrame(gCentroid(hmrc17, byid=TRUE), hmrc17@data, match.ID=FALSE)
hmrc17 <- hmrc17[-c(2,3,4,5,7)]



hospital <-  as.data.frame(hospitals17)
furthered <-  as.data.frame(furthered17)
hmrc <-  as.data.frame(hmrc17)
gps <-  as.data.frame(GPs17)
mental <-  as.data.frame(mental17)
station <-  as.data.frame(station17)
jobcent <-  as.data.frame(jobcent17)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital17.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered17.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc17.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental17.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station17.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps17.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent17.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])


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


mental18 <- over(uktowns, geometry(mental18sp), returnList = TRUE)
mental18 <- sapply(mental18, length)
mental18 <- spCbind(uktowns, mental18)
mental18 <- merge(mental18, popbua, by="BUA11CD")
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
station18 <- merge(station18, popbua, by="BUA11CD")
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
jobcent18 <- merge(jobcent18, popbua, by="BUA11CD")
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
furthered18 <- merge(furthered18, popbua, by="BUA11CD")
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
GPs18 <- merge(GPs18, popbua, by="BUA11CD")
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
hospitals18 <- merge(hospitals18, popbua, by="BUA11CD")
hospitals18 <- hospitals18[which(hospitals18$poplsoa<175000&hospitals18$poplsoa>10000), ]
hospitals18cent <- SpatialPointsDataFrame(gCentroid(hospitals18, byid=TRUE), hospitals18@data, match.ID=FALSE)
hospitals18 <- hospitals18[-c(2,3,4,5,7)]

hmrc18 <- gov_2018[which(gov_2018$code==06330417), ]
hmrc18 <- subset(hmrc18, "Feature Easting" != "" | "Feature Northing" != "")
hmrc18$poi_ID <- 1:nrow(hmrc18)
coords <- cbind(Easting = as.numeric(as.character(hmrc18$`Feature Easting`)),
                Northing = as.numeric(as.character(hmrc18$`Feature Northing`)))
hmrc18sp <- SpatialPointsDataFrame(coords, data = data.frame(hmrc18$Name,
                                                             hmrc18$poi_ID), proj4string = CRS("+init=epsg:27700"))


hmrc18 <- over(uktowns, geometry(hmrc18sp), returnList = TRUE)
hmrc18 <- sapply(hmrc18, length)
hmrc18 <- spCbind(uktowns, hmrc18)
hmrc18 <- merge(hmrc18, popbua, by="BUA11CD")
hmrc18 <- hmrc18[which(hmrc18$poplsoa<175000&hmrc18$poplsoa>10000), ]
hmrc18cent <- SpatialPointsDataFrame(gCentroid(hmrc18, byid=TRUE), hmrc18@data, match.ID=FALSE)
hmrc18 <- hmrc18[-c(2,3,4,5,7)]




hospital <-  as.data.frame(hospitals18)
furthered <-  as.data.frame(furthered18)
hmrc <-  as.data.frame(hmrc18)
gps <-  as.data.frame(GPs18)
mental <-  as.data.frame(mental18)
station <-  as.data.frame(station18)
jobcent <-  as.data.frame(jobcent18)




write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hospital18.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\furthered18.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\hmrc18.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\mental18.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\station18.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\gps18.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\jobcent18.csv")

rm(list=ls()[! ls() %in% c("buas", "buanames", "citybua", "popbua", "scotbua", "scotpop", "ukcrs", "uktowns", "leisure", "bet", "creative", "night", "sport", "cine", "hospital", "mental", "station", "gps", "furthered", "jobcent", "hmrc", "creatives", "nights", "artsreal", "sports", "cines", "bets"  )])

hospital <- merge(hospital, popbua, by="BUA11CD")
furthered <- merge(furthered, popbua, by="BUA11CD")
hmrc <- merge(hmrc, popbua, by="BUA11CD")
gps <- merge(gps, popbua, by="BUA11CD")
mental <- merge(mental, popbua, by="BUA11CD")
station <- merge(station, popbua, by="BUA11CD")
jobcent <- merge(jobcent, popbua, by="BUA11CD")




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

furthereds <- furthered[-c(2,3,4,5,6,7,8,9,10,11,12)]
stations <- station[-c(2,3,4,5,6,7,8,9,10,11,12)]
gpss <- gps[-c(2,3,4,5,6,7,8,9,10,11,12)]
hospitals <- hospital[-c(2,3,4,5,6,7,8,9,10,11,12)]
mentals <- mental[-c(2,3,4,5,6,7,8,9,10,11,12)]
jobcents <- jobcent[-c(2,3,4,5,6,7,8,9,10,11,12)]
hmrcs <- hmrc[-c(2,3,4,5,6,7,8,9,10,11,12)]


furthereds <- as.data.frame(furthereds)

furthereds <- reshape(furthereds, idvar = "BUA11CD", ids = furthereds$BUA11CD,
                        times = names(furthereds), timevar = "date",
                        varying = list(names(furthereds)),v.names="furthereds", new.row.names = 1:((dim(furthereds)[2])*(dim(furthereds)[1])),direction = "long")

furthereds$date <- as.Date(furthereds$date)

stations <- as.data.frame(stations)

stations <- reshape(stations, idvar = "BUA11CD", ids = stations$BUA11CD,
                      times = names(stations), timevar = "date",
                      varying = list(names(stations)),v.names="stations", new.row.names = 1:((dim(stations)[2])*(dim(stations)[1])),direction = "long")

stations$date <- as.Date(stations$date)

gpss <- as.data.frame(gpss)

gpss <- reshape(gpss, idvar = "BUA11CD", ids = gpss$BUA11CD,
                    times = names(gpss), timevar = "date",
                    varying = list(names(gpss)),v.names="gpss", new.row.names = 1:((dim(gpss)[2])*(dim(gpss)[1])),direction = "long")

gpss$date <- as.Date(gpss$date)

hospitals <- as.data.frame(hospitals)

hospitals <- reshape(hospitals, idvar = "BUA11CD", ids = hospitals$BUA11CD,
                times = names(hospitals), timevar = "date",
                varying = list(names(hospitals)),v.names="hospitals", new.row.names = 1:((dim(hospitals)[2])*(dim(hospitals)[1])),direction = "long")

hospitals$date <- as.Date(hospitals$date)

mentals <- as.data.frame(mentals)

mentals <- reshape(mentals, idvar = "BUA11CD", ids = mentals$BUA11CD,
                     times = names(mentals), timevar = "date",
                     varying = list(names(mentals)),v.names="mentals", new.row.names = 1:((dim(mentals)[2])*(dim(mentals)[1])),direction = "long")

mentals$date <- as.Date(mentals$date)

jobcents <- as.data.frame(jobcents)

jobcents <- reshape(jobcents, idvar = "BUA11CD", ids = jobcents$BUA11CD,
                     times = names(jobcents), timevar = "date",
                     varying = list(names(jobcents)),v.names="jobcents", new.row.names = 1:((dim(jobcents)[2])*(dim(jobcents)[1])),direction = "long")

jobcents$date <- as.Date(jobcents$date)

hmrcs <- as.data.frame(hmrcs)

hmrcs <- reshape(hmrcs, idvar = "BUA11CD", ids = hmrcs$BUA11CD,
                     times = names(hmrcs), timevar = "date",
                     varying = list(names(hmrcs)),v.names="hmrcs", new.row.names = 1:((dim(hmrcs)[2])*(dim(hmrcs)[1])),direction = "long")

hmrcs$date <- as.Date(hmrcs$date)




publicbinary <- merge(furthereds, stations, by=c("BUA11CD", "date"))
publicbinary <- merge(publicbinary, gpss, by=c("BUA11CD", "date"))
publicbinary <- merge(publicbinary, hospitals, by=c("BUA11CD", "date"))
publicbinary <- merge(publicbinary, mentals, by=c("BUA11CD", "date"))
publicbinary <- merge(publicbinary, hmrcs, by=c("BUA11CD", "date"))
publicbinary <- merge(publicbinary, jobcents, by=c("BUA11CD", "date"))

write.csv(publicbinary, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\apublicbinary.csv")
write.csv(hmrc, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\ahmrc.csv")
write.csv(jobcent, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\ajobcent.csv")
write.csv(hospital, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\ahospital.csv")
write.csv(mental, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\amental.csv")
write.csv(gps, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\agps.csv")
write.csv(station, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\astation.csv")
write.csv(furthered, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\afurthered.csv")
write.csv(hmrcs, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\ahmrcs.csv")
write.csv(jobcents, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\ajobcents.csv")
write.csv(hospitals, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\ahospitals.csv")
write.csv(mentals, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\amentals.csv")
write.csv(gpss, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\agpss.csv")
write.csv(stations, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\astations.csv")
write.csv(furthereds, "C:\\Users\\bjg55\\Documents\\Towns\\2019Ring\\afurthereds.csv")






