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


###master sets###

artspanel <- read_csv("artsmasterpanel.csv")
artscross <- read_csv("artsmastercross.csv")
OSmaster <- read_csv("osmaster.csv")


betlong <- read_csv("betlong.csv")
cinelong <- read_csv("cinelong.csv")
nightlong <- read_csv("nightlong.csv")
sportlong <- read_csv("sportlong.csv")
creativelong <- read_csv("creativelong.csv")

betlong <- betlong[-c(1)]
cinelong <- cinelong[-c(1)]
nightlong <- nightlong[-c(1)]
sportlong <- sportlong[-c(1)]
creativelong <- creativelong[-c(1)]

betlong <- merge(betlong, popbua, by="BUA11CD")
sportlong <- merge(sportlong, popbua, by="BUA11CD")
creativelong <- merge(creativelong, popbua, by="BUA11CD")
nightlong <- merge(nightlong, popbua, by="BUA11CD")
cinelong <- merge(cinelong, popbua, by="BUA11CD")

betlong$"2018-06-01" <- betlong$bet18/(betlong$poplsoa/1000)
betlong$"2017-06-01" <- betlong$bet17/(betlong$poplsoa/1000)
betlong$"2016-06-01" <- betlong$bet16/(betlong$poplsoa/1000)
betlong$"2015-06-01" <- betlong$bet15/(betlong$poplsoa/1000)
betlong$"2014-09-01" <- betlong$bet14/(betlong$poplsoa/1000)
betlong$"2013-06-01" <- betlong$bet13/(betlong$poplsoa/1000)
betlong$"2012-06-01" <- betlong$bet12/(betlong$poplsoa/1000)
betlong$"2011-06-01" <- betlong$bet11/(betlong$poplsoa/1000)
betlong$"2010-09-01" <- betlong$bet10/(betlong$poplsoa/1000)
betlong$"2009-06-01" <- betlong$bet09/(betlong$poplsoa/1000)
betlong$"2008-03-01" <- betlong$bet08/(betlong$poplsoa/1000)

creativelong$"2018-06-01" <- creativelong$creative18/(creativelong$poplsoa/1000)
creativelong$"2017-06-01" <- creativelong$creative17/(creativelong$poplsoa/1000)
creativelong$"2016-06-01" <- creativelong$creative16/(creativelong$poplsoa/1000)
creativelong$"2015-06-01" <- creativelong$creative15/(creativelong$poplsoa/1000)
creativelong$"2014-09-01" <- creativelong$creative14/(creativelong$poplsoa/1000)
creativelong$"2013-06-01" <- creativelong$creative13/(creativelong$poplsoa/1000)
creativelong$"2012-06-01" <- creativelong$creative12/(creativelong$poplsoa/1000)
creativelong$"2011-06-01" <- creativelong$creative11/(creativelong$poplsoa/1000)
creativelong$"2010-09-01" <- creativelong$creative10/(creativelong$poplsoa/1000)
creativelong$"2009-06-01" <- creativelong$creative09/(creativelong$poplsoa/1000)
creativelong$"2008-03-01" <- creativelong$creative08/(creativelong$poplsoa/1000)

sportlong$"2018-06-01" <- sportlong$sport18/(sportlong$poplsoa/1000)
sportlong$"2017-06-01" <- sportlong$sport17/(sportlong$poplsoa/1000)
sportlong$"2016-06-01" <- sportlong$sport16/(sportlong$poplsoa/1000)
sportlong$"2015-06-01" <- sportlong$sport15/(sportlong$poplsoa/1000)
sportlong$"2014-09-01" <- sportlong$sport14/(sportlong$poplsoa/1000)
sportlong$"2013-06-01" <- sportlong$sport13/(sportlong$poplsoa/1000)
sportlong$"2012-06-01" <- sportlong$sport12/(sportlong$poplsoa/1000)
sportlong$"2011-06-01" <- sportlong$sport11/(sportlong$poplsoa/1000)
sportlong$"2010-09-01" <- sportlong$sport10/(sportlong$poplsoa/1000)
sportlong$"2009-06-01" <- sportlong$sport09/(sportlong$poplsoa/1000)
sportlong$"2008-03-01" <- sportlong$sport08/(sportlong$poplsoa/1000)

nightlong$"2018-06-01" <- nightlong$night18/(nightlong$poplsoa/1000)
nightlong$"2017-06-01" <- nightlong$night17/(nightlong$poplsoa/1000)
nightlong$"2016-06-01" <- nightlong$night16/(nightlong$poplsoa/1000)
nightlong$"2015-06-01" <- nightlong$night15/(nightlong$poplsoa/1000)
nightlong$"2014-09-01" <- nightlong$night14/(nightlong$poplsoa/1000)
nightlong$"2013-06-01" <- nightlong$night13/(nightlong$poplsoa/1000)
nightlong$"2012-06-01" <- nightlong$night12/(nightlong$poplsoa/1000)
nightlong$"2011-06-01" <- nightlong$night11/(nightlong$poplsoa/1000)
nightlong$"2010-09-01" <- nightlong$night10/(nightlong$poplsoa/1000)
nightlong$"2009-06-01" <- nightlong$night09/(nightlong$poplsoa/1000)
nightlong$"2008-03-01" <- nightlong$night08/(nightlong$poplsoa/1000)

cinelong$"2018-06-01" <- cinelong$cine18/(cinelong$poplsoa/1000)
cinelong$"2017-06-01" <- cinelong$cine17/(cinelong$poplsoa/1000)
cinelong$"2016-06-01" <- cinelong$cine16/(cinelong$poplsoa/1000)
cinelong$"2015-06-01" <- cinelong$cine15/(cinelong$poplsoa/1000)
cinelong$"2014-09-01" <- cinelong$cine14/(cinelong$poplsoa/1000)
cinelong$"2013-06-01" <- cinelong$cine13/(cinelong$poplsoa/1000)
cinelong$"2012-06-01" <- cinelong$cine12/(cinelong$poplsoa/1000)
cinelong$"2011-06-01" <- cinelong$cine11/(cinelong$poplsoa/1000)
cinelong$"2010-09-01" <- cinelong$cine10/(cinelong$poplsoa/1000)
cinelong$"2009-06-01" <- cinelong$cine09/(cinelong$poplsoa/1000)
cinelong$"2008-03-01" <- cinelong$cine08/(cinelong$poplsoa/1000)


###findings###


bet <- plm(betz ~ jobz , data=artspanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
charity <- plm(charityz ~ jobz, data=artspanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
sport <- plm(sportz ~ jobz, data=artspanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
youth <- plm(youthz ~ jobz, data=artspanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
cine <- plm(cinez ~ jobz, data=artspanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
disco <- plm(discoz ~ jobz, data=artspanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
educat <- plm(educatz ~ jobz, data=artspanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")


stargazer(bet, charity, sport, youth, cine, disco, educat, style= "apsr", title="Results", align=TRUE, type="text")


engpanel <- artspanel[which(artspanel$country=="England"), ]
walespanel <- artspanel[which(artspanel$country=="Wales"), ]
scotpanel <- artspanel[which(artspanel$country=="Scotland"), ]


bet <- plm(betz ~ jobz , data=engpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
charity <- plm(charityz ~ jobz, data=engpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
sport <- plm(sportz ~ jobz, data=engpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
youth <- plm(youthz ~ jobz, data=engpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
cine <- plm(cinez ~ jobz, data=engpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
disco <- plm(discoz ~ jobz, data=engpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
educat <- plm(educatz ~ jobz, data=engpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")


stargazer(bet, charity, sport, youth, cine, disco, educat, style= "apsr", title="Results", align=TRUE, type="text")


bet <- plm(betz ~ jobz , data=scotpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
charity <- plm(charityz ~ jobz, data=scotpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
sport <- plm(sportz ~ jobz, data=scotpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
youth <- plm(youthz ~ jobz, data=scotpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
cine <- plm(cinez ~ jobz, data=scotpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
disco <- plm(discoz ~ jobz, data=scotpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
educat <- plm(educatz ~ jobz, data=scotpanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")


stargazer(bet, charity, sport, youth, cine, disco, educat, style= "apsr", title="Results", align=TRUE, type="text")

bet <- plm(betz ~ jobz , data=walespanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
charity <- plm(charityz ~ jobz, data=walespanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
sport <- plm(sportz ~ jobz, data=walespanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
youth <- plm(youthz ~ jobz, data=walespanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
cine <- plm(cinez ~ jobz, data=walespanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
disco <- plm(discoz ~ jobz, data=walespanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")
educat <- plm(educatz ~ jobz, data=walespanel, index=c( "BUA11CD", "date"), model="within", effect= "twoways")


stargazer(bet, charity, sport, youth, cine, disco, educat, style= "apsr", title="Results", align=TRUE, type="text")


artspanel <- artspanel %>% mutate(year= ifelse(artspanel$date=="2014-09-01"|artspanel$date=="2014-12-01", 2014, ifelse(artspanel$date=="2015-09-01"|artspanel$date=="2015-12-01"|artspanel$date=="2015-03-01"|artspanel$date=="2015-06-01", 2015, ifelse(artspanel$date=="2016-09-01"|artspanel$date=="2016-12-01"|artspanel$date=="2016-03-01"|artspanel$date=="2016-06-01", 2016, ifelse(artspanel$date=="2017-09-01"|artspanel$date=="2017-12-01"|artspanel$date=="2017-03-01"|artspanel$date=="2017-06-01", 2017, ifelse(artspanel$date=="2018-09-01"|artspanel$date=="2018-03-01"|artspanel$date=="2018-06-01", 2018, NA))))))



bet <- plm(betz ~ jobz +factor(year) , data=artspanel, index= "BUA11CD", model="within")
charity <- plm(charityz ~ jobz+factor(year), data=artspanel, index= "BUA11CD", model="within")
sport <- plm(sportz ~ jobz+factor(year), data=artspanel, index= "BUA11CD", model="within")
youth <- plm(youthz ~ jobz+factor(year), data=artspanel, index= "BUA11CD", model="within")
cine <- plm(cinez ~ jobz+factor(year), data=artspanel, index= "BUA11CD", model="within")
disco <- plm(discoz ~ jobz+factor(year), data=artspanel, index= "BUA11CD", model="within")
educat <- plm(educatz ~ jobz+factor(year), data=artspanel, index= "BUA11CD", model="within")


stargazer(bet, charity, sport, youth, cine, disco, educat, style= "apsr", title="Results", align=TRUE, type="text")



engmaster <- OSmaster[which(OSmaster$country=="England"), ]
walesmaster <- OSmaster[which(OSmaster$country=="Wales"), ]
scotmaster <- OSmaster[which(OSmaster$country=="Scotland"), ]


