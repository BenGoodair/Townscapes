
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


###creating the geographies and values###

buas <- readShapeSpatial("bua")
buanames <- buas[-c(1,4,5,6,7,8,9)]
names(buanames)[names(buanames)=="bua11cd"] <- "BUA11CD"
names(buas)[names(buas)=="bua11cd"] <- "BUA11CD"
perseq <- seq(-100, 100, by=20)
ukcrs <- CRS("+init=epsg:27700")
popbua <- read_csv("crosswalkpoplsoa.txt")
popbua <- popbua[-c(1)]
popbua <- aggregate(popbua[-1], popbua["BUA11CD"], sum)
proj4string(buas) <- crs(ukcrs)
collar <- (brewer.pal( 7, "RdBu"))
citybua <- merge(buas, popbua, by="BUA11CD")
citybua <- citybua[which(citybua$poplsoa>175000), ]
citylayer <- list( "sp.polygons",citybua, fill= "gray77", col="transparent")
theme_set(theme_minimal())

###Download my data###



businesscounts10 <- read_csv("businesscountsfeb.csv")
businesscounts18 <- read_csv("businesscountsfeb2.csv")
feb10 <- read_csv("feb10.csv")
feb18 <- read_csv("feb18.csv")
trtot10 <- read_csv("trtot10.csv")
trtot18 <- read_csv("trtot.csv")
tim11 <- read_csv("tottime11.csv")
tim16 <- read_csv("tottime16.csv")
pub10 <- read_csv("pub10.csv")
pub18 <- read_csv("pub18.csv")
pubtime <- read_csv("earlypubtime.csv")
pubtimes <- read.csv("pubtimeline.csv")
retimenew <- read_csv("retimenew.csv", col_types = cols(ret18 = col_integer()))
retail10 <- read_csv("retail10.csv")
retail18 <- read_csv("retail18.csv")
commuter <- read_csv("commute.csv")
manufacture <- read_csv("manufacture.csv")
unihalls <- read_csv("unihalls.csv")
retired <- read_csv("retired.csv")
popdensity <- read_csv("popdensity.csv")
rettimeold <- read_csv("rettimeold.csv")
brexit <- read_csv("brexit.csv", col_types = cols(Brexit = col_factor(levels = c("Remain", 
     "Leave")), Leave = col_number(), 
     `MP party` = col_factor(levels = c("Conservative", 
      "Labour", "Liberal Democrat", "Independent", 
     "Plaid Cymru"))))
msoa11 <- readShapeSpatial("20111")
msoa01 <- readShapeSpatial("2001")
deprivation <- read_csv("deprivation.csv")


###data cleaning and aggregating###

deprivation$dep <- deprivation$multipledeprivation*deprivation$pop
deprivation <- aggregate(deprivation[-1], deprivation["BUA11CD"], sum)
deprivation$dep <- deprivation$dep/deprivation$pop
deprivation <- deprivation[-c(2,3)]
retail10 <- aggregate(retail10[-1], retail10["BUA11CD"], sum)
retail18 <- aggregate(retail18[-1], retail18["BUA11CD"], sum)
feb18 <- aggregate(feb18[-1], feb18["BUA11CD"], sum)
feb10 <- aggregate(feb10[-1], feb10["BUA11CD"], sum)
pub10 <- aggregate(pub10[-1], pub10["BUA11CD"], sum)
pub18 <- aggregate(pub18[-1], pub18["BUA11CD"], sum)
tim11 <- aggregate(tim11[-1], tim11["BUA11CD"], sum)
tim16 <- aggregate(tim16[-1], tim16["BUA11CD"], sum)
trtot18 <- aggregate(trtot18[-1], trtot18["BUA11CD"], sum)
trtot10 <- aggregate(trtot10[-1], trtot10["BUA11CD"], sum)
pubtime <- aggregate(pubtime[-1], pubtime["BUA11CD"], sum)
pubtimes <- aggregate(pubtimes[-1], pubtimes["BUA11CD"], sum)
retimenew <- aggregate(retimenew[-1], retimenew["BUA11CD"], sum)
rettimeold <- aggregate(rettimeold[-1], rettimeold["BUA11CD"], sum)
businesscounts18 <- aggregate(businesscounts18[-1], businesscounts18["BUA11CD"], sum)
businesscounts10 <- aggregate(businesscounts10[-1], businesscounts10["BUA11CD"], sum)
businesscounts <- merge(businesscounts10, businesscounts18, by="BUA11CD")
businesscounts <- merge(businesscounts, buanames, by="BUA11CD")
businesscounts <- merge(businesscounts, popbua, by="BUA11CD")
businesscounts <- merge(businesscounts, feb10, by="BUA11CD")
businesscounts <- merge(businesscounts, feb18, by="BUA11CD")
businesscounts <- merge(businesscounts, tim16, by="BUA11CD")
businesscounts <- merge(businesscounts, tim11, by="BUA11CD")
businesscounts <- merge(businesscounts, trtot10, by="BUA11CD")
businesscounts <- merge(businesscounts, trtot18, by="BUA11CD")
businesscounts <- merge(businesscounts, retail10, by="BUA11CD")
businesscounts <- merge(businesscounts, retail18, by="BUA11CD")
businesscounts <- merge(businesscounts, pub18, by="BUA11CD")
businesscounts <- merge(businesscounts, pub10, by="BUA11CD")
businesscounts <- merge(businesscounts, pubtime, by="BUA11CD")
businesscounts <- merge(businesscounts, pubtimes, by="BUA11CD")
businesscounts <- merge(businesscounts, rettimeold, by="BUA11CD")
businesscounts <- merge(businesscounts, retimenew, by="BUA11CD")
businesscounts <- merge(businesscounts, unihalls, by="bua11nm")
businesscounts <- merge(businesscounts, manufacture, by="bua11nm")
businesscounts <- merge(businesscounts, retired, by="bua11nm")
businesscounts <- merge(businesscounts, popdensity, by="bua11nm")
businesscounts <- merge(businesscounts, commuter, by="bua11nm")
businesscounts <- merge(businesscounts, deprivation, by="BUA11CD", all.x=TRUE)

###gambling data insights###
poisbet <- read_csv("gamplingpoiclean.csv")
betper <- read_csv("betper.csv")
names(poisbet)[names(poisbet)=="PointX Classification Code"] <- "code"
poisbet <- poisbet[which(poisbet$code==04220279), ]
poisbet <- subset(poisbet, "Feature Easting" != "" | "Feature Northing" != "")
poisbet$poi_ID <- 1:nrow(poisbet)
coords <- cbind(Easting = as.numeric(as.character(poisbet$"Feature Easting")),
                Northing = as.numeric(as.character(poisbet$"Feature Northing")))
poisbetsp <- SpatialPointsDataFrame(coords, data = data.frame(poisbet$Name,
                                                              poisbet$poi_ID), proj4string = CRS("+init=epsg:27700"))


betlist <- over(buas, geometry(poisbetsp), returnList = TRUE)
betnum <- sapply(betlist, length)
betgeogs <- spCbind(buas, betnum)
betgeogs <- merge(betgeogs, popbua, by="BUA11CD")
betgeogs <- betgeogs[which(betgeogs$poplsoa<175000&betgeogs$poplsoa>10000), ]
betcent <- SpatialPointsDataFrame(gCentroid(betgeogs, byid=TRUE), betgeogs@data, match.ID=FALSE)
betgeogz <- betgeogs[-c(2,3,4,5,6,7,8,9,11)]
businesscounts <- merge(businesscounts, betgeogz, by="BUA11CD")
businesscounts <- merge(businesscounts, betper, by="BUA11CD")

buasize <- buas[-c(1,3,4,5,6,7,9)]
buasize$size <- (buasize$st_areasha/1000)/1000
buasize <- buasize[-c(2)]
businesscounts <- merge(businesscounts, buasize, by="BUA11CD")

rm(trtot10, trtot18, tim11, tim16, feb10, feb18, pub10, pub18, businesscounts10, businesscounts18, retail10, unihalls, commuter, manufacture, retired, popdensity, retimenew, rettimeold, retail18, pubtime, pubtimes)

###creating percentile changes###


businesscounts$foodper <- (businesscounts$food18-businesscounts$food10)/businesscounts$food10*100
businesscounts$totalper <- (businesscounts$total18-businesscounts$total10)/businesscounts$total10*100
businesscounts$sportper <- (businesscounts$sport18-businesscounts$sport10)/businesscounts$sport10*100
businesscounts$accomper <- (businesscounts$accom18-businesscounts$accom10)/businesscounts$accom10*100
businesscounts$edper <- (businesscounts$ed18-businesscounts$ed10)/businesscounts$ed10*100
businesscounts$healthper <- (businesscounts$health18-businesscounts$health10)/businesscounts$health10*100
businesscounts$betper <- (businesscounts$bet18-businesscounts$bet10)/businesscounts$bet10*100
businesscounts$vetper <- (businesscounts$vet18-businesscounts$vet10)/businesscounts$vet10*100
businesscounts$careper <- (businesscounts$care18-businesscounts$care10)/businesscounts$care10*100
businesscounts$artper <- (businesscounts$art18-businesscounts$art10)/businesscounts$art10*100
businesscounts$cultper <- (businesscounts$culture18-businesscounts$cult10)/businesscounts$cult10*100
businesscounts$nonspecper <- (businesscounts$nonspec18-businesscounts$nonspec10)/businesscounts$nonspec10*100
businesscounts$foodspecper <- (businesscounts$foodspec18-businesscounts$foodspec10)/businesscounts$foodspec10*100
businesscounts$techspecper <- (businesscounts$tech18-businesscounts$tech10)/businesscounts$tech10*100
businesscounts$otherspecper <- (businesscounts$otherspec18-businesscounts$otherspec10)/businesscounts$otherspec10*100
businesscounts$housespecper <- (businesscounts$housespec18-businesscounts$housespec10)/businesscounts$housespec10*100
businesscounts$recspecper <- (businesscounts$recspec18-businesscounts$recspec10)/businesscounts$recspec10*100
businesscounts$rettotalper <- (businesscounts$ret18-businesscounts$ret10)/businesscounts$ret10*100
businesscounts$pubedper <- (businesscounts$pubed18-businesscounts$pubed10)/businesscounts$pubed10*100
businesscounts$pubadper <- (businesscounts$pubad18-businesscounts$pubad10)/businesscounts$pubad10*100
businesscounts$pubhealthper <- (businesscounts$pubhealth18-businesscounts$pubhealth10)/businesscounts$pubhealth10*100
businesscounts$pubcoreper <- (businesscounts$pubcore18-businesscounts$pubcore10)/businesscounts$pubcore10*100
businesscounts$pubtotalper <- (businesscounts$pubtotal18.x-businesscounts$pubtotal10.x)/businesscounts$pubtotal10.x*100
businesscounts$retirementper <- businesscounts$aged65/businesscounts$poplsoa*100
businesscounts$commuteper <- businesscounts$`railway commuters`/businesscounts$poplsoa*100
businesscounts$manufactureper <- businesscounts$manufact/businesscounts$workinagepop*100
businesscounts$betpersqkm <- businesscounts$betnum/businesscounts$size

### getting rid of villages and cities###

businesscountstowns <- businesscounts[which(businesscounts$poplsoa<175000&businesscounts$poplsoa>10000), ]
townsample <- businesscountstowns[-c(3:117)]

towngeog <- merge(buas, townsample, by="BUA11CD")
towngeog <- sp.na.omit(towngeog)
towngeog <- merge(towngeog, businesscountstowns)

towngeogcent <- SpatialPointsDataFrame(gCentroid(towngeog, byid=TRUE), towngeog@data, match.ID=FALSE)


###categorising towns###

businesscounts <- businesscounts %>% mutate(retirement = ifelse(retirementper >25 , 1, ifelse(retirementper < 25 , 0, NA )))
businesscounts <- businesscounts %>% mutate(comm = ifelse(commuteper >5 , 1, ifelse(commuteper < 5 , 0, NA )))
businesscounts <- businesscounts %>% mutate(university = ifelse(universityhallsres>0, 1, ifelse(universityhallsres==0, 0, NA)))
businesscounts <- businesscounts %>% mutate(industrial = ifelse(manufactureper >10 , 1, ifelse(manufactureper < 10 , 0, NA )))



businesscounts$unif <- factor(businesscounts$university)


businesscounts <- businesscounts %>% mutate(category = ifelse(retirement==1 &industrial==0 & comm == 0 & unif==0, "Retirement", ifelse(unif==1 & comm==0 & industrial==0&retirement==0, "University", ifelse(industrial==1 &unif==0&retirement==0&comm==0, "Post-Industrial", ifelse(comm==1&industrial==0&unif==0&retirement==0, "Commuter", NA)))))
businesscounts$categoryf <- factor(businesscounts$category)
businesscountscatstowns <- businesscounts[which(businesscounts$poplsoa<175000&businesscounts$poplsoa>10000), ]                                                           
businesscountscategs <- businesscountscatstowns[which(businesscountscatstowns$categoryf != 0), ]
businesscountscategs <- merge(businesscountscategs, brexit, by="BUA11CD")

###making categorised geographies###





catsample <- businesscountscategs[-c(3:117)]

catsgeog <- merge(buas, catsample, by="BUA11CD")
catsgeog <- sp.na.omit(catsgeog)
catsgeog <- merge(catsgeog, businesscountscategs, by="BUA11CD")

catsgeogcent <- SpatialPointsDataFrame(gCentroid(catsgeog, byid=TRUE), catsgeog@data, match.ID=FALSE)






###category map ###

spplot(catsgeogcent, "categoryf.x", main=list(label="Town Categories",cex=1), key.space = "right", sp.layout=citylayer, pch=c(16,17,18,15))


##Brexit maps###

universitytowns <- businesscountscategs[which(businesscountscategs$categoryf=="University"), ]
industrialtowns <- businesscountscategs[which(businesscountscategs$categoryf=="Post-Industrial"), ]
commutertowns <- businesscountscategs[which(businesscountscategs$categoryf=="Commuter"), ]
retirementtowns <- businesscountscategs[which(businesscountscategs$categoryf=="Retirement"), ]
cols <- c("yellow2","blue")

universitysample <- universitytowns[-c(3:105)]
universitygeog <- merge(buas,universitysample, by="BUA11CD")
universitygeog <- sp.na.omit(universitygeog)

commutersample <- commutertowns[-c(3:105)]
commutergeog <- merge(buas,commutersample, by="BUA11CD")
commutergeog <- sp.na.omit(commutergeog)

retirementsample <- retirementtowns[-c(3:105)]
retirementgeog <- merge(buas,retirementsample, by="BUA11CD")
retirementgeog <- sp.na.omit(retirementgeog)

industrialsample <- industrialtowns[-c(3:105)]
industrialgeog <- merge(buas,industrialsample, by="BUA11CD")
industrialgeog <- sp.na.omit(industrialgeog)

universitytowns <- merge(universitygeog, universitytowns, by="BUA11CD")
retirementtowns <- merge(retirementgeog, retirementtowns, by="BUA11CD")
commutertowns <- merge(commutergeog, commutertowns, by="BUA11CD")
industrialtowns <- merge(industrialgeog, industrialtowns, by="BUA11CD")

universitytowns <- SpatialPointsDataFrame(gCentroid(universitytowns, byid=TRUE), universitytowns@data, match.ID=FALSE)
commutertowns <- SpatialPointsDataFrame(gCentroid(commutertowns, byid=TRUE), commutertowns@data, match.ID=FALSE)
industrialtowns <- SpatialPointsDataFrame(gCentroid(industrialtowns, byid=TRUE), industrialtowns@data, match.ID=FALSE)
retirementtowns <- SpatialPointsDataFrame(gCentroid(retirementtowns, byid=TRUE), retirementtowns@data, match.ID=FALSE)



spplot(universitytowns, "Brexit", col.regions = cols,   main=list(label="EU Referendum district result for University towns",cex=1), key.space = "right", sp.layout=citylayer, pch=c(15,17))
spplot(industrialtowns, "Brexit", col.regions = cols,   main=list(label="EU Referendum district result for Post-Indstrial towns",cex=1), key.space = "right", sp.layout=citylayer, pch=c(15,17))
spplot(commutertowns, "Brexit", col.regions = cols,   main=list(label="EU Referendum district result for Commuter towns",cex=1), key.space = "right", sp.layout=citylayer, pch=c(15,17))
spplot(retirementtowns, "Brexit", col.regions = cols,   main=list(label="EU Referendum district result for Retirement towns",cex=1), key.space = "right", sp.layout=citylayer, pch=c(15,17))

###stagnationmaps###

stagnanttowns <- businesscountstowns[ which(businesscountstowns$totalper<1 & businesscountstowns$totalper>-1), ]
successtowns <- businesscountstowns[ which(businesscountstowns$totalper>30 ), ]

stagnantsample <- stagnanttowns[-c(3:95)]
successsample <- successtowns[-c(3:95)]

stagnantgeog <- merge(buas, stagnantsample, by="BUA11CD")
successgeog <- merge(buas, successsample, by="BUA11CD")

stagnantgeog <- sp.na.omit(stagnantgeog)
successgeog <- sp.na.omit(successgeog)

stagnanttowns <- merge(stagnantgeog, stagnanttowns, by="BUA11CD")
successtowns <- merge(successgeog, successtowns, by="BUA11CD")

successtowns <- SpatialPointsDataFrame(gCentroid(successtowns, byid=TRUE), successtowns@data, match.ID=FALSE)
stagnanttowns <- SpatialPointsDataFrame(gCentroid(stagnanttowns, byid=TRUE), stagnanttowns@data, match.ID=FALSE)

spplot(stagnanttowns, "populationperhectare", sp.layout=citylayer, col.regions="red", main=list(label="Towns with less than 1% change in service businesses", colorkey=FALSE, legend=FALSE))
spplot(successtowns, "populationperhectare", sp.layout=citylayer, col.regions="red", main=list(label="Towns with more than 30% increase in service businesses", colorkey=FALSE, legend=FALSE))


rm(townsample, catsample, towngeog, catsgeog, brexit)

###all service maps###
spplot(towngeogcent, "totalper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of total service industry businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "betper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Gambling  businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "artper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Arts and Entertainment  businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "accomper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Accommodation  businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "healthper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Health  businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "edper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Education businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "sportper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Sports and Recreation businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "foodper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Food and Drink businesses between 2010-18 (%)",cex=0.9), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))



###all retail maps###
  
spplot(towngeogcent, "nonspecper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of 'non-specialised retail' businesses between 2010-18 (%)",cex=0.8), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "foodspecper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of 'food and drink specialised retail' businesses between 2010-18 (%)",cex=0.6), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "otherspecper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of 'other specialised retail' businesses between 2010-18 (%)",cex=0.8), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "recspecper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of 'cultural and recreation specialised retail' businesses between 2010-18 (%)",cex=0.6), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "housespecper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of 'household specialised retail' businesses between 2010-18 (%)",cex=0.7), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "techspecper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of 'technology specialised retail' businesses between 2010-18 (%)",cex=0.7), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))



###all public service maps###
spplot(towngeogcent, "pubedper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of State Education providers between 2010-18 (%)",cex=0.8), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "pubadper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of Public Administration providers between 2010-18 (%)",cex=0.8), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "pubhealthper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of public Health providers between 2010-18 (%)",cex=0.8), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "pubtotalper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of total public service providers between 2010-18 (%)",cex=0.8), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))
spplot(towngeogcent, "pubcoreper", col.regions = collar, cuts = perseq, col="transparent", main=list(label="Percentage change in number of 'Core Service' providers between 2010-18 (%)",cex=0.8), key.space = "right", sp.layout=citylayer, pch=c(15,15,15,16,16,16, 16, 17, 17, 17))

##city map###

citymap <- spplot(citybua, "poplsoa", col.regions="gray77", col="transparent", main=list(label="Cities"))
citynamelayer <- layer(sp.text(coordinates(citybua), txt=citybua$bua11nm, pos=1, cex=0.5))
citymap + citynamelayer


###all category plots###
  
publicplot <- xyplot(pubtotalper~categoryf, data=businesscountscategs, groups=categoryf, type=c("p"), xlab="Town category", ylab="Change in public service count (%)", main="Percentage change in public service providers between 2010 and 2018", pch=19)
update(publicplot, panel = function(...) {
       panel.abline(h = 0, v= 0, lty="dotted", col = "light gray")
       panel.xyplot(...)
   })

brexitcatplot <- xyplot(Leave~categoryf, data=businesscountscategs, groups=categoryf, type=c("p"), xlab="Town category", ylab="District Leave vote (%)", main="Percentage of Distict voting to Leave the EU", pch=19)
update(brexitcatplot, panel = function(...) {
  panel.abline(h = 50, v= 0, lty="dotted", col = "light gray")
  panel.xyplot(...)
   })


serviceplot <- xyplot(totalper~categoryf, data=businesscountscategs, groups=categoryf, type=c("p"), xlab="Town category", ylab="Change in business count (%)", main="Percentage change in service industry businesses between 2010 and 2018", pch=19)
update(serviceplot, panel = function(...) {
  panel.abline(h = 0, v= 0, lty="dotted", col = "light gray")
  panel.xyplot(...)
})

retailplot <- xyplot(rettotalper~categoryf, data=businesscountscategs, groups=categoryf, type=c("p"), xlab="Town category", ylab="Change in retail outlet business count (%)", main="Percentage change in retail outlets between 2010 and 2018", pch=19)
update(retailplot, panel = function(...) {
  panel.abline(h = 0, v= 0, lty="dotted", col = "light gray")
  panel.xyplot(...)
})


foodvrecretail <- xyplot(foodspecper~recspecper, data=businesscountstowns, type=c("p"),  xlab="Change in Recreation Specialised Retail (%)", ylab="Change in Food Specialised Retail (%)", main="Changes in different areas of retail", pch=19)
update(foodvrecretail, panel = function(...) {
  panel.abline(h = 0, v= 0, lty="dotted", col = "light gray")
  panel.xyplot(...)
})

housevnonspecretail <- xyplot(nonspecper~housespecper, data=businesscountstowns, type=c("p"),  xlab="Change in Household Specialised Retail (%)", ylab="Change in non-Specialised Retail (%)", main="Changes in different areas of retail", pch=19)
update(housevnonspecretail, panel = function(...) {
  panel.abline(h = 0, v= 0, lty="dotted", col = "light gray")
  panel.xyplot(...)
})

brexitplot <- ggplot(businesscountscategs, aes(x = Leave, y = totalper))
brexitplot + geom_point(aes(color = categoryf, shape = categoryf))+
            stat_ellipse(aes(color = categoryf), type = "t")+
            scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "black"))+
            scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "black"))+
            labs(x="Percent of District voting to Leave the EU (%)", y="Change in total service businesses between 2010 and 2018 (%)", title="2016 EU Referendum vote and change in business services")



### checking msoa boundary changes ###

output11 <- read_csv("pub18msoa.csv")
output11 <- output11[-c(2,4,5,6,7,8)]
names(output11)[names(output11)=="MSOA11CD"] <- "msoa11cd"
output11 <- merge(msoa11, output11, by="msoa11cd")
output11 <- sp.na.omit(output11)

output01 <- read_csv("pub10msoa.csv")
output01 <- output01[-c(2,3,5,6,7,8,9, 10)]
names(output01)[names(output01)=="MSOA01CD"] <- "msoa01cd"
output01 <- merge(msoa01, output01, by="msoa01cd")

output01 <- sp.na.omit(output01)

msoalayer <- list("sp.polygons", output01, fill="transparent", col="blue", lwd=0.45, lty=2)


spplot(output11, "BUA11CD", fill="transparent", sp.layout=msoalayer, main=list(label="Outlines of 2011 MSOAs and dotted blue 2001 MSOAs",cex=0.9))



### fancy map betting###
       
leafletry2 <- spTransform(towngeogcent, CRS("+init=epsg:4326"))
      
getColor <- function(leafletry2) {
              sapply(leafletry2$betpersqkm, function(betpersqkm) {
                  if(betpersqkm <= 0.75) {
                      "green"
                 } else if(betpersqkm <= 1.5) {
                      "orange"
                  } else {
                      "red"
                          } })
}


getColorz <- function(leafletry2) {
  sapply(leafletry2$totalper, function(totalper) {
    if(totalper >= 10) {
      "green"
    } else if(totalper >= 0) {
      "orange"
    } else {
      "red"
    } })
}

getColorr <- function(leafletry2) {
  sapply(leafletry2$rettotalper, function(rettotalper) {
    if(rettotalper >= 10) {
      "green"
    } else if(rettotalper >= 0) {
      "orange"
    } else {
      "red"
    } })
  
}

getColorp <- function(leafletry2) {
  sapply(leafletry2$pubtotalper, function(pubtotalper) {
    if(pubtotalper >= 10) {
      "green"
    } else if(pubtotalper >= 0) {
      "orange"
    } else {
      "red"
    } })
}


iconz <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorz(leafletry2))

iconp <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorp(leafletry2))

iconr <- awesomeIcons(
  icon = 'ios-close',
  iconColor = 'black',
  library = 'ion',
  markerColor = getColorr(leafletry2))

icons <- awesomeIcons(
              icon = 'ios-close',
              iconColor = 'black',
              library = 'ion',
              markerColor = getColor(leafletry2))
         


mymap <- leaflet(data=leafletry2) %>%
                 addTiles() %>%  
                 setView(lng=-1.3823, lat=53.0975, zoom=5) 
             addAwesomeMarkers(mymap, data=leafletry2, icon=icons,
                popup = ~as.character(towngeogcent$betname), label= ~as.character(towngeogcent$betname) layerId = "Betting Shops")
             
             
mymap <- leaflet(data=leafletry2) %>%
                 addTiles() %>%  
                  setView(lng=-1.3823, lat=53.0975, zoom=5) 
addAwesomeMarkers(mymap, data=leafletry2,icon=iconz,
                   popup = ~as.character(towngeogcent$totalper), label= ~as.character(towngeogcent$totalper), group = "Change in service industry businesses 10-18 (%)")%>%
  addAwesomeMarkers( data=leafletry2, icon=iconp,
                     popup = ~as.character(towngeogcent$pubtotalper), label= ~as.character(towngeogcent$pubtotalper), group = "Change in public service providers 10-18 (%)")%>%
  addAwesomeMarkers( data=leafletry2, icon=iconr,
                     popup = ~as.character(towngeogcent$rettotalper), label= ~as.character(towngeogcent$rettotalper), group = "Change in retail outlets 10-18 (%)")%>%
   addAwesomeMarkers( data=leafletry2, icon=icons,
                 popup = ~as.character(towngeogcent$betname), label= ~as.character(towngeogcent$betname), group = "Betting shops per square km")%>%
   addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("<0.75", "0.75-1.5", ">1.5"), opacity = 1,
    title = "Betting shops per square kilometer", 
    group = "Legends"
  ) %>%
  addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("<10%", "0-10%", "<0%"), opacity = 1,
    title = "Change in service industry businesses 10-18", 
    group = "Legends"
  ) %>%
   addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("<10%", "0-10%", "<0%"), opacity = 1,
    title = "Betting shops per square kilometer", 
    group = "Legends"
  ) %>%
   addLegend(
    position = "bottomright",
    colors = c("green", "orange", "red"),
    labels = c("<10%", "0-10%", "<0%"), opacity = 1,
    title = "Change in public service providers 10-18", 
    group = "Legends"
  ) %>%
              addLayersControl(
                baseGroups = c("Betting shops per square km", "Change in service industry businesses 10-18 (%)", "Change in public service providers 10-18 (%)", "Change in retail outlets 10-18 (%)"), data=leafletry2,
                 overlayGroups = "Legends", position= "topright",
                options = layersControlOptions(collapsed = FALSE)
              ) 
### Density Playing ###
cambgamb <- read_csv("cambgamb.csv", col_types = cols(`Feature Easting` = col_double(), 
                                                           `Feature Northing` = col_double()))
names(cambgamb)[names(cambgamb)=="PointX Classification Code"] <- "code"
cambgamb <- cambgamb[which(cambgamb$code==04220279), ]
cambgamb <- subset(cambgamb, "Feature Easting" != "" | "Feature Northing" != "")
cambgamb$poi_ID <- 1:nrow(cambgamb)
cambcoords <- cbind(Easting = as.numeric(as.character(cambgamb$"Feature Easting")),
                Northing = as.numeric(as.character(cambgamb$"Feature Northing")))
cambbetsp <- SpatialPointsDataFrame(cambcoords, data = data.frame(cambgamb$Name,
                                                              cambgamb$poi_ID), proj4string = CRS("+init=epsg:4326"))

cambetpoints <- as(SpatialPoints(cambbetsp), "ppp")
cambetdense <- density(cambetpoints, adjust = 0.2)

## chloropleth thyme ##

leafletry3 <- spTransform(towngeog, CRS("+init=epsg:4326"))
bins <- c(0,0.3,0.6,0.9,1.2,1.5,1.8,2.2)
pal <- colorBin("YlOrRd", domain = leafletry3$betpersqkm, bins = bins)
labels <- sprintf(
  "<strong>%s</strong><br/>%g Bookmakers / km <sup>2</sup>",
  leafletry3$name, leafletry3$betpersqkm
) %>% lapply(htmltools::HTML)

map2 <- leaflet(leafletry3) %>%
  setView(lng=-1.3823, lat=53.0975, zoom=8)%>%
  addTiles()%>%
  addPolygons(
    fillColor = ~pal(betpersqkm),
    weight = 1,
    opacity = 1,
    color = "white",
    dashArray = "3",
    fillOpacity = 0.7,
    highlight = highlightOptions(
      weight = 5,
      color = "#666",
      dashArray = "",
      fillOpacity = 0.7,
      bringToFront = TRUE),
    label = labels,
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",
      direction = "auto")) %>%
  addLegend(pal = pal, values = ~betpersqkm, opacity = 0.7, 
            position = "bottomright", title = "Betting shops per km <sup>2</sup>")


### panel data attempts###

jobseek <- read_csv("jobseekersallowance.csv")
jobz <- jobseek[-c(1)]
jobz <- reshape(jobseek, idvar = "name", ids = jobseek$bua11nm,
                times = names(jobseek), timevar = "date",
                varying = list(names(jobseek)),v.names="jobseekers", new.row.names = 1:((dim(jobseek)[2])*(dim(jobseek)[1])),direction = "long")