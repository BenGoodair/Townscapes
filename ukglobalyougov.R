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
###Constructing binary variables fo Culture vs Economic###


ukglobal <- read_csv("ukyougov.csv")
ukglobal$Glob_p02 <- factor(ukglobal$Glob_p02)
ukglobal$Glob_p05 <- factor(ukglobal$Glob_p05)
ukglobal <- ukglobal %>% mutate(populist=ifelse(ukglobal$Glob_p05=='Strongly agree' & ukglobal$Glob_p02=='Strongly agree', 1, 0))
ukglobal$populist <- factor(ukglobal$populist)
ukglobal$Glob_job <- factor(ukglobal$Glob_job)
ukglobal <- ukglobal %>% mutate(sacked=ifelse(ukglobal$Glob_job=='Very likely'| ukglobal$Glob_job=='Somewhat likely' , 1, 0))
ukglobal$sacked <- factor(ukglobal$sacked)
ukglobal$Glob_i10 <- factor(ukglobal$Glob_i10)
ukglobal <- ukglobal %>% mutate(reduceimmigration=ifelse(ukglobal$Glob_i10=='Should be reduced a lot'| ukglobal$Glob_i10=='Should be reduced a little' , 1, 0))
ukglobal$reduceimmigration <- factor(ukglobal$reduceimmigration)


cultureecon <- glm(populist~sacked+reduceimmigration,data=ukglobal, family = "binomial")
cultureeconmargins <- margins(cultureecon)
cultureeconmarginsplot <- summary(cultureeconmargins)
  
  
  
w <- ggplot(data = cultureeconmarginsplot) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Marginal effects of Economic/ Cultural on Populist")

###Now with only extremes in independent variables###

ukglobal <- read_csv("ukyougov.csv")
ukglobal$Glob_p02 <- factor(ukglobal$Glob_p02)
ukglobal$Glob_p05 <- factor(ukglobal$Glob_p05)
ukglobal <- ukglobal %>% mutate(populist=ifelse(ukglobal$Glob_p05=='Strongly agree' & ukglobal$Glob_p02=='Strongly agree', 1, 0))
ukglobal$populist <- factor(ukglobal$populist)
ukglobal$Glob_job <- factor(ukglobal$Glob_job)
ukglobal <- ukglobal %>% mutate(sacked=ifelse(ukglobal$Glob_job=='Very likely' , 1, 0))
ukglobal$sacked <- factor(ukglobal$sacked)
ukglobal$Glob_i10 <- factor(ukglobal$Glob_i10)
ukglobal <- ukglobal %>% mutate(reduceimmigration=ifelse(ukglobal$Glob_i10=='Should be reduced a lot' , 1, 0))
ukglobal$reduceimmigration <- factor(ukglobal$reduceimmigration)


cultureecon <- glm(populist~sacked+reduceimmigration,data=ukglobal, family = "binomial")
cultureeconmargins <- margins(cultureecon)
cultureeconmarginsplot <- summary(cultureeconmargins)



w <- ggplot(data = cultureeconmarginsplot) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Marginal effects of Economic/ Cultural on Populist")

#levels = c("under £5,000 per year","£5,000 to £9,999 per year","£10,000 to £14,999 per year","£15,000 to £19,999 per year","£20,000 to £24,999 per year","£25,000 to £29,999 per year","£30,000 to £34,999 per year","£35,000 to £39,999 per year","£40,000 to £44,999 per year","£45,000 to £49,999 per year","£50,000 to £59,999 per year","£60,000 to £69,999 per year","£70,000 to £99,999 per year","£100,000 to £149,999 per year","£150,000 and over","Don't know","Prefer not to answer"), labels = c("5k","10k","15k","20k","25k","30k","35k","40k","45k","50k","60k","70k","100k","155k","155+","Don't Know","Refuse"
###Place and control variables###
ukglobal$settlement <- factor(ukglobal$Glob_gov,  levels = c("Centre of a city/large town", "Suburb or part of a city/large town, which is outside its centre","Small town","Village","Settlement or isolated dwelling smaller than a village","Don't know","None of these" ), labels = c("City","Suburb","Town","Village","Hamlet","Don't know","Other"))
ukglobal$education <- factor(ukglobal$edu_grou, levels = c("Less than primary, lower secondary education:  (Low \u0096 GCSE and below)", "Upper secondary and post-secondary: (Medium \u0096 roughly completed A-levels)", "Tertiary: (High \u0096 tertiary education: advanced professional qualification/ degre"), labels = c("Lower","Secondary","Tertiary"))
ukglobal$age_grp <- factor(ukglobal$age_grp_, levels = c("18 - 24", "25 - 34", "35 - 44", "45 - 54", "55+"))
ukglobal$income <- factor(ukglobal$profile0)




ukglobal$Glob_gen <- factor(ukglobal$Glob_gen)
ukglobal$Glob_ge0 <- factor(ukglobal$Glob_ge0)
ukglobal$Glob_mu9 <- factor(ukglobal$Glob_mu9)
ukglobal$Glob_re0 <- factor(ukglobal$Glob_re0)


ukglobal <- ukglobal %>% mutate(globalisationoverallbad=ifelse(ukglobal$Glob_gen=='Fairly bad'| ukglobal$Glob_gen=='Very bad' , 1, 0))
ukglobal <- ukglobal %>% mutate(globalisationlocalbad=ifelse(ukglobal$Glob_ge0=='Fairly bad'| ukglobal$Glob_ge0=='Very bad' , 1, 0))
ukglobal <- ukglobal %>% mutate(antiislam=ifelse(ukglobal$Glob_re0=='Fairly unfavourable'| ukglobal$Glob_re0=='Very unfavourable' , 1, 0))
ukglobal <- ukglobal %>% mutate(EUbad=ifelse(ukglobal$Glob_mu9=='Fairly negative'| ukglobal$Glob_mu9=='Very negative' , 1, 0))

ukglobal$globalisationoverallbad <- factor(ukglobal$globalisationoverallbad)
ukglobal$globalisationlocalbad <- factor(ukglobal$globalisationlocalbad)
ukglobal$EUbad <- factor(ukglobal$EUbad)
ukglobal$antiislam <- factor(ukglobal$antiislam)





globoverbad <- glm(globalisationoverallbad~settlement+education+age_grp+income, data=ukglobal, family = "binomial")
globlocalbad <- glm(globalisationlocalbad~settlement+education+age_grp+income, data=ukglobal, family = "binomial")
EUbad <- glm(EUbad~settlement+education+age_grp+income, data=ukglobal, family = "binomial")
islambad <- glm(antiislam~settlement+education+age_grp+income, data=ukglobal, family = "binomial")

stargazer(globoverbad, globlocalbad, EUbad, islambad, style= "apsr", title="Results", align=TRUE, type="text")

globoverbadnocontrol <- glm(globalisationoverallbad~settlement, data=ukglobal, family = "binomial")
globlocalbadnocontrol <- glm(globalisationlocalbad~settlement, data=ukglobal, family = "binomial")
EUbadnocontrol <- glm(EUbad~settlement, data=ukglobal, family = "binomial")
islambadnocontrol <- glm(antiislam~settlement, data=ukglobal, family = "binomial")

stargazer(globoverbadnocontrol, globlocalbadnocontrol, EUbadnocontrol, islambadnocontrol, style= "apsr", title="Results", align=TRUE, type="text")

eubadmargins <- margins(EUbad)
islambadmargins <- margins(islambad)
globlocalbadmargins <- margins(globlocalbad)
globoverbadmargins <- margins(globoverbad)

eubadmarginsplot <- summary(eubadmargins)
islambadmarginsplot <- summary(islambadmargins)
globlocalbadmarginsplot <- summary(globlocalbadmargins)
globoverbadmarginsplot <- summary(globoverbadmargins)



w <- ggplot(data = globoverbadmarginsplot) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Globalisation has bad effect on UK economy")
x <- ggplot(data = globlocalbadmarginsplot) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Globalisation has bad effect on local economy")
y <- ggplot(data = islambadmarginsplot) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Negative view of Islam")
z <- ggplot(data = eubadmarginsplot) +
  geom_point(aes(factor, AME)) +
  geom_errorbar(aes(x = factor, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("EU has bad impact on world")

ggarrange(w, x, y, z,  ncol=2, nrow=2)

settlementmarginaleffects <- dydx(ukglobal, globoverbad, "settlement")

