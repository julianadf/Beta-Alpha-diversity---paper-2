# Species-habitat association analysis using the "indicspecies" package
rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(lattice)
library(ggpubr)
library(car)
library(psych)
library(lme4)
library(glmulti)
library(visreg)
library(glmmTMB)
library(blmeco)
library(DHARMa)
library(visreg)
library(lmerTest)
library(interactions)
library(MuMIn)
library(arm)
library(jtools)
library(emmeans)
library(effects)
library(indicspecies)
library(SpadeR)
library(BiodiversityR)

### butterflies----
bf <- read.csv2(here("data", "Butterfly.dataframe.csv"))
abund.bf <- bf[,4:ncol(bf)]
habitat.type.bf <- bf$Transect_type

bf.associations <- multipatt(abund.bf, habitat.type.bf, func = "r.g", control = how(nperm=9999))
summary(bf.associations)

# Species accumulation curves - total
bf.data <- bf[,4:ncol(bf)]
bf.plot <- specaccum(bf.data)
plot(bf.plot)

# Species accumulation / landscape
#subset each landscape into its own df
NoPL.LRD <- bf %>% filter(Category == "NoPL.LRD") 
NoPL.HRD <- bf %>% filter(Category == "NoPL.HRD")
PL.LRD <- bf %>% filter(Category == "PL.LRD")
PL.HRD <- bf %>% filter(Category == "PL.HRD")

#calc species accumulation curve for each landscape
curve_NoPL.LRD = specaccum(NoPL.LRD[, 4:56], method = "random")
curve_NoPL.HRD = specaccum(NoPL.HRD[, 4:56], method = "random")
curve_PL.LRD = specaccum(PL.LRD[, 4:56], method = "random")
curve_PL.HRD = specaccum(PL.HRD[, 4:56], method = "random")

#plot curve_all first
plot(curve_NoPL.LRD)
#then plot the rest
plot(curve_NoPL.HRD, add = TRUE, col = 2) #col is COLOUR setting, so change it to something else if you want
plot(curve_PL.LRD, add = TRUE, col = 3)
plot(curve_PL.HRD, add = TRUE, col = 4) 

# Species accumulation / habitat
#its own df
Big_road <- bf %>% filter(Transect_type == "Big road") 
Small_road <- bf %>% filter(Transect_type == "Small road")
Pasture <- bf %>% filter(Transect_type == "Pasture")
Between_fields <- bf %>% filter(Transect_type == "Between fields")
Powerline <- bf %>% filter(Transect_type == "Powerline")

#calc species accumulation curve for each habitat
Big.road = specaccum(Big_road[, 4:56], method = "random")
curve_Small.road = specaccum(Small_road[, 4:56], method = "random")
curve_Pasture = specaccum(Pasture[, 4:56], method = "random")
curve_Between.fields = specaccum(Between_fields[, 4:56], method = "random")
curve_Powerline= specaccum(Powerline[, 4:56], method = "random")

#plot curve_all first
plot(bf.plot)
plot(curve_Big.road, add = TRUE, col = 6)
plot(curve_Small.road, add = TRUE, col = 2) #col is COLOUR setting, so change it to something else if you want
plot(curve_Pasture, add = TRUE, col = 3)
plot(curve_Between.fields, add = TRUE, col = 4) 
plot(curve_Powerline, add = TRUE, col = 5) 

# bf.roadverges <- bf %>% 
#   mutate(rv.density =ifelse(Category %in% c("PL.LRD", "NoPL.LRD"), 0,1))
# 
# RV.landscape <- bf.roadverges$rv.density
# 
# bf.associations <- multipatt(abund.bf, RV.landscape, func = "r.g", control = how(nperm=9999))
# summary(bf.associations)

### bumblebees ----
bb <- read.csv2(here("data", "Bumblebee.dataframe.csv"))
abund.bb <- bb[,4:ncol(bb)]
habitat.type.bb <- bb$Transect_type

bb.associations <- multipatt(abund.bb, habitat.type.bb, func = "r.g", control = how(nperm=9999))
summary(bb.associations)

# Species accumulation curves
bb.data <- bb[,4:ncol(bb)]
bb.plot <- specaccum(bb.data)
plot(bb.plot)


# Species accumulation / landscape
#subset each landscape into its own df
NoPL.LRD <- bb %>% filter(Category == "NoPL.LRD") 
NoPL.HRD <- bb %>% filter(Category == "NoPL.HRD")
PL.LRD <- bb %>% filter(Category == "PL.LRD")
PL.HRD <- bb %>% filter(Category == "PL.HRD")

#calc species accumulation curve for each landscape
curve_NoPL.LRD = specaccum(NoPL.LRD[, 4:56], method = "random")
curve_NoPL.HRD = specaccum(NoPL.HRD[, 4:56], method = "random")
curve_PL.LRD = specaccum(PL.LRD[, 4:56], method = "random")
curve_PL.HRD = specaccum(PL.HRD[, 4:56], method = "random")

#plot curve_all first
plot(curve_NoPL.LRD)
#then plot the rest
plot(curve_NoPL.HRD, add = TRUE, col = 2) #col is COLOUR setting, so change it to something else if you want
plot(curve_PL.LRD, add = TRUE, col = 3)
plot(curve_PL.HRD, add = TRUE, col = 4) 

# Species accumulation / habitat
#its own df
Big_road <- bb %>% filter(Transect_type == "Big road") 
Small_road <- bb %>% filter(Transect_type == "Small road")
Pasture <- bb %>% filter(Transect_type == "Pasture")
Between_fields <- bb %>% filter(Transect_type == "Between fields")
Powerline <- bb %>% filter(Transect_type == "Powerline")

#calc species accumulation curve for each habitat
bf_Big.road = specaccum(Big_road[, 4:22], method = "random")
bf_Small.road = specaccum(Small_road[, 4:22], method = "random")
bf_Pasture = specaccum(Pasture[, 4:22], method = "random")
bf_Between.fields = specaccum(Between_fields[, 4:22], method = "random")
bf_Powerline= specaccum(Powerline[, 4:22], method = "random")

#plot curve_all first
plot(bb.plot)
plot(curve_Big.road, add = TRUE, col = 6)
plot(curve_Small.road, add = TRUE, col = 2) #col is COLOUR setting, so change it to something else if you want
plot(curve_Pasture)#, add = TRUE, col = 3)
plot(curve_Between.fields, add = TRUE, col = 4) 
plot(curve_Powerline, add = TRUE, col = 5) 

# bb.roadverges <- bb %>% 
#   mutate(rv.density =ifelse(Category %in% c("PL.LRD", "NoPL.LRD"), 0,1))
# 
# RV.landscapebb <- bb.roadverges$rv.density
# 
# bb.associations <- multipatt(abund.bb, RV.landscapebb, func = "r.g", control = how(nperm=9999))
# summary(bb.associations)

### plants ----
pp <- read.csv2(here("data", "Plant.dataframe.csv"))
abund.pp <- pp[,4:ncol(pp)]
habitat.type.pp <- pp$Transect_type

pp.associations <- multipatt(abund.pp, habitat.type.pp, func = "r.g", control = how(nperm=9999))
summary(pp.associations)

pp.roadverges <- pp %>%
  mutate(rv.density =ifelse(Category %in% c("PL.LRD", "NoPL.LRD"), 0,1))

RV.landscapepp <- pp.roadverges$rv.density

pp.associations <- multipatt(abund.pp, RV.landscapepp, func = "r.g", control = how(nperm=9999))
summary(pp.associations)

# Species accumulation curves
pp.data <- pp[,4:ncol(pp)]
pp.plot <- specaccum(pp.data)
plot(pp.plot)

# Species accumulation / landscape
#subset each landscape into its own df
NoPL.LRD <- pp %>% filter(Category == "NoPL.LRD") 
NoPL.HRD <- pp %>% filter(Category == "NoPL.HRD")
PL.LRD <- pp %>% filter(Category == "PL.LRD")
PL.HRD <- pp %>% filter(Category == "PL.HRD")

#calc species accumulation curve for each habitat
curve_NoPL.LRD = specaccum(NoPL.LRD[, 4:136], method = "random")
curve_NoPL.HRD = specaccum(NoPL.HRD[, 4:136], method = "random")
curve_PL.LRD = specaccum(PL.LRD[, 4:136], method = "random")
curve_PL.HRD = specaccum(PL.HRD[, 4:136], method = "random")

#plot curve_all first
plot(pp.plot)
plot(curve_NoPL.LRD, add = TRUE, col = 5) #blue
plot(curve_NoPL.HRD, add = TRUE, col = 2) #red
plot(curve_PL.LRD, add = TRUE, col = 3) #green
plot(curve_PL.HRD, add = TRUE, col = 6) #pink

# Species accumulation / habitat
#its own df
Big_road <- pp %>% filter(Transect_type == "Big road") 
Small_road <- pp %>% filter(Transect_type == "Small road")
Pasture <- pp %>% filter(Transect_type == "Pasture")
Between_fields <- pp %>% filter(Transect_type == "Between fields")
Powerline <- pp %>% filter(Transect_type == "Powerline")

#calc species accumulation curve for each habitat
curve_Big.road = specaccum(Big_road[, 4:136], method = "random")
curve_Small.road = specaccum(Small_road[, 4:136], method = "random")
curve_Pasture = specaccum(Pasture[, 4:136], method = "random")
curve_Between.fields = specaccum(Between_fields[, 4:136], method = "random")
curve_Powerline= specaccum(Powerline[, 4:136], method = "random")

#plot curve_all first
plot(pp.plot)
plot(curve_Big.road, add = TRUE, col = 5) #blue
plot(curve_Small.road, add = TRUE, col = 2) #red
plot(curve_Pasture, add = TRUE, col = 3) #green
plot(curve_Between.fields, add = TRUE, col = 4) #dark blue
plot(curve_Powerline, add = TRUE, col = 6) 
