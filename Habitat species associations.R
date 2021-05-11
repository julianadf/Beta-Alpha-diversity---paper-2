# 
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

### butterflies ###
bf <- read.csv2(here("data", "Butterfly.dataframe.csv"))
abund.bf <- bf[,4:ncol(bf)]
habitat.type.bf <- bf$Transect_type

bf.associations <- multipatt(abund.bf, habitat.type.bf, func = "r.g", control = how(nperm=9999))
summary(bf.associations)

### bumblebees ###
bb <- read.csv2(here("data", "Bumblebee.dataframe.csv"))
abund.bb <- bb[,4:ncol(bb)]
habitat.type.bb <- bb$Transect_type

bb.associations <- multipatt(abund.bb, habitat.type.bb, func = "r.g", control = how(nperm=9999))
summary(bb.associations)

### plants ###
pp <- read.csv2(here("data", "Plant.dataframe.csv"))
abund.pp <- pp[,4:ncol(pp)]
habitat.type.pp <- pp$Transect_type

pp.associations <- multipatt(abund.pp, habitat.type.pp, func = "r.g", control = how(nperm=9999))
summary(pp.associations)
