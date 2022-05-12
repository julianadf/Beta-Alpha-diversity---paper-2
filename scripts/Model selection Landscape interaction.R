# By: Juliana Dániel Ferreira
# Model selection species richness habitat-wise analyses

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
library(cowplot)
library(png)
library(magick)
library(ggeffects)
library(effects)

# Load data
data.1 <- read.csv2(here("data", "database.corrected.csv"))
str(data.1)

data <- data.1 %>% 
  mutate(PL=as.factor(PL)) %>% 
  mutate(RD=as.factor(RD)) %>% 
  mutate(Transect_type = ifelse(Transect_type=="Between fields", "Field border", Transect_type)) %>% 
  mutate(Transect_type = ifelse(Transect_type=="Powerline", "Power line", Transect_type)) %>% 
  mutate(Transect_type = as.factor(Transect_type))
data
str(data)
head(data)


# Model selection: Butterflies ------------
bf.glmm <- glmer(bf.rich ~ Transect_type * TUVA * PL * RD + (1 | Landscape),
                family= poisson,
                na.action= na.fail,
                data= data)
summary(bf.glmm)
plot(bf.glmm)
anova(bf.glmm, type=3)
options(contrasts = c("contr.sum","contr.poly"))
options('contrasts')
car::Anova(bf.glmm, type= "III")
aicc(bf.glmm)  # AICc = 741.5434
vif(bf.glmm)
dispersion_glmer(bf.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- bf.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Model averaging
bf.glmm <- standardize(bf.glmm) # yielded exact same output
MuMIn.bf <- dredge(bf.glmm)
MA.bf <- get.models(MuMIn.bf, subset=delta<4)
summary(model.avg(MA.bf))
importance(summary(model.avg(MA.bf)))

# Best model:
bf.glmm <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)
summary(bf.glmm)
plot(bf.glmm)
plot(allEffects(bf.glmm))
visreg(bf.glmm, scale="response")
tab_model(bf.glmm, show.=TRUE)
anova(bf.glmm)
Anova(bf.glmm, type="III")
options(contrasts = c("contr.sum","contr.poly"))
options('contrasts')
car::Anova(bf.glmm, type= "III")
aicc(bf.glmm)  # AICc = 741.5434
vif(bf.glmm)
dispersion_glmer(bf.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- bf.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()

# Compare least-square means of each habitat type
bf.em <- emmeans(bf.glmm, ~ Transect_type, type="response")
bf.emms<- emmeans(bf.glmm, "Transect_type")
pairs(bf.emms)
plot(bf.emms, comparisons=TRUE)

# Adam Flöhr:
em.bf <- as.data.frame(bf.em)
ggplot(em.bf, aes(Transect_type, rate)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  ylab("") +
  coord_cartesian(ylim = c(0,21)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

em <- multcomp::cld(bf.em, Letters = letters)
em <- as.data.frame(em)
butterflies <- ggplot(em, aes(Transect_type, rate)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +#, width = 0.1) +
  coord_cartesian(ylim = c(0,21)) +
  geom_text(aes(y = 21.1, label = .group)) +
  #xlab("Habitat type") +
  #ylab("Number of butterfly species")+
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        text = element_text(size=16))


# Nice boxplots but not the model:
jitter <- ggplot(data, aes(x=Transect_type , y=bf.rich, col=Transect_type))
jitter + scale_x_discrete() +
  geom_boxplot(outlier.colour = NULL, position = "dodge") +
  xlab("Habitat") + ylab("Number of butterfly species") +
  theme(axis.line = element_line(size = 0.5), legend.position = "none")

 # --
bf.glmm2 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL, 
                  family= poisson,
                  data= data)
summary(bf.glmm2)
anova(bf.glmm2)
Anova(bf.glmm2, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm2)  # AICc = 743.3977
vif(bf.glmm2)
dispersion_glmer(bf.glmm2) #to test if overdispersed, has to be > 1.4
mod_dharma2 <- bf.glmm2 %>% simulateResiduals(n=1000)
plot(mod_dharma2)
mod_dharma2 %>% testDispersion()

# --
bf.glmm3 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm3)
anova(bf.glmm3)
Anova(bf.glmm3, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm3)  # AICc = 742.9732
vif(bf.glmm3)
dispersion_glmer(bf.glmm3) #to test if overdispersed, has to be > 1.4
mod_dharma3 <- bf.glmm3 %>% simulateResiduals(n=1000)
plot(mod_dharma3)
mod_dharma3 %>% testDispersion()

# --
bf.glmm4 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm4)
anova(bf.glmm4)
Anova(bf.glmm4, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm4)  # AICc = 741.8303
vif(bf.glmm4)
dispersion_glmer(bf.glmm4) #to test if overdispersed, has to be > 1.4
mod_dharma4 <- bf.glmm4 %>% simulateResiduals(n=1000)
plot(mod_dharma4)
mod_dharma4 %>% testDispersion()


# Plot figure
visreg(bf.glmm4)

# --
bf.glmm5 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm5)
anova(bf.glmm5)
Anova(bf.glmm5, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm5)  # AICc = 745.3296
vif(bf.glmm5)
dispersion_glmer(bf.glmm5) #to test if overdispersed, has to be > 1.4
mod_dharma5 <- bf.glmm5 %>% simulateResiduals(n=1000)
plot(mod_dharma5)
mod_dharma5 %>% testDispersion()

# --
bf.glmm6 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm6)
anova(bf.glmm6)
car::Anova(bf.glmm6, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm6)  # AICc = 744.0299
vif(bf.glmm6)
dispersion_glmer(bf.glmm6) #to test if overdispersed, has to be > 1.4
mod_dharma6 <- bf.glmm6 %>% simulateResiduals(n=1000)
plot(mod_dharma6)
mod_dharma6 %>% testDispersion()

# --
bf.glmm7 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm7)
anova(bf.glmm7)
car::Anova(bf.glmm7, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm7)  # AICc = 743.7721
vif(bf.glmm7)
dispersion_glmer(bf.glmm7) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- bf.glmm7 %>% simulateResiduals(n=1000)
plot(mod_dharma7)
mod_dharma7 %>% testDispersion()

# --
bf.glmm8 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm8)
anova(bf.glmm8)
car::Anova(bf.glmm8, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm8)  # AICc = 746.1891
vif(bf.glmm8)
dispersion_glmer(bf.glmm8) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- bf.glmm8 %>% simulateResiduals(n=1000)
plot(mod_dharma8)
mod_dharma8 %>% testDispersion()

# --
bf.glmm9 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:TUVA, 
                  family= poisson,
                  data= data)
summary(bf.glmm9)
anova(bf.glmm9)
car::Anova(bf.glmm9, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm9)  # AICc = 744.1987
vif(bf.glmm9)
dispersion_glmer(bf.glmm9) #to test if overdispersed, has to be > 1.4
mod_dharma9 <- bf.glmm9 %>% simulateResiduals(n=1000)
plot(mod_dharma9)
mod_dharma9 %>% testDispersion()

# --
bf.glmm10 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:PL, 
                  family= poisson,
                  data= data)
summary(bf.glmm10)
anova(bf.glmm10)
car::Anova(bf.glmm10, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm10)  # AICc = 745.5532
vif(bf.glmm10)
dispersion_glmer(bf.glmm10) #to test if overdispersed, has to be > 1.4
mod_dharma10 <- bf.glmm10 %>% simulateResiduals(n=1000)
plot(mod_dharma10)
mod_dharma10 %>% testDispersion()

# --
bf.glmm11 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD, 
                   family= poisson,
                   data= data)
summary(bf.glmm11)
anova(bf.glmm11)
car::Anova(bf.glmm11, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm11)  # AICc = 750.0956
vif(bf.glmm11)
dispersion_glmer(bf.glmm11) #to test if overdispersed, has to be > 1.4
mod_dharma11 <- bf.glmm11 %>% simulateResiduals(n=1000)
plot(mod_dharma11)
mod_dharma11 %>% testDispersion()

# --
bf.glmm12 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + 
                      RD:PL  , 
                   family= poisson,
                   data= data)
summary(bf.glmm12)
anova(bf.glmm12)
car::Anova(bf.glmm12, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm12)  # AICc = 
vif(bf.glmm12)
dispersion_glmer(bf.glmm12) #to test if overdispersed, has to be > 1.4
mod_dharma12 <- bf.glmm12 %>% simulateResiduals(n=1000)
plot(mod_dharma12)
mod_dharma12 %>% testDispersion()
visreg(bf.glmm12, scale="response")

# --
bf.glmm0 <- glmer(bf.rich ~ 1 + (1 | Landscape), 
                  family= poisson,
                  data= data)
summary(bf.glmm0)
anova(bf.glmm0)
Anova(bf.glmm0, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm0)  # AICc = 792.2822
vif(bf.glmm0)
dispersion_glmer(bf.glmm0) #to test if overdispersed, has to be > 1.4
mod_dharma0 <- bf.glmm0 %>% simulateResiduals(n=1000)
plot(mod_dharma0)
mod_dharma0 %>% testDispersion()

# Model selection: Bumblebees ------------
bb.glmm <- glmer(bb.rich ~ Transect_type * TUVA * PL * RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)
summary(bb.glmm)
anova(bb.glmm)
plot(bb.glmm)
options(contrasts = c("contr.sum","contr.poly"))
car::Anova(bb.glmm, type="III")
aicc(bb.glmm)  # AICc = 599.883
vif(bb.glmm)
dispersion_glmer(bb.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- bb.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Model averaging:
bb.glmm <- standardize(bb.glmm)
MuMIn.bb <- dredge(bb.glmm)
MA.bb <- get.models(MuMIn.bb, subset=delta<4)
summary(model.avg(MA.bb))
importance(summary(model.avg(MA.bb)))

visreg(bb.glmm, scale = "response")

# Best model:
bb.glmm <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)
summary(bb.glmm)
anova(bb.glmm)
plot(allEffects(bb.glmm))
plot(bb.glmm)
options(contrasts = c("contr.sum","contr.poly"))
car::Anova(bb.glmm, type="III")
aicc(bb.glmm)  # AICc = 599.883
vif(bb.glmm)
dispersion_glmer(bb.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- bb.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Compare least-square means of each habitat type
bb.em <- emmeans(bb.glmm, ~ Transect_type, type="response")
bb.emms<- emmeans(bb.glmm, "Transect_type")
pairs(bb.emms)
plot(bb.emms, comparisons=TRUE)

# Adam Flöhr:
em.bb <- as.data.frame(bb.em)
ggplot(em.bb, aes(Transect_type, rate)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  ylab("") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

em.bb <- multcomp::cld(bb.em, Letters = letters)
em.bb <- as.data.frame(em.bb)
bumblebees <- ggplot(em.bb, aes(Transect_type, rate)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +#, width = 0.1) +
  geom_text(aes(y = 21.1, label = .group)) +
  coord_cartesian(ylim = c(0,21)) +
  #xlab("Habitat type") +
  ylab("Species richness")+
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 90), axis.title.x = element_blank(), 
        text = element_text(size=16))

# Better ggplot:

jitter.bb <- ggplot(data, aes(x=Transect_type , y=bb.rich, col=Transect_type)) 
jitter.bb + scale_x_discrete() +
  geom_boxplot(outlier.colour = NULL, position = "dodge") +  
  xlab("Habitat") + ylab("Number of bumblebee species") + 
  theme(axis.line = element_line(size = 0.5), legend.position = "none") 

# --
bb.glmm2 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL, 
                  family= poisson,
                  data= data)
summary(bb.glmm2)
anova(bb.glmm2)
car::Anova(bb.glmm2, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm2)  # AICc = 602.1567
vif(bb.glmm2)
dispersion_glmer(bb.glmm2) #to test if overdispersed, has to be > 1.4
mod_dharma2 <- bb.glmm2 %>% simulateResiduals(n=1000)
plot(mod_dharma2)
mod_dharma2 %>% testDispersion()

# --
bb.glmm3 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(bb.glmm3)
anova(bb.glmm3)
car::Anova(bb.glmm3, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm3)  # AICc = 601.3386
vif(bb.glmm3)
dispersion_glmer(bb.glmm3) #to test if overdispersed, has to be > 1.4
mod_dharma3 <- bb.glmm3 %>% simulateResiduals(n=1000)
plot(mod_dharma3)
mod_dharma3 %>% testDispersion()

# --
bb.glmm4 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + PL:RD, 
                  family= poisson,
                  data= data)
summary(bb.glmm4)
anova(bb.glmm4)
car::Anova(bb.glmm4, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm4)  # AICc = 601.9572
vif(bb.glmm4)
dispersion_glmer(bb.glmm4) #to test if overdispersed, has to be > 1.4
mod_dharma4 <- bb.glmm4 %>% simulateResiduals(n=1000)
plot(mod_dharma4)
mod_dharma4 %>% testDispersion()

# --
bb.glmm5 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(bb.glmm5)
anova(bb.glmm5)
car::Anova(bb.glmm5, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm5)  # AICc = 603.606
vif(bb.glmm5)
dispersion_glmer(bb.glmm5) #to test if overdispersed, has to be > 1.4
mod_dharma5 <- bb.glmm5 %>% simulateResiduals(n=1000)
plot(mod_dharma5)
mod_dharma5 %>% testDispersion()

# --
bb.glmm6 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + PL:RD, 
                  family= poisson,
                  data= data)
summary(bb.glmm6)
anova(bb.glmm6)
car::Anova(bb.glmm6, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm6)  # AICc = 604.3132
vif(bb.glmm6)
dispersion_glmer(bb.glmm6) #to test if overdispersed, has to be > 1.4
mod_dharma6 <- bb.glmm6 %>% simulateResiduals(n=1000)
plot(mod_dharma6)
mod_dharma6 %>% testDispersion()

# --
bb.glmm7 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(bb.glmm7)
anova(bb.glmm7)
car::Anova(bb.glmm7, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm7)  # AICc = 603.6334
vif(bb.glmm7)
dispersion_glmer(bb.glmm7) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- bb.glmm7 %>% simulateResiduals(n=1000)
plot(mod_dharma7)
mod_dharma7 %>% testDispersion()

# --
bb.glmm8 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(bb.glmm8)
anova(bb.glmm8)
car::Anova(bb.glmm8, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm8)  # AICc = 605.9239
vif(bb.glmm8)
dispersion_glmer(bb.glmm8) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- bb.glmm8 %>% simulateResiduals(n=1000)
plot(mod_dharma8)
mod_dharma8 %>% testDispersion()

# --
bb.glmm9 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:TUVA, 
                  family= poisson,
                  data= data)
summary(bb.glmm9)
anova(bb.glmm9)
car::Anova(bb.glmm9, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm9)  # AICc = 606.0112
vif(bb.glmm9)
dispersion_glmer(bb.glmm9) #to test if overdispersed, has to be > 1.4
mod_dharma9 <- bb.glmm9 %>% simulateResiduals(n=1000)
plot(mod_dharma9)
mod_dharma9 %>% testDispersion()

# --
bb.glmm10 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:PL, 
                   family= poisson,
                   data= data)
summary(bb.glmm10)
anova(bb.glmm10)
car::Anova(bb.glmm10, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm10)  # AICc = 605.8733
vif(bb.glmm10)
dispersion_glmer(bb.glmm10) #to test if overdispersed, has to be > 1.4
mod_dharma10 <- bb.glmm10 %>% simulateResiduals(n=1000)
plot(mod_dharma10)
mod_dharma10 %>% testDispersion()

# --
bb.glmm11 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD, 
                   family= poisson,
                   data= data)
summary(bb.glmm11)
anova(bb.glmm11)
car::Anova(bb.glmm11, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm11)  # AICc = 604.0637
vif(bb.glmm11)
dispersion_glmer(bb.glmm11) #to test if overdispersed, has to be > 1.4
mod_dharma11 <- bb.glmm11 %>% simulateResiduals(n=1000)
plot(mod_dharma11)
mod_dharma11 %>% testDispersion()

# --
bb.glmm12 <- glmer(bb.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD + Transect_type:TUVA, 
                   family= poisson,
                   data= data)
summary(bb.glmm12)
anova(bb.glmm12)
car::Anova(bb.glmm12, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm12)  # AICc = 
vif(bb.glmm12)
dispersion_glmer(bb.glmm12) #to test if overdispersed, has to be > 1.4
mod_dharma12 <- bb.glmm12 %>% simulateResiduals(n=1000)
plot(mod_dharma12)
mod_dharma12 %>% testDispersion()

# --
bb.glmm0 <- glmer(bb.rich ~ 1 + (1 | Landscape), 
                  family= poisson,
                  data= data)
summary(bb.glmm0)
anova(bb.glmm0)
car::Anova(bb.glmm0, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bb.glmm0)  # AICc = 595.3937
vif(bb.glmm0)
dispersion_glmer(bb.glmm0) #to test if overdispersed, has to be > 1.4
mod_dharma0 <- bb.glmm0 %>% simulateResiduals(n=1000)
plot(mod_dharma0)
mod_dharma0 %>% testDispersion()


# Model selection: Plants ------------
pp.spglmm <- glmer(pp.special ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)

# Need to remove missing sample to be able to run the model averaging
data.pp <- data %>% na.exclude(data)
pp.glmm <- glmer(pp.rich ~ Transect_type * TUVA * PL * RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data.pp)

summary(pp.glmm)
anova(pp.glmm)
car::Anova(pp.glmm, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm)  # AICc = 961.7039
vif(pp.glmm)
dispersion_glmer(pp.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- pp.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Model averaging:
MuMIn.pp <- dredge(pp.glmm)
MA.pp <- get.models(MuMIn.pp, subset=delta<4)
summary(model.avg(MA.pp))
importance(summary(model.avg(MA.pp)))

visreg(pp.glmm, scale="response")
# Best model:
# Need to remove missing sample to be able to run the model averaging
data.pp <- data[,1:24] %>% na.exclude(data)
pp.glmm <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data.pp)

summary(pp.glmm)
plot(allEffects(pp.glmm))
visreg(pp.glmm, scale="response")
pp.alfa <- effects::effect(term= "PL",  mod= pp.glmm)
summary(pp.alfa)
anova(pp.glmm)
car::Anova(pp.glmm, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm)  # AICc = 961.7039
vif(pp.glmm)
dispersion_glmer(pp.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- pp.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Compare least-square means of each habitat type
pp.em <- emmeans(pp.glmm, ~ Transect_type, type="response")
pp.emms<- emmeans(pp.glmm, "Transect_type")
pairs(pp.emms)
plot(pp.emms, comparisons=TRUE)


# Adam Flöhr:
em.pp <- as.data.frame(pp.em)
ggplot(em.pp, aes(Transect_type, rate)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  ylab("") +
  #coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

em.pp <- multcomp::cld(pp.em, Letters = letters)
em.pp <- as.data.frame(em.pp)
plants <- ggplot(em.pp, aes(Transect_type, rate)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +#, width = 0.1) +
  geom_text(aes(y = 21.1, label = .group)) +
  coord_cartesian(ylim = c(0,21)) +
  #xlab("Habitat type") +
  #ylab("Number of plant species")+
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 90), axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        text = element_text(size=16))

# Plot figure
pp.plot<-effect_plot(pp.glmm, pred = Transect_type, interval = TRUE, int.type = "confidence", int.width = .8,
                     outcome.scale = "response",
                     x.label = "Habitat", y.label = "Number of plant species",
                     line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
                     point.alpha = 0.3                    )
pp.plot + scale_x_discrete() +theme(axis.text=element_text(size = 15), legend.position = "none")

# --
pp.glmm2 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL, 
                  family= poisson,
                  data= data)
summary(pp.glmm2)
anova(pp.glmm2)
car::Anova(pp.glmm2, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm2)  # AICc = 963.9444
vif(pp.glmm2)
dispersion_glmer(pp.glmm2) #to test if overdispersed, has to be > 1.4
mod_dharma2 <- pp.glmm2 %>% simulateResiduals(n=1000)
plot(mod_dharma2)
mod_dharma2 %>% testDispersion()

# --
pp.glmm3 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm3)
anova(pp.glmm3)
car::Anova(pp.glmm3, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm3)  # AICc = 963.4168
vif(pp.glmm3)
dispersion_glmer(pp.glmm3) #to test if overdispersed, has to be > 1.4
mod_dharma3 <- pp.glmm3 %>% simulateResiduals(n=1000)
plot(mod_dharma3)
mod_dharma3 %>% testDispersion()

# --
pp.glmm4 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm4)
anova(pp.glmm4)
car::Anova(pp.glmm4, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm4)  # AICc = 963.8689
vif(pp.glmm4)
dispersion_glmer(pp.glmm4) #to test if overdispersed, has to be > 1.4
mod_dharma4 <- pp.glmm4 %>% simulateResiduals(n=1000)
plot(mod_dharma4)
mod_dharma4 %>% testDispersion()

# --
pp.glmm5 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm5)
anova(pp.glmm5)
car::Anova(pp.glmm5, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm5)  # AICc = 965.7573
vif(pp.glmm5)
dispersion_glmer(pp.glmm5) #to test if overdispersed, has to be > 1.4
mod_dharma5 <- pp.glmm5 %>% simulateResiduals(n=1000)
plot(mod_dharma5)
mod_dharma5 %>% testDispersion()

# --
pp.glmm6 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm6)
anova(pp.glmm6)
car::Anova(pp.glmm6, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm6)  # AICc = 966.1881
vif(pp.glmm6)
dispersion_glmer(pp.glmm6) #to test if overdispersed, has to be > 1.4
mod_dharma6 <- pp.glmm6 %>% simulateResiduals(n=1000)
plot(mod_dharma6)
mod_dharma6 %>% testDispersion()

# --
pp.glmm7 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm7)
anova(pp.glmm7)
car::Anova(pp.glmm7, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm7)  # AICc = 965.7231
vif(pp.glmm7)
dispersion_glmer(pp.glmm7) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- pp.glmm7 %>% simulateResiduals(n=1000)
plot(mod_dharma7)
mod_dharma7 %>% testDispersion()

# --
pp.glmm8 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm8)
anova(pp.glmm8)
car::Anova(pp.glmm8, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm8)  # AICc = 968.0945
vif(pp.glmm8)
dispersion_glmer(pp.glmm8) #to test if overdispersed, has to be > 1.4
mod_dharma8 <- pp.glmm8 %>% simulateResiduals(n=1000)
plot(mod_dharma8)
mod_dharma8 %>% testDispersion()

# --
pp.glmm9 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:TUVA, 
                  family= poisson,
                  data= data)
summary(pp.glmm9)
anova(pp.glmm9)
car::Anova(pp.glmm9, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm9)  # AICc = 965.715
vif(pp.glmm9)
dispersion_glmer(pp.glmm9) #to test if overdispersed, has to be > 1.4
mod_dharma9 <- pp.glmm9 %>% simulateResiduals(n=1000)
plot(mod_dharma9)
mod_dharma9 %>% testDispersion()

# --
pp.glmm10 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:PL, 
                   family= poisson,
                   data= data)
summary(pp.glmm10)
anova(pp.glmm10)
car::Anova(pp.glmm10, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm10)  # AICc = 964.1856
vif(pp.glmm10)
dispersion_glmer(pp.glmm10) #to test if overdispersed, has to be > 1.4
mod_dharma10 <- pp.glmm10 %>% simulateResiduals(n=1000)
plot(mod_dharma10)
mod_dharma10 %>% testDispersion()

# --
pp.glmm11 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD, 
                   family= poisson,
                   data= data)
summary(pp.glmm11)
anova(pp.glmm11)
car::Anova(pp.glmm11, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm11)  # AICc = 966.4854
vif(pp.glmm11)
dispersion_glmer(pp.glmm11) #to test if overdispersed, has to be > 1.4
mod_dharma11 <- pp.glmm11 %>% simulateResiduals(n=1000)
plot(mod_dharma11)
mod_dharma11 %>% testDispersion()

# --
pp.glmm12 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD + Transect_type:TUVA, 
                   family= poisson,
                   data= data)
summary(pp.glmm12)
anova(pp.glmm12)
car::Anova(pp.glmm12, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm12)  # AICc = 971.8976
vif(pp.glmm12)
dispersion_glmer(pp.glmm12) #to test if overdispersed, has to be > 1.4
mod_dharma12 <- pp.glmm12 %>% simulateResiduals(n=1000)
plot(mod_dharma12)
mod_dharma12 %>% testDispersion()

# --
pp.glmm0 <- glmer(pp.rich ~ 1 + (1 | Landscape), 
                  family= poisson,
                  data= data)
summary(pp.glmm0)
anova(pp.glmm0)
car::Anova(pp.glmm0, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm0)  # AICc = 995.1401
vif(pp.glmm0)
dispersion_glmer(pp.glmm0) #to test if overdispersed, has to be > 1.4
mod_dharma0 <- pp.glmm0 %>% simulateResiduals(n=1000)
plot(mod_dharma0)
mod_dharma0 %>% testDispersion()

# Plot butterflies, bumblebees and plants together ----

Figure.2 <- ggarrange(humla,maniola, plant, bumblebees, butterflies, plants, ncol=3, nrow=2, widths = c(1.15,1,1), heights = c(1,6))
Figure.2
# add silhouettes

humla <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/bumble.png")
maniola <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/butterfly.png")
plant <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/Fragaria.png")
ggdraw(Figure.2) + draw_image(humla, x = 0.27, y = 0.97, hjust = 1, vjust = 1, width = 0.12, height = 0.12) +
  draw_image(maniola, x = 0.593, y = 0.982, hjust = 1, vjust = 1, width = 0.13, height = 0.13) + 
  draw_image(plant, x = 0.9, y = 0.99, hjust = 1, vjust = 1, width = 0.13, height = 0.13)
  


# Model selection: Indicator plants (positive only) ------------
pi.glmm <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 data= data)
# Need to remove missing sample to be able to run the model averaging
data.pi <- data %>% na.exclude(data)
str(data.pi)
pi.glmm <- glmer(pi.rich ~ Transect_type * TUVA * PL * RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data.pi)

summary(pi.glmm)
anova(pi.glmm)
car::Anova(pi.glmm, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm)  # AICc = 
vif(pi.glmm)
dispersion_glmer(pi.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- pi.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Model averaging:
MuMIn.pi <- dredge(pi.glmm)
MA.pi <- get.models(MuMIn.pi, subset=delta<4)
summary(model.avg(MA.pi))
importance(summary(model.avg(MA.pi)))

# Best model
# Need to remove missing sample to be able to run the model averaging
data.pi <- data %>% na.exclude(data)
pi.glmm <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),
                 data= data.pi)

summary(pi.glmm)
anova(pi.glmm)
car::Anova(pi.glmm, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm)  # AICc = 466.8095
vif(pi.glmm)
dispersion_glmer(pi.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- pi.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()

# Compare least-square means of each habitat type
pi.emms<- emmeans(pi.glmm, "Transect_type")
pairs(pi.emms)
plot(pi.emms, comparisons=TRUE)

# Plot figure
# Transect type:
pi.plot<-cat_plot(pi.glmm, pred = Transect_type, interval = TRUE, int.type = "confidence", int.width = .8,
                     outcome.scale = "response",
                     x.label = "Habitat", y.label = "Number of indicator plant species",
                     line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
                     point.alpha = 0.3                    )
pi.plot + scale_x_discrete() + theme(axis.text=element_text(size = 15), legend.position = "none")
# TUVA:
visreg(pi.glmm, scale="response", "TUVA", line=list(col="black"), rug=FALSE, bty="n")
# OR:
pi.plot<-effect_plot(pi.glmm, pred = TUVA, interval = TRUE, int.type = "confidence", int.width = .8,
                     outcome.scale = "response",
                     x.label = "Area TUVA", y.label = "Number of indicator plant species",
                     line.thickness = 1, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
                     point.alpha = 0.3                    )
pi.plot + theme(axis.text=element_text(size = 15), legend.position = "none")

# --
pi.glmm2 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL, 
                  family= poisson,
                  data= data)
summary(pi.glmm2)
anova(pi.glmm2)
car::Anova(pi.glmm2, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm2)  # AICc = 468.902
vif(pi.glmm2)
dispersion_glmer(pi.glmm2) #to test if overdispersed, has to be > 1.4
mod_dharma2 <- pi.glmm2 %>% simulateResiduals(n=1000)
plot(mod_dharma2)
mod_dharma2 %>% testDispersion()

# --
pi.glmm3 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(pi.glmm3)
anova(pi.glmm3)
car::Anova(pi.glmm3, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm3)  # AICc = 468.902
vif(pi.glmm3)
dispersion_glmer(pi.glmm3) #to test if overdispersed, has to be > 1.4
mod_dharma3 <- pi.glmm3 %>% simulateResiduals(n=1000)
plot(mod_dharma3)
mod_dharma3 %>% testDispersion()

# --
pi.glmm4 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + PL:RD, 
                  family= poisson,
                  data= data)
summary(pi.glmm4)
anova(pi.glmm4)
car::Anova(pi.glmm4, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm4)  # AICc = 467.532
vif(pi.glmm4)
dispersion_glmer(pi.glmm4) #to test if overdispersed, has to be > 1.4
mod_dharma4 <- pi.glmm4 %>% simulateResiduals(n=1000)
plot(mod_dharma4)
mod_dharma4 %>% testDispersion()

# --
pi.glmm5 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(pi.glmm5)
anova(pi.glmm5)
car::Anova(pi.glmm5, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm5)  # AICc = 469.7707
vif(pi.glmm5)
dispersion_glmer(pi.glmm5) #to test if overdispersed, has to be > 1.4
mod_dharma5 <- pi.glmm5 %>% simulateResiduals(n=1000)
plot(mod_dharma5)
mod_dharma5 %>% testDispersion()

# --
pi.glmm6 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + PL:RD, 
                  family= poisson,
                  data= data)
summary(pi.glmm6)
anova(pi.glmm6)
car::Anova(pi.glmm6, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm6)  # AICc = 469.8772
vif(pi.glmm6)
dispersion_glmer(pi.glmm6) #to test if overdispersed, has to be > 1.4
mod_dharma6 <- pi.glmm6 %>% simulateResiduals(n=1000)
plot(mod_dharma6)
mod_dharma6 %>% testDispersion()

# --
pi.glmm7 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(pi.glmm7)
anova(pi.glmm7)
car::Anova(pi.glmm7, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm7)  # AICc = 468.8535
vif(pi.glmm7)
dispersion_glmer(pi.glmm7) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- pi.glmm7 %>% simulateResiduals(n=1000)
plot(mod_dharma7)
mod_dharma7 %>% testDispersion()

# --
pi.glmm8 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(pi.glmm8)
anova(pi.glmm8)
car::Anova(pi.glmm8, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm8)  # AICc = 470.2239
vif(pi.glmm8)
dispersion_glmer(pi.glmm8) #to test if overdispersed, has to be > 1.4
mod_dharma8 <- pi.glmm8 %>% simulateResiduals(n=1000)
plot(mod_dharma8)
mod_dharma8 %>% testDispersion()

# --
pi.glmm9 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:TUVA, 
                  family= poisson,
                  data= data)
summary(pi.glmm9)
anova(pi.glmm9)
car::Anova(pi.glmm9, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm9)  # AICc = 469.4743
vif(pi.glmm9)
dispersion_glmer(pi.glmm9) #to test if overdispersed, has to be > 1.4
mod_dharma9 <- pi.glmm9 %>% simulateResiduals(n=1000)
plot(mod_dharma9)
mod_dharma9 %>% testDispersion()

# --
pi.glmm10 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:PL, 
                   family= poisson,
                   data= data)
summary(pi.glmm10)
anova(pi.glmm10)
car::Anova(pi.glmm10, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm10)  # AICc = 472.3182
vif(pi.glmm10)
dispersion_glmer(pi.glmm10) #to test if overdispersed, has to be > 1.4
mod_dharma10 <- pi.glmm10 %>% simulateResiduals(n=1000)
plot(mod_dharma10)
mod_dharma10 %>% testDispersion()

# --
pi.glmm11 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD, 
                   family= poisson,
                   data= data)
summary(pi.glmm11)
anova(pi.glmm11)
car::Anova(pi.glmm11, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm11)  # AICc = 473.3584
vif(pi.glmm11)
dispersion_glmer(pi.glmm11) #to test if overdispersed, has to be > 1.4
mod_dharma11 <- pi.glmm11 %>% simulateResiduals(n=1000)
plot(mod_dharma11)
mod_dharma11 %>% testDispersion()

# --
pi.glmm12 <- glmer(pi.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD + Transect_type:TUVA, 
                   family= poisson,
                   data= data)
summary(pi.glmm12)
anova(pi.glmm12)
car::Anova(pi.glmm12, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pi.glmm12)  # AICc = 477.4417
vif(pi.glmm12)
dispersion_glmer(pi.glmm12) #to test if overdispersed, has to be > 1.4
mod_dharma12 <- pi.glmm12 %>% simulateResiduals(n=1000)
plot(mod_dharma12)
mod_dharma12 %>% testDispersion()


# --
pp.glmm0 <- glmer(pp.rich ~ 1 + (1 | Landscape), 
                  family= poisson,
                  data= data)
summary(pp.glmm0)
anova(pp.glmm0)
Anova(pp.glmm0, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm0)  # AICc = 
vif(pp.glmm0)
dispersion_glmer(pp.glmm0) #to test if overdispersed, has to be > 1.4
mod_dharma0 <- pp.glmm0 %>% simulateResiduals(n=1000)
plot(mod_dharma0)
mod_dharma0 %>% testDispersion()
# Model selection: grassland specialist butterflies -----
bf.spglmm <- glmer(bf.special ~ Transect_type * TUVA * PL * RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)
summary(bf.spglmm)
plot(bf.spglmm)
anova(bf.spglmm)
options(contrasts = c("contr.sum","contr.poly"))
options('contrasts')
car::Anova(bf.spglmm, type= "III")
aicc(bf.spglmm)  # AICc = 741.5434
vif(bf.spglmm)
dispersion_glmer(bf.spglmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- bf.spglmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Model averaging
bf.spglmm <- standardize(bf.glmm) # yielded exact same output
MuMIn.spbf <- dredge(bf.spglmm)
MA.spbf <- get.models(MuMIn.spbf, subset=delta<4)
summary(model.avg(MA.spbf))
importance(summary(model.avg(MA.spbf)))

# Best model:
spbf.glmm <- glmer(bf.special ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)
summary(spbf.glmm)
plot(spbf.glmm)
anova(spbf.glmm)
options(contrasts = c("contr.sum","contr.poly"))
options('contrasts')
car::Anova(spbf.glmm, type= "III")
aicc(spbf.glmm)  # AICc = 551.5809
vif(spbf.glmm)
dispersion_glmer(spbf.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- spbf.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Compare least-square means of each habitat type
spbf.emms<- emmeans(spbf.glmm, "Transect_type")
pairs(spbf.emms)
plot(spbf.emms, comparisons=TRUE)
# Plot figure
spbf.plot<-effect_plot(spbf.glmm, pred = Transect_type, interval = TRUE, int.type = "confidence", int.width = .8,
                     outcome.scale = "response",
                     x.label = "Habitat", y.label = "Number of butterfly species",
                     line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
                     point.alpha = 0.3                    )
spbf.plot + scale_x_discrete() +theme(axis.text=element_text(size =15), legend.position = "none")
# theme_update(text = element_text(size=20))
visreg(bf.glmm, "Transect_type", scale="response", band=TRUE, partial=FALSE, rug=FALSE, type="contrast")#, gg=TRUE)

# --
spbf.glmm2 <- glmer(bf.special ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL, 
                  family= poisson,
                  data= data)
summary(spbf.glmm2)
anova(bf.glmm2)
Anova(bf.glmm2, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm2)  # AICc = 743.3977
vif(bf.glmm2)
dispersion_glmer(bf.glmm2) #to test if overdispersed, has to be > 1.4
mod_dharma2 <- bf.glmm2 %>% simulateResiduals(n=1000)
plot(mod_dharma2)
mod_dharma2 %>% testDispersion()

# --
bf.glmm3 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm3)
anova(bf.glmm3)
Anova(bf.glmm3, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm3)  # AICc = 742.9732
vif(bf.glmm3)
dispersion_glmer(bf.glmm3) #to test if overdispersed, has to be > 1.4
mod_dharma3 <- bf.glmm3 %>% simulateResiduals(n=1000)
plot(mod_dharma3)
mod_dharma3 %>% testDispersion()

# --
bf.glmm4 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm4)
anova(bf.glmm4)
Anova(bf.glmm4, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm4)  # AICc = 741.8303
vif(bf.glmm4)
dispersion_glmer(bf.glmm4) #to test if overdispersed, has to be > 1.4
mod_dharma4 <- bf.glmm4 %>% simulateResiduals(n=1000)
plot(mod_dharma4)
mod_dharma4 %>% testDispersion()

# Plot figure
visreg(bf.glmm4)

# --
bf.glmm5 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm5)
anova(bf.glmm5)
Anova(bf.glmm5, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm5)  # AICc = 745.3296
vif(bf.glmm5)
dispersion_glmer(bf.glmm5) #to test if overdispersed, has to be > 1.4
mod_dharma5 <- bf.glmm5 %>% simulateResiduals(n=1000)
plot(mod_dharma5)
mod_dharma5 %>% testDispersion()

# --
bf.glmm6 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm6)
anova(bf.glmm6)
car::Anova(bf.glmm6, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm6)  # AICc = 744.0299
vif(bf.glmm6)
dispersion_glmer(bf.glmm6) #to test if overdispersed, has to be > 1.4
mod_dharma6 <- bf.glmm6 %>% simulateResiduals(n=1000)
plot(mod_dharma6)
mod_dharma6 %>% testDispersion()

# --
bf.glmm7 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm7)
anova(bf.glmm7)
car::Anova(bf.glmm7, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm7)  # AICc = 743.7721
vif(bf.glmm7)
dispersion_glmer(bf.glmm7) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- bf.glmm7 %>% simulateResiduals(n=1000)
plot(mod_dharma7)
mod_dharma7 %>% testDispersion()

# --
bf.glmm8 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(bf.glmm8)
anova(bf.glmm8)
car::Anova(bf.glmm8, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm8)  # AICc = 746.1891
vif(bf.glmm8)
dispersion_glmer(bf.glmm8) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- bf.glmm8 %>% simulateResiduals(n=1000)
plot(mod_dharma8)
mod_dharma8 %>% testDispersion()

# --
bf.glmm9 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:TUVA, 
                  family= poisson,
                  data= data)
summary(bf.glmm9)
anova(bf.glmm9)
car::Anova(bf.glmm9, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm9)  # AICc = 744.1987
vif(bf.glmm9)
dispersion_glmer(bf.glmm9) #to test if overdispersed, has to be > 1.4
mod_dharma9 <- bf.glmm9 %>% simulateResiduals(n=1000)
plot(mod_dharma9)
mod_dharma9 %>% testDispersion()

# --
bf.glmm10 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:PL, 
                   family= poisson,
                   data= data)
summary(bf.glmm10)
anova(bf.glmm10)
car::Anova(bf.glmm10, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm10)  # AICc = 745.5532
vif(bf.glmm10)
dispersion_glmer(bf.glmm10) #to test if overdispersed, has to be > 1.4
mod_dharma10 <- bf.glmm10 %>% simulateResiduals(n=1000)
plot(mod_dharma10)
mod_dharma10 %>% testDispersion()

# --
bf.glmm11 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD, 
                   family= poisson,
                   data= data)
summary(bf.glmm11)
anova(bf.glmm11)
car::Anova(bf.glmm11, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm11)  # AICc = 750.0956
vif(bf.glmm11)
dispersion_glmer(bf.glmm11) #to test if overdispersed, has to be > 1.4
mod_dharma11 <- bf.glmm11 %>% simulateResiduals(n=1000)
plot(mod_dharma11)
mod_dharma11 %>% testDispersion()

# --
bf.glmm12 <- glmer(bf.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + 
                     RD:PL  , 
                   family= poisson,
                   data= data)
summary(bf.glmm12)
anova(bf.glmm12)
car::Anova(bf.glmm12, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm12)  # AICc = 
vif(bf.glmm12)
dispersion_glmer(bf.glmm12) #to test if overdispersed, has to be > 1.4
mod_dharma12 <- bf.glmm12 %>% simulateResiduals(n=1000)
plot(mod_dharma12)
mod_dharma12 %>% testDispersion()
visreg(bf.glmm12, scale="response")

# --
bf.glmm0 <- glmer(bf.rich ~ 1 + (1 | Landscape), 
                  family= poisson,
                  data= data)
summary(bf.glmm0)
anova(bf.glmm0)
Anova(bf.glmm0, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(bf.glmm0)  # AICc = 792.2822
vif(bf.glmm0)
dispersion_glmer(bf.glmm0) #to test if overdispersed, has to be > 1.4
mod_dharma0 <- bf.glmm0 %>% simulateResiduals(n=1000)
plot(mod_dharma0)
mod_dharma0 %>% testDispersion()
# Model selection: grassland specialist plants ----
ppsp.glmm <- glmer(pp.special ~ Transect_type * TUVA * PL * RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)

summary(ppsp.glmm)
anova(ppsp.glmm)
car::Anova(ppsp.glmm, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(ppsp.glmm)  # AICc = 
vif(ppsp.glmm)
dispersion_glmer(ppsp.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- ppsp.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Model averaging:
MuMIn.ppsp <- dredge(ppsp.glmm)
MA.ppsp <- get.models(MuMIn.ppsp, subset=delta<4)
summary(model.avg(MA.ppsp))
importance(summary(model.avg(MA.ppsp)))

visreg(ppsp.glmm, scale="response")
# Best model:
ppsp.glmm <- glmer(pp.special ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                 family= poisson,
                 na.action= na.fail,
                 data= data)

summary(ppsp.glmm)
anova(ppsp.glmm)
car::Anova(ppsp.glmm, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(ppsp.glmm)  # AICc = 
vif(ppsp.glmm)
dispersion_glmer(ppsp.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- ppsp.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()

# Compare least-square means of each habitat type
ppsp.emms<- emmeans(ppsp.glmm, "Transect_type")
pairs(ppsp.emms)
plot(ppsp.emms, comparisons=TRUE)

# Plot figure
ppsp.plot<-effect_plot(ppsp.glmm, pred = Transect_type, interval = TRUE, int.type = "confidence", int.width = .8,
                       outcome.scale = "response",
                       x.label = "Habitat ", y.label = "Number of specialist plant species",
                       line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
                       point.alpha = 0.3                    )
ppsp.plot + scale_x_discrete() +theme(axis.text=element_text(size = 15), legend.position = "none")

ppsp.plot<-effect_plot(ppsp.glmm, pred = PL, interval = TRUE, int.type = "confidence", int.width = .8,
                     outcome.scale = "response",
                     x.label = "Power line ", y.label = "Number of specialist plant species",
                     line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
                     point.alpha = 0.3                    )
ppsp.plot + scale_x_discrete() +theme(axis.text=element_text(size = 15), legend.position = "none")



# --
pp.glmm2 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL, 
                  family= poisson,
                  data= data)
summary(pp.glmm2)
anova(pp.glmm2)
car::Anova(pp.glmm2, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm2)  # AICc = 963.9444
vif(pp.glmm2)
dispersion_glmer(pp.glmm2) #to test if overdispersed, has to be > 1.4
mod_dharma2 <- pp.glmm2 %>% simulateResiduals(n=1000)
plot(mod_dharma2)
mod_dharma2 %>% testDispersion()

# --
pp.glmm3 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm3)
anova(pp.glmm3)
car::Anova(pp.glmm3, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm3)  # AICc = 963.4168
vif(pp.glmm3)
dispersion_glmer(pp.glmm3) #to test if overdispersed, has to be > 1.4
mod_dharma3 <- pp.glmm3 %>% simulateResiduals(n=1000)
plot(mod_dharma3)
mod_dharma3 %>% testDispersion()

# --
pp.glmm4 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm4)
anova(pp.glmm4)
car::Anova(pp.glmm4, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm4)  # AICc = 963.8689
vif(pp.glmm4)
dispersion_glmer(pp.glmm4) #to test if overdispersed, has to be > 1.4
mod_dharma4 <- pp.glmm4 %>% simulateResiduals(n=1000)
plot(mod_dharma4)
mod_dharma4 %>% testDispersion()

# --
pp.glmm5 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm5)
anova(pp.glmm5)
car::Anova(pp.glmm5, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm5)  # AICc = 965.7573
vif(pp.glmm5)
dispersion_glmer(pp.glmm5) #to test if overdispersed, has to be > 1.4
mod_dharma5 <- pp.glmm5 %>% simulateResiduals(n=1000)
plot(mod_dharma5)
mod_dharma5 %>% testDispersion()

# --
pp.glmm6 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm6)
anova(pp.glmm6)
car::Anova(pp.glmm6, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm6)  # AICc = 966.1881
vif(pp.glmm6)
dispersion_glmer(pp.glmm6) #to test if overdispersed, has to be > 1.4
mod_dharma6 <- pp.glmm6 %>% simulateResiduals(n=1000)
plot(mod_dharma6)
mod_dharma6 %>% testDispersion()

# --
pp.glmm7 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm7)
anova(pp.glmm7)
car::Anova(pp.glmm7, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm7)  # AICc = 965.7231
vif(pp.glmm7)
dispersion_glmer(pp.glmm7) #to test if overdispersed, has to be > 1.4
mod_dharma7 <- pp.glmm7 %>% simulateResiduals(n=1000)
plot(mod_dharma7)
mod_dharma7 %>% testDispersion()

# --
pp.glmm8 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + TUVA:PL + TUVA:RD + PL:RD, 
                  family= poisson,
                  data= data)
summary(pp.glmm8)
anova(pp.glmm8)
car::Anova(pp.glmm8, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm8)  # AICc = 968.0945
vif(pp.glmm8)
dispersion_glmer(pp.glmm8) #to test if overdispersed, has to be > 1.4
mod_dharma8 <- pp.glmm8 %>% simulateResiduals(n=1000)
plot(mod_dharma8)
mod_dharma8 %>% testDispersion()

# --
pp.glmm9 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:TUVA, 
                  family= poisson,
                  data= data)
summary(pp.glmm9)
anova(pp.glmm9)
car::Anova(pp.glmm9, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm9)  # AICc = 965.715
vif(pp.glmm9)
dispersion_glmer(pp.glmm9) #to test if overdispersed, has to be > 1.4
mod_dharma9 <- pp.glmm9 %>% simulateResiduals(n=1000)
plot(mod_dharma9)
mod_dharma9 %>% testDispersion()

# --
pp.glmm10 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:PL, 
                   family= poisson,
                   data= data)
summary(pp.glmm10)
anova(pp.glmm10)
car::Anova(pp.glmm10, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm10)  # AICc = 964.1856
vif(pp.glmm10)
dispersion_glmer(pp.glmm10) #to test if overdispersed, has to be > 1.4
mod_dharma10 <- pp.glmm10 %>% simulateResiduals(n=1000)
plot(mod_dharma10)
mod_dharma10 %>% testDispersion()

# --
pp.glmm11 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD, 
                   family= poisson,
                   data= data)
summary(pp.glmm11)
anova(pp.glmm11)
car::Anova(pp.glmm11, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm11)  # AICc = 966.4854
vif(pp.glmm11)
dispersion_glmer(pp.glmm11) #to test if overdispersed, has to be > 1.4
mod_dharma11 <- pp.glmm11 %>% simulateResiduals(n=1000)
plot(mod_dharma11)
mod_dharma11 %>% testDispersion()

# --
pp.glmm12 <- glmer(pp.rich ~ Transect_type + TUVA + PL + RD + (1 | Landscape) + Transect_type:RD + Transect_type:TUVA, 
                   family= poisson,
                   data= data)
summary(pp.glmm12)
anova(pp.glmm12)
car::Anova(pp.glmm12, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm12)  # AICc = 971.8976
vif(pp.glmm12)
dispersion_glmer(pp.glmm12) #to test if overdispersed, has to be > 1.4
mod_dharma12 <- pp.glmm12 %>% simulateResiduals(n=1000)
plot(mod_dharma12)
mod_dharma12 %>% testDispersion()

# --
pp.glmm0 <- glmer(pp.rich ~ 1 + (1 | Landscape), 
                  family= poisson,
                  data= data)
summary(pp.glmm0)
anova(pp.glmm0)
car::Anova(pp.glmm0, type="III")
options(contrasts = c("contr.sum","contr.poly"))
aicc(pp.glmm0)  # AICc = 995.1401
vif(pp.glmm0)
dispersion_glmer(pp.glmm0) #to test if overdispersed, has to be > 1.4
mod_dharma0 <- pp.glmm0 %>% simulateResiduals(n=1000)
plot(mod_dharma0)
mod_dharma0 %>% testDispersion()

# Model selection: red listed butterfly species (binomial presence absence) -----
redbf.glmm <- glmer(bf.redpa ~ Transect_type * TUVA * PL * RD + (1 | Landscape),
                   family= binomial,
                   na.action= na.fail,
                   glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),
                   data= data)
summary(redbf.glmm)
plot(redbf.glmm)
anova(redbf.glmm)
options(contrasts = c("contr.sum","contr.poly"))
options('contrasts')
car::Anova(redbf.glmm, type= "III")
aicc(redbf.glmm)  # AICc = 741.5434
vif(redbf.glmm)
dispersion_glmer(redbf.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- redbf.spglmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()
# Model averaging
redbf.glmm <- standardize(redbf.glmm) # yielded exact same output
MuMIn.redbf <- dredge(redbf.glmm)
MA.redbf <- get.models(MuMIn.redbf, subset=delta<4)
summary(model.avg(MA.redbf))
importance(summary(model.avg(MA.redbf)))

# Best model:
redbf.glmm <- glmer(bf.redpa ~ Transect_type + TUVA + PL + RD + (1 | Landscape),
                   family= binomial,
                   na.action= na.fail,
                   glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),
                   data= data)
summary(redbf.glmm)
plot(redbf.glmm)
anova(redbf.glmm)
options(contrasts = c("contr.sum","contr.poly"))
options('contrasts')
car::Anova(redbf.glmm, type= "III")
aicc(redbf.glmm)  # AICc = 551.5809
vif(redbf.glmm)
dispersion_glmer(redbf.glmm) #to test if overdispersed, has to be > 1.4
mod_dharma <- redbf.glmm %>% simulateResiduals(n=1000)
plot(mod_dharma)
mod_dharma %>% testDispersion()

# Plot figure
spbf.plot<-effect_plot(spbf.glmm, pred = Transect_type, interval = TRUE, int.type = "confidence", int.width = .8,
                       outcome.scale = "response",
                       x.label = "Habitat", y.label = "Number of butterfly species",
                       line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
                       point.alpha = 0.3                    )
spbf.plot + scale_x_discrete() +theme(axis.text=element_text(size =15), legend.position = "none")
# theme_update(text = element_text(size=20))
visreg(bf.glmm, "Transect_type", scale="response", band=TRUE, partial=FALSE, rug=FALSE, type="contrast")#, gg=TRUE)

# Plot figure old:
# Butterflies:
# # Plot figure
# bf.plot<-effect_plot(bf.glmm, pred = Transect_type, interval = TRUE, int.type = "confidence", int.width = .8,
#                     outcome.scale = "response",
#                     x.label = "Habitat", y.label = "Number of butterfly species",
#                     line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
#                     point.alpha = 0.3 )
# bf.plot + scale_x_discrete() +theme(axis.text=element_text(size =15), legend.position = "none")
# levels <- c("Between fields", "Big road", "Small road", "Pasture", "Powerline")
# bf.plot + ggplot(bf.plot, aes(x=factor(Transect_type, level=levels), y=bf.sp ))
# # theme_update(text = element_text(size=20))
# visreg(bf.glmm, "Transect_type", scale="response", band=TRUE, partial=FALSE, rug=FALSE, type="contrast")#, gg=TRUE)



# Plot figure - bumblebees
# bb.plot<-effect_plot(bb.glmm, pred = Transect_type, interval = TRUE, int.type = "confidence", int.width = .8,
#                      outcome.scale = "response",
#                      x.label = "Habitat", y.label = "Number of bumblebee species",
#                      line.thickness = 0.7, plot.points=FALSE, partial.residuals= FALSE, point.size =0.1, jitter=0.5,
#                      point.alpha = 0.3                    )
# bb.plot + scale_x_discrete() +theme(axis.text=element_text(size = 15), legend.position = "none")
# # theme_update(text = element_text(size=20))