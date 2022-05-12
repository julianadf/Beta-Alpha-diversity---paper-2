rm(list=ls())
library(here)
library(tidyverse)
library(emmeans)
library(vegan)
library(betapart)
library(lmerTest)
library(cowplot)
library(knitr)
library(pairwiseAdonis)
library(png)
library(magick)
library(ggpubr)
library(zetadiv)


### butterflies ###
bf <- read.csv2(here("data", "Butterfly.dataframe.csv"))
abund.bf <- bf[,4:ncol(bf)]
habitat.type.bf <- bf$Transect_type
habitat <- bf$Transect_type

NMDS.bf <- metaMDS(abund.bf, try = 100, trymax = 1000, maxit = 20000)

# Create figure with ggplot
bf.scores <- as.data.frame(scores(NMDS.bf))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
bf.scores$landscape <- rownames(bf.scores)  # create a column of site names, from the rownames of data.scores
bf.scores$category <- category  #  add the grp variable created earlier
head(bf.scores)  #look at the data


bf.scores <- as.data.frame(scores(NMDS.bf))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
bf.scores$habitat <- rownames(bf.scores)  # create a column of site names, from the rownames of data.scores
bf.scores$habitat <- habitat  #  add the grp variable created earlier
head(bf.scores)  #look at the data

bfsp.scores <- as.data.frame(scores(NMDS.bf, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
bfsp.scores$species <- rownames(bfsp.scores)  # create a column of species, from the rownames of species.scores
head(bfsp.scores)  #look at the data

bf.pasture <- bf.scores[bf.scores$habitat == "Pasture", ][chull(bf.scores[bf.scores$habitat =="Pasture", c("NMDS1", "NMDS2")]), ]  # hull values 
bf.powerline <- bf.scores[bf.scores$habitat == "Powerline", ][chull(bf.scores[bf.scores$habitat =="Powerline", c("NMDS1", "NMDS2")]), ]  # hull 


hull.bf <- rbind(bf.pasture, bf.powerline)  #combine landscapes
hull.bf

# main text
gg.bf<- ggplot() + 
  geom_polygon(data=hull.bf,aes(x=NMDS1,y=NMDS2,fill=habitat),alpha=0.25) + # add the convex hulls
  #geom_text(data=bfsp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3)   # add the species labels
  #geom_point(data=bf.scores,aes(x=NMDS1,y=NMDS2,shape=habitat,colour=habitat),size=2) +  # add the point markers
  #geom_text(data=bf.scores,aes(x=NMDS1,y=NMDS2,label=landscape),size=2,vjust=0) +  # add the site labels
  scale_color_manual("habitat", values = c("#2c7bb6",  "#fdae61", "#e6f598", "#abd9e9", "#d7191c")) +
  scale_fill_manual("habitat", values = c("#2c7bb6",  "#fdae61", "#e6f598", "#abd9e9", "#d7191c")) +
  scale_shape_discrete("habitat") +
  coord_equal() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
       axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
       axis.text.x = element_text(angle = 90, vjust=0.5),
       axis.title.y = element_text(size=28), axis.title.x = element_text(size=28))
gg.bf

# add silhouette
maniola <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/butterfly.png")
a <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/a.png")
ggdraw(gg.bf) + draw_image(maniola, x = 0.87, y = 0.99, hjust = 1, vjust = 1, width = 0.15, height = 0.15) +
  draw_image(a, x = 0.3, y = 1, hjust = 1, vjust = 1, width = 0.1, height = 0.1)

### bumblebees ###
bb <- read.csv2(here("data", "Bumblebee.dataframe.csv"))
abund.bb <- bb[,4:ncol(bb)]
habitat.type.bb <- bb$Transect_type
habitat <- bb$Transect_type

NMDS.bb <- metaMDS(abund.bb, try = 100, trymax = 1000, maxit = 20000)

# Create figure with ggplot
bb.scores <- as.data.frame(scores(NMDS.bb))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
bb.scores$habitat <- rownames(bb.scores)  # create a column of site names, from the rownames of data.scores
bb.scores$habitat <- habitat  #  add the grp variable created earlier
head(bb.scores)  #look at the data

bbsp.scores <- as.data.frame(scores(NMDS.bb, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
bbsp.scores$species <- rownames(bbsp.scores)  # create a column of species, from the rownames of species.scores
head(bbsp.scores)  #look at the data

bb.pasture <- bb.scores[bb.scores$habitat == "Pasture", ][chull(bb.scores[bb.scores$habitat =="Pasture", c("NMDS1", "NMDS2")]), ]  # hull values 
bb.powerline <- bb.scores[bb.scores$habitat == "Powerline", ][chull(bb.scores[bb.scores$habitat =="Powerline", c("NMDS1", "NMDS2")]), ]  # hull 


hull.bb <- rbind(bb.pasture, bb.powerline)  #combine landscapes
hull.bb

# main text
gg.bb<- ggplot() + 
  geom_polygon(data=hull.bb,aes(x=NMDS1,y=NMDS2,fill=habitat),alpha=0.25) + # add the convex hulls
  #geom_text(data=bfsp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3)   # add the species labels
  #geom_point(data=bf.scores,aes(x=NMDS1,y=NMDS2,shape=habitat,colour=habitat),size=2) +  # add the point markers
  #geom_text(data=bf.scores,aes(x=NMDS1,y=NMDS2,label=landscape),size=2,vjust=0) +  # add the site labels
  scale_color_manual("habitat", values = c("#2c7bb6",  "#fdae61", "#e6f598", "#abd9e9", "#d7191c")) +
  scale_fill_manual("habitat", values = c("#2c7bb6",  "#fdae61", "#e6f598", "#abd9e9", "#d7191c")) +
  scale_shape_discrete("habitat") +
  coord_equal() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28))
gg.bb

# add silhouette
terrestris <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/bumble.png")
b <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/b.png")
ggdraw(gg.bb) + draw_image(terrestris, x = 0.84, y = 0.98, hjust = 1, vjust = 1, width = 0.15, height = 0.15) +
  draw_image(b, x = 0.325, y = 1, hjust = 1, vjust = 1, width = 0.1, height = 0.1) 

### plants ###
pp <- read.csv2(here("data", "Plant.dataframe.csv"))
abund.pp <- pp[,4:ncol(pp)]
habitat.type.pp <- pp$Transect_type
habitat <- pp$Transect_type

NMDS.pp <- metaMDS(abund.pp, try = 100, trymax = 1000, maxit = 20000)

# Create figure with ggplot
pp.scores <- as.data.frame(scores(NMDS.pp))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
pp.scores$habitat <- rownames(pp.scores)  # create a column of site names, from the rownames of data.scores
pp.scores$habitat <- habitat  #  add the grp variable created earlier
head(pp.scores)  #look at the data

ppsp.scores <- as.data.frame(scores(NMDS.pp, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
ppsp.scores$species <- rownames(ppsp.scores)  # create a column of species, from the rownames of species.scores
head(ppsp.scores)  #look at the data

pp.pasture <- pp.scores[pp.scores$habitat == "Pasture", ][chull(pp.scores[pp.scores$habitat =="Pasture", c("NMDS1", "NMDS2")]), ]  # hull values 
pp.powerline <- pp.scores[pp.scores$habitat == "Powerline", ][chull(pp.scores[pp.scores$habitat =="Powerline", c("NMDS1", "NMDS2")]), ]  # hull 


hull.pp <- rbind(pp.pasture, pp.powerline)  #combine landscapes
hull.pp

# main text
gg.pp<- ggplot() + 
  geom_polygon(data=hull.pp,aes(x=NMDS1,y=NMDS2,fill=habitat),alpha=0.25) + # add the convex hulls
  #geom_text(data=bfsp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3)   # add the species labels
  #geom_point(data=bf.scores,aes(x=NMDS1,y=NMDS2,shape=habitat,colour=habitat),size=2) +  # add the point markers
  #geom_text(data=bf.scores,aes(x=NMDS1,y=NMDS2,label=landscape),size=2,vjust=0) +  # add the site labels
  scale_color_manual("habitat", values = c("#2c7bb6",  "#fdae61", "#e6f598", "#abd9e9", "#d7191c")) +
  scale_fill_manual("habitat", values = c("#2c7bb6",  "#fdae61", "#e6f598", "#abd9e9", "#d7191c")) +
  scale_shape_discrete("habitat") +
  coord_equal() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28))
gg.pp


# add silhouette
fragaria <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/Fragaria.png")
c <- readPNG("//storage.slu.se/Home$/jada0002/My Documents/My Pictures/Illustrations/c.png")
ggdraw(gg.pp) + draw_image(fragaria, x = 0.95, y = 0.98, hjust = 1, vjust = 1, width = 0.26, height = 0.26) +
  draw_image(c, x = 0.27, y = 1, hjust = 1, vjust = 1, width = 0.1, height = 0.1)


