setwd("C:/Users/b0040/Documents/abelhas")


library(raster)
library(GGally)
library(statmod)
library(fitdistrplus)
library(betareg)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(sjPlot)
library(sjmisc)
library(coefplot)
library(ggeffects)
library(multcompView)
library(emmeans)
library(effects)




area.change <- read.csv("bees_traits_final.csv", header = T)
area.change<-area.change[,c(1:4,7:10)]
area.change$size.cat<-"medium"
area.change$size.cat[which(area.change$ITDmeasured<2.2)]<-"small"
area.change$size.cat[which(area.change$ITDmeasured>3.9)]<-"large"



#############################################################################################################################
#############################################################################################################################

listagem_ensemble<-list.files(path = "models", pattern = ".tif$", full.names = T, recursive = T)
head(listagem_ensemble)


db_nomes<-data.frame(strsplit(listagem_ensemble, "/|_|.tif$"))
db_nomes<-db_nomes[c(4:6),]
colnames(db_nomes)<-apply(db_nomes, 2, function(x) paste(x, collapse="_"))


ensembles<-stack(listagem_ensemble)
names(ensembles)<-colnames(db_nomes)
head(names(ensembles))
plot(ensembles[[1:5]])


dir.create("bin")
#i <- "Lestrimelitta_glabrata"
for(i in names(ensembles)){
  ensemblesbin <- (ensemble[[grep(i, names(ensemble), value = T)]])
  ensemblesbin[ensemblesbin < 60] <- 0
  ensemblesbin[ensemblesbin >= 60] <- 1
  #unique(ensemblesbin)
  #plot(ensemblesbin)
  
  
  writeRaster(ensemblesbin, paste("bin/", i, sep=""), format="HFA")
  cat(i, "pronto..\n")
  
}

rm(list= ls()[!(ls() %in% c("area.change"))])
#############################################################################################################################
#############################################################################################################################


bin.list <-list.files(path = "bin", pattern = ".img$", full.names = T, recursive = T)
head(bin.list)

bin<-stack(bin.list)
head(names(bin))
plot(bin[[1:5]])

PA<-shapefile("C:\\Users\\b0040\\Documents\\shapes\\base\\Politico\\DivisaÌfo PoliÌtica BR\\PA.shp")
plot(PA, add=T)

sp.list <- read.csv("abelhas_modeladas.csv")
area.change <- area.change[area.change$apelido %in% sp.list$apelido,]


#spp <- "Lestrimelitta_glabrata"
for(spp in sp.list$apelido){
  multi_cenarios<-(bin[[grep(spp, names(bin), value = T)]]) 
  
  Cmulti_cenarios <- crop(multi_cenarios, extent(PA), snap="out")
  PAraster <- rasterize(PA, Cmulti_cenarios) 
  CFmulti_cenarios <- mask(x=Cmulti_cenarios, mask=PAraster)
  
  
  
  atual<-CFmulti_cenarios[[1]]
  #atual
  futuro1a<-CFmulti_cenarios[[2]]
  #futuro1a
  futuro1b<-CFmulti_cenarios[[3]]
  #futuro1b
  futuro2a<-CFmulti_cenarios[[4]]
  #futuro2a
  futuro2b<-CFmulti_cenarios[[5]]
  #futuro2b
  
  area.change$area_current[which(area.change$apelido==spp)] <- sum(atual[atual==1])*10^2
  area.change$area_future1A[which(area.change$apelido==spp)] <- sum(futuro1a[futuro1a==1])*10^2
  area.change$area_future1B[which(area.change$apelido==spp)] <- sum(futuro1b[futuro1b==1])*10^2
  area.change$area_future2A[which(area.change$apelido==spp)] <- sum(futuro2a[futuro2a==1])*10^2
  area.change$area_future2B[which(area.change$apelido==spp)] <- sum(futuro2b[futuro2b==1])*10^2
  
  }


area.change$respF1A <- (area.change$area_future1A - area.change$area_current)/area.change$area_current
area.change$respF1B <- (area.change$area_future1B - area.change$area_current)/area.change$area_current
area.change$respF2A <- (area.change$area_future2A - area.change$area_current)/area.change$area_current
area.change$respF2B <- (area.change$area_future2B - area.change$area_current)/area.change$area_current


area.change$area.cat<-"medium"
area.change$area.cat[which(area.change$area_current<=441000)]<-"small"
area.change$area.cat[which(area.change$area_current>882000 )]<-"large"




str(area.change)



write.csv(area.change, "bees_traits_final[a].csv")



rm(list= ls()[!(ls() %in% c("area.change"))])
#############################################################################################################################
#############################################################################################################################

#area.change <- read.csv("bees_traits_final[a].csv", header = T)
area.change <- area.change[complete.cases(area.change),]

write.csv(area.change, "bees_traits_final[b].csv")


###################
### area loss #####
###################

area.loss <- area.change[which(area.change$respF1A<0),]
area.loss$respF1A <- area.loss$respF1A*(-1)
area.loss$respF1B <- area.loss$respF1B*(-1)
area.loss$respF2A <- area.loss$respF2A*(-1)
area.loss$respF2B <- area.loss$respF2B*(-1)

#if y also assumes the extremes 0 and 1, 
#

n<-length(area.loss$apelido)
area.loss$respF1A <- (area.loss$respF1A * (n -1) +0.5)/n
area.loss$respF1B <- (area.loss$respF1B * (n -1) +0.5)/n
area.loss$respF2A <- (area.loss$respF2A * (n -1) +0.5)/n
area.loss$respF2B <- (area.loss$respF2B * (n -1) +0.5)/n

hist(area.loss$respF1A)

write.csv(area.change, "bees_traits_final[c].csv")



ft1A <- betareg(respF1A ~ area.cat + size.cat + sociality + nest.location + crop.pollination, data = area.loss)
summary(ft1A)
AIC(ft1A)
plot(ft1A)

ft1B <- betareg(respF1B ~ area.cat + size.cat + sociality + nest.location + crop.pollination, data = area.loss)
summary(ft1B)
plot(ft1B)

ft2A <- betareg(respF2A ~ area.cat + size.cat + sociality + nest.location + crop.pollination, data = area.loss)
summary(ft2A)
plot(ft2A)

ft2B <- betareg(respF2B ~ area.cat + size.cat + sociality + nest.location + crop.pollination, data = area.loss)
summary(ft2B)
plot(ft2B)




###################
###   area    #####
###################

marg.area1A <- emmeans(ft1A, ~ area.cat)
pairs(marg.area1A, adjust="tukey")
sum.area1A <- cld(marg.area1A, alpha = .05, Letters = letters, adjust = "tukey")
sum.area1A
sum.area1A = sum.area1A[order(factor(sum.area1A$area.cat, levels=c("restrict", "medium", "wide"))),]
sum.area1A$year <- "2050"
sum.area1A$scenario <- "rcp 4.5"
sum.area1A



marg.area1B <- emmeans(ft1B, ~ area.cat)
pairs(marg.area1B, adjust="tukey")
sum.area1B <- cld(marg.area1B, alpha = .05, Letters = letters, adjust = "tukey")
sum.area1B
sum.area1B = sum.area1B[order(factor(sum.area1B$area.cat, levels=c("restrict", "medium", "wide"))),]
sum.area1B$year <- "2050"
sum.area1B$scenario <- "rcp 8.5"
sum.area1B



marg.area2A <- emmeans(ft2A, ~ area.cat)
pairs(marg.area2A, adjust="tukey")
sum.area2A <- cld(marg.area2A, alpha = .05, Letters = letters, adjust = "tukey")
sum.area2A
sum.area2A = sum.area2A[order(factor(sum.area2A$area.cat, levels=c("restrict", "medium", "wide"))),]
sum.area2A$year <- "2070"
sum.area2A$scenario <- "rcp 4.5"
sum.area2A



marg.area2B <- emmeans(ft2B, ~ area.cat)
pairs(marg.area2B, adjust="tukey")
sum.area2B <- cld(marg.area2B, alpha = .05, Letters = letters, adjust = "tukey")
sum.area2B
sum.area2B = sum.area2B[order(factor(sum.area2B$area.cat, levels=c("restrict", "medium", "wide"))),]
sum.area2B$year <- "2070"
sum.area2B$scenario <- "rcp 8.5"
sum.area2B


sum.area <- rbind(sum.area1A, sum.area1B, sum.area2A, sum.area2B)
write.csv(sum.area, "area_effects.csv", row.names = F)

grph.area <- ggplot(sum.area, aes(x = year, y = emmean, color=scenario, shape = area.cat))+
  geom_point(size=8, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .15, size=1, position = position_dodge(0.5))+
  scale_y_continuous(name = "", breaks = seq(0,1,.1))+ 
  scale_x_discrete(limits=c("2050", "2070"), expand = expand_scale(add = .4))+
  scale_color_manual(values = c("#1b9e77", "#d95f02"))+
  scale_shape_manual(limits=c("restrict", "medium", "wide"), values = c(19,17,15))+
  coord_cartesian(ylim=c(0.5,1))+
  #annotate("text", x = 1:12, y = sum.area$asymp.UCL+0.5, label = gsub(" ", "", sum.area$.group), area=8, position = position_dodge(0.5))+
  labs(x ="" , y = "", title = NULL)+
  theme_bw()+
  theme(axis.title = element_text(lineheight=.8, hjust=0.5, size = 36, face = "bold"),
        text = element_text(family = "serif", size = 28),
        panel.background = element_rect(linetype = "blank", fill = NA),
        panel.grid = element_blank())

grph.area

png(filename = "area_effects.png",
    width = 30, height = 20, units = "cm", 
    bg = "white", res = 600)

grph.area
dev.off()




rm(list= ls()[!(ls() %in% c("area.change", "area.loss", "ft1A", "ft1B", "ft2A", "ft2B", "sum.area"))])
###################
###   size    #####
###################

marg.size1A <- emmeans(ft1A, ~ size.cat)
pairs(marg.size1A, adjust="tukey")
sum.size1A <- cld(marg.size1A, alpha = .05, Letters = letters, adjust = "tukey")
sum.size1A
sum.size1A = sum.size1A[order(factor(sum.size1A$size.cat, levels=c("small", "medium", "large"))),]
sum.size1A$year <- "2050"
sum.size1A$scenario <- "rcp 4.5"
sum.size1A



marg.size1B <- emmeans(ft1B, ~ size.cat)
pairs(marg.size1B, adjust="tukey")
sum.size1B <- cld(marg.size1B, alpha = .05, Letters = letters, adjust = "tukey")
sum.size1B
sum.size1B = sum.size1B[order(factor(sum.size1B$size.cat, levels=c("small", "medium", "large"))),]
sum.size1B$year <- "2050"
sum.size1B$scenario <- "rcp 8.5"
sum.size1B



marg.size2A <- emmeans(ft2A, ~ size.cat)
pairs(marg.size2A, adjust="tukey")
sum.size2A <- cld(marg.size2A, alpha = .05, Letters = letters, adjust = "tukey")
sum.size2A
sum.size2A = sum.size2A[order(factor(sum.size2A$size.cat, levels=c("small", "medium", "large"))),]
sum.size2A$year <- "2070"
sum.size2A$scenario <- "rcp 4.5"
sum.size2A



marg.size2B <- emmeans(ft2B, ~ size.cat)
pairs(marg.size2B, adjust="tukey")
sum.size2B <- cld(marg.size2B, alpha = .05, Letters = letters, adjust = "tukey")
sum.size2B
sum.size2B = sum.size2B[order(factor(sum.size2B$size.cat, levels=c("small", "medium", "large"))),]
sum.size2B$year <- "2070"
sum.size2B$scenario <- "rcp 8.5"
sum.size2B


sum.size <- rbind(sum.size1A, sum.size1B, sum.size2A, sum.size2B)
write.csv(sum.size, "size_effects.csv", row.names = F)

grph.size <- ggplot(sum.size, aes(x = year, y = emmean, color=scenario, shape = size.cat))+
  geom_point(size=8, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .15, size=1, position = position_dodge(0.5))+
  scale_y_continuous(name = "", breaks = seq(0,1,.1))+ 
  scale_x_discrete(limits=c("2050", "2070"), expand = expand_scale(add = .4))+
  scale_color_manual(values = c("#1b9e77", "#d95f02"))+
  coord_cartesian(ylim=c(0.5,1))+
  #annotate("text", x = 1:12, y = sum.size$asymp.UCL+0.5, label = gsub(" ", "", sum.size$.group), size=8, position = position_dodge(0.5))+
  labs(x ="" , y = "", title = NULL)+
  theme_bw()+
  theme(axis.title = element_text(lineheight=.8, hjust=0.5, size = 36, face = "bold"),
        text = element_text(family = "serif", size = 28),
        panel.background = element_rect(linetype = "blank", fill = NA),
        panel.grid = element_blank())

grph.size

png(filename = "size_effects.png",
    width = 30, height = 20, units = "cm", 
    bg = "white", res = 600)

grph.size
dev.off()




rm(list= ls()[!(ls() %in% c("area.change", "area.loss", "ft1A", "ft1B", "ft2A", "ft2B", "sum.size"))])
###################
###  sociality  ###
###################

marg.sociality1A <- emmeans(ft1A, ~ sociality)
pairs(marg.sociality1A, adjust="tukey")
sum.sociality1A <- cld(marg.sociality1A, alpha = .05, Letters = letters, adjust = "tukey")
sum.sociality1A
sum.sociality1A = sum.sociality1A[order(factor(sum.sociality1A$sociality, levels=c("solitary", "social", "cleptoparasitic"))),]
sum.sociality1A$year <- "2050"
sum.sociality1A$scenario <- "rcp 4.5"
sum.sociality1A



marg.sociality1B <- emmeans(ft1B, ~ sociality)
pairs(marg.sociality1B, adjust="tukey")
sum.sociality1B <- cld(marg.sociality1B, alpha = .05, Letters = letters, adjust = "tukey")
sum.sociality1B
sum.sociality1B = sum.sociality1B[order(factor(sum.sociality1B$sociality, levels=c("solitary", "social", "cleptoparasitic"))),]
sum.sociality1B$year <- "2050"
sum.sociality1B$scenario <- "rcp 8.5"
sum.sociality1B



marg.sociality2A <- emmeans(ft2A, ~ sociality)
pairs(marg.sociality2A, adjust="tukey")
sum.sociality2A <- cld(marg.sociality2A, alpha = .05, Letters = letters, adjust = "tukey")
sum.sociality2A
sum.sociality2A = sum.sociality2A[order(factor(sum.sociality2A$sociality, levels=c("solitary", "social", "cleptoparasitic"))),]
sum.sociality2A$year <- "2070"
sum.sociality2A$scenario <- "rcp 4.5"
sum.sociality2A



marg.sociality2B <- emmeans(ft2B, ~ sociality)
pairs(marg.sociality2B, adjust="tukey")
sum.sociality2B <- cld(marg.sociality2B, alpha = .05, Letters = letters, adjust = "tukey")
sum.sociality2B
sum.sociality2B = sum.sociality2B[order(factor(sum.sociality2B$sociality, levels=c("solitary", "social", "cleptoparasitic"))),]
sum.sociality2B$year <- "2070"
sum.sociality2B$scenario <- "rcp 8.5"
sum.sociality2B


sum.sociality <- rbind(sum.sociality1A, sum.sociality1B, sum.sociality2A, sum.sociality2B)
write.csv(sum.sociality, "sociality_effects.csv", row.names = F)

grph.sociality <- ggplot(sum.sociality, aes(x = year, y = emmean, color=scenario, shape = sociality))+
  geom_point(size=8, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .15, size=1, position = position_dodge(0.5))+
  scale_y_continuous(name = "", breaks = seq(0,1,.1))+ 
  scale_x_discrete(limits=c("2050", "2070"), expand = expand_scale(add = .4))+
  scale_color_manual(values = c("#1b9e77", "#d95f02"))+
  coord_cartesian(ylim=c(0.5,1))+
  #annotate("text", x = 1:12, y = sum.sociality$asymp.UCL+0.5, label = gsub(" ", "", sum.sociality$.group), sociality=8, position = position_dodge(0.5))+
  labs(x ="" , y = "", title = NULL)+
  theme_bw()+
  theme(axis.title = element_text(lineheight=.8, hjust=0.5, size = 36, face = "bold"),
        text = element_text(family = "serif", size = 28),
        panel.background = element_rect(linetype = "blank", fill = NA),
        panel.grid = element_blank())

grph.sociality

png(filename = "sociality_effects.png",
    width = 30, height = 20, units = "cm", 
    bg = "white", res = 600)

grph.sociality
dev.off()




rm(list= ls()[!(ls() %in% c("area.change", "area.loss", "ft1A", "ft1B", "ft2A", "ft2B", "sum.size", "sum.sociality"))])
###################
###   nest     ###
###################

marg.nest1A <- emmeans(ft1A, ~ nest.location)
pairs(marg.nest1A, adjust="tukey")
sum.nest1A <- cld(marg.nest1A, alpha = .05, Letters = letters, adjust = "tukey")
sum.nest1A
sum.nest1A = sum.nest1A[order(factor(sum.nest1A$nest.location, levels=c("exposed", "soil", "cavity", "termite", "multi"))),]
sum.nest1A$year <- "2050"
sum.nest1A$scenario <- "rcp 4.5"
sum.nest1A



marg.nest1B <- emmeans(ft1B, ~ nest.location)
pairs(marg.nest1B, adjust="tukey")
sum.nest1B <- cld(marg.nest1B, alpha = .05, Letters = letters, adjust = "tukey")
sum.nest1B
sum.nest1B = sum.nest1B[order(factor(sum.nest1B$nest.location, levels=c("exposed", "soil", "cavity", "termite", "multi"))),]
sum.nest1B$year <- "2050"
sum.nest1B$scenario <- "rcp 8.5"
sum.nest1B



marg.nest2A <- emmeans(ft2A, ~ nest.location)
pairs(marg.nest2A, adjust="tukey")
sum.nest2A <- cld(marg.nest2A, alpha = .05, Letters = letters, adjust = "tukey")
sum.nest2A
sum.nest2A = sum.nest2A[order(factor(sum.nest2A$nest.location, levels=c("exposed", "soil", "cavity", "termite", "multi"))),]
sum.nest2A$year <- "2070"
sum.nest2A$scenario <- "rcp 4.5"
sum.nest2A



marg.nest2B <- emmeans(ft2B, ~ nest.location)
pairs(marg.nest2B, adjust="tukey")
sum.nest2B <- cld(marg.nest2B, alpha = .05, Letters = letters, adjust = "tukey")
sum.nest2B
sum.nest2B = sum.nest2B[order(factor(sum.nest2B$nest.location, levels=c("exposed", "soil", "cavity", "termite", "multi"))),]
sum.nest2B$year <- "2070"
sum.nest2B$scenario <- "rcp 8.5"
sum.nest2B


sum.nest <- rbind(sum.nest1A, sum.nest1B, sum.nest2A, sum.nest2B)
write.csv(sum.nest, "nest_effects.csv", row.names = F)

grph.nest <- ggplot(sum.nest, aes(x = year, y = emmean, color=scenario, shape = nest.location))+
  geom_point(size=8, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .15, size=1, position = position_dodge(0.5))+
  scale_y_continuous(name = "", breaks = seq(0,1,.1))+ 
  scale_x_discrete(limits=c("2050", "2070"), expand = expand_scale(add = .4))+
  scale_color_manual(values = c("#1b9e77", "#d95f02"))+
  scale_shape_manual(values = c(16,17,15,18,8))+
  coord_cartesian(ylim=c(0.5,1))+
  #annotate("text", x = 1:12, y = sum.nest$asymp.UCL+0.5, label = gsub(" ", "", sum.nest$.group), nest=8, position = position_dodge(0.5))+
  labs(x ="" , y = "", title = NULL)+
  theme_bw()+
  theme(axis.title = element_text(lineheight=.8, hjust=0.5, size = 36, face = "bold"),
        text = element_text(family = "serif", size = 28),
        panel.background = element_rect(linetype = "blank", fill = NA),
        panel.grid = element_blank())

grph.nest

png(filename = "nest_effects.png",
    width = 30, height = 20, units = "cm", 
    bg = "white", res = 600)

grph.nest
dev.off()




rm(list= ls()[!(ls() %in% c("area.change", "area.loss", "ft1A", "ft1B", "ft2A", "ft2B", "sum.size", "sum.sociality", "sum.nest"))])
###################
###   crop     ###
###################

marg.crop1A <- emmeans(ft1A, ~ crop.pollination)
pairs(marg.crop1A, adjust="tukey")
sum.crop1A <- cld(marg.crop1A, alpha = .05, Letters = letters, adjust = "tukey")
sum.crop1A
sum.crop1A = sum.crop1A[order(factor(sum.crop1A$crop.pollination, levels=c("no", "yes"))),]
sum.crop1A$year <- "2050"
sum.crop1A$scenario <- "rcp 4.5"
sum.crop1A



marg.crop1B <- emmeans(ft1B, ~ crop.pollination)
pairs(marg.crop1B, adjust="tukey")
sum.crop1B <- cld(marg.crop1B, alpha = .05, Letters = letters, adjust = "tukey")
sum.crop1B
sum.crop1B = sum.crop1B[order(factor(sum.crop1B$crop.pollination, levels=c("no", "yes"))),]
sum.crop1B$year <- "2050"
sum.crop1B$scenario <- "rcp 8.5"
sum.crop1B



marg.crop2A <- emmeans(ft2A, ~ crop.pollination)
pairs(marg.crop2A, adjust="tukey")
sum.crop2A <- cld(marg.crop2A, alpha = .05, Letters = letters, adjust = "tukey")
sum.crop2A
sum.crop2A = sum.crop2A[order(factor(sum.crop2A$crop.pollination, levels=c("no", "yes"))),]
sum.crop2A$year <- "2070"
sum.crop2A$scenario <- "rcp 4.5"
sum.crop2A



marg.crop2B <- emmeans(ft2B, ~ crop.pollination)
pairs(marg.crop2B, adjust="tukey")
sum.crop2B <- cld(marg.crop2B, alpha = .05, Letters = letters, adjust = "tukey")
sum.crop2B
sum.crop2B = sum.crop2B[order(factor(sum.crop2B$crop.pollination, levels=c("no", "yes"))),]
sum.crop2B$year <- "2070"
sum.crop2B$scenario <- "rcp 8.5"
sum.crop2B


sum.crop <- rbind(sum.crop1A, sum.crop1B, sum.crop2A, sum.crop2B)
write.csv(sum.crop, "crop_effects.csv", row.names = F)

grph.crop <- ggplot(sum.crop, aes(x = year, y = emmean, color=scenario, shape = crop.pollination))+
  geom_point(size=8, position = position_dodge(0.5))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .15, size=1, position = position_dodge(0.5))+
  scale_y_continuous(name = "", breaks = seq(0,1,.1))+ 
  scale_x_discrete(limits=c("2050", "2070"), expand = expand_scale(add = .4))+
  scale_color_manual(values = c("#1b9e77", "#d95f02"))+
  coord_cartesian(ylim=c(0.5,1))+
  #annotate("text", x = 1:12, y = sum.crop$asymp.UCL+0.5, label = gsub(" ", "", sum.crop$.group), crop=8, position = position_dodge(0.5))+
  labs(x ="" , y = "", title = NULL)+
  theme_bw()+
  theme(axis.title = element_text(lineheight=.8, hjust=0.5, size = 36, face = "bold"),
        text = element_text(family = "serif", size = 28),
        panel.background = element_rect(linetype = "blank", fill = NA),
        panel.grid = element_blank())

grph.crop

png(filename = "crop_effects.png",
    width = 30, height = 20, units = "cm", 
    bg = "white", res = 600)

grph.crop
dev.off()



rm(list= ls()[!(ls() %in% c("area.change", "area.loss", "ft1A", "ft1B", "ft2A", "ft2B", "sum.size", "sum.sociality", "sum.nest", "sum.crop"))])
#############################################################################################################################
#############################################################################################################################

listagem_ensemble<-list.files(path = "bin", pattern = ".img$", full.names = T, recursive = T)
head(listagem_ensemble)


PA<-shapefile("C:\\Users\\b0040\\Documents\\shapes\\base\\Politico\\DivisaÌfo PoliÌtica BR\\PA.shp")
plot(PA)


###################
###   total   #####
###################

ensemble<-stack(listagem_ensemble)
head(names(ensemble))
plot(ensemble[[1:5]])



scen_list<-c("current", "Futuro1A")

dir.create("sprich")

#scen <- "current"
for(scen in scen_list){
  scenarioX <- ensemble[[grep(scen, names(ensemble), value = T)]]
  multispecies<-sum(scenarioX)
  #plot(multispecies)
  
  multispecies_croped <- crop(multispecies, extent(PA), snap="out")
  PA_raster <- rasterize(PA, multispecies_croped) 
  multispecies_final <- mask(x=multispecies_croped, mask=PA_raster)
  #plot(multispecies_final)
  
  writeRaster(multispecies_final, paste("sprich/Total", scen, sep="_"), format="HFA")
  cat(scen, "pronto..\n")
}


###################
###   area    #####
###################
sp.list1 <- as.vector(as.character(area.loss$apelido[which(area.loss$area.cat!="wide")]))

k <-1; vet <- NULL
for(i in 1:length(listagem_ensemble)){
  split <- strsplit(listagem_ensemble[i],"/")[[1]]
  temp <- split[length(split)]
  term <- paste(strsplit(temp,"_")[[1]][1], strsplit(temp,"_")[[1]][2], sep = "_")
  cat(term,"\n")
  if((term %in% sp.list1) == TRUE){
    vet[k] <- i
    k <- k +1
  }
}


area_listagem_ensemble <- listagem_ensemble[vet]

area_ensemble<-stack(area_listagem_ensemble)
head(names(area_ensemble))
plot(area_ensemble[[1:5]])



#scen <- "current"
for(scen in scen_list){
  scenarioX <- area_ensemble[[grep(scen, names(area_ensemble), value = T)]]
  multispecies<-sum(scenarioX)
  #plot(multispecies)
  
  multispecies_croped <- crop(multispecies, extent(PA), snap="out")
  PA_raster <- rasterize(PA, multispecies_croped) 
  multispecies_final <- mask(x=multispecies_croped, mask=PA_raster)
  #plot(multispecies_final)
  
  writeRaster(multispecies_final, paste("sprich/RestricArea", scen, sep="_"), format="HFA")
  cat(scen, "pronto..\n")
}


###################
###   crop     ###
###################
sp.list2 <- as.vector(as.character(area.loss$apelido[which(area.loss$crop.pollination=="yes")]))

k <-1; vet <- NULL
for(i in 1:length(listagem_ensemble)){
  split <- strsplit(listagem_ensemble[i],"/")[[1]]
  temp <- split[length(split)]
  term <- paste(strsplit(temp,"_")[[1]][1], strsplit(temp,"_")[[1]][2], sep = "_")
  cat(term,"\n")
  if((term %in% sp.list2) == TRUE){
    vet[k] <- i
    k <- k +1
  }
}


crop_listagem_ensemble <- listagem_ensemble[vet]

crop_ensemble<-stack(crop_listagem_ensemble)
head(names(crop_ensemble))
plot(crop_ensemble[[1:5]])



#scen <- "current"
for(scen in scen_list){
  scenarioX <- crop_ensemble[[grep(scen, names(crop_ensemble), value = T)]]
  multispecies<-sum(scenarioX)
  #plot(multispecies)
  
  multispecies_croped <- crop(multispecies, extent(PA), snap="out")
  PA_raster <- rasterize(PA, multispecies_croped) 
  multispecies_final <- mask(x=multispecies_croped, mask=PA_raster)
  #plot(multispecies_final)
  
  writeRaster(multispecies_final, paste("sprich/CropPollinat", scen, sep="_"), format="HFA")
  cat(scen, "pronto..\n")
}
#############################################################################################################################
#############################################################################################################################

persist <- data.frame(do.call(rbind, strsplit(names(ensemble),"_")))
persist[,4]<-"stay"
colnames(persist) <- c("Genus", "specie", "Scen", "Persist")
persist$apelido <- paste(persist$Genus, persist$specie, sep = "_")



espaco <- shapefile("sprich/carajas.shp")
crop.persist <- crop(ensemble, extent(espaco), snap="out")
ras.espaco <- rasterize(espaco, crop.persist)
crop.persist.final <-mask(x=crop.persist, mask = ras.espaco)



j=1
for(i in names(crop.persist.final)){
  
  specieX<-crop.persist.final[[grep(i, names(crop.persist.final), value = T)]]
  
  if (maxValue(specieX) == F)  { persist$Persist[j] <- "out" }
  j=j+1
  
  cat(i, "pronto..\n") 
  
}

write.csv(persist, "sprich/persist_carajas.csv", row.names = F) 
