#######correlation between bacterial host and virus within every treatment
#####organize file first
bac.mag.ra <- read.csv("virus_host_cor/bac.mag.csv")
head(bac.mag.ra);dim(bac.mag.ra)
virus.otu.ra <- read.csv("virus_host_cor/t.votu.tab.csv")
head(virus.otu.ra);dim(virus.otu.ra)
colnames(virus.otu.ra)
virus.otu.ra.sub <- virus.otu.ra[,c(1,182,15,153,4,20,49,64,94,154,204,82,116,164,253,162,107,149)]
head(virus.otu.ra.sub)
#c("vOTU_181","vOTU_14","vOTU_152","vOTU_3","vOTU_19","vOTU_48","vOTU_63","vOTU_93","vOTU_153","vOTU_203","vOTU_81","vOTU_115","vOTU_163","vOTU_252","vOTU_161","vOTU_106","vOTU_148"))
##combine two dataset

virus.bac.tab <- left_join(virus.otu.ra.sub,bac.mag.ra,by=c("X"="X"))
head(virus.bac.tab);dim(virus.bac.tab)
virus.bac.df <- mutate(virus.bac.tab, Treat=c(rep("NCNTN0",3),rep("NCNTN60",3),rep("VNTN0",3),rep("VNTN60",3)), 
                       Nitrogen=c(rep("N0",3),rep("N60",3),rep("N0",3),rep("N60",3)),
                       Cover=c(rep("NC",6),rep("V",6)),
                       block=c(rep(c("b1","b2","b4"),4)))
head(virus.bac.df)
virus.bac.df.1 <- virus.bac.df[,c(1,32:35,2:31)];head(virus.bac.df)
colnames(virus.bac.df.1)[1] <- c("Sample")
head(virus.bac.df.1)

#Pearson's correlation richness between bacteria and viruses 
#let's test the treatment effect first on bac and virus, seperately
library(car)
library(fBasics)
library(MASS)
library(lme4)
#bac
bac.112.7.mod <- lmer(sqrt(bin112.7) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
bac.112.7.mod <- lm(sqrt(bin112.7) ~ Cover* Nitrogen, data = virus.bac.df.1)
Anova(bac.112.7.mod)
summary(bac.112.7.mod)#nitrogen
coef(bac.112.7.mod)

bac.37.7.mod <- lmer(sqrt(bin37.7) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.37.7.mod)
summary(bac.37.7.mod)#nitrogen
coef(bac.37.7.mod)

bac.37.3.mod <- lmer(sqrt(bin37.3) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
#bac.37.3.mod <- lm((bin37.3) ~ Cover* Nitrogen, data = virus.bac.df.1)

Anova(bac.37.3.mod)
summary(bac.37.3.mod)
coef(bac.37.3.mod)

bac.16.2.mod <- lmer(sqrt(bin16.2) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.16.2.mod)
summary(bac.16.2.mod)#nitrogen
coef(bac.16.2.mod)
###
bac.62.2.mod <- lmer(sqrt(bin62.2) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.62.2.mod)
summary(bac.62.2.mod)
coef(bac.62.2.mod)

bac.16.1.mod <- lmer(sqrt(bin16.1) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.16.1.mod)
summary(bac.16.1.mod)#significant on nitrogen
coef(bac.16.1.mod)

bac.37.1.mod <- lmer(sqrt(bin37.1)~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.37.1.mod)
summary(bac.37.1.mod)
coef(bac.37.1.mod)

bac.35.6.mod <- lmer(sqrt(bin35.6) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.35.6.mod)
summary(bac.35.6.mod)
coef(bac.35.6.mod)

bac.62.1.mod <- lmer(sqrt(bin62.1) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.62.1.mod)
summary(bac.62.1.mod)
coef(bac.62.1.mod)

bac.9.3.mod <- lmer(sqrt(bin9.3)~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.9.3.mod)
summary(bac.9.3.mod)###nitrogen
coef(bac.9.3.mod)

bac.51.3.mod <- lmer(sqrt(bin51.3) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.51.3.mod)
summary(bac.51.3.mod)##nitrogen
coef(bac.51.3.mod)

bac.107.3.mod <- lmer(sqrt(bin107.3) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.107.3.mod)
summary(bac.107.3.mod)
coef(bac.107.3.mod)

bac.35.9.mod <- lmer(sqrt(bin35.9) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.35.9.mod)
summary(bac.35.9.mod)
coef(bac.35.9.mod)

###virus
# virus.181.mod <- lmer(sqrt(vOTU_181) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.181.mod)
# summary(virus.181.mod)
# coef(virus.181.mod)
# 
# virus.14.mod <- lmer(sqrt(vOTU_14) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.14.mod)
# summary(virus.14.mod)
# coef(virus.14.mod)
# 
# virus.152.mod <- lmer(vOTU_152 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.152.mod)
# summary(virus.152.mod)
# coef(virus.152.mod)
# 
# virus.3.mod <- lmer(vOTU_3 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.3.mod)
# summary(virus.3.mod)
# coef(virus.3.mod)

virus.19.mod <- lmer(sqrt(vOTU_19) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.19.mod)
summary(virus.19.mod)#nitrogen
coef(virus.19.mod)

# virus.48.mod <- lmer(vOTU_48 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.48.mod)
# summary(virus.48.mod)
# coef(virus.48.mod)
# 
# 
# virus.63.mod <- lmer(vOTU_63 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.63.mod)
# summary(virus.63.mod)
# coef(virus.63.mod)
# 
# 
# virus.93.mod <- lmer(vOTU_93 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.93.mod)
# summary(virus.93.mod)
coef(virus.93.mod)


virus.153.mod <- lmer(sqrt(vOTU_153) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.153.mod)
summary(virus.153.mod)
coef(virus.153.mod)

virus.203.mod <- lmer(sqrt(vOTU_203) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.203.mod)
summary(virus.203.mod)#nitrogen
coef(virus.203.mod)

virus.81.mod <- lmer(sqrt(vOTU_81) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.81.mod)
summary(virus.81.mod)#nitrogen
coef(virus.81.mod)

virus.115.mod <- lmer(sqrt(vOTU_115) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.115.mod)
summary(virus.115.mod)#nitrogen
coef(virus.115.mod)

virus.163.mod <- lmer(sqrt(vOTU_163) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.163.mod)
summary(virus.163.mod)#nitrogen
coef(virus.163.mod)

virus.252.mod <- lmer(sqrt(vOTU_252) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.252.mod)
anova(virus.252.mod)
summary(virus.252.mod)#nitrogen
coef(virus.252.mod)

virus.161.mod <- lmer(sqrt(vOTU_161) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
a <- lm(vOTU_161 ~ Cover* Nitrogen, data = virus.bac.df.1)
Anova(a)
stepAIC(a)
Anova(virus.161.mod)
summary(virus.161.mod)
coef(virus.161.mod)

# virus.106.mod <- lmer(vOTU_106 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.106.mod)
# summary(virus.106.mod)
# coef(virus.106.mod)
# 
# virus.148.mod <- lmer(vOTU_148 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
# Anova(virus.148.mod)
# summary(virus.148.mod)
# coef(virus.148.mod)




normalTest((virus.bac.df.1$vOTU_181),"sw")#not normal
normalTest(sqrt(virus.bac.df.1$bin112.7),"sw")

normalTest(sqrt(virus.bac.df.1$vOTU_14),"sw")#one value
normalTest(sqrt(virus.bac.df.1$vOTU_19),"da")#not normal
normalTest(sqrt(virus.bac.df.1$vOTU_161),"sw")#not normal
normalTest((virus.bac.df.1$bin35.9),"sw")#not normal

normalTest((virus.bac.df.1$vOTU_153),"sw")
normalTest((virus.bac.df.1$vOTU_203),"sw")#good
normalTest((virus.bac.df.1$vOTU_81),"sw")#good
normalTest((virus.bac.df.1$vOTU_115),"sw")#good
normalTest((virus.bac.df.1$vOTU_163),"sw")#good
normalTest((virus.bac.df.1$vOTU_252),"sw")#good

#bac
normalTest((virus.bac.df.1$bin16.1),"sw")#not bad
normalTest((virus.bac.df.1$bin9.3),"sw")#good
# a <- lm(bin16.1~Nitrogen*Cover,data = virus.bac.df.1)
# Anova(a)

# Calculating Pearson's product-moment correlation
# cor.test(virus.bac.df.1$vOTU_181,virus.bac.df.1$bin112.7, method = "pearson", conf.level = 0.95)
# a <- cor.test(virus.bac.df.1$vOTU_181,virus.bac.df.1$bin112.7, method = "pearson", conf.level = 0.95)
# a
########I am going to make bin9.3 figures with multiple different vOTUs with good normality
#subset data first and try to figure out if there are correlation bewtween virus and host on N0 and N60 seperately
#if you want to conduct a pearson correlation test you must fullfill all!! the preassumption test
#1.random sample 2. variables are continous data 3. data contains paired samples 4. Independence of observations 5. variables are approxmately normally distributed 6. a linear association exists 7. absence of outliers
virus.bac_N0 <- subset(virus.bac.df.1, Nitrogen=="N0");dim(virus.bac_N0)
virus.bac_N60 <- subset(virus.bac.df.1, Nitrogen=="N60");dim(virus.bac_N60)
N0_bin9.3_153 <-cor.test(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_153) 
N0_bin9.3_153#no sig
N60_bin9.3_153 <-cor.test(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_153) 
N60_bin9.3_153#no
#r=0.9694543, R2=r^2 which is coefficient of determination, R2 percentage of the variability in host is explaned by the variability in viral abundance, reverse also the same. 1-R2 of variance is explained by unknown factors
N0_bin9.3_203 <-cor.test(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_203) 
N0_bin9.3_203#sig
N60_bin9.3_203 <-cor.test(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_203) 
N60_bin9.3_203#sig


N0_bin9.3_81 <-cor.test(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_81) 
N0_bin9.3_81#sig
N60_bin9.3_81 <-cor.test(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_81) 
N60_bin9.3_81#sig

#all viral
N0_bin9.3_115 <-cor.test(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_115) 
N0_bin9.3_115#sig
N60_bin9.3_115 <-cor.test(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_115) 
N60_bin9.3_115#sig


N0_bin9.3_163 <-cor.test(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_163) 
N0_bin9.3_163#sig
N60_bin9.3_163 <-cor.test(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_163) 
N60_bin9.3_163#sig

N0_bin9.3_252 <-cor.test(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_252) 
N0_bin9.3_252#sig
N60_bin9.3_252 <-cor.test(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_252) 
N60_bin9.3_252#sig


library(ggplot2)
library(grid)

annot_N0_bin9.3_153 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_153), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
#annot_N60_bin9.3_153 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_153), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N0_bin9.3_153.p = grobTree(textGrob(paste("P < 0.05"), x = 0.05, y = 0.97, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))

  bin9.3_153 <- ggplot(virus.bac.df.1, aes(Nitrogen,vOTU_153))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_153)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_153)+
  #  annotation_custom(annot_N60_bin9.3_153)+
  annotation_custom(annot_N0_bin9.3_153.p)+
  theme_bw()


annot_N0_bin9.3_203 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_203), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N60_bin9.3_203 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_203), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N0_bin9.3_203.p = grobTree(textGrob(paste("P < 0.01"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))

bin9.3_203 <- ggplot(virus.bac.df.1, aes(vOTU_203,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_203)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_203)+
  annotation_custom(annot_N60_bin9.3_203)+
  annotation_custom(annot_N0_bin9.3_203.p)+
  theme_bw()
#  scale_x_continuous(name = "Host abundance", limits = c(0, 30), breaks = seq(0, 30, 5)) + 
#  scale_y_continuous(name = "Viral abundance", limits = c(3000,5000), breaks = seq(3000,5000,500)) + 
#  annotation_custom(annot_N0) + 
# theme(plot.title = element_text(hjust = 0.5),
#       panel.background = element_blank(),
#       axis.line = element_line(color="black"),
#       axis.line.x = element_line(color="black"))

annot_N0_bin9.3_81 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_81), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N60_bin9.3_81 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_81), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N0_bin9.3_81.p = grobTree(textGrob(paste("P < 0.01"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))

bin9.3_81 <- ggplot(virus.bac.df.1, aes(vOTU_81,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_81)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_81)+
  annotation_custom(annot_N60_bin9.3_81)+
  annotation_custom(annot_N0_bin9.3_81.p)+
  theme_bw()

annot_N0_bin9.3_115 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_115), 3) ), x = 0.2, y = 0.55, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N60_bin9.3_115 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_115), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N0_bin9.3_115.p = grobTree(textGrob(paste("P < 0.05"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))

bin9.3_115 <- ggplot(virus.bac.df.1, aes(vOTU_115,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_115)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_115)+
  annotation_custom(annot_N60_bin9.3_115)+
  annotation_custom(annot_N0_bin9.3_115.p)+
  theme_bw()

annot_N0_bin9.3_163 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_163), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N60_bin9.3_163 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_163), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N0_bin9.3_163.p = grobTree(textGrob(paste("P < 0.01"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))

bin9.3_163 <- ggplot(virus.bac.df.1, aes(vOTU_163,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_163)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_163)+
  annotation_custom(annot_N60_bin9.3_163)+
  annotation_custom(annot_N0_bin9.3_163.p)+
  theme_bw()

annot_N0_bin9.3_252 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_252), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N60_bin9.3_252 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_252), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))
annot_N0_bin9.3_252.p = grobTree(textGrob(paste("P < 0.05"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 11, fontface = "plain")))

bin9.3_252 <- ggplot(virus.bac.df.1, aes(vOTU_252,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_252)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_252)+
  annotation_custom(annot_N60_bin9.3_252)+
  annotation_custom(annot_N0_bin9.3_252.p)+
  theme_bw()
