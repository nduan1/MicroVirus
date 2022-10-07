library(tidyverse)
bin112.3 <- read.csv("dataset/bin.112.3.contigs.tsv", header = F)
bin16.10 <- read.csv("dataset/bin.16.10.contigs.tsv", header = F)
bin16.4<- read.csv("dataset/bin.16.4.contigs.tsv", header = F)
bin16.5 <- read.csv("dataset/bin.16.5.contigs.tsv", header = F)
bin35.3 <- read.csv("dataset/bin.35.3.contigs.tsv", header = F)
bin63.4 <- read.csv("dataset/bin.63.4.contigs.tsv", header = F)
bin79.8 <- read.csv("dataset/bin.79.8.contigs.tsv", header = F)

bins_names <- c("bin112.3","bin16.10","bin16.4","bin16.5","bin35.3",
                "bin63.4","bin79.8")

bin112.3.df <-mutate(bin112.3,bin=c(rep("bin112.3",dim(bin112.3)[1])))
bin112.3.df <- bin112.3.df[,c(2,1)]

bin16.10.df <-mutate(bin16.10,bin=c(rep("bin16.10",dim(bin16.10)[1])))
bin16.10.df <- bin16.10.df[,c(2,1)]

bin16.4.df <-mutate(bin16.4,bin=c(rep("bin16.4",dim(bin16.4)[1])))
bin16.4.df <- bin16.4.df[,c(2,1)]

bin16.5.df <-mutate(bin16.5,bin=c(rep("bin16.5",dim(bin16.5)[1])))
bin16.5.df <- bin16.5.df[,c(2,1)]

bin35.3.df <-mutate(bin35.3,bin=c(rep("bin35.3",dim(bin35.3)[1])))
bin35.3.df <- bin35.3.df[,c(2,1)]

bin63.4.df <-mutate(bin63.4,bin=c(rep("bin63.4",dim(bin63.4)[1])))
bin63.4.df <- bin63.4.df[,c(2,1)]

bin79.8.df <-mutate(bin79.8,bin=c(rep("bin79.8",dim(bin79.8)[1])))
bin79.8.df <- bin79.8.df[,c(2,1)]

all.mag.contigs <- rbind(bin112.3.df,bin16.10.df,bin16.4.df,bin16.5.df,bin35.3.df,
                         bin63.4.df,bin79.8.df)

colnames(all.mag.contigs) <- c("bin","contigs")
head(all.mag.contigs)
dim(all.mag.contigs)
contig.cover.tab <- read.csv("dataset/mag.coverage.table.tsv",sep = "\t",)
dim(contig.cover.tab)
head(contig.cover.tab)
old.new.header <- read.csv("dataset/clean.header.tsv",sep = "\t",header = F)
dim(old.new.header)
head(old.new.header)
old.new.header.tab <- separate(old.new.header,1,sep=" ",into = c("old","flag","multi","length"))# seperate by white space
dim(old.new.header.tab)
head(old.new.header.tab)
new.header <- old.new.header.tab[,-c(2,3)]
colnames(new.header) <- c("old","length","new")
cov.name.list <- c(contig.cover.tab$X.contig)
cov.name.sub <- left_join(as.data.frame(all.mag.contigs),as.data.frame(new.header), by=c("contigs"="old"))
dim(cov.name.sub)
head(cov.name.sub)
bin.cov.tab <- right_join(as.data.frame(cov.name.sub),as.data.frame(contig.cover.tab),by=c("new"="X.contig"))
dim(bin.cov.tab)
head(bin.cov.tab)
tail(bin.cov.tab)
bin.cov.tab.clean <- subset(bin.cov.tab,contigs!="NA")
dim(bin.cov.tab.clean)
tail(bin.cov.tab.clean)
bin.cov.tab.clean.df <-bin.cov.tab.clean[,-3] 
head(bin.cov.tab.clean.df)

########now we have clean bins with contigs length and coverage table in every sample
# calculate the normalization factor for the bins (= contig_len / average_bin_len)
#I made a little change manually on genomeInfo.csv
genome.info <- read.csv("dataset/genomeInfo.csv",sep=",")
head(genome.info)
bin.cov.clean.binlength <- left_join(as.data.frame(bin.cov.tab.clean.df),
                                     as.data.frame(genome.info),by=c("bin"="bins"))
head(bin.cov.clean.binlength)
colnames(bin.cov.clean.binlength)
bin.cov.clean.binlength.df <-bin.cov.clean.binlength[,-c(17,18,19,21)]
colnames(bin.cov.clean.binlength.df) <- c("bin","contigs","new_contig_name","contig_len",
                                          "NCNTN0_4","VNTN0_4","VNTN60_2","NCNTN0_1","VNTN0_1",
                                          "VNTN60_1","NCNTN60_1","VNTN0_2","NCNTN0_2","NCNTN60_4",
                                          "VNTN60_4","NCNTN60_2","bin_len")
head(bin.cov.clean.binlength.df)
colnames(bin.cov.clean.binlength.df)
bin.cov.clean.binlength.df <- bin.cov.clean.binlength.df[,c(1,17,2:16)]
head(bin.cov.clean.binlength.df)

###############first time normalization coverage across different contigs length cover number * contiglen/avergae contig len for the bin
total_contigs_len <- aggregate(bin.cov.clean.binlength.df$contig_len, list(bin.cov.clean.binlength.df$bin), FUN=mean)
total_contigs_len
bin.cov.clean.df <- inner_join(bin.cov.clean.binlength.df,total_contigs_len,by=c("bin"="Group.1"))
head(bin.cov.clean.df)
bin.cov.clean.nor1 <- mutate(bin.cov.clean.df,nor1=contig_len/x)
head(bin.cov.clean.nor1)
nor.matrix.bac.bin <- matrix(nrow=0,ncol =12)
for (i in 1:dim(bin.cov.clean.nor1)[1]){
  df <- bin.cov.clean.nor1[i,19]*bin.cov.clean.nor1[i,c(6:17)]
  nor.matrix.bac.bin <- rbind(nor.matrix.bac.bin,df)
}
dim(nor.matrix.bac.bin)
head(nor.matrix.bac.bin)
normalied.bin.cov <-cbind(bin.cov.clean.nor1[,c(1:5,18,19)], nor.matrix.bac.bin)
head(normalied.bin.cov)
colnames(normalied.bin.cov) <- c("bin","bin_len","contigs","new_contig_name","contig_len","avg_contig.len","nor1",
                                 "NCNTN0_4","VNTN0_4","VNTN60_2","NCNTN0_1","VNTN0_1",
                                 "VNTN60_1","NCNTN60_1","VNTN0_2","NCNTN0_2","NCNTN60_4",
                                 "VNTN60_4","NCNTN60_2")

head(normalied.bin.cov)
#view(normalied.bin.cov)
seq.depth <- read.csv("dataset/seq.depth.csv", row.names = 1)
head(seq.depth)
normalied.bin.cov.longer <- pivot_longer(normalied.bin.cov,cols =colnames(normalied.bin.cov)[8:19])
head(normalied.bin.cov.longer)
normalied.bin.cov.longer.reads.nor <- inner_join(normalied.bin.cov.longer,seq.depth,by=c("name"="sample"))
head(normalied.bin.cov.longer.reads.nor)
#####second time normalization for the different library size of reads
nor2.matrix.bac.bin <- mutate(normalied.bin.cov.longer.reads.nor,nor2.cov=value*mutiplier)
head(nor2.matrix.bac.bin)
nor2.matrix.bac.bin$nor2.cov
######Then avergae the coverage for each bin
####clean our table and do the average of each contigs and finally get the coverage table of each bin
###clean
colnames(nor2.matrix.bac.bin)
clean.bin.coverage <- nor2.matrix.bac.bin[,c(1,4,6,8,13)]
head(clean.bin.coverage)
clean.bin.coverage.short <-pivot_wider(clean.bin.coverage,names_from  = name, values_from = nor2.cov) 
head(clean.bin.coverage.short)
#view(clean.bin.coverage.short)
clean.bin.coverage.short.mean <- aggregate(clean.bin.coverage.short, list(clean.bin.coverage.short$bin), mean)
head(clean.bin.coverage.short.mean)
clean.bin.otu <-clean.bin.coverage.short.mean[,-c(2,3,4)] 
head(clean.bin.otu)
#####then filter the coverage less than 0.25*, and turn it into 0
# clean.bin.otu.1$NCNTN0_4[clean.bin.otu.1$NCNTN0_4<0.25] <- 0
# clean.bin.otu.1$NCNTN0_4
replace.df <- function(file){
  for (i in 2:dim(file)[2]){
    for (j in 1:dim(file)[1]){
      a <- file[j,i] < 0.25
      if (a==TRUE){
        file[j,i] <-  0
      }
    }
  }
  return(file)
}

clean.bin.otu.filtered <- replace.df(clean.bin.otu)
View(clean.bin.otu.filtered)
######now you have a clean archaea_mag_otu table
######now we have normalized average number of reads per base-pair per contigs across samples
#next we plan to generate the phyloseq
library(phyloseq)
library(ggplot2)
library(tidyverse)
meta.data <- read.csv("dataset/meta_data.csv",header = T, row.names = 1)
meta.data.F <- read.csv("dataset/meta_data.csv",header = T)
head(meta.data)
arc.mag.sample.data <- sample_data(meta.data)
head(clean.bin.otu.filtered)
row.names(clean.bin.otu.filtered) <- clean.bin.otu.filtered[,1]
head(clean.bin.otu.filtered)
arc.mag.otu <- otu_table(clean.bin.otu.filtered[,-1], taxa_are_rows = TRUE)
arc.mag.phyloseq <-phyloseq(arc.mag.otu, arc.mag.sample.data)
arc.mag.phyloseq
r_arc.mag.phyloseq<-transform_sample_counts(arc.mag.phyloseq,function(x) x/sum(x))
head(otu_table(r_arc.mag.phyloseq))
write.csv(t(otu_table(r_arc.mag.phyloseq)),"virus_host_cor/arc.mag.csv")
arc.mag.hp.df <- mutate(as.data.frame(otu_table(r_arc.mag.phyloseq)),bins=rownames(clean.bin.otu.filtered)) %>% pivot_longer(cols = colnames(otu_table(r_arc.mag.phyloseq)))
##reorder the messy rows
write.csv(arc.mag.hp.df,"dataset/arc.mag.hp.df.csv")
write.csv(t(otu_table(arc.mag.phyloseq)),"dataset/arc.mag.phylotu.csv")#I reorder the files
arc.mag.hp.df.1 <- read.csv("dataset/arc.mag.hp.df.1.csv")
tail(arc.mag.hp.df.1)
arc.mag.hp <- inner_join(arc.mag.hp.df.1,meta.data.F,by=c("name"="X"))
head(arc.mag.hp)
arc.mag.hp.mean <- arc.mag.hp %>% group_by(bins,Treat) %>% summarise_at(vars(value),list(mean.value=mean))
arc.mag.hp.tab <- left_join(arc.mag.hp,arc.mag.hp.mean,by=c("bins"="bins","Treat"="Treat"))
head(arc.mag.hp.tab)

bins.order <- c("bin35.3","bin63.4",
                "bin16.4","bin16.5","bin16.10","bin79.8","bin112.3")#make sure the order is the same to the arc.mag.tree
arc.mag.hp.tab$bins <- factor(arc.mag.hp.tab$bins, levels = as.factor(bins.order))

tiff('figures/arc.mag.hp.tiff', units="in", width=5, height=10, res=300)
arc.hp <- ggplot(arc.mag.hp.tab, aes(x =Treat, y = bins,fill=mean.value)) + 
  geom_tile(aes(fill = round(mean.value,3)))+
  geom_text(aes(label = round(mean.value, 3)))+
  #  geom_hline(yintercept = 3.5, size=1.5)+
  #geom_vline(xintercept = c(3.5,6.5,9.5), size=1)+
  #geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,11.5), size=0.2)+
  #  geom_hline(yintercept = c(1.5,2.5,4.5,5.5,6.5), size=0.2)+
  xlab(label = "Treat")+
  ylab(label = "MAGs")+
  #  scale_fill_gradient(name="Relative abundance",low="#3794bf", high="#a34c1d")+#low = "violetred", high = "aquamarine"
  scale_fill_gradient(name="MAGs",low="#f7f9fa", high="#174d80")+#low = "violetred", high = "aquamarine"
  
  #  facet_grid(~ Depth)+
  theme_bw(base_size=15,base_family = "",) 
dev.off()
head(arc.mag.hp.tab)

virus.otu.relative.ab <- read.csv("dataset/t.votu.tab.csv")

head(virus.otu.relative.ab);dim(virus.otu.relative.ab)
virus.otu.relative.sub <- virus.otu.relative.ab[,c(1,138,200,139,167,57,211,244,238,192,235,131,190,130,233,121,160,143,129,252)]
virus.otu.relative.sub
virus.otu.relative <- mutate(virus.otu.relative.sub, Treat=c(rep("NCNTN0",3),rep("NCNTN60",3),rep("VNTN0",3),rep("VNTN60",3)), 
                             Nitrogen=c(rep("N0",3),rep("N60",3),rep("N0",3),rep("N60",3)),
                             Cover=c(rep("NC",6),rep("V",6)),
                             block=c(rep(c("b1","b2","b4"),4)))
colnames(virus.otu.relative)
virus.otu.relative.1 <- virus.otu.relative[,c(1,21:24,2:20)]
colnames(virus.otu.relative.1)[1] <- "Sample" 
head(virus.otu.relative.1)
virus.otu.relative.df <- pivot_longer(virus.otu.relative.1, cols = c("vOTU_137","vOTU_199","vOTU_138","vOTU_166","vOTU_56",
                                                                     "vOTU_210","vOTU_243","vOTU_237","vOTU_191","vOTU_234",
                                                                     "vOTU_130","vOTU_189","vOTU_129","vOTU_232","vOTU_120",
                                                                     "vOTU_159","vOTU_142","vOTU_128","vOTU_251"))
head(virus.otu.relative.df);dim(virus.otu.relative.df)

virus.hp.mean <- virus.otu.relative.df %>% group_by(name,Treat) %>% 
  summarise_at(vars(value),list(mean.value=mean))
virus.hp.df <- inner_join(virus.otu.relative.df,virus.hp.mean,by=c("name"="name","Treat"="Treat"))
head(virus.hp.df)

virus.order <- c("vOTU_251","vOTU_189","vOTU_130","vOTU_129","vOTU_128","vOTU_191",
                 "vOTU_237","vOTU_243","vOTU_210","vOTU_56","vOTU_234","vOTU_166",
                 "vOTU_138","vOTU_232","vOTU_159","vOTU_120","vOTU_142","vOTU_199","vOTU_137")
#make sure the order is the same to the bac.mag.tree
# 
virus.hp.df$name <- factor(virus.hp.df$name, levels = as.factor(virus.order))

virus.hp <- ggplot(virus.hp.df, aes(x =Treat, y = name,fill=round(mean.value,3))) + 
  geom_tile(aes(fill = mean.value))+
  geom_text(aes(label = round(mean.value, 3)))+
  #  geom_hline(yintercept = 3.5, size=1.5)+
  #geom_vline(xintercept = c(3.5,6.5,9.5), size=1)+
  #geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,15.5), size=0.2)+
  #  geom_hline(yintercept = c(1.5,2.5,4.5,5.5,6.5), size=0.2)+
  xlab(label = "Treat")+
  ylab(label = "vOTUs")+
  #  scale_fill_gradient(name="Relative abundance",low="#3794bf", high="#a34c1d")+#low = "violetred", high = "aquamarine"
  scale_fill_gradient(name="vOTU",low="#fafafa", high="#b06690")+#low = "violetred", high = "aquamarine"
  
  #  facet_grid(~ Depth)+
  theme_bw(base_size=15,base_family = "",)
# theme(legend.position="none",panel.background = element_rect(fill = "white",color="black"),
#       panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#       plot.title=(element_text(size=15,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
#       text=element_text(family = "Arial",face="plain",size=13),
#       panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#       axis.text.x=element_text(size=15,family="Arial",face = "plain"),
#       axis.text.y=element_text(size=15,family="Arial",face = "bold"),
#       axis.title=element_text(size = 15,face="plain",family = "Arial",vjust = 1),
#       legend.title = element_text(size=15,face="plain",family = "Arial"),
#       legend.text = (element_text(size=15,family = "Arial")))
library("ggpubr")
tiff('figures/arc.virus.hp.tiff', units="in", width=9, height=7, res=300)
ggarrange(arc.hp,virus.hp,labels = c("", ""),ncol = 2, nrow = 1,common.legend = F,legend = "bottom") 
dev.off()

#######correlation between bacterial host and virus within every treatment
#####organize file first
arc.mag.ra <- read.csv("dataset/arc.mag.phylotu.csv")
head(arc.mag.ra);dim(arc.mag.ra)
virus.otu.ra <- read.csv("dataset/virus.phyloseq.cov.csv")
head(virus.otu.ra);dim(virus.otu.ra)
colnames(virus.otu.ra)
virus.otu.ra.sub <- virus.otu.ra[,c(1,138,200,139,167,57,211,244,238,192,235,131,190,130,233,121,160,143,129,252)]
head(virus.otu.ra.sub)
# vOTU_137
# vOTU_199
# vOTU_138
# vOTU_166
# vOTU_56
# vOTU_210
# vOTU_243
# vOTU_237
# vOTU_191
# vOTU_234
# vOTU_130
# vOTU_189
# vOTU_129
# vOTU_232
# vOTU_120
# vOTU_159
# vOTU_142
# vOTU_128
vOTU_251##combine two dataset

virus.arc.tab <- left_join(virus.otu.ra.sub,arc.mag.ra,by=c("X"="X"))#remove the bin16.4 and 16.5 due to there is no connection between these two
head(virus.arc.tab);dim(virus.arc.tab)
virus.arc.df <- mutate(virus.arc.tab, Treat=c(rep("NCNTN0",3),rep("NCNTN60",3),rep("VNTN0",3),rep("VNTN60",3)), 
                       Nitrogen=c(rep("N0",3),rep("N60",3),rep("N0",3),rep("N60",3)),
                       Cover=c(rep("NC",6),rep("V",6)),
                       block=c(rep(c("b1","b2","b4"),4)))
head(virus.arc.df)
virus.arc.df.1 <- virus.arc.df[,c(1,28:31,2:27)];head(virus.arc.df.1)
colnames(virus.arc.df.1)[1] <- c("Sample")
head(virus.arc.df.1)
#Pearson's correlation richness between bacteria and viruses 
#let's test the treatment effect first on arc and virus, seperately
library(car)
library(fBasics)
library(MASS)
library(lme4)
colnames(virus.arc.df.1)
#arc
arc.112.3.mod <- lmer(sqrt(bin112.3) ~ Cover* Nitrogen +(1|block), data = virus.arc.df.1)
Anova(arc.112.3.mod)#nitrogen
summary(arc.112.3.mod)
coef(arc.112.3.mod)

arc.16.10.mod <- lmer(sqrt(bin16.10) ~ Cover +(1|block), data = subset(virus.arc.df.1,Nitrogen=="N0"))
Anova(arc.16.10.mod)# cover and nitrogen
summary(arc.16.10.mod)
coef(arc.16.10.mod)

arc.16.4.mod <- lmer(sqrt(bin16.4) ~ Cover* Nitrogen +(1|block), data = virus.arc.df.1)
Anova(arc.16.4.mod)#
summary(arc.16.4.mod)
coef(arc.16.4.mod)

arc.16.5.mod <- lmer(sqrt(bin16.5)~ Cover*Nitrogen +(1|block), data = virus.arc.df.1)
Anova(arc.16.5.mod)#Cover
summary(arc.16.5.mod)
coef(arc.16.5.mod)

arc.35.3.mod <- lmer(sqrt(bin35.3) ~ Cover* Nitrogen +(1|block), data = virus.arc.df.1)
Anova(arc.35.3.mod)#nitrogen
summary(arc.35.3.mod)
coef(arc.35.3.mod)

arc.63.4.mod <- lmer(sqrt(bin63.4) ~ Cover* Nitrogen +(1|block), data = virus.arc.df.1)
Anova(arc.63.4.mod)#
summary(arc.63.4.mod)#nitrogen
coef(arc.63.4.mod)

arc.79.8.mod <- lmer(sqrt(bin79.8) ~ Cover* Nitrogen +(1|block), data = virus.arc.df.1)
Anova(arc.79.8.mod)#nitrogen
summary(arc.79.8.mod)
coef(arc.79.8.mod)

###virus
viral.name.list <- colnames(virus.arc.df.1)[6:24]
for (i in 1:19){
  viral.mod <- lmer(sqrt(virus.arc.df.1[,i+5]) ~ Cover* Nitrogen +(1|block), data = virus.arc.df.1)
  ANOVA <- Anova(viral.mod)
  print(viral.name.list[i])
  print(ANOVA)
}
viral.name.list 
# [1] "vOTU_137"Nnitrogen "vOTU_199" nitrogen"vOTU_138"nitrogen*cover "vOTU_166"nitrogen*cover
# [5] "vOTU_56"nitrogen*cover  "vOTU_210"nitrogen*cover "vOTU_243"nitrogen*cover "vOTU_237"nitrogen*cover
# [9] "vOTU_191" "vOTU_234"nitrogen*cover "vOTU_130" cover or nitrogen "vOTU_189" cover or nitrogen
# [13] "vOTU_129"nitrogen "vOTU_232" cover*nitrogen"vOTU_120"cover*nitrogen "vOTU_159"
# [17] "vOTU_142" "vOTU_128"cover*nitrogen "vOTU_251 cover or nitrogen
#######normality test for host and virus#######
all.name.list <- colnames(virus.arc.df.1)[6:31]
for (i in 1:26){
  nol <- normalTest(sqrt(virus.arc.df.1[,i+5]),"sw")#sqrt transformed
  print(all.name.list[i])
  print(nol)
}
all.name.list

#bin112.3 
virus.arc_N0 <- subset(virus.arc.df.1, Nitrogen=="N0");dim(virus.arc_N0)
virus.arc_N60 <- subset(virus.arc.df.1, Nitrogen=="N60");dim(virus.arc_N60)
N0_bin112.3_137 <-cor.test(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_137)) 
N0_bin112.3_137#sig
N60_bin112.3_137 <-cor.test(sqrt(virus.arc_N60$bin112.3),sqrt(virus.arc_N60$vOTU_137))
N60_bin112.3_137#sig

N0_bin112.3_199 <-cor.test(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_199)) 
N0_bin112.3_199#sig
N60_bin112.3_199 <-cor.test(sqrt(virus.arc_N60$bin112.3),sqrt(virus.arc_N60$vOTU_199))
N60_bin112.3_199#sig

N0_bin112.3_138 <-cor.test(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_138)) 
N0_bin112.3_138
N60_bin112.3_138 <-cor.test(sqrt(virus.arc_N60$bin112.3),sqrt(virus.arc_N60$vOTU_138))
N60_bin112.3_138

N0_bin112.3_120 <-cor.test(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_120)) 
N0_bin112.3_120
N60_bin112.3_120 <-cor.test(sqrt(virus.arc_N60$bin112.3),sqrt(virus.arc_N60$vOTU_120))
N60_bin112.3_120


N0_bin112.3_159 <-cor.test(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_159)) 
N0_bin112.3_159
N60_bin112.3_159 <-cor.test(sqrt(virus.arc_N60$bin112.3),sqrt(virus.arc_N60$vOTU_159))
N60_bin112.3_159

N0_bin112.3_142 <-cor.test(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_142)) 
N0_bin112.3_142
N60_bin112.3_142 <-cor.test(sqrt(virus.arc_N60$bin112.3),sqrt(virus.arc_N60$vOTU_142))
N60_bin112.3_142

N0_bin112.3_234 <-cor.test(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_234)) 
N0_bin112.3_234
N60_bin112.3_234 <-cor.test(sqrt(virus.arc_N60$bin112.3),sqrt(virus.arc_N60$vOTU_234))
N60_bin112.3_234

#bin16.10
N0_bin16.10_234 <-cor.test(sqrt(virus.arc_N0$bin16.10),sqrt(virus.arc_N0$vOTU_234)) 
N0_bin16.10_234#sig
N60_bin16.10_234 <-cor.test(sqrt(virus.arc_N60$bin16.10),sqrt(virus.arc_N60$vOTU_234))
N60_bin16.10_234#sig

N0_bin16.10_159 <-cor.test(sqrt(virus.arc_N0$bin16.10),sqrt(virus.arc_N0$vOTU_159)) 
N0_bin16.10_159#sig
N60_bin16.10_159 <-cor.test(sqrt(virus.arc_N60$bin16.10),sqrt(virus.arc_N60$vOTU_159))
N60_bin16.10_159

N0_bin16.10_191 <-cor.test(sqrt(virus.arc_N0$bin16.10),sqrt(virus.arc_N0$vOTU_191)) 
N0_bin16.10_191#sig
N60_bin16.10_191 <-cor.test(sqrt(virus.arc_N60$bin16.10),sqrt(virus.arc_N60$vOTU_191))
N60_bin16.10_191

N0_bin16.10_210 <-cor.test(sqrt(virus.arc_N0$bin16.10),sqrt(virus.arc_N0$vOTU_210)) 
N0_bin16.10_210#sig
N60_bin16.10_210 <-cor.test(sqrt(virus.arc_N60$bin16.10),sqrt(virus.arc_N60$vOTU_210))
N60_bin16.10_210#sig

N0_bin16.10_138 <-cor.test(sqrt(virus.arc_N0$bin16.10),sqrt(virus.arc_N0$vOTU_138)) 
N0_bin16.10_138#sig
N60_bin16.10_138 <-cor.test(sqrt(virus.arc_N60$bin16.10),sqrt(virus.arc_N60$vOTU_138))
N60_bin16.10_138#sig

N0_bin16.10_166 <-cor.test(sqrt(virus.arc_N0$bin16.10),sqrt(virus.arc_N0$vOTU_166)) 
N0_bin16.10_166#sig
N60_bin16.10_166 <-cor.test(sqrt(virus.arc_N60$bin16.10),sqrt(virus.arc_N60$vOTU_166))
N60_bin16.10_166#sig



#for votu or bins have cover and fertlization interaction effect
virus.arc_ncntn0 <- subset(virus.arc_N0, Cover=="NC");dim(virus.arc_ncntn0)
virus.arc_vntn0 <- subset(virus.arc_N0, Cover=="V");dim(virus.arc_vntn0)

virus.arc_ncntn60 <- subset(virus.arc_N60, Cover=="NC");dim(virus.arc_ncntn60)
virus.arc_vntn60 <- subset(virus.arc_N60, Cover=="V");dim(virus.arc_vntn60)

ncntn0_bin16.10_159 <-cor.test(sqrt(virus.arc_ncntn0$bin16.10),sqrt(virus.arc_ncntn0$vOTU_159)) 
ncntn0_bin16.10_159
vntn0_bin16.10_159 <-cor.test(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_159))
vntn0_bin16.10_159
ncntn60_bin16.10_159 <-cor.test(sqrt(virus.arc_ncntn60$bin16.10),sqrt(virus.arc_ncntn60$vOTU_159)) 
ncntn60_bin16.10_159
vntn60_bin16.10_159 <-cor.test(sqrt(virus.arc_vntn60$bin16.10),sqrt(virus.arc_vntn60$vOTU_159))
vntn60_bin16.10_159

ncntn0_bin112.3_166 <-cor.test(sqrt(virus.arc_ncntn0$bin112.3),sqrt(virus.arc_ncntn0$vOTU_166)) 
ncntn0_bin112.3_166
vntn0_bin112.3_166 <-cor.test(sqrt(virus.arc_vntn0$bin112.3),sqrt(virus.arc_vntn0$vOTU_166))
vntn0_bin112.3_166
ncntn60_bin112.3_166 <-cor.test(sqrt(virus.arc_ncntn60$bin112.3),sqrt(virus.arc_ncntn60$vOTU_166)) 
ncntn60_bin112.3_166
vntn60_bin112.3_166 <-cor.test(sqrt(virus.arc_vntn60$bin112.3),sqrt(virus.arc_vntn60$vOTU_166))
vntn60_bin112.3_166


ncntn0_bin16.10_56 <-cor.test(sqrt(virus.arc_ncntn0$bin16.10),sqrt(virus.arc_ncntn0$vOTU_56)) 
ncntn0_bin16.10_56#sig
vntn0_bin16.10_56 <-cor.test(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_56))
vntn0_bin16.10_56#sig
ncntn60_bin16.10_56 <-cor.test(sqrt(virus.arc_ncntn60$bin16.10),sqrt(virus.arc_ncntn60$vOTU_56)) 
ncntn60_bin16.10_56
vntn60_bin16.10_56 <-cor.test(sqrt(virus.arc_vntn60$bin16.10),sqrt(virus.arc_vntn60$vOTU_56))
vntn60_bin16.10_56

ncntn0_bin16.10_243 <-cor.test(sqrt(virus.arc_ncntn0$bin16.10),sqrt(virus.arc_ncntn0$vOTU_243)) 
ncntn0_bin16.10_243
vntn0_bin16.10_243 <-cor.test(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_243))
vntn0_bin16.10_243#sig
ncntn60_bin16.10_243 <-cor.test(sqrt(virus.arc_ncntn60$bin16.10),sqrt(virus.arc_ncntn60$vOTU_243)) 
ncntn60_bin16.10_243
vntn60_bin16.10_243 <-cor.test(sqrt(virus.arc_vntn60$bin16.10),sqrt(virus.arc_vntn60$vOTU_243))
vntn60_bin16.10_243

ncntn0_bin16.10_237 <-cor.test(sqrt(virus.arc_ncntn0$bin16.10),sqrt(virus.arc_ncntn0$vOTU_237)) 
ncntn0_bin16.10_237#sig
vntn0_bin16.10_237 <-cor.test(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_237))
vntn0_bin16.10_237#sig
ncntn60_bin16.10_237 <-cor.test(sqrt(virus.arc_ncntn60$bin16.10),sqrt(virus.arc_ncntn60$vOTU_237)) 
ncntn60_bin16.10_237
vntn60_bin16.10_237 <-cor.test(sqrt(virus.arc_vntn60$bin16.10),sqrt(virus.arc_vntn60$vOTU_237))
vntn60_bin16.10_237

ncntn0_bin16.10_128 <-cor.test(sqrt(virus.arc_ncntn0$bin16.10),sqrt(virus.arc_ncntn0$vOTU_128)) 
ncntn0_bin16.10_128
vntn0_bin16.10_128 <-cor.test(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_128))
vntn0_bin16.10_128#sig
ncntn60_bin16.10_128 <-cor.test(sqrt(virus.arc_ncntn60$bin16.10),sqrt(virus.arc_ncntn60$vOTU_128)) 
ncntn60_bin16.10_128
vntn60_bin16.10_128 <-cor.test(sqrt(virus.arc_vntn60$bin16.10),sqrt(virus.arc_vntn60$vOTU_128))
vntn60_bin16.10_128#sig

#bin35.3
virus.arc_N0 <- subset(virus.arc.df.1, Nitrogen=="N0");dim(virus.arc_N0)
virus.arc_N60 <- subset(virus.arc.df.1, Nitrogen=="N60");dim(virus.arc_N60)
virus.arc_NC <- subset(virus.arc.df.1, Cover=="NC");dim(virus.arc_NC)
virus.arc_V <- subset(virus.arc.df.1, Cover=="V");dim(virus.arc_V)

N0_bin35.3_130 <-cor.test(sqrt(virus.arc_N0$bin35.3),sqrt(virus.arc_N0$vOTU_130)) 
N0_bin35.3_130#sig
N60_bin35.3_130 <-cor.test(sqrt(virus.arc_N60$bin35.3),sqrt(virus.arc_N60$vOTU_130))
N60_bin35.3_130#sig

NC_bin35.3_130 <-cor.test(sqrt(virus.arc_NC$bin35.3),sqrt(virus.arc_NC$vOTU_130)) 
NC_bin35.3_130#sig
V_bin35.3_130 <-cor.test(sqrt(virus.arc_V$bin35.3),sqrt(virus.arc_V$vOTU_130))
V_bin35.3_130#sig

ncntn0_bin35.3_130 <-cor.test(sqrt(virus.arc_ncntn0$bin35.3),sqrt(virus.arc_ncntn0$vOTU_130)) 
ncntn0_bin35.3_130#sig
vntn0_bin35.3_130 <-cor.test(sqrt(virus.arc_vntn0$bin35.3),sqrt(virus.arc_vntn0$vOTU_130))
vntn0_bin35.3_130#sig
ncntn60_bin35.3_130 <-cor.test(sqrt(virus.arc_ncntn60$bin35.3),sqrt(virus.arc_ncntn60$vOTU_130)) 
ncntn60_bin35.3_130#sig
vntn60_bin35.3_130 <-cor.test(sqrt(virus.arc_vntn60$bin35.3),sqrt(virus.arc_vntn60$vOTU_130))
vntn60_bin35.3_130#sig


N0_bin35.3_189 <-cor.test(sqrt(virus.arc_N0$bin35.3),sqrt(virus.arc_N0$vOTU_189)) 
N0_bin35.3_189#sig
N60_bin35.3_189 <-cor.test(sqrt(virus.arc_N60$bin35.3),sqrt(virus.arc_N60$vOTU_189))
N60_bin35.3_189#sig

NC_bin35.3_189 <-cor.test(sqrt(virus.arc_NC$bin35.3),sqrt(virus.arc_NC$vOTU_189)) 
NC_bin35.3_189#sig
V_bin35.3_189 <-cor.test(sqrt(virus.arc_V$bin35.3),sqrt(virus.arc_V$vOTU_189))
V_bin35.3_189#sig

ncntn0_bin35.3_189 <-cor.test(sqrt(virus.arc_ncntn0$bin35.3),sqrt(virus.arc_ncntn0$vOTU_189)) 
ncntn0_bin35.3_189#sig
vntn0_bin35.3_189 <-cor.test(sqrt(virus.arc_vntn0$bin35.3),sqrt(virus.arc_vntn0$vOTU_189))
vntn0_bin35.3_189#sig
ncntn60_bin35.3_189 <-cor.test(sqrt(virus.arc_ncntn60$bin35.3),sqrt(virus.arc_ncntn60$vOTU_189)) 
ncntn60_bin35.3_189
vntn60_bin35.3_189 <-cor.test(sqrt(virus.arc_vntn60$bin35.3),sqrt(virus.arc_vntn60$vOTU_189))
vntn60_bin35.3_189#sig


N0_bin35.3_251 <-cor.test(sqrt(virus.arc_N0$bin35.3),sqrt(virus.arc_N0$vOTU_251)) 
N0_bin35.3_251#sig
N60_bin35.3_251 <-cor.test(sqrt(virus.arc_N60$bin35.3),sqrt(virus.arc_N60$vOTU_251))
N60_bin35.3_251#sig

NC_bin35.3_251 <-cor.test(sqrt(virus.arc_NC$bin35.3),sqrt(virus.arc_NC$vOTU_251)) 
NC_bin35.3_251#sig
V_bin35.3_251 <-cor.test(sqrt(virus.arc_V$bin35.3),sqrt(virus.arc_V$vOTU_251))
V_bin35.3_251#sig

ncntn0_bin35.3_251 <-cor.test(sqrt(virus.arc_ncntn0$bin35.3),sqrt(virus.arc_ncntn0$vOTU_251)) 
ncntn0_bin35.3_251#sig
vntn0_bin35.3_251 <-cor.test(sqrt(virus.arc_vntn0$bin35.3),sqrt(virus.arc_vntn0$vOTU_251))
vntn0_bin35.3_251
ncntn60_bin35.3_251 <-cor.test(sqrt(virus.arc_ncntn60$bin35.3),sqrt(virus.arc_ncntn60$vOTU_251)) 
ncntn60_bin35.3_251
vntn60_bin35.3_251 <-cor.test(sqrt(virus.arc_vntn60$bin35.3),sqrt(virus.arc_vntn60$vOTU_251))#many 0s
vntn60_bin35.3_251

#bin63.4
N0_bin63.4_129 <-cor.test(sqrt(virus.arc_N0$bin63.4),sqrt(virus.arc_N0$vOTU_129)) 
N0_bin63.4_129#sig
N60_bin63.4_129 <-cor.test(sqrt(virus.arc_N60$bin63.4),sqrt(virus.arc_N60$vOTU_129))
N60_bin63.4_129#sig

#bin79.8
N0_bin79.8_232 <-cor.test(sqrt(virus.arc_N0$bin79.8),sqrt(virus.arc_N0$vOTU_232)) 
N0_bin79.8_232#sig
N60_bin79.8_232 <-cor.test(sqrt(virus.arc_N60$bin79.8),sqrt(virus.arc_N60$vOTU_232))
N60_bin79.8_232#sig

ncntn0_bin79.8_120 <-cor.test(sqrt(virus.arc_ncntn0$bin79.8),sqrt(virus.arc_ncntn0$vOTU_120)) 
ncntn0_bin79.8_120#sig
vntn0_bin79.8_120 <-cor.test(sqrt(virus.arc_vntn0$bin79.8),sqrt(virus.arc_vntn0$vOTU_120))
vntn0_bin79.8_120
ncntn60_bin79.8_120 <-cor.test(sqrt(virus.arc_ncntn60$bin79.8),sqrt(virus.arc_ncntn60$vOTU_120)) 
ncntn60_bin79.8_120#sig
vntn60_bin79.8_120 <-cor.test(sqrt(virus.arc_vntn60$bin79.8),sqrt(virus.arc_vntn60$vOTU_120))
vntn60_bin79.8_120#sig

N0_bin79.8_159 <-cor.test(sqrt(virus.arc_N0$bin79.8),sqrt(virus.arc_N0$vOTU_159)) 
N0_bin79.8_159#sig
N60_bin79.8_159 <-cor.test(sqrt(virus.arc_N60$bin79.8),sqrt(virus.arc_N60$vOTU_159))
N60_bin79.8_159#sig

#correlation plot
library(grid)

#bin112.3 not significant on any votu after squre transform
# annot_N0_bin112.3_166 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_N0$bin112.3),sqrt(virus.arc_N0$vOTU_166)), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
# #annot_N60_bin112.3_166 = grobTree(textGrob(paste("r =", round(cor(virus.arc_N60$bin112.3,virus.arc_N60$vOTU_166), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
# annot_N0_bin112.3_166.p = grobTree(textGrob(paste("p"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold.italic")))
# annot_N0_bin112.3_166.p2 = grobTree(textGrob(paste("< 0.05"), x = 0.08, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))
# 
# bin112.3_166 <- ggplot(virus.arc.df.1, aes(sqrt(vOTU_166),sqrt(bin112.3),color=Nitrogen))+
#   geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_166)")+
#   geom_smooth(method="lm",se=TRUE)+
#   scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
#   ggtitle("Nitrososphaeraceae-UBA10452")+
#   annotation_custom(annot_N0_bin112.3_166)+
# #  annotation_custom(annot_N60_bin112.3_166)+
#   annotation_custom(annot_N0_bin112.3_166.p)+
#   annotation_custom(annot_N0_bin112.3_166.p2)+
#   theme_bw()

#bin16.10 and only and 128
# annot_ncntn0_bin16.10_159 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_ncntn0$bin16.10),sqrt(virus.arc_ncntn0$vOTU_159)), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
# annot_ncntn60_bin16.10_159 = grobTree(textGrob(paste("r =", round(cor(virus.arc_ncntn60$bin16.10,virus.arc_ncntn60$vOTU_159), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
# annot_vntn0_bin16.10_159 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_159)), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
# annot_vntn60_bin16.10_159 = grobTree(textGrob(paste("r =", round(cor(virus.arc_vntn60$bin16.10,virus.arc_vntn60$vOTU_159), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
# 
# 
# annot_N0_bin16.10_159.p = grobTree(textGrob(paste("p"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold.italic")))
# annot_N0_bin16.10_159.p2 = grobTree(textGrob(paste("< 0.05"), x = 0.08, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))
# 
# bin16.10_159 <- ggplot(virus.arc.df.1, aes(sqrt(vOTU_159),sqrt(bin16.10),color=Nitrogen))+
#   geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_159)")+
#   geom_smooth(method="lm",se=TRUE)+
#   scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
#   ggtitle("Nitrososphaeraceae-UBA10452")+
#   annotation_custom(annot_N0_bin16.10_159)+
#   annotation_custom(annot_N60_bin16.10_159)+
#   annotation_custom(annot_N0_bin16.10_159.p)+
#   annotation_custom(annot_N0_bin16.10_159.p2)+
#   theme_bw()


ncntn0_bin16.10_128 <-cor.test(sqrt(virus.arc_ncntn0$bin16.10),sqrt(virus.arc_ncntn0$vOTU_128)) 
ncntn0_bin16.10_128
vntn0_bin16.10_128 <-cor.test(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_128))
vntn0_bin16.10_128#sig
ncntn60_bin16.10_128 <-cor.test(sqrt(virus.arc_ncntn60$bin16.10),sqrt(virus.arc_ncntn60$vOTU_128)) 
ncntn60_bin16.10_128
vntn60_bin16.10_128 <-cor.test(sqrt(virus.arc_vntn60$bin16.10),sqrt(virus.arc_vntn60$vOTU_128))
vntn60_bin16.10_128#sig
annot_vntn0_bin16.10_128 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_vntn0$bin16.10),sqrt(virus.arc_vntn0$vOTU_128)), 3) ), x = 0.7, y = 0.92, hjust = 0, gp = gpar(col = "#29588a", fontsize = 15, fontface = "bold")))
annot_vntn60_bin16.10_128 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_vntn60$bin16.10),sqrt(virus.arc_vntn60$vOTU_128)), 3) ), x = 0.4, y = 0.8, hjust = 0, gp = gpar(col = "#7a103c", fontsize = 15, fontface = "bold")))
annot_bin16.10_128.p = grobTree(textGrob(paste("p"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold.italic")))
annot_bin16.10_128.p2 = grobTree(textGrob(paste("< 0.05"), x = 0.08, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin16.10_128 <- ggplot(virus.arc.df.1, aes(sqrt(vOTU_128),sqrt(bin16.10),color=Treat))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_128)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a","#29588a","#7a103c"),"Treatment")+
  ggtitle("Nitrososphaeraceae-UBA10452-bin16.10")+
  annotation_custom(annot_vntn0_bin16.10_128)+
annotation_custom(annot_vntn60_bin16.10_128)+
  annotation_custom(annot_bin16.10_128.p)+
  annotation_custom(annot_bin16.10_128.p2)+
  theme_bw()

#bin 35.3
annot_ncntn0_bin35.3_251 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_ncntn0$bin35.3),sqrt(virus.arc_ncntn0$vOTU_251)), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#1f8772", fontsize = 15, fontface = "bold")))
#annot_vntn60_bin35.3_251 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_vntn60$bin35.3),sqrt(virus.arc_vntn60$vOTU_251)), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#7a103c", fontsize = 15, fontface = "bold")))
annot_bin35.3_251.p = grobTree(textGrob(paste("p"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold.italic")))
annot_bin35.3_251.p2 = grobTree(textGrob(paste("< 0.05"), x = 0.08, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin35.3_251 <- ggplot(virus.arc.df.1, aes(sqrt(vOTU_251),sqrt(bin35.3),color=Treat))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_251)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a","#29588a","#7a103c"),"Treatment")+
  ggtitle("Nitrososphaeraceae-bin.35.3")+
  annotation_custom(annot_ncntn0_bin35.3_251)+
#  annotation_custom(annot_vntn60_bin35.3_251)+
  annotation_custom(annot_bin35.3_251.p)+
  annotation_custom(annot_bin35.3_251.p2)+
  theme_bw()

#bin79.8 sig with120 159
#annot_ncntn0_bin79.8_120 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_ncntn0$bin79.8),sqrt(virus.arc_ncntn0$vOTU_120)), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "#1f8772", fontsize = 15, fontface = "bold")))
annot_vntn60_bin79.8_120 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_vntn60$bin79.8),sqrt(virus.arc_vntn60$vOTU_120)), 3) ), x = 0.9, y = 0.7, hjust = 0, gp = gpar(col = "#7a103c", fontsize = 15, fontface = "bold")))
annot_ncntn60_bin79.8_120 = grobTree(textGrob(paste("r =", round(cor(sqrt(virus.arc_ncntn60$bin79.8),sqrt(virus.arc_ncntn60$vOTU_120)), 3) ), x = 0.9, y = 0.55, hjust = 0, gp = gpar(col = "#b07f2a", fontsize = 15, fontface = "bold")))

annot_bin79.8_120.p = grobTree(textGrob(paste("p"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold.italic")))
annot_bin79.8_120.p2 = grobTree(textGrob(paste("< 0.05"), x = 0.08, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin79.8_120 <- ggplot(virus.arc.df.1, aes(sqrt(vOTU_120),sqrt(bin79.8),color=Treat))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_120)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a","#29588a","#7a103c"),"Treatment")+
  ggtitle("Nitrososphaeraceae-UBA10452-bin79.8")+
#  annotation_custom(annot_ncntn0_bin79.8_120)+#120 in ncntn0 don't have enough data for correlation
  annotation_custom(annot_vntn60_bin79.8_120)+
  annotation_custom(annot_ncntn60_bin79.8_120)+
  annotation_custom(annot_bin79.8_120.p)+
  annotation_custom(annot_bin79.8_120.p2)+
  theme_bw()

annot_N0_bin79.8_159 = grobTree(textGrob(paste("r =", round(cor(virus.arc_N0$bin79.8,virus.arc_N0$vOTU_159), 3) ), x = 0.6, y = 0.3, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
annot_N60_bin79.8_159 = grobTree(textGrob(paste("r =", round(cor(virus.arc_N60$bin79.8,virus.arc_N60$vOTU_159), 3) ), x = 0.5, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
annot_bin79.8_159.p = grobTree(textGrob(paste("p"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold.italic")))
annot_bin79.8_159.p2 = grobTree(textGrob(paste("< 0.05"), x = 0.08, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin79.8_159 <- ggplot(virus.arc.df.1, aes(sqrt(vOTU_159),sqrt(bin79.8),color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_159)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Nitrososphaeraceae-UBA10452-bin79.8")+
  annotation_custom(annot_N0_bin79.8_159)+
  annotation_custom(annot_N60_bin79.8_159)+
  annotation_custom(annot_bin79.8_159.p)+
  annotation_custom(annot_bin79.8_159.p2)+
  theme_bw()

#######correlation figures#########
library(ggpubr)
tiff('figures/arc.virus.cor.tiff', units="in", width=12, height=12, res=180)
ggarrange(bin16.10_128,bin35.3_251,bin79.8_120,bin79.8_159,
          labels = c("(A)", "(B)","(C)","(D)"),ncol = 2, nrow = 2,common.legend = F,legend = "top",
          hjust = -0.4) 
dev.off()
#stacked bar plot
library(dplyr)
cols_phylum <- c("#994e03", "#9B110E", "#688787" ,"#d49783", "#550307", "#446455", "#FDD262", "#D3DDDC", "#C7B19C","#899DA4", "#FAEFD1", "#DC863B","#798E87", "#C27D38","#CCC591","#29211F")
arc.mag.hp.tab.avg <- unique(arc.mag.hp.tab[,-c(1,3,4,5)])
arc.mag.hp.tab.avg.host <- arc.mag.hp.tab.avg[which(arc.mag.hp.tab.avg$bins=="bin112.3"|arc.mag.hp.tab.avg$bins=="bin16.10"|
                                                      arc.mag.hp.tab.avg$bins=="bin35.3"|arc.mag.hp.tab.avg$bins=="bin63.4"|
                                                      arc.mag.hp.tab.avg$bins=="bin79.8"),]

tiff('figures/stacked_arc_host.tiff', units="in", width=10, height=8, res=300)
#arc.hp <- 
  ggplot(arc.mag.hp.tab.avg.host, aes(x =Treat, y = mean.value,fill=bins)) + 
  geom_bar(stat='identity',colour="grey",size=0.1, width = 0.6,alpha=1)+
    ggtitle("Archaeal MAGs")+
  #  geom_hline(yintercept = 3.5, size=1.5)+
  #geom_vline(xintercept = c(3.5,6.5,9.5), size=1)+
  #geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,11.5), size=0.2)+
  #  geom_hline(yintercept = c(1.5,2.5,4.5,5.5,6.5), size=0.2)+
  xlab(label = "Treat")+
  ylab(label = "Relative abundance")+
  #  scale_fill_gradient(name="Relative abundance",low="#3794bf", high="#a34c1d")+#low = "violetred", high = "aquamarine"
  scale_fill_manual(name=NULL,values=cols_phylum)+#low = "violetred", high = "aquamarine"
  
  #  facet_grid(~ Depth)+
  theme_bw(base_size=15,base_family = "",) 
dev.off()
