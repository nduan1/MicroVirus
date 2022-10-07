
#The bin.9 10.50602416 value means that on average, the contigs in bin.9 have a read coverage of ~10. So in your example, the total estimated number of reads mapping to your bin is 1073567*10.50602416/100 = 113,755 reads (I used 100bp for read length because I assumed it was a paired-end library). So of your total library size of 23,887,914 reads, at least 0.5% came from this organism. However, if you goal is to know the bin's "abundance", using the original read depth as a representation of the abundance (10.50) is probably more accurate. Because you may have only a part of the whole genome, read depth is a better metric. Just dont forget to standardize them if you start comparing samples. Hope that makes sense.
#Ok, so here's the issue. You calculated that 0.5% of the reads map back to bin.9. However you dont have the entire genome - just 78%. That means that in theory even more reads map to the genome, you just dont have it. On the other hand, the original bin read coverage (i.e. 10.50602416) should still be roughly the same even if you recovered 50% of the genome or 100% of the genome. So why not just use that as your abundance - this is the estimate read depth of the bin, and corresponds to its abundance in the sample.
#To put it simply, if your goal is to relate the relative abundance of the bins to each other (which one is more abundant in each sample), then use the 10.5 value (if you have many samples, you will need to standardize this). If you need to determine the relative abundance of the bin in the whole community, use the 0.5% of total reads value (roughly what % of the reads came from the organism).

#I did this because I made a coverage table per contig, and a bin (or MAG) contains multiple contigs, this is why I did it that way.
#So indeed, I started with coverage numbers x (contig length/total length of all contigs in bin).
#And then those coverage numbers x (average sequence depth for all samples / average sequence depth per sample), so the other way around as what you are doing I think?
  
#  Then after those two steps you filter out the abundance < 0.25 coverage. 




library(tidyverse)
bin107.3 <- read.csv("dataset/bin.107.3.contigs.tsv", header = F)
bin112.7 <- read.csv("dataset/bin.112.7.contigs.tsv", header = F)
bin16.1<- read.csv("dataset/bin.16.1.contigs.tsv", header = F)
bin16.2 <- read.csv("dataset/bin.16.2.contigs.tsv", header = F)
bin35.6 <- read.csv("dataset/bin.35.6.contigs.tsv", header = F)
bin35.9 <- read.csv("dataset/bin.35.9.contigs.tsv", header = F)
bin37.1 <- read.csv("dataset/bin.37.1.contigs.tsv", header = F)
bin37.3 <- read.csv("dataset/bin.37.3.contigs.tsv", header = F)
bin37.7 <- read.csv("dataset/bin.37.7.contigs.tsv", header = F)
bin51.3 <- read.csv("dataset/bin.51.3.contigs.tsv", header = F)
bin62.1 <- read.csv("dataset/bin.62.1.contigs.tsv", header = F)
bin62.2 <- read.csv("dataset/bin.62.2.contigs.tsv", header = F)
bin9.3 <- read.csv("dataset/bin.9.3.contigs.tsv", header = F)
head(bin107.3)
bins_names <- c("bin107.3","bin112.7","bin16.1","bin16.2","bin35.6",
                "bin35.9","bin37.1","bin37.3","bin37.7","bin51.3",
                "bin62.1","bin62.2","bin9.3")
# add.bin <- function(file){
#   bins.df <-mutate(file, bin=c(rep(bins_names[i],dim(file)[1])))
#   bins.df <-bins.df[,c(2,1)]
#   }
#   return(bins.df)
# }

# for (i in length(bins_names)){
#   eval(as.symbol(paste((bins_names)[i],sep = ".","df"))) <- add.bin(as.data.frame(eval(as.symbol(paste((bins_names)[1])))))
#   return()
# }
bins_names <- c("bin107.3","bin112.7","bin16.1","bin16.2","bin35.6",
                "bin35.9","bin37.1","bin37.3","bin37.7","bin51.3",
                "bin62.1","bin62.2","bin9.3")
bin107.3.df <-mutate(bin107.3,bin=c(rep("bin107.3",dim(bin107.3)[1])))
bin107.3.df <- bin107.3.df[,c(2,1)]

bin112.7.df <-mutate(bin112.7,bin=c(rep("bin112.7",dim(bin112.7)[1])))
bin112.7.df <- bin112.7.df[,c(2,1)]

bin16.1.df <-mutate(bin16.1,bin=c(rep("bin16.1",dim(bin16.1)[1])))
bin16.1.df <- bin16.1.df[,c(2,1)]

bin16.2.df <-mutate(bin16.2,bin=c(rep("bin16.2",dim(bin16.2)[1])))
bin16.2.df <- bin16.2.df[,c(2,1)]

bin35.6.df <-mutate(bin35.6,bin=c(rep("bin35.6",dim(bin35.6)[1])))
bin35.6.df <- bin35.6.df[,c(2,1)]

bin35.9.df <-mutate(bin35.9,bin=c(rep("bin35.9",dim(bin35.9)[1])))
bin35.9.df <- bin35.9.df[,c(2,1)]

bin37.1.df <-mutate(bin37.1,bin=c(rep("bin37.1",dim(bin37.1)[1])))
bin37.1.df <- bin37.1.df[,c(2,1)]

bin37.3.df <-mutate(bin37.3,bin=c(rep("bin37.3",dim(bin37.3)[1])))
bin37.3.df <- bin37.3.df[,c(2,1)]

bin37.7.df <-mutate(bin37.7,bin=c(rep("bin37.7",dim(bin37.7)[1])))
bin37.7.df <- bin37.7.df[,c(2,1)]

bin51.3.df <-mutate(bin51.3,bin=c(rep("bin51.3",dim(bin51.3)[1])))
bin51.3.df <- bin51.3.df[,c(2,1)]

bin62.1.df <-mutate(bin62.1,bin=c(rep("bin62.1",dim(bin62.1)[1])))
bin62.1.df <- bin62.1.df[,c(2,1)]

bin62.2.df <-mutate(bin62.2,bin=c(rep("bin62.2",dim(bin62.2)[1])))
bin62.2.df <- bin62.2.df[,c(2,1)]

bin9.3.df <-mutate(bin9.3,bin=c(rep("bin9.3",dim(bin9.3)[1])))
bin9.3.df <- bin9.3.df[,c(2,1)]
head(bin9.3.df)
all.mag.contigs <- rbind(bin107.3.df,bin112.7.df,bin16.1.df,bin16.2.df,bin35.6.df,
                         bin35.9.df,bin37.1.df,bin37.3.df,bin37.7.df,bin51.3.df,
                         bin62.1.df,bin62.2.df,bin9.3.df)

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
head(bin.cov.tab.clean.df )
########now we have clean bins with contigs length and coverage table in every sample
# calculate the normalization factor for the bins (= contig_len / total_bin_len)
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
view(clean.bin.otu.filtered)
######now you have a clean bacteria_mag_otu table
######now we have normalized average number of reads per base-pair per contigs across samples
#next we plan to generate the phyloseq
library(phyloseq)
library(ggplot2)
library(tidyverse)
meta.data <- read.csv("dataset/meta_data.csv",header = T, row.names = 1)
meta.data.F <- read.csv("dataset/meta_data.csv",header = T)
head(meta.data)
bac.mag.sample.data <- sample_data(meta.data)
head(clean.bin.otu.filtered)
row.names(clean.bin.otu.filtered) <- clean.bin.otu.filtered[,1]
head(clean.bin.otu.filtered)
bac.mag.otu <- otu_table(clean.bin.otu.filtered[,-1], taxa_are_rows = TRUE)
bac.mag.phyloseq <-phyloseq(bac.mag.otu, bac.mag.sample.data)
bac.mag.phyloseq
r_bac.mag.phyloseq<-transform_sample_counts(bac.mag.phyloseq,function(x) x/sum(x))
head(otu_table(r_bac.mag.phyloseq))
write.csv(t(otu_table(r_bac.mag.phyloseq)),"virus_host_cor/bac.mag.csv")
bac.mag.hp.df <- mutate(as.data.frame(otu_table(r_bac.mag.phyloseq)),bins=rownames(clean.bin.otu.filtered)) %>% pivot_longer(cols = colnames(otu_table(r_bac.mag.phyloseq)))
##reorder the messy rows
write.csv(bac.mag.hp.df,"dataset/bac.mag.hp.df.csv")
write.csv(t(otu_table(bac.mag.phyloseq)),"dataset/bac.mag.phylotu.csv")#I reorder the files
bac.mag.hp.df.1 <- read.csv("dataset/bac.mag.hp.df.1.csv")
head(bac.mag.hp.df.1)
tail(bac.mag.hp.df.1)
# treat.order <- c(rep(c("NCNTN0_1","NCNTN0_2","NCNTN0_4","NCNTN60_1","NCNTN60_2","NCNTN60_4",
#                  "VNTN0_1","VNTN0_2","VNTN0_4","VNTN60_1","VNTN60_2","VNTN60_4"),13))
# bac.mag.hp.df.order <- bac.mag.hp.df %>% slice(match(treat.order,name))
# head(bac.mag.hp.df.order)
# bac.mag.hp.df
bac.mag.hp <- inner_join(bac.mag.hp.df.1,meta.data.F,by=c("name"="X"))
head(bac.mag.hp)
bac.mag.hp.mean <- bac.mag.hp %>% group_by(bins,Treat) %>% summarise_at(vars(value),list(mean.value=mean))
bac.mag.hp.tab <- left_join(bac.mag.hp,bac.mag.hp.mean,by=c("bins"="bins","Treat"="Treat"))
head(bac.mag.hp.tab)

bins.order <- c("bin35.9","bin107.3","bin51.3","bin9.3","bin62.1","bin35.6",
                "bin37.1","bin16.1","bin62.2","bin16.2","bin37.3","bin37.7",
                "bin112.7")#make sure the order is the same to the bac.mag.tree
bac.mag.hp.tab$bins <- factor(bac.mag.hp.tab$bins, levels = as.factor(bins.order))
tiff('figures/bac.mag.hp.tiff', units="in", width=5, height=10, res=300)
ggplot(bac.mag.hp.tab, aes(x =Treat, y = bins,fill=mean.value)) + 
  geom_tile(aes(fill = mean.value))+
  geom_text(aes(label = round(mean.value, 2)))+
  #  geom_hline(yintercept = 3.5, size=1.5)+
  #geom_vline(xintercept = c(3.5,6.5,9.5), size=1)+
  #geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,11.5), size=0.2)+
  #  geom_hline(yintercept = c(1.5,2.5,4.5,5.5,6.5), size=0.2)+
  xlab(label = "Treat")+
  ylab(label = "MAGs")+
#  scale_fill_gradient(name="Relative abundance",low="#3794bf", high="#a34c1d")+#low = "violetred", high = "aquamarine"
  scale_fill_gradient(name="Relative abundance",low="#edf0ed", high="#4c784c")+#low = "violetred", high = "aquamarine"
  
  #  facet_grid(~ Depth)+
  theme_bw()
dev.off()
head(bac.mag.hp.tab)

#######correlation between bacterial host and virus within every treatment
#####organize file first
bac.mag.ra <- read.csv("virus_host_cor/bac.mag.phylotu.csv")
head(bac.mag.ra);dim(bac.mag.ra)
virus.otu.ra <- read.csv("virus_host_cor/virus.phyloseq.cov.csv")
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
bac.112.7.mod <- lmer(bin112.7^(1/4) ~  Cover +(1|block), 
                      data = subset(virus.bac.df.1,Nitrogen=="N60"))
Anova(bac.112.7.mod)
summary(bac.112.7.mod)
coef(bac.112.7.mod)

bac.37.7.mod <- lmer(bin37.7^(1/3) ~ Nitrogen+(1|block), 
                     data = subset(virus.bac.df.1,Cover=="V"))
Anova(bac.37.7.mod)
summary(bac.37.7.mod)
coef(bac.37.7.mod)

bac.37.3.mod <- lmer(bin37.3^(1/2) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.37.3.mod)
summary(bac.37.3.mod)
coef(bac.37.3.mod)

bac.16.2.mod <- lmer(bin16.2 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.16.2.mod)
summary(bac.16.2.mod)
coef(bac.16.2.mod)
###
bac.62.2.mod <- lmer(bin62.2 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.62.2.mod)
summary(bac.62.2.mod)
coef(bac.62.2.mod)

bac.16.1.mod <- lmer(bin16.1 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.16.1.mod)
summary(bac.16.1.mod)#significant on nitrogen
coef(bac.16.1.mod)

bac.37.1.mod <- lmer(sqrt(bin37.1) ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.37.1.mod)
summary(bac.37.1.mod)
coef(bac.37.1.mod)

bac.35.6.mod <- lmer(bin35.6 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.35.6.mod)
summary(bac.35.6.mod)
coef(bac.35.6.mod)

bac.62.1.mod <- lmer(bin62.1 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.62.1.mod)
summary(bac.62.1.mod)
coef(bac.62.1.mod)

bac.9.3.mod <- lmer(bin9.3 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.9.3.mod)
summary(bac.9.3.mod)###nitrogen
coef(bac.9.3.mod)

bac.51.3.mod <- lmer(bin51.3 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.51.3.mod)
summary(bac.51.3.mod)
coef(bac.51.3.mod)

bac.107.3.mod <- lmer(bin107.3 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.107.3.mod)
summary(bac.107.3.mod)
coef(bac.107.3.mod)

bac.35.9.mod <- lmer(bin35.9 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(bac.35.9.mod)
summary(bac.35.9.mod)
coef(bac.35.9.mod)

###virus
virus.181.mod <- lmer(vOTU_181 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.181.mod)
summary(virus.181.mod)
coef(virus.181.mod)

virus.14.mod <- lmer(vOTU_14 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.14.mod)
summary(virus.14.mod)
coef(virus.14.mod)

virus.152.mod <- lmer(vOTU_152 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.152.mod)
summary(virus.152.mod)
coef(virus.152.mod)

virus.3.mod <- lmer(vOTU_3 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.3.mod)
summary(virus.3.mod)
coef(virus.3.mod)

virus.19.mod <- lmer(vOTU_19 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.19.mod)
summary(virus.19.mod)#nitrogen
coef(virus.19.mod)

virus.48.mod <- lmer(vOTU_48 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.48.mod)
summary(virus.48.mod)
coef(virus.48.mod)


virus.63.mod <- lmer(vOTU_63 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.63.mod)
summary(virus.63.mod)
coef(virus.63.mod)


virus.93.mod <- lmer(vOTU_93 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.93.mod)
summary(virus.93.mod)
coef(virus.93.mod)


virus.153.mod <- lmer(vOTU_153 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.153.mod)
summary(virus.153.mod)
coef(virus.153.mod)

virus.203.mod <- lmer(vOTU_203 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.203.mod)
summary(virus.203.mod)#nitrogen
coef(virus.203.mod)

virus.81.mod <- lmer(vOTU_81 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.81.mod)
summary(virus.81.mod)#nitrogen
coef(virus.81.mod)

virus.115.mod <- lmer(vOTU_115 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.115.mod)
summary(virus.115.mod)#nitrogen
coef(virus.115.mod)

virus.163.mod <- lmer(vOTU_163 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.163.mod)
summary(virus.163.mod)#nitrogen
coef(virus.163.mod)

virus.252.mod <- lmer(vOTU_252 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.252.mod)
anova(virus.252.mod)
summary(virus.252.mod)#nitrogen
coef(virus.252.mod)

virus.161.mod <- lmer(vOTU_161 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
a <- lm(vOTU_161 ~ Cover* Nitrogen, data = virus.bac.df.1)
Anova(a)
stepAIC(a)
Anova(virus.161.mod)
summary(virus.161.mod)
coef(virus.161.mod)

virus.106.mod <- lmer(vOTU_106 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.106.mod)
summary(virus.106.mod)
coef(virus.106.mod)

virus.148.mod <- lmer(vOTU_148 ~ Cover* Nitrogen +(1|block), data = virus.bac.df.1)
Anova(virus.148.mod)
summary(virus.148.mod)
coef(virus.148.mod)




normalTest(sqrt(virus.bac.df.1$vOTU_181),"sw")#not normal
normalTest(sqrt(virus.bac.df.1$bin112.7),"sw")

normalTest((virus.bac.df.1$vOTU_14),"sw")#one value
normalTest(sqrt(virus.bac.df.1$vOTU_19),"da")#not normal
normalTest(sqrt(virus.bac.df.1$vOTU_161),"sw")#not normal
normalTest((virus.bac.df.1$bin35.9),"sw")#not normal

normalTest((virus.bac.df.1$vOTU_153),"sw")
normalTest((virus.bac.df.1$vOTU_203),"sw")#good
normalTest((virus.bac.df.1$vOTU_81),"sw")#good
normalTest((virus.bac.df.1$vOTU_115),"sw")#good
normalTest((virus.bac.df.1$vOTU_163),"sw")#good
normalTest((virus.bac.df.1$vOTU_252),"sw")#not bad

#bac
normalTest((virus.bac.df.1$bin16.1),"sw")#not bad
normalTest((virus.bac.df.1$bin9.3),"sw")#good

#the rest of bacteria and viruses don't have significant influence by two treatment
#therefore, we could put all things together to do a test
cor.test(virus.bac.df.1$vOTU_181,virus.bac.df.1$bin112.7, method = "pearson", conf.level = 0.95)
t.test(virus.bac.df.1$vOTU_14,virus.bac.df.1$bin112.7)#sig but votu only have two values, so don't consider it
cor.test(virus.bac.df.1$vOTU_181,virus.bac.df.1$bin37.7, method = "pearson", conf.level = 0.95)
cor.test(virus.bac.df.1$vOTU_181,virus.bac.df.1$bin37.7, method = "pearson", conf.level = 0.95)
#14 and 112.7 14 and 37.7 14 and 37.3 , votu 3 votu48 will not be considered due to the lack of number of information
normalTest(sqrt(virus.bac.df.1$bin37.1),"sw")
normalTest(sqrt(virus.bac.df.1$vOTU_19),"sw")
cor.test(sqrt(virus.bac.df.1$vOTU_19),sqrt(virus.bac.df.1$bin37.1), method = "pearson", conf.level = 0.95)
#bin35.6 lack info, therefore 63 and 93 won't considered
#bin35.9 lack info, therefore 161 106 148 wont consider


# Calculating Pearson's product-moment correlation
cor.test(virus.bac.df.1$vOTU_181,virus.bac.df.1$bin112.7, method = "pearson", conf.level = 0.95)
a <- cor.test(virus.bac.df.1$vOTU_181,virus.bac.df.1$bin112.7, method = "pearson", conf.level = 0.95)
a
########I am going to make bin9.3 figures with multiple different vOTUs with good normality
#subset data first and try to figure out if there are correlation bewtween virus and host on N0 and N60 seperately
#if you want to conduct a pearson correlation test you must fullfill all!! the preassumption test
#1.random sample 2. variables are continous data 3. data contains paired samples 4. Independence of observations 5. variables are approxmately normally distributed 6. a linear association exists 7. absence of outliers
virus.bac_N0 <- subset(virus.bac.df.1, Nitrogen=="N0");dim(virus.bac_N0)
virus.bac_N60 <- subset(virus.bac.df.1, Nitrogen=="N60");dim(virus.bac_N60)
N0_bin9.3_153 <-cor.test(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_153) 
N0_bin9.3_153#sig
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

annot_N0_bin9.3_153 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_153), 3) ), x = 0.7, y = 0.15, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
#annot_N60_bin9.3_153 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_153), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))
annot_N0_bin9.3_153.p = grobTree(textGrob(paste("P < 0.05"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin9.3_153 <- ggplot(virus.bac.df.1, aes(vOTU_153,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_153)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_153)+
#  annotation_custom(annot_N60_bin9.3_153)+
  annotation_custom(annot_N0_bin9.3_153.p)+
  theme_bw()


annot_N0_bin9.3_203 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_203), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
annot_N60_bin9.3_203 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_203), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
annot_N0_bin9.3_203.p = grobTree(textGrob(paste("P < 0.01"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

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

annot_N0_bin9.3_81 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_81), 3) ), x = 0.1, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
annot_N60_bin9.3_81 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_81), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
annot_N0_bin9.3_81.p = grobTree(textGrob(paste("P < 0.01"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin9.3_81 <- ggplot(virus.bac.df.1, aes(vOTU_81,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_81)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_81)+
  annotation_custom(annot_N60_bin9.3_81)+
  annotation_custom(annot_N0_bin9.3_81.p)+
  theme_bw()

annot_N0_bin9.3_115 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_115), 3) ), x = 0.05, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
annot_N60_bin9.3_115 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_115), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
annot_N0_bin9.3_115.p = grobTree(textGrob(paste("P < 0.05"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin9.3_115 <- ggplot(virus.bac.df.1, aes(vOTU_115,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_115)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_115)+
  annotation_custom(annot_N60_bin9.3_115)+
  annotation_custom(annot_N0_bin9.3_115.p)+
  theme_bw()

annot_N0_bin9.3_163 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_163), 3) ), x = 0.1, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
annot_N60_bin9.3_163 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_163), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col ="#4f3152", fontsize = 15, fontface = "bold")))
annot_N0_bin9.3_163.p = grobTree(textGrob(paste("P < 0.01"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin9.3_163 <- ggplot(virus.bac.df.1, aes(vOTU_163,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_163)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_163)+
  annotation_custom(annot_N60_bin9.3_163)+
  annotation_custom(annot_N0_bin9.3_163.p)+
  theme_bw()

annot_N0_bin9.3_252 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N0$bin9.3,virus.bac_N0$vOTU_252), 3) ), x = 0.3, y = 0.55, hjust = 0, gp = gpar(col = "#145448", fontsize = 15, fontface = "bold")))
annot_N60_bin9.3_252 = grobTree(textGrob(paste("r =", round(cor(virus.bac_N60$bin9.3,virus.bac_N60$vOTU_252), 3) ), x = 0.7, y = 0.9, hjust = 0, gp = gpar(col = "#4f3152", fontsize = 15, fontface = "bold")))
annot_N0_bin9.3_252.p = grobTree(textGrob(paste("P < 0.05"), x = 0.05, y = 0.95, hjust = 0, gp = gpar(col = "black", fontsize = 15, fontface = "bold")))

bin9.3_252 <- ggplot(virus.bac.df.1, aes(vOTU_252,bin9.3,color=Nitrogen))+
  geom_point(size=2)+labs(x = "Host abundance", y = "Viral abundance (vOTU_252)")+
  geom_smooth(method="lm",se=TRUE)+
  scale_color_manual(values = c("#1f8772","#b07f2a"),"Fertilization")+
  ggtitle("Verrucomicrobiota-AV55")+
  annotation_custom(annot_N0_bin9.3_252)+
  annotation_custom(annot_N60_bin9.3_252)+
  annotation_custom(annot_N0_bin9.3_252.p)+
  theme_bw()

#########heatmap figures for both bac and viruses###########
tiff('figures/bac.mag.hp.tiff', units="in", width=5, height=10, res=300)
bac.hp <- ggplot(bac.mag.hp.tab, aes(x =Treat, y = bins,fill=mean.value)) + 
  geom_tile(aes(fill = mean.value))+
  geom_text(aes(label = round(mean.value, 3)))+
  #  geom_hline(yintercept = 3.5, size=1.5)+
  #geom_vline(xintercept = c(3.5,6.5,9.5), size=1)+
  #geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,15.5), size=0.2)+
  #  geom_hline(yintercept = c(1.5,2.5,4.5,5.5,6.5), size=0.2)+
  xlab(label = "Treat")+
  ylab(label = "MAGs")+
  #  scale_fill_gradient(name="Relative abundance",low="#3794bf", high="#a34c1d")+#low = "violetred", high = "aquamarine"
  scale_fill_gradient(name="MAGs",low="#edf0ed", high="#4c784c")+#low = "violetred", high = "aquamarine"
  #scale_fill_gradient(name="MAGs",low="#fafbfc", high="#083778")+
  #  facet_grid(~ Depth)+
  theme_bw(base_size=15)
dev.off()

virus.otu.relative.ab <- read.csv("virus_host_cor/t.votu.tab.csv")

head(virus.otu.relative.ab);dim(virus.otu.relative.ab)
virus.otu.relative.sub <- virus.otu.relative.ab[,c(1,182,15,153,4,20,49,64,94,154,204,82,156,164,253,162,107,149)]
virus.otu.relative.sub
virus.otu.relative <- mutate(virus.otu.relative.sub, Treat=c(rep("NCNTN0",3),rep("NCNTN60",3),rep("VNTN0",3),rep("VNTN60",3)), 
                       Nitrogen=c(rep("N0",3),rep("N60",3),rep("N0",3),rep("N60",3)),
                       Cover=c(rep("NC",6),rep("V",6)),
                       block=c(rep(c("b1","b2","b4"),4)))
colnames(virus.otu.relative)
virus.otu.relative.1 <- virus.otu.relative[,c(1,19:22,2:18)]
colnames(virus.otu.relative.1)[1] <- "Sample" 
head(virus.otu.relative.1)
virus.otu.relative.df <- pivot_longer(virus.otu.relative.1, cols = c("vOTU_181","vOTU_14","vOTU_152","vOTU_3","vOTU_19",
                                                                     "vOTU_48","vOTU_63","vOTU_93","vOTU_153","vOTU_203",
                                                                     "vOTU_81","vOTU_155","vOTU_163","vOTU_252","vOTU_161",
                                                                     "vOTU_106","vOTU_148"))
head(virus.otu.relative.df);dim(virus.otu.relative.df)

virus.hp.mean <- virus.otu.relative.df %>% group_by(name,Treat) %>% 
  summarise_at(vars(value),list(mean.value=mean))
virus.hp.df <- left_join(virus.otu.relative.df,virus.hp.mean,by=c("name"="name","Treat"="Treat"))
head(virus.hp.df)

virus.order <- c("vOTU_148","vOTU_106","vOTU_161","vOTU_252","vOTU_163","vOTU_155",
                 "vOTU_81","vOTU_203","vOTU_153","vOTU_93","vOTU_63","vOTU_48","vOTU_19",
                 "vOTU_3","vOTU_152","vOTU_14","vOTU_181")#make sure the order is the same to the bac.mag.tree

virus.hp.df$name <- factor(virus.hp.df$name, levels = as.factor(virus.order))

virus.hp <- ggplot(virus.hp.df, aes(x =Treat, y = name,fill=mean.value)) + 
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
tiff('figures/bac.virus.hp.tiff', units="in", width=9, height=9, res=300)
ggarrange(bac.hp,virus.hp,labels = c("", ""),ncol = 2, nrow = 1,common.legend = F,legend = "bottom") 
dev.off()

#######correlation figures#########
tiff('figures/bac.virus.cor.tiff', units="in", width=14.4, height=9, res=190)
ggarrange(bin9.3_203,bin9.3_81,bin9.3_115,bin9.3_163,bin9.3_252,
          labels = c("(A)", "(B)","(C)","(D)","(E)"),ncol = 3, nrow = 2,common.legend = T,legend = "top") 
dev.off()

#stacked bar plot
library(dplyr)
cols_phylum <- c("#994e03", "#9B110E", "#688787" ,"#d49783", "#550307", "#446455", "#FDD262", "#D3DDDC", "#C7B19C","#899DA4", "#FAEFD1", "#DC863B","#798E87", "#C27D38","#CCC591","#29211F")
bac.mag.hp.tab.avg <- unique(bac.mag.hp.tab[,-c(2,3,4)])
bac.mag.hp.tab.avg.host <- bac.mag.hp.tab.avg[which(bac.mag.hp.tab.avg$bins=="bin112.7"|bac.mag.hp.tab.avg$bins=="bin16.2"|
                                                      bac.mag.hp.tab.avg$bins=="bin35.9"|bac.mag.hp.tab.avg$bins=="bin37.7"|
                                                      bac.mag.hp.tab.avg$bins=="bin9.3"|bac.mag.hp.tab.avg$bins=="bin35.6"|
                                                      bac.mag.hp.tab.avg$bins=="bin37.1"|bac.mag.hp.tab.avg$bins=="bin37.3"),]
tiff('figures/stacked_bac_host.tiff', units="in", width=10, height=8, res=300)
#arc.hp <- 
ggplot(bac.mag.hp.tab.avg.host , aes(x =Treat, y = mean.value,fill=bins)) + 
  geom_bar(stat='identity',colour="grey",size=0.1, width = 0.6,alpha=1)+
  ggtitle("Bacterial MAGs")+
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
