library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
#################normalize your abundance table################
votu_metagenome_rev <- read.csv("datataset/filtered_otu_table.csv", header = T, row.names = 1)
dim(votu_metagenome_rev)
head(votu_metagenome_rev)
head(metag.quality)
read.depth.12 <- unique(metag.quality$Sequencing_depth)
sample_list <- unique(metag.quality$Samples)
avg.read.depth <- sum(read.depth.12)/12##average number of reads across 12 samples
seq.depth <- data.frame(sample=sample_list, depth=read.depth.12) %>% 
  mutate(avg.depth=rep(avg.read.depth,12),
         mutiplier=avg.read.depth/depth)
write.csv(seq.depth,"seq.depth.csv")##got a normalization coeffcient
head(seq.depth)
head(votu_metagenome_rev)
seq.depth$mutiplier[2]
votu.nor <- function(file){
  nor.matrix <- matrix(nrow=260,ncol = 0)#if you don't add 0 here, you will get an extra useless columns
  for(a in 1:12){
    newcol <- round(votu_metagenome_rev[,a]*(seq.depth$mutiplier[a]),3)
    nor.matrix <- cbind(nor.matrix,newcol)
  }
  return(nor.matrix)
}

votu.nor.tab <- votu.nor(votu_metagenome_rev)
colnames(votu.nor.tab) <- colnames(votu_metagenome_rev)
rownames(votu.nor.tab) <- rownames(votu_metagenome_rev)
head(votu.nor.tab)
dim(votu.nor.tab)
######now we have normalized average number of reads per base-pair per contigs across samples
#next we plan to generate the phyloseq
meta.data <- read.csv("datataset/meta_data.csv",header = T, row.names = 1)
head(meta.data)
vsample.data <- sample_data(meta.data)
votu.table <- otu_table(votu.nor.tab, taxa_are_rows = TRUE)
virus.phyloseq <-phyloseq(votu.table, vsample.data)
virus.phyloseq
#############create venn figure among four treatment
install.packages("VennDiagram")
library(VennDiagram)
votu_venn <- read.csv("datataset/votu_venn.csv",header = T, row.names = 1)
head(votu_venn)
#I process and use countif function in excel and find out how many shared otus across four treatment
a1=114
a2=112
a3=107
a4=135
vn12=53
vn13=61
vn14=51
vn23=63
vn24=75
vn34=62
vn123=51
vn134=50
vn234=56
vn124=48
vn1234=48
virus.venn.plot <- draw.quad.venn(area1 = a1,area2 = a2, area3 = a3,area4 = a4,
                                  n12 = vn12,n13 = vn13,n14 = vn14,n23 = vn23,n24 = vn24,n34 = vn34,
                                  n123 = vn123,n134 = vn134,n124 = vn124,n234 = vn234,n1234 = vn1234,
                                  category = c("NCNTN0","NCNTN60(2)","VNTN0","VNTN60(4)"),
                                  fill = c("#4D897C", "#C6B7EC", "#ba892f","#c76d70"),lty = "blank",
                                  cex = 2,col = c("#242323","#242323","#242323","#242323"),
                                  cat.cex = 1.5,alpha = 0.65,
                                  euler.d = TRUE,
                                  scaled = TRUE)
ggsave('/Users/Ning/Desktop/R/New project/metagenome_2020/figures/virus_venn_4way_treatment.pdf',dpi = 300, virus.venn.plot) 
dev.off()

# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")
# library("ggVennDiagram")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")
new_venn_data <- read.csv("datataset/new_venn_data.csv")
venn.data <- list(
  A = new_venn_data[,1], 
  B = new_venn_data[,2], 
  C = new_venn_data[,3],
  D = new_venn_data[,4]
)
names(venn.data) <- c("NCNTN0","NCNTN60","VNTN0","VNTN60")
venn <- ggvenn(venn.data,fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       stroke_size = 0.3, set_name_size = 4.5,fill_alpha = 0.5,text_size = 3.8)

# ggVennDiagram(venn.data, label_alpha = 0,color = NULL,
#   category.names = c("NCNTN0","NCNTN60","VNTN0","VNTN60")
# ) +scale_color_manual(values =c("#4D897C", "#C6B7EC", "#ba892f","#c76d70"))



############box plot for number of votus in each sample#############
virus.phyloseq
sample_sums(otu_table(virus.phyloseq))
head(meta.data)
tiff('box_votus.diversity.tiff', units="in", width=5, height=5, res=300)
votu_num <- ggplot(meta.data,aes(x=Treat,y=vOTU,fill=Treat))+
  geom_boxplot(size=0.5,alpha=0.8)+
  geom_jitter(shape=16, position=position_jitter(0),aes(fill=Treat))+
  scale_fill_manual(values = c("#4D897C", "#C6B7EC", "#ba892f","#c76d70"),name="Treatment")+
  labs(title="",x="Treatment",y="Number of vOTUs")+
  theme_bw() 
dev.off()

library(lmerTest)
library(car)
library(fBasics)
votu.mod <- lmer(vOTU ~ Cover* Nitrogen +(1|block), data = meta.data)
#votu.mod.lm <- lm(vOTU ~ Cover* Nitrogen, data = meta.data)
#normalTest(votu.mod.lm$residuals,"sw")
Anova(votu.mod)
normalTest(meta.data$vOTU,"sw")#votu normality is good
votu.mod
############heatmap barplot###########
library(tidyverse)
##heatmap for samples
r_virus.phyloseq<-transform_sample_counts(virus.phyloseq,function(x) x/sum(x))
head(otu_table(r_virus.phyloseq))
t.votu.tab <-t(otu_table(r_virus.phyloseq))

write.csv(t.votu.tab,"t.votu.tab.csv")
write.csv(t(otu_table(virus.phyloseq)),"datataset/virus.phyloseq.cov.csv")
votu.hp.df <- pivot_longer(as.data.frame(t.votu.tab),cols = colnames(t.votu.tab)) %>% 
  mutate(sample=c(rep("NCNTN0_1",260),rep("NCNTN0_2",260),rep("NCNTN0_4",260),
                  rep("NCNTN60_1",260),rep("NCNTN60_2",260),rep("NCNTN60_4",260),
                  rep("VNTN0_1",260),rep("VNTN0_2",260),rep("VNTN0_4",260),
                  rep("VNTN60_1",260),rep("VNTN60_2",260),rep("VNTN60_4",260)))
head(votu.hp.df)
tiff('votu.hp.tiff', units="in", width=13, height=30, res=300)
ggplot(votu.hp.df, aes(x =sample, y = name,fill=value)) + 
  geom_tile()+
#  geom_hline(yintercept = 3.5, size=1.5)+
   geom_vline(xintercept = c(3.5,6.5,9.5), size=1)+
   geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,11.5), size=0.2)+
#  geom_hline(yintercept = c(1.5,2.5,4.5,5.5,6.5), size=0.2)+
  xlab(label = "Sample")+
  ylab(label = "vOTU")+
  scale_fill_gradient(name="Relative abundance",low="#3794bf", high="#a34c1d")+
#  facet_grid(~ Depth)+
  theme_bw()
dev.off()

###heatmap for treatment
merged_virus.phyloseq<-merge_samples(virus.phyloseq,"Treat")#merge
merged_virus.phyloseq
r_merged_virus.phyloseq<-transform_sample_counts(merged_virus.phyloseq,function(x) x/sum(x))
r_merged_virus.phyloseq
head(otu_table(r_merged_virus.phyloseq))
t.merged.votu.tab <-(otu_table(r_merged_virus.phyloseq))
merged.votu.hp.df <- pivot_longer(as.data.frame(t.merged.votu.tab),cols = colnames(t.merged.votu.tab)) %>% 
  mutate(sample=c(rep("NCNTN0",260),rep("NCNTN60",260),rep("VNTN0",260),rep("VNTN60",260)))
head(merged.votu.hp.df)
tiff('merged.votu.hp.tiff', units="in", width=13, height=30, res=300)
ggplot(merged.votu.hp.df, aes(x =sample, y = name,fill=value)) + 
  geom_tile()+
  #  geom_hline(yintercept = 3.5, size=1.5)+
  geom_vline(xintercept = c(1.5,2.5,3.5), size=1)+
#  geom_vline(xintercept = c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,11.5), size=0.2)+
  #  geom_hline(yintercept = c(1.5,2.5,4.5,5.5,6.5), size=0.2)+
  xlab(label = "Sample")+
  ylab(label = "vOTU")+
#  labs(title = "Relative abundance")+
  scale_fill_gradient(name="Relative abundance",low="#3794bf", high="#a34c1d")+
  #  facet_grid(~ Depth)+
  theme_bw()
dev.off()


############diversity############

#alpha diversity
votu.diversity<-data.frame(estimate_richness(virus.phyloseq,split=TRUE,measures = NULL),sample_data(virus.phyloseq))
votu.diversity
tiff('boxvotu_Chao1 diversity.tiff', units="in", width=5, height=5, res=300)
votu_Chao1 <- ggplot(votu.diversity,aes(x=Treat,y=Chao1,fill=Treat))+
  geom_boxplot(size=0.5,alpha=0.8)+
  geom_jitter(shape=16, position=position_jitter(0))+
  scale_fill_manual(values = c("#4D897C", "#C6B7EC", "#ba892f","#c76d70"),name="Treatment")+
  labs(title="",x="Treatment",y="Chao1")+
  theme_bw()
dev.off()

votu.mod.Shannon <- lmer(Shannon ~ Cover* Nitrogen +(1|block), data = votu.diversity)
votu.mod.lm <- lmer(vOTU ~ Cover* Nitrogen+(1|block), data = votu.diversity)
Anova(votu.mod.Shannon)
normalTest(votu.mod.lm$residuals,"sw")
Anova(votu.mod.chao1)
normalTest(votu.diversity$Shannon,"sw")#votu normality is good
normalTest(votu.diversity$Chao1,"sw")

tiff('boxvotu_Shannon diversity.tiff', units="in", width=5, height=5, res=300)
votu_Shannon <- ggplot(votu.diversity,aes(x=Treat,y=Shannon,fill=Treat))+
  geom_boxplot(size=0.5,alpha=0.8,width=0.3)+
  geom_jitter(shape=16, position=position_jitter(0))+
  scale_fill_manual(values = c("#106494", "#c7a446", "#b3b2af", "#873207"),name="Treatment")+
  labs(title="",x="Treatment",y="Shannon")+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=13,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=13,family="Arial",face = "plain"),
        axis.title=element_text(size =13,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="plain",family = "Arial"),
        legend.text = (element_text(size=13,family = "Arial")),
        legend.box.margin=margin(8,8,8,8))
dev.off()

votu.mod.Shannon <- lmer(Shannon ~ Cover* Nitrogen +(1|block), data = votu.diversity)
Anova(votu.mod.Shannon)
normalTest(votu.diversity$Shannon,"sw")#votu normality is good

tiff('boxvotu_Simpson diversity.tiff', units="in", width=7, height=5, res=300)
votu_Simpson <- ggplot(votu.diversity,aes(x=Treat,y=Simpson,fill=Treat))+
  geom_boxplot(size=0.5,alpha=0.8)+
  geom_jitter(shape=16, position=position_jitter(0))+
  scale_fill_manual(values = c("#4D897C", "#C6B7EC", "#ba892f","#c76d70"),name="Treatment")+
  labs(title="",x="Treatment",y="Simpson")+
  theme_bw()
dev.off()

votu.mod.Simpson <- lmer(Simpson ~ Cover* Nitrogen +(1|block), data = votu.diversity)
Anova(votu.mod.Simpson)
normalTest(votu.diversity$Simpson,"sw")

#######PcoA##########
library(vegan)
#install.packages("ggforce")
library(ggforce)
virus.phyloseq_bray<-vegdist(t(otu_table(virus.phyloseq)),method="bray",binary = FALSE)
virus.phyloseq_PCoA<-ordinate(virus.phyloseq,method = "PCoA", virus.phyloseq_bray)
virus.phyloseq_nmds<-ordinate(virus.phyloseq,method = "NMDS", virus.phyloseq_bray)
tiff('PCoA_votu_nitrogen.tiff', units="in", width=7, height=5, res=300)
pcoa <- plot_ordination(virus.phyloseq,virus.phyloseq_PCoA,
                type="samples",color="Nitrogen", shape = "Cover")+
  geom_point(size=3)+
  geom_text(aes(label=sample_names(otu_table(virus.phyloseq)),vjust=-1),size=2.5)+
  #ggtitle("Principal Coordinates Analysis")+
  scale_color_manual(values = c("#106494", "#bf6a0a"))+
  xlab(label = "PCoA1 [31.9%]")+
  ylab(label = "PCoA2 [14.3%]")+
  geom_mark_ellipse(aes(color = Nitrogen, group=Nitrogen),expand=0)+
  geom_vline(xintercept = c(0), size=0.2,linetype=2)+
  geom_hline(yintercept = c(0), size=0.2,linetype=2)+
#  geom_polygon(aes(fill=c(rep("#1a7856",6), rep("#bf6a0a",6)),color = Nitrogen,group=Nitrogen))+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=13,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=13,family="Arial",face = "plain"),
        axis.title=element_text(size =13,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="plain",family = "Arial"),
        legend.text = (element_text(size=13,family = "Arial")))
dev.off()

tiff('PCoA_votu_cover.tiff', units="in", width=7, height=5, res=300)
plot_ordination(virus.phyloseq,virus.phyloseq_PCoA,
                type="samples",color="Cover", shape = "Nitrogen")+
  geom_point(size=3)+
  geom_text(aes(label=sample_names(otu_table(virus.phyloseq)),vjust=-1),size=2.5)+
  #ggtitle("Principal Coordinates Analysis")+
  scale_color_manual(values = c("#1a7856", "#bf6a0a"))+
  xlab(label = "PCoA1 [31.9%]")+
  ylab(label = "PCoA2 [14.3%]")+
  theme_bw()
dev.off()


tiff('NMDS_votu_nitrogen.tiff', units="in", width=7, height=5, res=300) #+ ggtitle("Agriculture soil ~ Treatment \n nmds plot ")#save as tiff
plot_ordination(virus.phyloseq,virus.phyloseq_nmds,color ="Nitrogen")+
  geom_point(size=3)+
  #stat_ellipse(geom = "polygon", type = "norm",aes(fill=Treat))+
  #geom_path(data=as.data.frame(sample_data(its.phyloseq.v20_30_sub)), aes(x = NMDS1, y = NMDS2, group = Treat))+
  geom_text(aes(x=NMDS1,y=NMDS2,label=sample_names(otu_table(virus.phyloseq))),size=1.5,vjust=-2)+
  scale_color_manual(values = c("#1a7856", "#bf6a0a")) +
  annotate(geom="text", x=0.4, y=0.34, label=paste("Stress:",round(virus.phyloseq_nmds$stress,digits = 2)),color="black")+
  annotate(geom="text", x=0.4, y=0.3, label=paste("p < 0.01"),color="black")+
theme_bw()
dev.off()

tiff('NMDS_votu_cover.tiff', units="in", width=7, height=5, res=300) #+ ggtitle("Agriculture soil ~ Treatment \n nmds plot ")#save as tiff
plot_ordination(virus.phyloseq,virus.phyloseq_nmds,color ="Cover")+
  geom_point(size=3)+
  #stat_ellipse(geom = "polygon", type = "norm",aes(fill=Treat))+
  #geom_path(data=as.data.frame(sample_data(its.phyloseq.v20_30_sub)), aes(x = NMDS1, y = NMDS2, group = Treat))+
  geom_text(aes(x=NMDS1,y=NMDS2,label=sample_names(otu_table(virus.phyloseq))),size=1.5,vjust=-2)+
  scale_color_manual(values = c("#1a7856", "#bf6a0a")) +
  annotate(geom="text", x=0.4, y=0.34, label=paste("Stress:",round(virus.phyloseq_nmds$stress,digits = 2)),color="black")+
  theme_bw()
dev.off()

set.seed(100)
may_pcoa_ord <- ordinate(physeq = virus.phyloseq, method = "CCA", distance =virus.phyloseq_bray, 
                        formula = ~Nitrogen+Cover+POXC+pH+GWC+NH4+NO3)
step.backward <-
  ordistep(may_pcoa_ord,
           permutations = how(nperm =9999)
  )
RsquareAdj(step.backward)

adonis(t(otu_table(virus.phyloseq))~Nitrogen+Cover+POXC+pH+GWC+NH4+NO3,data=data.frame(sample_data(virus.phyloseq)),permutations=9999,by="margin")#at least 1000
adonis(virus.phyloseq_bray~Nitrogen*Cover+POXC+pH+GWC+NH4+NO3,data=data.frame(sample_data(virus.phyloseq)),permutations=9999,by="margin")#at least 1000

#######stacked barplot#############
dim(votu.nor.tab)
head(votu.nor.tab)

rownames(votu.nor.tab)
votu.nor.tab.classified <- votu.nor.tab[c(55,44,41,27,53,103,30,229,63,93,69,100,149,9,145),]
dim(votu.nor.tab.classified)
head(votu.nor.tab.classified)
rownames(votu.nor.tab.classified) <-c("Podoviridae_uncl_1","Podoviridae_uncl_1",
                                       "Siphoviridae_uncl_1","Siphoviridae_uncl_1",
                                       "Podoviridae_uncl_2","Podoviridae_uncl_2",
                                       "Myoviridae_uncl_1 \n (g_Cr3virus)",
                                       "Siphoviridae_uncl_2",
                                       "g_Peduovirus","g_Peduovirus",
                                       "Podoviridae_uncl_3",
                                       "g_Bpp1virus",
                                       "Myoviridae_uncl_2",
                                       "g_Cecivirus",
                                       "Siphoviridae_uncl_3")
#write.csv(votu.nor.tab.classified.taxa,"votu.nor.tab.classified.taxa.csv")
votu.nor.tab.classified.taxa.unique <- rowsum(as.data.frame(votu.nor.tab.classified),rownames(votu.nor.tab.classified))
votu.nor.tab.classified.taxa.unique$NCNTN0 <- rowSums(votu.nor.tab.classified.taxa.unique[,1:3])
votu.nor.tab.classified.taxa.unique$NCNTN60 <- rowSums(votu.nor.tab.classified.taxa.unique[,4:6])
votu.nor.tab.classified.taxa.unique$VNTN0 <- rowSums(votu.nor.tab.classified.taxa.unique[,7:9])
votu.nor.tab.classified.taxa.unique$VNTN60 <- rowSums(votu.nor.tab.classified.taxa.unique[,10:12])
colnames(votu.nor.tab.classified.taxa.unique)
votu.nor.tab.classified.taxa.unique.new <- votu.nor.tab.classified.taxa.unique[,-c(1:12)]
r.votu.nor.tab.classified.taxa.unique.new <- apply(votu.nor.tab.classified.taxa.unique.new,2,function(x) x/sum(x))
r.votu.nor.tab <- mutate(as.data.frame(t(r.votu.nor.tab.classified.taxa.unique.new)),sample=colnames(r.votu.nor.tab.classified.taxa.unique.new)) 
r.votu.nor.trans <- pivot_longer(r.votu.nor.tab,cols = rownames(r.votu.nor.tab.classified.taxa.unique.new))

meta.data.new <- unique(meta.data[,-c(1,6)]) 
votu_stacked_df <- inner_join(r.votu.nor.trans,meta.data.new,by=c("sample"="Treat"))

head(votu_stacked_df)
votu_stacked_df$name


cols_phylum <- c("#994e03", "#9B110E", "#688787" ,"#d49783", "#550307", "#446455", "#FDD262", "#D3DDDC", "#C7B19C","#899DA4", "#FAEFD1", "#DC863B","#798E87", "#C27D38","#CCC591","#29211F")
tiff('figures/stacked_votu.treat.tiff', units="in", width=7, height=5, res=300)
stacked_bar <- ggplot(votu_stacked_df, aes(x =sample, y = value, fill=name))+
  geom_bar(stat='identity',colour="grey",size=0.1, width = 0.6,alpha=1)+
  labs(title="",x="Treatment",y="Relative abundance \n within classified viruses (%)")+
  scale_fill_manual(name=NULL,values=cols_phylum)+
  theme(legend.position="bottom",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor = element_blank(),
        plot.title=(element_text(size=13,family = "Arial",face="plain",vjust = 3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(size=13,family="Arial",face = "plain"),
        axis.title=element_text(size =13,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="plain",family = "Arial"),
        legend.text = (element_text(size=13,family = "Arial")))
dev.off()

###put together
install.packages("wesanderson")
library(wesanderson)
names(wes_palettes)
#https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
library(ggpubr)
tiff('figures/bar_venn.tiff', units="in", width=12, height=5, res=300)
ggarrange(venn,stacked_bar,
          labels = c("A","B"),widths = c(1,1.5),ncol = 2, nrow = 1,legend = "right",common.legend = F)
dev.off()
tiff('figures/shannon_pcoa.tiff', units="in", width=12, height=5, res=300)
ggarrange(votu_Shannon,pcoa,
          labels = c("C","D"),widths = c(1,1),ncol = 2, nrow = 1,legend = "bottom",common.legend = F)
dev.off()

tiff('figures/all.tiff', units="in", width=15, height=10, res=180)
ggarrange(venn,stacked_bar,votu_Shannon,pcoa,
          labels = c("(A)","(B)","(C)","(D)"),heights  = c(1,1),widths = c(1,1), ncol = 2, nrow = 2,legend = "right",common.legend = F,
          hjust=-1.2)
dev.off()

###########pie chart

# votu.total=260
# match.PIEGON=28
# unknow=260-28
# votu.network=107
# votu.network.outilier=30 #unknown=outliers+unknown,
# votu.network.vc=77
# votu.network.classified=15#vOTU_55,vOTU_44,vOTU_41,vOTU_27,vOTU_53,vOTU_103,vOTU_30,vOTU_229,vOTU_63,vOTU_93,vOTU_69,vOTU_100,vOTU_149,vOTU_9,vOTU_145 
# votu.network.unclass.vc=55









