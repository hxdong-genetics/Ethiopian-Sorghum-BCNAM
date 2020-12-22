setwd("~/Dropbox/Postdoc Projects/Ethiopia_BCNAM/2 association mapping in R")
load("Figure 5.RData")
#save.image("Figure 5.RData")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(grid)
library(gridExtra)

# import QTLs from JL analysis:
JL.QTL <- read.csv("QTL visualization_effect_blup.csv", header = TRUE)

# plot QTL effects by trait
tiff("Figure S4.2 JL QTL effect.tiff", units = "in", width = 6, height = 6, res = 300)
JL.QTL.1 <- JL.QTL[!(JL.QTL$Trait=="DF" | JL.QTL$Trait=="DM" | JL.QTL$Trait=="PH"),]
ggplot(data = JL.QTL, aes(x=Effect)) +
#ggplot(data = JL.QTL, aes(x=Effect)) +
  facet_wrap(.~Trait, scales = "free") +
  geom_histogram()+
  ylab("") +
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 9), 
    axis.text.y = element_text(color = "black", size = 9))
dev.off()

#import map, genotypes, and phenotypes (BLUPs)
map=read.csv("map_4395 snps.csv",header=T)
map <- map[,2:4]
geno=read.csv("1203 x 4395 SNPs.csv",header=T)
pheno=read.csv("BLUP_Meiso.csv",header=T)

# Remove IS32234:
pheno=pheno[!pheno$Population=="IS32234",]

# subset data of the common set of taxa:
geno <- geno[geno$Taxa%in%pheno$Taxa,]
pheno <- pheno[pheno$Taxa%in%geno$Taxa,]
geno=as.matrix(apply(geno[,-1], 2, as.numeric)) #change marker data to class numeric
rownames(geno) <- pheno$Taxa
#geno <- geno[,colnames(geno)%in%map$rs.]
dim(geno) #1031 x 4395 snps in Kobo; 1155 x 4395 snps in Meiso; 1046 x 4395 snps in Sheraro; 


#STEP 1. conduct single marker tests: regress phenotype on all markers, one at a time 
# lm(my_pheno~family + my_geno) regresses phenotype on genotype, and returns the intercept and slope of the regression line
# anova()[2,5] performs ANOVA on the linear model and returns the value from row 2, col 5 of the resulting table (p-value for my_geno)
pvals=c()
for (i in 1:ncol(geno)){
  pvals[i]= anova(lm(pheno$HE ~ pheno$Population + geno[,i]))[2,5]
  print(i)
}
plot(-log10(pvals),pch=19, cex=0.5, main = "")

#which marker is most significant?
#which(pvals==min(pvals)) 
#colnames(geno)[which.min(pvals)]

# plot results with qqman
gwasResult <- cbind(map, pvals)
gwasResult[c(is.na(gwasResult$pvals)),4] <- 1 # convert NA to 1.

#write.csv(gwasResult[gwasResult$pvals <= 0.001,], "Sheraro PGY GWAS hits.csv")

Meiso.DF.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Meiso.DM.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Meiso.PH.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Meiso.NL.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Meiso.HE.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Meiso.CPP.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Meiso.PGY.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Kobo.DF.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Kobo.DM.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Kobo.PH.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Kobo.NL.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Kobo.PGY.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Kobo.CPP.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Sheraro.DF.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Sheraro.DM.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Sheraro.PH.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Sheraro.NL.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Sheraro.CPP.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

Sheraro.PGY.don1 <- gwasResult %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResult, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum=pos+tot)

axisdf = Kobo.DF.don1 %>% group_by(chrom) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
vline <- Kobo.DF.don1 %>% select(chrom, BPcum) %>% group_by(chrom) %>% summarise(vline.max=max(BPcum))




tiff("Figure 5.tiff", units = "in", width = 8.5, height = 11, res = 300)

Kobo.DF.p.GWAS <- ggplot(Kobo.DF.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =5.5, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 6745069, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=23145069, y =5.5, label="Ma5", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 60852774, lty=1, lwd=1, alpha=0.4, color="green") + # Ma3
  annotate("text", x=76852774, y =5.5, label="Ma3", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 478434065, lty=1, lwd=1, alpha=0.4, color="green") + # DREB1A
  annotate("text", x=478434065, y =5.5, label="DREB1A", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 475364000, lty=1, lwd=1, alpha=0.4, color="green") + # Dw3
  #annotate("text", x=495364000, y =5.5, label="Dw3", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 154284883, lty=1, lwd=1, alpha=0.4, color="green") + # SbGI
  #geom_vline(data = NULL, mapping = NULL, xintercept = 213228688, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN12
  #annotate("text", x=233228688, y =6.5, label="SbCN12") +
  #geom_vline(data = NULL, xintercept = 591237756, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN8
  #annotate("text", x=571237756, y =6.5, label="SbCN8") +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 9), 
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))



Kobo.DF.p.JL <- ggplot(Kobo.DF.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 6745069, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=23145069, y =8.5, label="Ma5", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 60852774, lty=1, lwd=1, alpha=0.4, color="green") + # Ma3
  annotate("text", x=76852774, y =8.5, label="Ma3", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 478434065, lty=1, lwd=1, alpha=0.4, color="green") + # DREB1A
  annotate("text", x=478434065, y =8.5, label="DREB1A", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 475364000, lty=1, lwd=1, alpha=0.4, color="green") + # Dw3
  #annotate("text", x=495364000, y =8.5, label="Dw3", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 421521725, lty=1, lwd=1, alpha=0.4, color="green") + # LHY
  #annotate("text", x=441521725, y =8.5, label="LHY", size=3.5, fontface="italic") +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  
  geom_point(data=JL.QTL[JL.QTL$Environment=="Kobo" & JL.QTL$Trait=="DF",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  scale_color_gradient(low = "orange", high = "blue") +
  
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Kobo.DF.p <- ggarrange(Kobo.DF.p.JL, Kobo.DF.p.GWAS, ncol=1, nrow=2, align = "hv", heights = c(1.2,1))




# Meiso DF
Meiso.DF.p.GWAS <- ggplot(Meiso.DF.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,8), expand = c(0,0)) +
  
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =7, label="Ma6", size=3.5) +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 154284883, lty=1, lwd=1, alpha=0.4, color="green") + # SbGI
  #geom_vline(data = NULL, mapping = NULL, xintercept = 213228688, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN12
  #annotate("text", x=233228688, y =7, label="SbCN12", size=3.5, fontface="italic") +
  geom_vline(data = NULL, xintercept = 591237756, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN8
  annotate("text", x=591237756, y =7, label="SbCN8", size=3.5, fontface="italic") +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 9), 
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))



Meiso.DF.p.JL <- ggplot(Meiso.DF.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, xintercept = 591237756, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN8
  annotate("text", x=591237756, y =8.5, label="SbCN8", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 478434065, lty=1, lwd=1, alpha=0.4, color="green") + # DREB1A
  annotate("text", x=478434065, y =8.5, label="DREB1A", size=3.5, fontface="italic") +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  geom_point(data=JL.QTL[JL.QTL$Environment=="Meiso" & JL.QTL$Trait=="DF",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Meiso.DF.p <- ggarrange(Meiso.DF.p.JL, Meiso.DF.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))


# Sheraro DF
Sheraro.DF.p.GWAS <- ggplot(Sheraro.DF.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  
  labs(x="Chromosome", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 6745069, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=23145069, y =8.5, label="Ma5", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 421521725, lty=1, lwd=1, alpha=0.4, color="green") + # LHY
  annotate("text", x=441521725, y =8.5, label="LHY", size=3.5, fontface="italic") +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 9), 
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))



Sheraro.DF.p.JL <- ggplot(Sheraro.DF.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 6745069, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=23145069, y =8.5, label="Ma5", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 421521725, lty=1, lwd=1, alpha=0.4, color="green") + # LHY
  annotate("text", x=441521725, y =8.5, label="LHY", size=3.5, fontface="italic") +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  
  geom_point(data=JL.QTL[JL.QTL$Environment=="Sheraro" & JL.QTL$Trait=="DF",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Sheraro.DF.p <- ggarrange(Sheraro.DF.p.JL, Sheraro.DF.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))



# Kobo NL
Kobo.NL.p.GWAS <- ggplot(Kobo.NL.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,7), expand = c(0,0)) +
  
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =6, label="Ma6", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 218863379, lty=1, lwd=1, alpha=0.4, color="green") + # P5CS2
  #annotate("text", x=238863379, y =6, label="P5CS2", size=3.5) +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 475364000, lty=1, lwd=1, alpha=0.4, color="green") + # LHY
  #annotate("text", x=495364000, y =6, label="Dw3", size=3.5) +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 154284883, lty=1, lwd=1, alpha=0.4, color="green") + # SbGI
  #geom_vline(data = NULL, mapping = NULL, xintercept = 213228688, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN12
  #annotate("text", x=233228688, y =5.5, label="SbCN12", size=3.5) +
  #geom_vline(data = NULL, xintercept = 591237756, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN8
  #annotate("text", x=571237756, y =6.5, label="SbCN8") +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 9), 
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))



Kobo.NL.p.JL <- ggplot(Kobo.NL.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, xintercept = 218863379, lty=1, lwd=1, alpha=0.4, color="green") + # P5CS2
  annotate("text", x=218863379, y =8.5, label="P5CS2", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 295802230, lty=1, lwd=1, alpha=0.4, color="green") + # LEA
  annotate("text", x=305802230, y =8.5, label="LEA", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 6745069, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  #annotate("text", x=23145069, y =8.5, label="Ma5", size=3.5) +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 213228688, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN12
  #annotate("text", x=243228688, y =8.5, label="SbCN12", size=3.5) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  geom_point(data=JL.QTL[JL.QTL$Environment=="Kobo" & JL.QTL$Trait=="NL",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Kobo.NL.p <- ggarrange(Kobo.NL.p.JL, Kobo.NL.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))


# Meiso NL
Meiso.NL.p.GWAS <- ggplot(Meiso.NL.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =13, label="Ma6", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 6745069, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  #annotate("text", x=23145069, y =13, label="Ma5", size=3.5) +
  geom_vline(data = NULL, mapping = NULL, xintercept = 478434065, lty=1, lwd=1, alpha=0.4, color="green") + # DREB1A
  annotate("text", x=478434065, y =13, label="DREB1A", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 295802230, lty=1, lwd=1, alpha=0.4, color="green") + # LEA
  annotate("text", x=305802230, y =13, label="LEA", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, mapping = NULL, xintercept = 154284883, lty=1, lwd=1, alpha=0.4, color="green") + # SbGI
  #geom_vline(data = NULL, mapping = NULL, xintercept = 213228688, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN12
  #annotate("text", x=233228688, y =6.5, label="SbCN12") +
  #geom_vline(data = NULL, xintercept = 591237756, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN8
  #annotate("text", x=571237756, y =6.5, label="SbCN8") +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 9), 
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))



Meiso.NL.p.JL <- ggplot(Meiso.NL.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, xintercept = 218863379, lty=1, lwd=1, alpha=0.4, color="green") + # P5CS2
  #annotate("text", x=238863379, y =8.5, label="P5CS2", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 295802230, lty=1, lwd=1, alpha=0.4, color="green") + # LEA
  annotate("text", x=305802230, y =8.5, label="LEA", size=3.5, fontface="italic") +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  geom_point(data=JL.QTL[JL.QTL$Environment=="Meiso" & JL.QTL$Trait=="NL",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Meiso.NL.p <- ggarrange(Meiso.NL.p.JL, Meiso.NL.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))






# Sheraro NL
Sheraro.NL.p.GWAS <- ggplot(Sheraro.NL.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,7), expand = c(0,0)) +
  
  labs(x="Chromosome", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  geom_vline(data = NULL, mapping = NULL, xintercept = 60852774, lty=1, lwd=1, alpha=0.4, color="green") + # Ma3
  annotate("text", x=72852774, y =5.5, label="Ma3", size=3.5) +

  #geom_vline(data = NULL, mapping = NULL, xintercept = 154284883, lty=1, lwd=1, alpha=0.4, color="green") + # SbGI
  #geom_vline(data = NULL, mapping = NULL, xintercept = 213228688, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN12
  #annotate("text", x=233228688, y =6.5, label="SbCN12") +
  #geom_vline(data = NULL, xintercept = 591237756, lty=1, lwd=1, alpha=0.4, color="green") + # SbCN8
  #annotate("text", x=571237756, y =6.5, label="SbCN8") +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 9), 
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(margin = margin(0,-40,0,0)))



Sheraro.NL.p.JL <- ggplot(Sheraro.NL.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 60852774, lty=1, lwd=1, alpha=0.4, color="green") + # Ma3
  annotate("text", x=72852774, y =8.5, label="Ma3", size=3.5) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  geom_point(data=JL.QTL[JL.QTL$Environment=="Sheraro" & JL.QTL$Trait=="NL",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Sheraro.NL.p <- ggarrange(Sheraro.NL.p.JL, Sheraro.NL.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))



# Meiso HE
Meiso.HE.p.GWAS <- ggplot(Meiso.HE.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,8), expand = c(0,0)) +

  labs(x="Chromosome", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =7, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, xintercept = 69302601, lty=1, lwd=1, alpha=0.4, color="green") + # KCS20
  annotate("text", x=121302601, y =7, label="KCS20", size=3.5, fontface="italic") +
  #geom_vline(data = NULL, xintercept = 218863379, lty=1, lwd=1, alpha=0.4, color="green") + # P5CS2
  #annotate("text", x=218863379, y =7, label="P5CS2", size=3.5, fontface="italic") +
  
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.y = element_text(margin = margin(0,-40,0,0)))



Meiso.HE.p.JL <- ggplot(Meiso.HE.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  geom_vline(data = NULL, xintercept = 69302601, lty=1, lwd=1, alpha=0.4, color="green") + # KCS20
  annotate("text", x=121302601, y =8.5, label="KCS20", size=3.5, fontface="italic") +
  geom_vline(data = NULL, xintercept = 218863379, lty=1, lwd=1, alpha=0.4, color="green") + # P5CS2
  annotate("text", x=218863379, y =8.5, label="P5CS2", size=3.5, fontface="italic") +
  geom_vline(data = NULL, mapping = NULL, xintercept = 295802230, lty=1, lwd=1, alpha=0.4, color="green") + # LEA
  annotate("text", x=305802230, y =8.5, label="LEA", size=3.5, fontface="italic") +

  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  geom_point(data=JL.QTL[JL.QTL$Environment=="Meiso" & JL.QTL$Trait=="HE",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Meiso.HE.p <- ggarrange(Meiso.HE.p.JL, Meiso.HE.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))


ggarrange(Kobo.DF.p, Kobo.NL.p, Meiso.DF.p, Meiso.NL.p, Sheraro.DF.p, Meiso.HE.p, ncol=2, nrow=3, align = "v",
          labels = c("A","B","C","D","E","F"))

dev.off()

save.image("Figure 5.RData")

tiff("Figure 5.2.tiff", units = "in", width = 6, height = 9, res = 300)
ggarrange(Kobo.DF.p, Meiso.DF.p, Sheraro.DF.p, ncol=1, nrow=3, align = "v",
          labels = c("A","B","C"))

dev.off()











# Kobo PH
Kobo.PH.p.GWAS <- ggplot(Kobo.PH.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +

  geom_vline(data = NULL, mapping = NULL, xintercept = 403299966, lty=1, lwd=1, alpha=0.4, color="green") + # Dw2
  annotate("text", x=403299966, y =5.5, label="Dw2", size=3.5, fontface="italic") +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.y = element_text(margin = margin(0,-40,0,0)))



Kobo.PH.p.JL <- ggplot(Kobo.PH.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  #geom_vline(data = NULL, mapping = NULL, xintercept = 475364000, lty=1, lwd=1, alpha=0.4, color="green") + # Dw3
  #annotate("text", x=495364000, y =8.5, label="Dw3", size=3.5, fontface="italic") +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  
  geom_point(data=JL.QTL[JL.QTL$Environment=="Kobo" & JL.QTL$Trait=="PH",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  scale_color_gradient(low = "orange", high = "blue") +
  
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Kobo.PH.p <- ggarrange(Kobo.PH.p.JL, Kobo.PH.p.GWAS, ncol=1, nrow=2, align = "hv", heights = c(1.2,1))




# Meiso PH
Meiso.PH.p.GWAS <- ggplot(Meiso.PH.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,8), expand = c(0,0)) +
  
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  #geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  #annotate("text", x=375515140, y =7, label="Ma6", size=3.5) +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.y = element_text(margin = margin(0,-40,0,0)))



Meiso.PH.p.JL <- ggplot(Meiso.PH.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  #geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  #annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +

  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  geom_point(data=JL.QTL[JL.QTL$Environment=="Meiso" & JL.QTL$Trait=="PH",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Meiso.PH.p <- ggarrange(Meiso.PH.p.JL, Meiso.PH.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))


# Sheraro PH
Sheraro.PH.p.GWAS <- ggplot(Sheraro.PH.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  
  labs(x="Chromosome", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  geom_vline(data = NULL, xintercept = 288340176, lty=1, lwd=1, alpha=0.4, color="green") + # Dw4
  annotate("text", x=288340176, y =7, label="Dw4", size=3.5) +
  geom_vline(data = NULL, xintercept = 585896198, lty=1, lwd=1, alpha=0.4, color="green") + # Dw1
  annotate("text", x=585896198, y =7, label="Dw1", size=3.5) +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.y = element_text(margin = margin(0,-40,0,0)))



Sheraro.PH.p.JL <- ggplot(Sheraro.PH.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  
  geom_vline(data = NULL, xintercept = 288340176, lty=1, lwd=1, alpha=0.4, color="green") + # Dw4
  annotate("text", x=288340176, y =8.5, label="Dw4", size=3.5) +

  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  
  geom_point(data=JL.QTL[JL.QTL$Environment=="Sheraro" & JL.QTL$Trait=="PH",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Sheraro.PH.p <- ggarrange(Sheraro.PH.p.JL, Sheraro.PH.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))

tiff("Figure 5.3.tiff", units = "in", width = 6, height = 9, res = 300)
ggarrange(Kobo.PH.p, Meiso.PH.p, Sheraro.PH.p, ncol=1, nrow=3, align = "v",
          labels = c("A","B","C"))

dev.off()












# Kobo DM
Kobo.DM.p.GWAS <- ggplot(Kobo.DM.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) +
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.y = element_text(margin = margin(0,-40,0,0)))



Kobo.DM.p.JL <- ggplot(Kobo.DM.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  #geom_vline(data = NULL, mapping = NULL, xintercept = 475364000, lty=1, lwd=1, alpha=0.4, color="green") + # Dw3
  #annotate("text", x=495364000, y =8.5, label="Dw3", size=3.5, fontface="italic") +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  
  geom_point(data=JL.QTL[JL.QTL$Environment=="Kobo" & JL.QTL$Trait=="DM",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  scale_color_gradient(low = "orange", high = "blue") +
  
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Kobo.DM.p <- ggarrange(Kobo.DM.p.JL, Kobo.DM.p.GWAS, ncol=1, nrow=2, align = "hv", heights = c(1.2,1))




# Meiso DM
Meiso.DM.p.GWAS <- ggplot(Meiso.DM.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,8), expand = c(0,0)) +
  
  labs(x="", y=expression(paste("-log"[10],"(",italic("P"),")"))) +
  
  #geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  #annotate("text", x=375515140, y =7, label="Ma6", size=3.5) +
  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.y = element_text(margin = margin(0,-40,0,0)))



Meiso.DM.p.JL <- ggplot(Meiso.DM.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  #geom_vline(data = NULL, xintercept = 355515140, lty=1, lwd=1, alpha=0.4, color="green") + # Ma6
  #annotate("text", x=375515140, y =8.5, label="Ma6", size=3.5, fontface="italic") +
  
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  geom_point(data=JL.QTL[JL.QTL$Environment=="Meiso" & JL.QTL$Trait=="DM",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Meiso.DM.p <- ggarrange(Meiso.DM.p.JL, Meiso.DM.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))


# Sheraro DM
Sheraro.DM.p.GWAS <- ggplot(Sheraro.DM.don1, aes(x=BPcum, y=-log10(pvals), color=factor(chrom))) + 
  #scale_color_manual(values = c("red","blue","green","orange","grey","darkgreen","purple","gold","black","cyan")) +
  scale_color_manual(values = c(rep(c("cadetblue", "darkslateblue"),5))) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) +
  
  labs(x="Chromosome", y=expression(paste("-log"[10],"(",italic("P"),")"))) +

  geom_hline(yintercept = 3, lty=2) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="gray") +
  geom_point(pch=19, cex=1) +
  theme_bw()+
  # Custom the theme:
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.y = element_text(margin = margin(0,-40,0,0)))



Sheraro.DM.p.JL <- ggplot(Sheraro.DM.don1, aes(x=BPcum, y=rep(0,4395))) + geom_point(pch="|", cex=0.01) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,9.5), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,1)) +
  labs(x="", y="") +
  # Custom the theme:
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"pt"),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9)) +
  
  scale_color_gradient(low = "orange", high = "blue") +
  
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2, color="grey") +
  
  geom_point(data=JL.QTL[JL.QTL$Environment=="Sheraro" & JL.QTL$Trait=="DM",], 
             mapping = aes(x=BPcum, y=Population, colour=sign(Effect))) +
  theme_bw()+
  guides(size=FALSE) + # remove legend for size
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9), 
        axis.text.y = element_text(color = "black", size = 9))

Sheraro.DM.p <- ggarrange(Sheraro.DM.p.JL, Sheraro.DM.p.GWAS, ncol=1, nrow=2, align = "v", heights = c(1.2,1))

tiff("Figure 5.5.tiff", units = "in", width = 6, height = 9, res = 300)
ggarrange(Kobo.NL.p, Meiso.NL.p, Meiso.HE.p, ncol=1, nrow=3, align = "v",
          labels = c("A","B","C"))

dev.off()
