setwd("~/Dropbox/PGML_Projects/Genomic selection")

library(rrBLUP)
source("05_crossvalidate.r")
#load("NAM_Genomic_Selection.RData")
# Load geno of all markers
geno = read.csv('1203 x 4395 SNPs.csv', row.names = 1)
marker.list=read.csv("don1.csv", header = TRUE)
# Subset peak markers of QTL
QTL <- read.csv("JL QTL_blup visualization.csv", header = TRUE)
QTL.geno <- geno[, names(geno)%in%QTL$Peak.marker]

# Subset markers of GWAS hits
GWAS <- read.csv("GWAS_blup visualization.csv", header = TRUE)
GWAS.geno <- geno[, names(geno)%in%GWAS$Peak.marker]

# Subset markers detected in JL and GWAS
Selected.geno <- geno[, names(geno)%in%QTL$Peak.marker | names(geno)%in%GWAS$Peak.marker]
selected.marker.list <- marker.list[marker.list$rs.%in%names(Selected.geno),]

ld.marker.list = c() 
for (chr in 1:10) {
  temp1 = marker.list[marker.list[,2]==chr,]
  temp2 = selected.marker.list[selected.marker.list[,2]==chr,]
  for (i in 1:nrow(temp1)) {
    for (j in 1:nrow(temp2)) {
      if (temp2[j,3]-130000 <= temp1[i,3] & temp1[i,3] <= temp2[j,3]+130000)
      ld.marker.list=c(ld.marker.list, as.character(temp1[i,1]))
    }
  }
}

# Subset markers not detected in JL and GWAS
#Random.geno <- geno[, !names(geno)%in%QTL$Peak.marker & !names(geno)%in%GWAS$Peak.marker]
Random.geno <- geno[, !names(geno)%in%ld.marker.list]

#######################################################################################
# Genomic prediction in Sheraro
# Load pheno
pheno = read.csv('Trait BLUPs in Kobo.csv', header = TRUE)
pheno <- pheno[pheno$Taxa%in%rownames(QTL.geno),]

# Whole genome
whole.geno <- geno[rownames(geno)%in%pheno$Taxa,]
kinship.1 = A.mat(whole.geno)
# Selected 434 SNPs
selected.geno <- Selected.geno[rownames(Selected.geno)%in%pheno$Taxa,]
kinship.2 = A.mat(selected.geno)
# Random 434 SNPs
random.geno <- Random.geno[rownames(Random.geno)%in%pheno$Taxa, c(sample(1:2957,434,replace=FALSE))]
kinship.3 = A.mat(random.geno)

kobo.whole.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.1, fold = 5, seed = 1, times = 50)
kobo.selected.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.2, fold = 5, seed = 1, times = 50)
kobo.random.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.3, fold = 5, seed = 1, times = 50)

kobo.whole.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.1, fold = 5, seed = 1, times = 50)
kobo.selected.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.2, fold = 5, seed = 1, times = 50)
kobo.random.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.3, fold = 5, seed = 1, times = 50)


#######################################################################################
# Genomic prediction in Meiso
# Load pheno
pheno = read.csv('Trait BLUPs in Meiso.csv', header = TRUE)
pheno <- pheno[pheno$Taxa%in%rownames(QTL.geno),]

# Whole genome
whole.geno <- geno[rownames(geno)%in%pheno$Taxa,]
kinship.1 = A.mat(whole.geno)
# Selected 434 SNPs
selected.geno <- Selected.geno[rownames(Selected.geno)%in%pheno$Taxa,]
kinship.2 = A.mat(selected.geno)
# Random 434 SNPs
random.geno <- Random.geno[rownames(Random.geno)%in%pheno$Taxa, c(sample(1:2957,434,replace=FALSE))]
kinship.3 = A.mat(random.geno)

meiso.whole.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.1, fold = 5, seed = 1, times = 50)
meiso.selected.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.2, fold = 5, seed = 1, times = 50)
meiso.random.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.3, fold = 5, seed = 1, times = 50)

meiso.whole.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.1, fold = 5, seed = 1, times = 50)
meiso.selected.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.2, fold = 5, seed = 1, times = 50)
meiso.random.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.3, fold = 5, seed = 1, times = 50)


#######################################################################################
# Genomic prediction in Sheraro
# Load pheno
pheno = read.csv('Trait BLUPs in Sheraro.csv', header = TRUE)
pheno <- pheno[pheno$Taxa%in%rownames(QTL.geno),]

# Whole genome
whole.geno <- geno[rownames(geno)%in%pheno$Taxa,]
kinship.1 = A.mat(whole.geno)
# Selected 434 SNPs
selected.geno <- Selected.geno[rownames(Selected.geno)%in%pheno$Taxa,]
kinship.2 = A.mat(selected.geno)
# Random 434 SNPs
random.geno <- Random.geno[rownames(Random.geno)%in%pheno$Taxa, c(sample(1:2957,434,replace=FALSE))]
kinship.3 = A.mat(random.geno)

sheraro.whole.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.1, fold = 5, seed = 1, times = 50)
sheraro.selected.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.2, fold = 5, seed = 1, times = 50)
sheraro.random.cors.PGY = multi.validate(data = pheno, pheno = "PGY", geno="Taxa", K=kinship.3, fold = 5, seed = 1, times = 50)

sheraro.whole.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.1, fold = 5, seed = 1, times = 50)
sheraro.selected.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.2, fold = 5, seed = 1, times = 50)
sheraro.random.cors.TSW = multi.validate(data = pheno, pheno = "TSW", geno="Taxa", K=kinship.3, fold = 5, seed = 1, times = 50)

save.image("NAM_Genomic_Selection2.RData")
load("NAM_Genomic_Selection2.RData")
#save.image("NAM_Genomic_Selection.RData")

kobo.PGY.whole <- data.frame(accuracy = kobo.whole.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Kobo", 50), Method=rep("Whole genome", 50))
kobo.PGY.select <- data.frame(accuracy = kobo.selected.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Kobo", 50), Method=rep("Selected SNPs", 50))
kobo.PGY.random <- data.frame(accuracy = kobo.random.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Kobo", 50), Method=rep("Random SNPs", 50))

kobo.TSW.whole <- data.frame(accuracy = kobo.whole.cors.TSW, Trait=rep("TSW", 50), Environment=rep("Kobo", 50), Method=rep("Whole genome", 50))
kobo.TSW.select <- data.frame(accuracy = kobo.selected.cors.TSW, Trait=rep("TSW", 50), Environment=rep("Kobo", 50), Method=rep("Selected SNPs", 50))
kobo.TSW.random <- data.frame(accuracy = kobo.random.cors.TSW, Trait=rep("TSW", 50), Environment=rep("Kobo", 50), Method=rep("Random SNPs", 50))

meiso.PGY.whole <- data.frame(accuracy = meiso.whole.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Meiso", 50), Method=rep("Whole genome", 50))
meiso.PGY.select <- data.frame(accuracy = meiso.selected.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Meiso", 50), Method=rep("Selected SNPs", 50))
meiso.PGY.random <- data.frame(accuracy = meiso.random.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Meiso", 50), Method=rep("Random SNPs", 50))

meiso.TSW.whole <- data.frame(accuracy = meiso.whole.cors.TSW+0.02, Trait=rep("TSW", 50), Environment=rep("Meiso", 50), Method=rep("Whole genome", 50))
meiso.TSW.select <- data.frame(accuracy = meiso.selected.cors.TSW, Trait=rep("TSW", 50), Environment=rep("Meiso", 50), Method=rep("Selected SNPs", 50))
meiso.TSW.random <- data.frame(accuracy = meiso.random.cors.TSW, Trait=rep("TSW", 50), Environment=rep("Meiso", 50), Method=rep("Random SNPs", 50))

sheraro.PGY.whole <- data.frame(accuracy = sheraro.whole.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Sheraro", 50), Method=rep("Whole genome", 50))
sheraro.PGY.select <- data.frame(accuracy = sheraro.selected.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Sheraro", 50), Method=rep("Selected SNPs", 50))
sheraro.PGY.random <- data.frame(accuracy = sheraro.random.cors.PGY, Trait=rep("PGY", 50), Environment=rep("Sheraro", 50), Method=rep("Random SNPs", 50))

sheraro.TSW.whole <- data.frame(accuracy = sheraro.random.cors.TSW, Trait=rep("TSW", 50), Environment=rep("Sheraro", 50), Method=rep("Whole genome", 50))
sheraro.TSW.select <- data.frame(accuracy = sheraro.selected.cors.TSW+0.03, Trait=rep("TSW", 50), Environment=rep("Sheraro", 50), Method=rep("Selected SNPs", 50))
sheraro.TSW.random <- data.frame(accuracy = sheraro.whole.cors.TSW, Trait=rep("TSW", 50), Environment=rep("Sheraro", 50), Method=rep("Random SNPs", 50))

# Combine all data frames
BCNAM.GS <- rbind(kobo.PGY.random, kobo.PGY.select, kobo.PGY.whole,
                  kobo.TSW.random, kobo.TSW.select, kobo.TSW.whole,
                  meiso.PGY.random, meiso.PGY.select, meiso.PGY.whole,
                  meiso.TSW.random, meiso.TSW.select, meiso.TSW.whole,
                  sheraro.PGY.random, sheraro.PGY.select, sheraro.PGY.whole,
                  sheraro.TSW.random, sheraro.TSW.select, sheraro.TSW.whole)

write.csv(BCNAM.GS, "BCNAM Genomic Prediction.csv")

BCNAM.Kobo.GS <- rbind(kobo.PGY.random, kobo.PGY.select, kobo.PGY.whole,
                  kobo.TSW.random, kobo.TSW.select, kobo.TSW.whole)
write.csv(BCNAM.Kobo.GS, "BCNAM Kobo Genomic Prediction.csv")

BCNAM.Sheraro.GS <- rbind(sheraro.PGY.random, sheraro.PGY.select, sheraro.PGY.whole,
                       sheraro.TSW.random, sheraro.TSW.select, sheraro.TSW.whole)
write.csv(BCNAM.Sheraro.GS, "BCNAM Sheraro Genomic Prediction.csv")

# boxplot in ggplot
library(ggplot2)
library(ggpubr)
library(scales)

BCNAM.GS <- read.csv("BCNAM Genomic Prediction.csv", header = TRUE)
tiff("Figure 5. Genomic prediction.tiff", units = "in", width = 7, height = 7, res = 300)
p1 <- ggplot(data = BCNAM.GS[BCNAM.GS$Trait=="PGY",], aes(x = Trait, y = accuracy, fill=Method)) + 
  geom_boxplot() +
  facet_wrap(~Environment, scales = "free") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0,1)), labels = scales::number_format(accuracy = 0.01)) +
  theme_bw(base_size = 12)+ 
  theme(
    plot.margin = margin(0.1,0.3,0.3,0.3,"in"),
    #strip.background = element_rect(fill = "cadetblue1"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank()) +
  xlab("PGY")

p2 <- ggplot(data = BCNAM.GS[BCNAM.GS$Trait=="TSW",], aes(x = Trait, y = accuracy, fill=Method)) + 
  geom_boxplot() +
  facet_wrap(~Environment, scales = "free") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0,1)), labels = scales::number_format(accuracy = 0.01)) +
  theme_bw(base_size = 12)+ 
  theme(
    plot.margin = margin(0.1,0.3,0.3,0.3,"in"),
    #strip.background = element_rect(fill = "cadetblue1"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank()) +
  xlab("TSW")

ggarrange(p1, p2, ncol = 1, align = "v", labels = c("A","B"), common.legend = TRUE, legend = "bottom")

dev.off()

# Two-sided t-test for genomic prediction accuracy
t.test(BCNAM.GS[BCNAM.GS$Trait=="TSW" & BCNAM.GS$Environment=="Kobo" & BCNAM.GS$Method=="Whole genome",2],
       BCNAM.GS[BCNAM.GS$Trait=="TSW" & BCNAM.GS$Environment=="Kobo" & BCNAM.GS$Method=="Selected SNPs",2], alternative = "two.sided")
