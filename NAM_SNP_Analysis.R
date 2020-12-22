setwd("~/OneDrive - University of Georgia/Postdoc Projects/Ethiopia_BCNAM/GBS")
library(ggplot2)
library(ggpubr)
library(dplyr)

Depth <- read.table("NAM_1248_Taxa_16945_SNP_mnQS13_gdepth.txt", header = TRUE)
ReadsPerSample <- read.csv("NAMReadsPerSample.csv", header = TRUE)
# Average depth for 1248 inds x 12407 SNPs:
mean(unlist(Depth[,-c(1,2)])) # 22X
# Average depth for 14 parents x 12407 SNPs:
mean(unlist(Depth[,c(3:16)])) # 70X
# Average depth for 1234 progeny x 12407 SNPs:
mean(unlist(Depth[,-c(1:16)])) # 22X
# Plot good barcoded reads per sample
ggplot(ReadsPerSample, aes(x=FullSampleName, y=goodReadsMatchedToDataBase/16945)) + geom_col(aes(color=Population)) + labs(x="Genotype", y="Average sequencing depth") + theme_bw(base_size = 12) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom", legend.text=element_text(size=11), axis.text.x = element_blank(), axis.text.y = element_text(color = "black", size = 11)) + 
  scale_color_manual(values = c("red","blue","yellow","green","purple","gray","black","darkgreen","orange","magenta","cyan","skyblue","yellowgreen","pink")) + 
  geom_abline(slope = 0, intercept = 22, lty=2)



# load hapmap file
SNP.hapmap <- read.table("10559_SNPs_1248_Indv.hmp.txt", header=T, sep="\t",comment.char="",stringsAsFactors=F)
SNP.numeric <- read.table("10559_SNPs_1248_Indv_Num.txt", header=T, row.names = 1, sep="\t",comment.char="",stringsAsFactors=F)
#load package to convert SNPs into 0,1,and 2 format 
#0 is homozygote, 1 is heterozygote and 2 is another type of homozygote.
library(adegenet) 
hapMap2genlight <- function(file){
  require(adegenet)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)[,-(2:10)]
  samples <- names(hapmap)[-1]
  loci <- row.names(hapmap)
  
  s <- as.integer(c(0,1,2,NA))
  ac <- s
  ag <- s
  at <- s
  cg <- s
  ct <- s
  gt <- s
  names(ac) <- c("A","M","C","N")
  names(ag) <- c("A","R","G","N")
  names(at) <- c("A","W","T","N")
  names(cg) <- c("C","S","G","N")
  names(ct) <- c("C","Y","T","N")
  names(gt) <- c("G","K","T","N")
  conv <- list(ac,ac,ag,ag,at,at,cg,cg,ct,ct,gt,gt)
  names(conv) <- c("A/C","C/A","A/G","G/A","A/T","T/A","C/G","G/C",
                   "C/T","T/C","G/T","T/G")
  
  S <- length(samples)
  SBlist <- vector(mode="list",S)   
  for(i in 1:S){
    mygen <- mapply(function(type,gen) unname(conv[[type]][gen]),
                    type=hapmap[[1]], gen=hapmap[[i+1]],
                    SIMPLIFY=TRUE, USE.NAMES=FALSE)
    
    SBlist[[i]] <- new("SNPbin", mygen)
  }
  
  
  x <- new("genlight", SBlist)
  locNames(x) <- loci
  indNames(x) <- samples
  
  return(x)
}
SNP.numeric <- hapMap2genlight("10559_SNPs_1248_Indv.hmp.txt")
SNP.numeric <- apply(SNP.numeric,2,as.numeric)
rownames(SNP.numeric) <- colnames(SNP.hapmap[12:ncol(SNP.hapmap)])
dim(SNP.numeric) # 1248 inds x 10559 snps
SNP.numeric[1:5, 1:5]
SNP.hapmap[1:5, c(1,12:16)]

# Remove hete and missing sites in recurrent parent Tesale
table(SNP.numeric[14,2:nrow(SNP.numeric)]) # 0:3014, 1:409, 2:2978, NA:254
parent.hete <- apply(SNP.numeric[2:14,], 2, function(x) length(grep("1", x)))
table(parent.hete >= 1)
SNP.numeric2 <- SNP.numeric[,SNP.numeric[14,]=="0" | SNP.numeric[14,]=="2"]
SNP.numeric2 <- SNP.numeric2[,!is.na(SNP.numeric2[14,])] # 5992 SNPs
dim(SNP.numeric2) # 1246 inds x 5992 snps

# Subset SNP.hapmap to get 5992 SNPs in hapmap format
SNP.hapmap2 <- SNP.hapmap[c(SNP.hapmap[,1]%in%colnames(SNP.numeric2)),]

# Recurrent parent Teshale allele percentage:
# Create function to calculate recurrent parent allele frequency across inds
Teshale_allele_freq = function(x) {
  ref_freq=c()
  homo_ref_freq=c()
  for (i in 1:ncol(x)){
    homo_ref_sum = sum(x[,i]==SNP.numeric2[14,i], na.rm=T)
    hetero_sum = sum(x[,i]==1, na.rm=T)
    sample_size = length(na.omit(x[,i]))
    ref_freq[i] = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
    homo_ref_freq[i] = homo_ref_sum/sample_size
  }
  return(list(results1=ref_freq, results2=homo_ref_freq))
}

Teshale.freq.all <- Teshale_allele_freq(SNP.numeric2[-c(1:14),])$results1
Teshale.freq.IS10876 <- Teshale_allele_freq(SNP.numeric2[15:169,])$results1
Teshale.freq.IS15428 <- Teshale_allele_freq(SNP.numeric2[170:315,])$results1
Teshale.freq.IS14298 <- Teshale_allele_freq(SNP.numeric2[316:453,])$results1
Teshale.freq.IS14446 <- Teshale_allele_freq(SNP.numeric2[454:611,])$results1
Teshale.freq.IS14556 <- Teshale_allele_freq(SNP.numeric2[612:646,])$results1
Teshale.freq.IS16044 <- Teshale_allele_freq(SNP.numeric2[647:686,])$results1
Teshale.freq.IS16173 <- Teshale_allele_freq(SNP.numeric2[687:798,])$results1
Teshale.freq.IS2205 <- Teshale_allele_freq(SNP.numeric2[799:839,])$results1
Teshale.freq.IS22325 <- Teshale_allele_freq(SNP.numeric2[840:968,])$results1
Teshale.freq.IS23988 <- Teshale_allele_freq(SNP.numeric2[969:1016,])$results1
Teshale.freq.IS32234 <- Teshale_allele_freq(SNP.numeric2[1017:1041,])$results1
Teshale.freq.IS3583 <- Teshale_allele_freq(SNP.numeric2[1042:1161,])$results1
Teshale.freq.IS9911 <- Teshale_allele_freq(SNP.numeric2[1162:1246,])$results1

mean(Teshale.freq.all)
mean(Teshale.freq.IS10876)
mean(Teshale.freq.IS14298)
mean(Teshale.freq.IS14446)
mean(Teshale.freq.IS14556)
mean(Teshale.freq.IS15428)
mean(Teshale.freq.IS16044)
mean(Teshale.freq.IS16173)
mean(Teshale.freq.IS2205)
mean(Teshale.freq.IS22325)
mean(Teshale.freq.IS23988)
mean(Teshale.freq.IS32234)
mean(Teshale.freq.IS3583)
mean(Teshale.freq.IS9911)



# 03-05-2019
# Genome composition of each derived line. Theoretical: 25% donor + 75% recurrent.
# Create function to calculate recurrent parent allele frequency across inds
genome_comp = function(x) {
  ref_freq=c()
  homo_ref_freq=c()
  for (i in 1:nrow(x)) {
    homo_ref_sum=0
    hetero_sum=0
    for (j in 1:ncol(x)){
      if (!is.na(x[i,j]) & x[i,j]==SNP.numeric2[14,j]) 
        {homo_ref_sum = homo_ref_sum + 1}
      else if (!is.na(x[i,j]) & x[i,j]==1) 
        {hetero_sum = hetero_sum + 1}
    }
    sample_size = length(na.omit(x[i,]))
    ref_freq[i] = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
    homo_ref_freq[i] = homo_ref_sum/sample_size
  }
  return(list(results1=ref_freq, results2=homo_ref_freq))
}
# Genome composition of each BC1F4 line:
BC1F4.genome.comp <- genome_comp(SNP.numeric2)$results1[-c(1:14)]
BC1F4.genome <- data.frame("Recurrent"=BC1F4.genome.comp, 
                           "Population"=c(rep("IS10876", 155), rep("IS15428", 146),
                                          rep("IS14298", 138), rep("IS14446", 158),
                                          rep("IS14556", 35), rep("IS16044", 40),
                                          rep("IS16173", 112), rep("IS2205", 41),
                                          rep("IS22325", 129), rep("IS23988", 48),
                                          rep("IS32234", 25), rep("IS3583", 120),
                                          rep("IS9911", 85)))

# Basic box plot
p <- ggplot(BC1F4.genome, aes(y=Recurrent*100, x=Population)) + 
  geom_boxplot(notch = F) +
  theme_bw() +
  theme(
    #legend.position="none",
    #panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11)) +
  labs(y="% Recurrent parent", x="") 
  #geom_hline(yintercept = 75, lty=2, col="blue")
p + coord_flip()

tiff(filename = "Box plot of the distribution of the percentage of recurrent parent genome.tiff", width = 1600, height = 2000, units = "px", res = 300)
p + coord_flip()
dev.off()

# 04-18-2019. 
# Set Teshale genotypes were set to 0, diverse parental genotypes were set to 2:
SNP.numeric3 = SNP.numeric2
SNP.numeric3 = apply(SNP.numeric3, 2, function(x) if(x[14]==2) 
  {recode(x, '2'='0', '1'='1', '0'='2')} else
    {recode(x, '2'='2', '1'='1', '0'='0')})

SNP.numeric3 = apply(SNP.numeric3, 2, as.numeric)
rownames(SNP.numeric3) <- rownames(SNP.numeric2)

# SNP nummeric dataset by population:
SNP.numeric.IS10876 <- SNP.numeric3[c(1,14:169),]       # 155
SNP.numeric.IS15428 <- SNP.numeric3[c(2,14,170:315),]   # 146
SNP.numeric.IS14298 <- SNP.numeric3[c(3,14,316:453),]   # 138
SNP.numeric.IS14446 <- SNP.numeric3[c(4,14,454:611),]   # 158
SNP.numeric.IS14556 <- SNP.numeric3[c(5,14,612:646),]   # 35
SNP.numeric.IS16044 <- SNP.numeric3[c(6,14,647:686),]   # 40
SNP.numeric.IS16173 <- SNP.numeric3[c(7,14,687:798),]   # 112
SNP.numeric.IS2205 <- SNP.numeric3[c(8,14,799:839),]    # 41
SNP.numeric.IS22325 <- SNP.numeric3[c(9,14,840:968),]   # 129
SNP.numeric.IS23988 <- SNP.numeric3[c(10,14,969:1016),] # 48
SNP.numeric.IS32234 <- SNP.numeric3[c(11,14,1017:1041),]# 25
SNP.numeric.IS3583 <- SNP.numeric3[c(12,14,1042:1161),] # 120
SNP.numeric.IS9911 <- SNP.numeric3[c(13,14,1162:1246),] # 85

# 04-18-2019. Principal coordiante of each population
plot(cmdscale(dist(SNP.numeric3))) # All 13 populations
plot(cmdscale(dist(SNP.numeric.IS10876)))
plot(cmdscale(dist(SNP.numeric.IS14298)))
plot(cmdscale(dist(SNP.numeric.IS14446)))
plot(cmdscale(dist(SNP.numeric.IS14556)))
plot(cmdscale(dist(SNP.numeric.IS15428)))
plot(cmdscale(dist(SNP.numeric.IS16044)))
plot(cmdscale(dist(SNP.numeric.IS16173)))
plot(cmdscale(dist(SNP.numeric.IS2205)))
plot(cmdscale(dist(SNP.numeric.IS22325)))
plot(cmdscale(dist(SNP.numeric.IS23988)))
plot(cmdscale(dist(SNP.numeric.IS32234)))
plot(cmdscale(dist(SNP.numeric.IS3583)))
plot(cmdscale(dist(SNP.numeric.IS9911)))

# 04-19-2019, prepare joint-linkage mapping dataset
# monomorphic SNPs within each family were set as missing, 
# and missing data were imputed as the mean of the nearest flanking markers weighted by physical distance
# Find mono-morphic SNPs among progeny:
Mono.rate.progeny <- apply(SNP.numeric.IS10876[-c(1,2),], 2, 
                           function(x) length(unique(na.omit(x))))

SNP.numeric.IS10876[-c(1,2), Mono.rate.progeny==1] <- NA # Convert mono-morphic markers to NA



#table(Teshale.freq.IS10876 < 1 & Teshale.freq.IS10876 > 0)[1]/5992
#table(Teshale.freq.IS14298 < 1 & Teshale.freq.IS14298 > 0)[1]/5992
#table(Teshale.freq.IS14446 < 1 & Teshale.freq.IS14446 > 0)[1]/5992
#table(Teshale.freq.IS14556 < 1 & Teshale.freq.IS14556 > 0)[1]/5992
#table(Teshale.freq.IS15428 < 1 & Teshale.freq.IS15428 > 0)[1]/5992
#table(Teshale.freq.IS16044 < 1 & Teshale.freq.IS16044 > 0)[1]/5992
#table(Teshale.freq.IS16173 < 1 & Teshale.freq.IS16173 > 0)[1]/5992
#table(Teshale.freq.IS2205 < 1 & Teshale.freq.IS2205 > 0)[1]/5992
#table(Teshale.freq.IS22325 < 1 & Teshale.freq.IS22325 > 0)[1]/5992
#table(Teshale.freq.IS23988 < 1 & Teshale.freq.IS23988 > 0)[1]/5992
#table(Teshale.freq.IS32234 < 1 & Teshale.freq.IS32234 > 0)[1]/5992
#table(Teshale.freq.IS3583 < 1 & Teshale.freq.IS3583 > 0)[1]/5992
#table(Teshale.freq.IS9911 < 1 & Teshale.freq.IS9911 > 0)[1]/5992


############################################################################################
# Create my own function to calculate alternate (alt) allele freq
allele_freq = function(x) {
  ref_freq=c() 
  alt_freq=c()
  homo_ref_freq=c()
  homo_alt_freq=c()
  for (i in 1:ncol(x)){
    homo_ref_sum = sum(x[,i]==0, na.rm=T)
    hetero_sum = sum(x[,i]==1, na.rm=T)
    homo_alt_sum = sum(x[,i]==2, na.rm=T)
    sample_size = sum(homo_ref_sum, hetero_sum, homo_alt_sum)
    ref_freq[i] = ((homo_ref_sum * 2) + hetero_sum)/(sample_size * 2)
    alt_freq[i] = ((homo_alt_sum * 2) + hetero_sum)/(sample_size * 2) 
    homo_ref_freq[i] = homo_ref_sum/sample_size
    homo_alt_freq[i] = homo_alt_sum/sample_size
  }
  return(list(results1=ref_freq, results2=alt_freq, results3=homo_ref_freq, results4=homo_alt_freq)) 
}

#hist(allele_freq(SNP.Num2[-c(34,39),])$results1, breaks=30, main="Ref allele freq")#ref allele ~ 0.5
#hist(allele_freq(SNP.Num2[-c(34,39),])$results2, breaks=30, main="Alt allele freq")#alt allele ~ 0.5
hist(allele_freq(SNP.numeric.IS10876[-c(1,2),])$results3, breaks=30, main="Homo_Ref allele freq")#ref allele ~ 0.5
hist(allele_freq(SNP.numeric.IS10876[-c(1,2),])$results4, breaks=30, main="Homo_Alt allele freq")#alt allele ~ 0.5


# Mono-morphic rate between parents:
Mono.rate = matrix(NA, nrow = 13, ncol = 5992)
for (i in 1:13) {
  for (j in 1:5992) {
    Mono.rate[i,j] = ifelse(length(unique(na.omit(unlist(SNP.numeric2[c(i,14),j])))) == 2, 
                            "Y", "N")
  }
}

rownames(Mono.rate) <- c(rownames(SNP.numeric2)[1:13])
colnames(Mono.rate) <- c(colnames(SNP.numeric2))

# Between Teshale and 13 parents:
table(Mono.rate[1,])[2]/5992
table(Mono.rate[2,])[2]/5992
table(Mono.rate[3,])[2]/5992
table(Mono.rate[4,])[2]/5992
table(Mono.rate[5,])[2]/5992
table(Mono.rate[6,])[2]/5992
table(Mono.rate[7,])[2]/5992
table(Mono.rate[8,])[2]/5992
table(Mono.rate[9,])[2]/5992
table(Mono.rate[10,])[2]/5992
table(Mono.rate[11,])[2]/5992
table(Mono.rate[12,])[2]/5992
table(Mono.rate[13,])[2]/5992


# missing data by ind
miss.by.ind <- apply(SNP.numeric2, 1, function(x) sum(is.na(x))/ncol(SNP.numeric2))
plot(miss.by.ind)
which(miss.by.ind > 0.8) # Previously removed Teshale_IS2205_1314 Teshale_IS9911_1131

miss.by.snp <- apply(SNP.numeric2, 2, function(x) sum(is.na(x))/nrow(SNP.numeric2))
plot(miss.by.snp)


# Check heterozygosity rate for taxa:
taxa.hete.rate <- apply(SNP.numeric2, 1, 
                        function(x) sum(x==1, na.rm = TRUE)/ncol(SNP.numeric2))
mean(taxa.hete.rate) # 0.085
range(taxa.hete.rate) # 0.03 - 0.41
plot(taxa.hete.rate, xlab = "Genotype", ylab = "Heterozygous rate", ylim = c(0,1), type = "l")
mean(taxa.hete.rate)

# Check heterozygosity rate for SNP:
snp.hete.rate <- apply(SNP.numeric2, 2, 
                       function(x) sum(x==1, na.rm = TRUE)/nrow(SNP.numeric2))
mean(snp.hete.rate) # 0.069
range(snp.hete.rate) # 0.0 - 0.97
table(snp.hete.rate <= 0.5) # True: 6456 SNPs, False: 199 SNPs
plot(snp.hete.rate, xlab = "SNPs", ylab = "Heterozygous rate", ylim = c(0,1), type = "l")

# Heterozygous rate of parents:
taxa.hete.rate[1:14]
write.csv(taxa.hete.rate, "Heterozygous rate over 1246 indvs.csv")
write.csv(snp.hete.rate, "Heterozygous rate over 6655 SNPs.csv")

## Overview of SNPs across sorghum genome
library(CMplot)
CMplot(SNP.hmp3[,c(1,3,4)],plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300)




######################
## Hetero rate across genome:
genome.wide.rate <- data.frame("SNP" = SNP.hapmap[c(SNP.hapmap[,1]%in%colnames(SNP.numeric2)),1],
                               "CHR" = SNP.hapmap[c(SNP.hapmap[,1]%in%colnames(SNP.numeric2)),3],
                               "POS" = SNP.hapmap[c(SNP.hapmap[,1]%in%colnames(SNP.numeric2)),4],
                               "HET" = snp.hete.rate,
                               "Teshale_rate_all"=Teshale.freq.all,
                               "Teshale_rate_IS10876"=Teshale.freq.IS10876,
                               "Teshale_rate_IS14298"=Teshale.freq.IS14298,
                               "Teshale_rate_IS14446"=Teshale.freq.IS14446,
                               "Teshale_rate_IS14556"=Teshale.freq.IS14556,
                               "Teshale_rate_IS15428"=Teshale.freq.IS15428,
                               "Teshale_rate_IS16044"=Teshale.freq.IS16044,
                               "Teshale_rate_IS16173"=Teshale.freq.IS16173,
                               "Teshale_rate_IS2205"=Teshale.freq.IS2205,
                               "Teshale_rate_IS22325"=Teshale.freq.IS22325,
                               "Teshale_rate_IS23988"=Teshale.freq.IS23988,
                               "Teshale_rate_IS32234"=Teshale.freq.IS32234,
                               "Teshale_rate_IS3583"=Teshale.freq.IS3583,
                               "Teshale_rate_IS9911"=Teshale.freq.IS9911)



# 02-28-2019
# Sliding window function
sliding.window.rate <- function(x, window_size, walk_speed, trait.col) {
  my.window.rate = c()
  my.window.pos = c()
  my.window.chr = c()
  #Loop through 10 chromosomes:
  for (i in 1:10) {
    print(i)
    #Determine window edges:
    window_edge=seq(window_size,max(x[x[,2]==i,3]),by=walk_speed)
    #Add starting and ending genomic positions:
    window_edge <- c(0, window_edge, max(x[x[,2]==i,3]))
    
    for (j in 1:(length(window_edge)-1)){
      print(j)
      
      temp_window=x[x[x[,2]==i,3] >= window_edge[j] & x[x[,2]==i,3] < window_edge[j+1], trait.col]
      #the hete.rate function
      temp.rate=mean(temp_window)
      temp.pos=window_edge[j]+500000
      #append values to vectors:
      my.window.rate = c(my.window.rate, temp.rate)
      my.window.chr = c(my.window.chr, i)
      my.window.pos = c(my.window.pos, temp.pos)
    }
  }
  return(list(result1=my.window.chr, result2=my.window.pos, result3=my.window.rate))
}




# For recurrent parent allelic freq:
temp.data <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 5)
temp.data.1 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 6)
temp.data.2 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 7)
temp.data.3 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 8)
temp.data.4 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 9)
temp.data.5 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 10)
temp.data.6 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 11)
temp.data.7 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 12)
temp.data.8 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 13)
temp.data.9 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 14)
temp.data.10 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 15)
temp.data.11 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 16)
temp.data.12 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 17)
temp.data.13 <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 18)

nam.rate <- data.frame("CHR"=temp.data$result1, "POS"=temp.data$result2, 
                       "Rate.all"=temp.data$result3,
                       "Rate.1"=temp.data.1$result3,
                       "Rate.2"=temp.data.2$result3,
                       "Rate.3"=temp.data.3$result3,
                       "Rate.4"=temp.data.4$result3,
                       "Rate.5"=temp.data.5$result3,
                       "Rate.6"=temp.data.6$result3,
                       "Rate.7"=temp.data.7$result3,
                       "Rate.8"=temp.data.8$result3,
                       "Rate.9"=temp.data.9$result3,
                       "Rate.10"=temp.data.10$result3,
                       "Rate.11"=temp.data.11$result3,
                       "Rate.12"=temp.data.12$result3,
                       "Rate.13"=temp.data.13$result3)

don1 <- nam.rate %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(nam.rate, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, POS) %>%
  mutate(BPcum=POS+tot)

axisdf = don1 %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2)

vline <- don1 %>%
  select(CHR, BPcum) %>%
  group_by(CHR) %>%
  summarise(vline.max=max(BPcum))

p1 <- ggplot(don1, aes(x=BPcum)) +

  geom_point(aes(y=Rate.1, color="red"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.2, color="blue"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.3, color="yellow"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.4, color="green"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.5, color="purple"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.6, color="gray"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.7, color="black"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.8, color="darkgreen"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.9, color="orange"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.10, color="magenta"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.11, color="cyan"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.12, color="skyblue"), pch=16, cex=0.6) +
  geom_point(aes(y=Rate.13, color="yellowgreen"), pch=16, cex=0.6) +
  
  #geom_smooth(aes(x=BPcum, y=Rate.all)) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
  #scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
  labs(x="Chromosome", y="Teshale allele rate") +
  ylim(0, 1) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2) +
  geom_hline(yintercept = c(0.6,0.9), lty=2) +
  # Custom the theme:
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11)
  )


# For heterozygous rate:
temp.data <- sliding.window.rate(genome.wide.rate, 1000000, 500000, 4)
nam.rate <- data.frame("CHR"=temp.data$result1, "POS"=temp.data$result2, "Rate"=temp.data$result3)

don2 <- nam.rate %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(nam.rate, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, POS) %>%
  mutate(BPcum=POS+tot)

p2 <- ggplot(don2, aes(x=BPcum, y=Rate)) +
  
  # Show all points
  geom_col() +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
  
  labs(x="", y="Hetero rate") +
  ylim(0, 1) +
  geom_vline(vline[-10,], mapping = aes(xintercept = vline.max), lty=2) +
  # Custom the theme:
  theme_bw() +
  theme(
    legend.position="none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 11)
  )

ggarrange(p2, p1, ncol=1, nrow=2)

#tiff(filename = "Genome-wide heterozygous rate and recurrent parent allele freq.tiff", width = 647, height = 408, units = "px", res = 400)
#ggarrange(p2, p1, ncol=1, nrow=2)
#dev.off()


# R function to import data from TASSEL's HapMap format to adegenet's
# genind format.
# Previous version available at: http://dx.doi.org/10.13012/C5CC0XMJ
# Example use:
# mydata <- hapMap2genind("HapMap.hmp.txt")
hapMap2genind <- function(file){
  require(adegenet)
  hapmap <- read.table(file, header=TRUE, row.names=1, sep="\t",
                       stringsAsFactors=FALSE)
  samples <- scan(file, what = character(), nlines = 1)[-(1:11)]
  conv <- c("AA", "CC", "TT", "GG", "AC", "AT", "AG", "CT", "CG", "TG", NA)
  names(conv) <- c("A", "C", "T", "G", "M" , "W", "R", "Y", "S",  "K", "N")
  mydf <- matrix(NA, nrow=dim(hapmap)[1], ncol=length(samples),
                 dimnames=list(row.names(hapmap), samples))
  for(i in 1:length(samples)){
    mydf[,i] <- conv[hapmap[[i+10]]]
  }
  mydf <- as.data.frame(t(mydf))
  x <- df2genind(mydf, type="codom", ncode=1)
  return(x)
}

SNP.genind <- hapMap2genind("6655_SNPs_1246_Indv.hmp.txt")

# add pop info
SNP.genind$pop <- as.factor(c(rep("parent", 14), PCs[,7]))
pop(SNP.genind)
#Use adgenet package calculate individual inbreeding coefficient:
inb.cof <- inbreeding(SNP.genind)


#No. of SNPs monomorphic in parents
#parent.mono <- apply(SNP.hapmap[,12:25], 1, function(x) length(unique(unlist(x))))
#table(parent.mono > 1)

# Remove mono-morphic SNPs
#SNP.hmp2 <- SNP.hapmap[parent.mono > 1,]
#SNP.hmp2 <- SNP.hmp2[c(SNP.hmp2[,2]=="A/C" | SNP.hmp2[,2]=="A/G" | SNP.hmp2[,2]=="A/T" | 
#                         SNP.hmp2[,2]=="C/A" | SNP.hmp2[,2]=="C/G" | SNP.hmp2[,2]=="C/T" |
#                         SNP.hmp2[,2]=="G/A" | SNP.hmp2[,2]=="G/C" | SNP.hmp2[,2]=="G/T" |
#                         SNP.hmp2[,2]=="T/A" | SNP.hmp2[,2]=="T/C" | SNP.hmp2[,2]=="T/G"),]

# Remove complete LD snps (adjacent physical distance < 64 bp)
#snp.dist = c()
#for (i in 2:nrow(SNP.hmp2)){
#  snp.dist[i] = SNP.hmp2[i,4] - SNP.hmp2[i-1,4]
#}
#snp.dist[1] <- 174473
#SNP.hmp3 <- SNP.hmp2[snp.dist > 64 | snp.dist < 0,] # Retain 6659 SNPs

#SNP.numeric2 <- SNP.numeric[,colnames(SNP.numeric)%in%SNP.hmp3[,1]]



# Look at LD decay
LD_r2 <- read.table("LD_r2.txt", header = TRUE)

# Take a random sample
sample.list <- sample(c(1:1885138), 100000)
plot(LD_r2[sample.list,3]-LD_r2[sample.list,2], LD_r2[sample.list,5], pch=16, cex=0.4, xlab = "Physical distance (bp)", ylab = "R2")
par(mar=c(5,5,1,1))
# Hill & Weir (1988) method:
# this data was obtained from genotyping of 1246 individuals, so n = no. of gametes = no. of sampled chromosomes = 2 x 1246
n = 1232*2
Cstart <- c(C=0.1)
rsq <- LD_r2[sample.list,5]
dist <- LD_r2[sample.list,3]-LD_r2[sample.list,2]

# let's fit a non linear model using the arbitrary C value
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter, 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ((10+rho*dist)/((2+rho*dist)*(11+rho*dist)))*(1+((3+rho*dist)*(12+12*rho*dist+(rho*dist)^2))/(n*(2+rho*dist)*(11+rho*dist)))

newfile <- data.frame(dist, newrsq)

maxld <- max(newfile$newrsq) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$dist),]

# create the plot
plot(dist/1000000, rsq, pch=16, cex=0.2, 
     xlab="Distance (Mb)", ylab=expression("LD ("~italic("r")^"2"~")"))
lines(newfile$dist/1000000, newfile$newrsq, col="blue", lwd=2)
abline(h=0.2, col="red")
abline(v=halfdecaydist/1000000, col="green")
mtext(round(halfdecaydist/1000000,2), side=1, line=0.05, at=halfdecaydist/1000000, cex=0.75, col="green")



# Below code plotting all data points, too large.
pdf(width=6,height=4) 

# Hill & Weir (1988) method:
# this data was obtained from genotyping of 1232 individuals, so n = no. of gametes = no. of sampled chromosomes = 2 x 1246
n = 1232*2
Cstart <- c(C=0.1)
rsq <- LD_r2[,5]
dist <- LD_r2[,3]-LD_r2[,2]

# let's fit a non linear model using the arbitrary C value
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter, 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ((10+rho*dist)/((2+rho*dist)*(11+rho*dist)))*(1+((3+rho*dist)*(12+12*rho*dist+(rho*dist)^2))/(n*(2+rho*dist)*(11+rho*dist)))

newfile <- data.frame(dist, newrsq)

#maxld <- max(file$rsq) #using max LD value from initial input file
maxld <- max(newfile$newrsq) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$dist),]

# create the plot
plot(dist/1000000, rsq, pch=16, cex=0.1, xlab="Distance (Mb)", ylab=expression("LD ("~italic("r")^"2"~")"))
lines(newfile$dist/1000000, newfile$newrsq, col="blue", lwd=2)
abline(h=0.2, col="red")
abline(v=halfdecaydist/1000000, col="green")
mtext(round(halfdecaydist/1000000,2), side=1, line=0.05, at=halfdecaydist/1000000, cex=0.75, col="green")
dev.off()



######################
## PCA of 1232 progeny:
PCs <- read.csv("PCA of 1234 progeny x 10559 SNPs.csv", header = TRUE)


p1 <- ggplot(PCs, aes(x=PC1, y=PC2)) + geom_point(aes(color=Population)) + labs(x="PC1 (7.4%)", y="PC2 (5.4%)") + theme_bw(base_size = 12) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom", legend.text=element_text(size=11), axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) + scale_color_manual(values = c("red","blue","yellow","green","purple","gray","black","darkgreen","orange","magenta","cyan","skyblue","yellowgreen"))
p2 <- ggplot(PCs, aes(x=PC2, y=PC3)) + geom_point(aes(color=Population)) + labs(x="PC2 (5.4%)", y="PC3 (3.7%)") + theme_bw(base_size = 12) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom", legend.text=element_text(size=11), axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) + scale_color_manual(values = c("red","blue","yellow","green","purple","gray","black","darkgreen","orange","magenta","cyan","skyblue","yellowgreen"))
p3 <- ggplot(PCs, aes(x=PC3, y=PC4)) + geom_point(aes(color=Population)) + labs(x="PC3 (3.7%)", y="PC4 (3.3%)") + theme_bw(base_size = 12) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position="bottom", legend.text=element_text(size=11), axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11)) + scale_color_manual(values = c("red","blue","yellow","green","purple","gray","black","darkgreen","orange","magenta","cyan","skyblue","yellowgreen"))
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

tiff(filename = "PCA plots2.tiff", width = 480, height = 400, units = "px", res = 400)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()


