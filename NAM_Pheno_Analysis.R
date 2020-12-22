setwd("~/Dropbox/PGML_Projects/Raw Pheno")

## load raw data
Kobo.raw = read.csv("Raw Data Kobo.csv", header = TRUE)
Meiso.raw = read.csv("Raw Data Meiso.csv", header = TRUE)
Sheraro.raw = read.csv("Raw Data Sheraro.csv", header = TRUE)

# remove NA ID rows
# DNA of some BCNAM lines were no collected due to missing leaf samples, those lines were eliminated from analysis
Kobo = Kobo.raw[!is.na(Kobo.raw$ID),]
Meiso = Meiso.raw[!is.na(Meiso.raw$ID),]
Sheraro = Sheraro.raw[!is.na(Sheraro.raw$ID),]

# remove raw dataset to save space
rm(Kobo.raw)
rm(Meiso.raw)
rm(Sheraro.raw)

length(unique(Kobo[!Kobo$Population=="Parent",]$ID)) # 1077 BCNAM genotypes + 12 parents
length(unique(Meiso[!Meiso$Population=="Parent",]$ID)) # 1212 BCNAM genotypes + 12 parents
length(unique(Sheraro$ID)) # 1074 BCNAM genotypes, parents were not included in Sheraro
# ----------------------------------------------------------------------------

# Look number of reps per env:
table(table(Kobo[!Kobo$Population=="Parent",]$ID)) # 362 Taxa with 1 rep, 715 Taxa with 2 reps
table(table(Meiso[!Meiso$Population=="Parent",]$ID)) # 195 Taxa with 1 rep, 1017 Taxa with 2 reps
table(table(Sheraro[!Sheraro$Population=="Parent",]$ID)) # 2 Taxa with 1 rep, 1072 Taxa with 2 reps

####################################################################################
####################################################################################
library(lme4)
library(MASS)
library(rcompanion)
library(ggplot2)
library(plyr)
library(grid)
library(ggpubr)
library(reshape2)
####################################################################################
####################################################################################
lm <- lmer(GYP~1+(1|Population) + (1|ID), 
           data = Kobo[!Kobo$Population=="IS32234" & !Kobo$Population=="Parent",])
summary(lm)

# BLUPs in Kobo:
# Days to 50% flowering
Replication <- as.factor(Kobo$Replication)
Block <- as.factor(Kobo$Block)
Kobo.DF.lm <- lmer(DF ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.DF.lm))
qqline(resid(Kobo.DF.lm), col= "red") # add a perfect fit line
summary(Kobo.DF.lm)
Kobo.DF.BLUPs <- ranef(Kobo.DF.lm) # Extract BLUPs of random terms
hist(Kobo.DF.BLUPs$ID[,1]+78.163, xlab = "Days to 50% flowering", main = "")

# Days to maturity
Kobo.DM.lm <- lmer(DM ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.DM.lm))
qqline(resid(Kobo.DM.lm), col= "red") # add a perfect fit line
summary(Kobo.DM.lm)
Kobo.DM.BLUPs <- ranef(Kobo.DM.lm) # Extract BLUPs of random terms
hist(Kobo.DM.BLUPs$ID[,1]+127.04, xlab = "Days to maturity", main = "")

# Plant height
Kobo.PH.lm <- lmer(PH ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.PH.lm))
qqline(resid(Kobo.PH.lm), col= "red") # add a perfect fit line
summary(Kobo.PH.lm)
Kobo.PH.BLUPs <- ranef(Kobo.PH.lm) # Extract BLUPs of random terms
hist(Kobo.PH.BLUPs$ID[,1]+175.75, xlab = "Plant height", main = "")

# Leaf senescence
Kobo.LS.lm <- lmer(LS ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.LS.lm))
qqline(resid(Kobo.LS.lm), col= "red") # add a perfect fit line
summary(Kobo.LS.lm)
Kobo.LS.BLUPs <- ranef(Kobo.LS.lm) # Extract BLUPs of random terms
hist(Kobo.LS.BLUPs$ID[,1]+2.4891, xlab = "Leaf senescence", main = "")

# Number of panicle per plant
Kobo.NPP.lm <- lmer(NPP ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.NPP.lm)) 
qqline(resid(Kobo.NPP.lm), col= "red") # add a perfect fit line
bc <- boxcox(NPP ~ 1, lambda = seq(-10,10,0.1), data = Kobo) # box-cox transformation
lambda <- bc$x[which.max(bc$y)]
plotNormalHistogram(Kobo$NPP)
plotNormalHistogram(((Kobo$NPP)^lambda-1)/lambda)
# Redo lmer with transformed data:
Kobo.NPP.lm <- lmer(((Kobo$NPP)^lambda-1)/lambda ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
summary(Kobo.NPP.lm)
Kobo.NPP.BLUPs <- ranef(Kobo.NPP.lm) # Extract BLUPs of random terms
hist(Kobo.NPP.BLUPs$ID[,1]+1.41812, xlab = "Number of panicles per plant", main = "")

# Number of leaves
Kobo.NL.lm <- lmer(NL ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.NL.lm)) 
qqline(resid(Kobo.NL.lm), col= "red") # add a perfect fit line
summary(Kobo.NL.lm)
Kobo.NL.BLUPs <- ranef(Kobo.NL.lm) # Extract BLUPs of random terms
hist(Kobo.NL.BLUPs$ID[,1]+8.7612, xlab = "Number of leaves", main = "")
range(Kobo.NL.BLUPs$ID[,1]+8.7612)

# Head exsertion, left-skewed distribution
Kobo.HE.lm <- lmer(HE ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.HE.lm)) 
qqline(resid(Kobo.HE.lm), col= "red") # add a perfect fit line
#bc <- boxcox(HE+1 ~ 1, lambda = seq(-10,10,0.1), data = Kobo) # box-cox transformation
#lambda <- bc$x[which.max(bc$y)]
#plotNormalHistogram(Kobo$HE)
#plotNormalHistogram(((Kobo$HE+1)^lambda-1)/lambda)
# Redo lmer with transformed data:
#Kobo.HE.lm <- lmer(((Kobo$HE+1)^lambda-1)/lambda ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
summary(Kobo.HE.lm)
Kobo.HE.BLUPs <- ranef(Kobo.HE.lm) # Extract BLUPs of random terms
hist(Kobo.HE.BLUPs$ID[,1]+3.8879, xlab = "Head exsertion", main = "")

# Number of tillers, left-skewed distribution
Kobo.NT.lm <- lmer(NT ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.NT.lm)) 
qqline(resid(Kobo.NT.lm), col= "red") # add a perfect fit line
#bc <- boxcox(NT ~ 1, lambda = seq(-10,10,0.1), data = Kobo) # box-cox transformation
#lambda <- bc$x[which.max(bc$y)]
#plotNormalHistogram(Kobo$NT)
#plotNormalHistogram((Kobo$NT^lambda-1)/lambda)
# Redo lmer with transformed data:
#Kobo.NT.lm <- lmer((Kobo$NT^lambda-1)/lambda ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
summary(Kobo.NT.lm)
Kobo.NT.BLUPs <- ranef(Kobo.NT.lm) # Extract BLUPs of random terms
hist(Kobo.NT.BLUPs$ID[,1]+2.5066, xlab = "Number of tillers", main = "")

# Grain yield per primary panicle
Kobo.GYP.lm <- lmer(GYP ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.GYP.lm)) 
qqline(resid(Kobo.GYP.lm), col= "red") # add a perfect fit line
summary(Kobo.GYP.lm)
Kobo.GYP.BLUPs <- ranef(Kobo.GYP.lm) # Extract BLUPs of random terms
hist(Kobo.GYP.BLUPs$ID[,1]+62.753, xlab = "Grain yield per primary panicle", main = "")
range(Kobo.GYP.BLUPs$ID[,1]+62.753)

# Plant per plot
Kobo.PPP.lm <- lmer(plant.per.plot ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.PPP.lm)) 
qqline(resid(Kobo.PPP.lm), col= "red") # add a perfect fit line
summary(Kobo.PPP.lm)
Kobo.PPP.BLUPs <- ranef(Kobo.PPP.lm) # Extract BLUPs of random terms
hist(Kobo.PPP.BLUPs$ID[,1]+8.45678, xlab = "Avg. plant per plot", main = "")
range(Kobo.PPP.BLUPs$ID[,1]+8.45678)

# Grain yield per primary panicle, adjusted for plant per plot
Kobo.GYP.adj.lm <- lmer(GYP ~ 1 + plant.per.plot + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.GYP.adj.lm)) 
qqline(resid(Kobo.GYP.adj.lm), col= "red") # add a perfect fit line
summary(Kobo.GYP.adj.lm)
Kobo.GYP.adj.BLUPs <- ranef(Kobo.GYP.adj.lm) # Extract BLUPs of random terms
hist(Kobo.GYP.adj.BLUPs$ID[,1]+40.1797, xlab = "Grain yield per primary panicle", main = "")
range(Kobo.GYP.adj.BLUPs$ID[,1]+40.1797)

cor.test(Kobo.GYP.BLUPs$ID[,1], Kobo.GYP.adj.BLUPs$ID[,1])

# Thousand seed weight
Kobo.TSW.lm <- lmer(TSW ~ 1 + (1|ID) + (1|Replication/Block), data = Kobo)
qqnorm(resid(Kobo.TSW.lm)) 
qqline(resid(Kobo.TSW.lm), col= "red") # add a perfect fit line
summary(Kobo.TSW.lm)
Kobo.TSW.BLUPs <- ranef(Kobo.TSW.lm) # Extract BLUPs of random terms
hist(Kobo.TSW.BLUPs$ID[,1]+29.4754, xlab = "Thousand seed weight", main = "")

Trait.BLUPs.Kobo <- cbind(Kobo.DF.BLUPs$ID, Kobo.DM.BLUPs$ID, Kobo.HE.BLUPs$ID, Kobo.LS.BLUPs$ID,
                          Kobo.NL.BLUPs$ID, Kobo.NPP.BLUPs$ID, Kobo.NT.BLUPs$ID, Kobo.PH.BLUPs$ID,
                          Kobo.GYP.BLUPs$ID, Kobo.GYP.adj.BLUPs$ID, Kobo.TSW.BLUPs$ID, Kobo.PPP.BLUPs$ID)
colnames(Trait.BLUPs.Kobo) <- c("DF", "DM", "HE", "LS", "NL", "NPP", "NT", "PH", "GYP", "GYP.adj", "TSW", "PPP")

Trait.BLUPs.Kobo$DF <- Trait.BLUPs.Kobo$DF + 78.145
Trait.BLUPs.Kobo$DM <- Trait.BLUPs.Kobo$DM + 127.04
Trait.BLUPs.Kobo$HE <- Trait.BLUPs.Kobo$HE + 1.33985
Trait.BLUPs.Kobo$LS <- Trait.BLUPs.Kobo$LS + 2.4891
Trait.BLUPs.Kobo$NL <- Trait.BLUPs.Kobo$NL + 8.7612
Trait.BLUPs.Kobo$NPP <- Trait.BLUPs.Kobo$NPP + 0.1693
Trait.BLUPs.Kobo$NT <- Trait.BLUPs.Kobo$NT + 0.66576
Trait.BLUPs.Kobo$PH <- Trait.BLUPs.Kobo$PH + 175.75
Trait.BLUPs.Kobo$GYP <- Trait.BLUPs.Kobo$GYP + 62.753
Trait.BLUPs.Kobo$GYP.adj <- Trait.BLUPs.Kobo$GYP.adj + 40.1797
Trait.BLUPs.Kobo$TSW <- Trait.BLUPs.Kobo$TSW + 29.4754
Trait.BLUPs.Kobo$PPP <- Trait.BLUPs.Kobo$PPP + 8.45678

library(psych)
pairs.panels(Trait.BLUPs.Kobo, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             cex.cor = 1.5)

##############################################################################
# Meiso
# Days to 50% flowering
Replication <- as.factor(Meiso$Replication)
Block <- as.factor(Meiso$Block)
Meiso.DF.lm <- lmer(DF ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.DF.lm)) 
qqline(resid(Meiso.DF.lm), col= "red") # add a perfect fit line
summary(Meiso.DF.lm)
Meiso.DF.BLUPs <- ranef(Meiso.DF.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.DF.BLUPs$ID[,1]+85.03227, xlab = "Days to 50% flowering", main = "")


# Days to maturity
Meiso.DM.lm <- lmer(DM ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.DM.lm)) 
qqline(resid(Meiso.DM.lm), col= "red") # add a perfect fit line
summary(Meiso.DM.lm)
Meiso.DM.BLUPs <- ranef(Meiso.DM.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.DM.BLUPs$ID[,1]+122.6602, xlab = "Days to maturity", main = "")

# Plant height
Meiso.PH.lm <- lmer(PH ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.PH.lm)) 
qqline(resid(Meiso.PH.lm), col= "red") # add a perfect fit line
summary(Meiso.PH.lm)
Meiso.PH.BLUPs <- ranef(Meiso.PH.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.PH.BLUPs$ID[,1]+177.797, xlab = "Plant height", main = "")

# Leaf senescence
Meiso.LS.lm <- lmer(LS ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.LS.lm)) 
qqline(resid(Meiso.LS.lm), col= "red") # add a perfect fit line
summary(Meiso.LS.lm)
Meiso.LS.BLUPs <- ranef(Meiso.LS.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.LS.BLUPs$ID[,1]+2.8095, xlab = "Leaf senescence", main = "")

# Number of panicle per plant, no variation!
Meiso.NPP.lm <- lmer(NPP ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.NPP.lm)) 
qqline(resid(Meiso.NPP.lm), col= "red") # add a perfect fit line
#bc <- boxcox(NPP+1 ~ 1, lambda = seq(-10,10,0.1), data = Meiso) # box-cox transformation
#lambda <- bc$x[which.max(bc$y)]
#plotNormalHistogram(Meiso$NPP+1)
#plotNormalHistogram(((Meiso$NPP+1)^lambda-1)/lambda)
# Redo lmer with transformed data:
#Meiso.NPP.lm <- lmer(((Meiso$NPP+1)^lambda-1)/lambda ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
summary(Meiso.NPP.lm)
Meiso.NPP.BLUPs <- ranef(Meiso.NPP.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.NPP.BLUPs$ID[,1]+1.3153, xlab = "Number of panicles per plant", main = "")

# Number of leaves, need correction for some outliers.
Meiso.NL.lm <- lmer(NL ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.NL.lm)) 
qqline(resid(Meiso.NL.lm), col= "red") # add a perfect fit line
summary(Meiso.NL.lm)
Meiso.NL.BLUPs <- ranef(Meiso.NL.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.NL.BLUPs$ID[,1]+9.9454, xlab = "Number of leaves", main = "")

# Head exsertion, left-skewed distribution
Meiso.HE.lm <- lmer(HE ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.HE.lm)) 
qqline(resid(Meiso.HE.lm), col= "red") # add a perfect fit line
bc <- boxcox(HE+1 ~ 1, lambda = seq(-10,10,0.1), data = Meiso) # box-cox transformation
lambda <- bc$x[which.max(bc$y)] # 0.4
plotNormalHistogram(Meiso$HE, xlab="Original head exsertion")
plotNormalHistogram(((Meiso$HE+1)^lambda-1)/lambda, xlab="Box-Cox transformed head exsertion")
# Redo lmer with transformed data:
Meiso.HE.lm <- lmer(((Meiso$HE+1)^lambda-1)/lambda ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
summary(Meiso.HE.lm)
Meiso.HE.BLUPs <- ranef(Meiso.HE.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.HE.BLUPs$ID[,1]+2.9275, xlab = "Head exsertion", main = "")

par(mfrow=c(1,2))
plotNormalHistogram(Meiso$HE, xlab="Original head exsertion")
plotNormalHistogram(((Meiso$HE+1)^lambda-1)/lambda, xlab="Box-Cox transformed head exsertion")


# Number of tillers, left-skewed distribution
Meiso.NT.lm <- lmer(NT ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.NT.lm)) 
qqline(resid(Meiso.NT.lm), col= "red") # add a perfect fit line
bc <- boxcox(NT ~ 1, lambda = seq(-10,10,0.1), data = Meiso) # box-cox transformation
lambda <- bc$x[which.max(bc$y)] # -1.5
plotNormalHistogram(Meiso$NT)
plotNormalHistogram((Meiso$NT^lambda-1)/lambda)
# Redo lmer with transformed data:
Meiso.NT.lm <- lmer((Meiso$NT^lambda-1)/lambda ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
summary(Meiso.NT.lm)
Meiso.NT.BLUPs <- ranef(Meiso.NT.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.NT.BLUPs$ID[,1]+0.21, xlab = "Number of tillers", main = "")

# Grain yield per primary panicle
Meiso.GYP.lm <- lmer(GYP ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.GYP.lm)) 
qqline(resid(Meiso.GYP.lm), col= "red") # add a perfect fit line
summary(Meiso.GYP.lm)
Meiso.GYP.BLUPs <- ranef(Meiso.GYP.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.GYP.BLUPs$ID[,1]+54.5603, xlab = "Grain yield per primary panicle", main = "")

# Plant per plot
Meiso.PPP.lm <- lmer(plant.per.plot ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.PPP.lm)) 
qqline(resid(Meiso.PPP.lm), col= "red") # add a perfect fit line
summary(Meiso.PPP.lm)
Meiso.PPP.BLUPs <- ranef(Meiso.PPP.lm) # Extract BLUPs of random terms
hist(Meiso.PPP.BLUPs$ID[,1]+9.08, xlab = "Avg. plant per plot", main = "")
range(Meiso.PPP.BLUPs$ID[,1]+9.08)

# Grain yield per primary panicle, adjusted for plant per plot
Meiso.GYP.adj.lm <- lmer(GYP ~ 1 + plant.per.plot + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.GYP.adj.lm)) 
qqline(resid(Meiso.GYP.adj.lm), col= "red") # add a perfect fit line
summary(Meiso.GYP.adj.lm)
Meiso.GYP.adj.BLUPs <- ranef(Meiso.GYP.adj.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.GYP.adj.BLUPs$ID[,1]+57.4977, xlab = "Grain yield per primary panicle", main = "")

# Thousand seed weight
Meiso.TSW.lm <- lmer(TSW ~ 1 + (1|ID) + (1|Replication/Block), data = Meiso)
qqnorm(resid(Meiso.TSW.lm)) 
qqline(resid(Meiso.TSW.lm), col= "red") # add a perfect fit line
summary(Meiso.TSW.lm)
Meiso.TSW.BLUPs <- ranef(Meiso.TSW.lm) # Extract BLUPs of random terms
plotNormalHistogram(Meiso.TSW.BLUPs$ID[,1]+20.814, xlab = "Thousand seed weight", main = "")

Trait.BLUPs.Meiso <- cbind(Meiso.DF.BLUPs$ID, Meiso.DM.BLUPs$ID, Meiso.HE.BLUPs$ID, Meiso.LS.BLUPs$ID,
                          Meiso.NL.BLUPs$ID, Meiso.NPP.BLUPs$ID, Meiso.NT.BLUPs$ID, Meiso.PH.BLUPs$ID,
                          Meiso.GYP.BLUPs$ID, Meiso.GYP.adj.BLUPs$ID,Meiso.TSW.BLUPs$ID, Meiso.PPP.BLUPs$ID)
colnames(Trait.BLUPs.Meiso) <- c("DF", "DM", "HE", "LS", "NL", "NPP", "NT", "PH", "GYP", "GYP.adj", "TSW", "PPP")

Trait.BLUPs.Meiso$DF <- Trait.BLUPs.Meiso$DF + 85.032
Trait.BLUPs.Meiso$DM <- Trait.BLUPs.Meiso$DM + 122.66
Trait.BLUPs.Meiso$HE <- Trait.BLUPs.Meiso$HE + 2.9275
Trait.BLUPs.Meiso$LS <- Trait.BLUPs.Meiso$LS + 2.8095
Trait.BLUPs.Meiso$NL <- Trait.BLUPs.Meiso$NL + 9.9454
Trait.BLUPs.Meiso$NPP <- Trait.BLUPs.Meiso$NPP + 1.3153
Trait.BLUPs.Meiso$NT <- Trait.BLUPs.Meiso$NT + 0.21
Trait.BLUPs.Meiso$PH <- Trait.BLUPs.Meiso$PH + 177.797
Trait.BLUPs.Meiso$GYP <- Trait.BLUPs.Meiso$GYP + 54.5603
Trait.BLUPs.Meiso$GYP.adj <- Trait.BLUPs.Meiso$GYP.adj + 57.4977
Trait.BLUPs.Meiso$TSW <- Trait.BLUPs.Meiso$TSW + 20.814
Trait.BLUPs.Meiso$PPP <- Trait.BLUPs.Meiso$PPP + 9.08

#library(psych)
pairs.panels(Trait.BLUPs.Meiso, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             cex.cor = 1.5)


##############################################################################
# Sheraro
Replication <- as.factor(Sheraro$Replication)
Block <- as.factor(Sheraro$Block)

Sheraro.DF.lm <- lmer(DF ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.DF.lm)) 
qqline(resid(Sheraro.DF.lm), col= "red") # add a perfect fit line
summary(Sheraro.DF.lm)
Sheraro.DF.BLUPs <- ranef(Sheraro.DF.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.DF.BLUPs$ID[,1]+64.99, xlab = "Days to 50% flowering", main = "")

# Days to maturity, DO NOT USE!
Sheraro.DM.lm <- lmer(DM ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.DM.lm)) 
qqline(resid(Sheraro.DM.lm), col= "red") # add a perfect fit line
summary(Sheraro.DM.lm)
Sheraro.DM.BLUPs <- ranef(Sheraro.DM.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.DM.BLUPs$ID[,1]+90.2196, xlab = "Days to maturity", main = "")

# Plant height, need correction for outliers.
Sheraro.PH.lm <- lmer(PH ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.PH.lm)) 
qqline(resid(Sheraro.PH.lm), col= "red") # add a perfect fit line
summary(Sheraro.PH.lm)
Sheraro.PH.BLUPs <- ranef(Sheraro.PH.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.PH.BLUPs$ID[,1]+279.5986, xlab = "Plant height", main = "")

# Number of panicle per plant
Sheraro.NPP.lm <- lmer(NPP ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.NPP.lm)) 
qqline(resid(Sheraro.NPP.lm), col= "red") # add a perfect fit line
bc <- boxcox(NPP ~ 1, lambda = seq(-10,10,0.1), data = Sheraro) # box-cox transformation
lambda <- bc$x[which.max(bc$y)]
plotNormalHistogram(Sheraro$NPP)
plotNormalHistogram(((Sheraro$NPP)^lambda-1)/lambda)
# Redo lmer with transformed data:
Sheraro.NPP.lm <- lmer(((Sheraro$NPP)^lambda-1)/lambda ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
summary(Sheraro.NPP.lm)
Sheraro.NPP.BLUPs <- ranef(Sheraro.NPP.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.NPP.BLUPs$ID[,1]+0.063, xlab = "Number of panicles per plant", main = "")

# Number of leaves
Sheraro.NL.lm <- lmer(NL ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.NL.lm)) 
qqline(resid(Sheraro.NL.lm), col= "red") # add a perfect fit line
summary(Sheraro.NL.lm)
Sheraro.NL.BLUPs <- ranef(Sheraro.NL.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.NL.BLUPs$ID[,1]+13.5781, xlab = "Number of leaves", main = "")

# Grain yield per primary panicle
Sheraro.GYP.lm <- lmer(GYP ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.GYP.lm)) 
qqline(resid(Sheraro.GYP.lm), col= "red") # add a perfect fit line
summary(Sheraro.GYP.lm)
Sheraro.GYP.BLUPs <- ranef(Sheraro.GYP.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.GYP.BLUPs$ID[,1]+38.204, xlab = "Grain yield per primary panicle", main = "")

# Plant per plot
Sheraro.PPP.lm <- lmer(plant.per.plot ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.PPP.lm)) 
qqline(resid(Sheraro.PPP.lm), col= "red") # add a perfect fit line
summary(Sheraro.PPP.lm)
Sheraro.PPP.BLUPs <- ranef(Sheraro.PPP.lm) # Extract BLUPs of random terms
hist(Sheraro.PPP.BLUPs$ID[,1]+15.672, xlab = "Avg. plant per plot", main = "")
range(Sheraro.PPP.BLUPs$ID[,1]+15.672)

# Grain yield per primary panicle, adjusted for plant per plot
Sheraro.GYP.adj.lm <- lmer(GYP ~ 1 + plant.per.plot + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.GYP.adj.lm)) 
qqline(resid(Sheraro.GYP.adj.lm), col= "red") # add a perfect fit line
summary(Sheraro.GYP.adj.lm)
Sheraro.GYP.adj.BLUPs <- ranef(Sheraro.GYP.adj.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.GYP.adj.BLUPs$ID[,1]+38.4083, xlab = "Grain yield per primary panicle", main = "")

# Thousand seed weight
Sheraro.TSW.lm <- lmer(TSW ~ 1 + (1|ID) + (1|Replication/Block), data = Sheraro)
qqnorm(resid(Sheraro.TSW.lm)) 
qqline(resid(Sheraro.TSW.lm), col= "red") # add a perfect fit line
summary(Sheraro.TSW.lm)
Sheraro.TSW.BLUPs <- ranef(Sheraro.TSW.lm) # Extract BLUPs of random terms
plotNormalHistogram(Sheraro.TSW.BLUPs$ID[,1]+19.6736, xlab = "Thousand seed weight", main = "")

Trait.BLUPs.Sheraro <- cbind(Sheraro.DF.BLUPs$ID, Sheraro.DM.BLUPs$ID,
                           Sheraro.NL.BLUPs$ID, Sheraro.NPP.BLUPs$ID, Sheraro.PH.BLUPs$ID,
                           Sheraro.GYP.BLUPs$ID, Sheraro.GYP.adj.BLUPs$ID, Sheraro.TSW.BLUPs$ID, Sheraro.PPP.BLUPs$ID)
colnames(Trait.BLUPs.Sheraro) <- c("DF", "DM", "NL", "NPP", "PH", "GYP", "GYP.adj", "TSW", "PPP")

Trait.BLUPs.Sheraro$DF <- Trait.BLUPs.Sheraro$DF + 64.99
Trait.BLUPs.Sheraro$DM <- Trait.BLUPs.Sheraro$DM + 90.22
Trait.BLUPs.Sheraro$NL <- Trait.BLUPs.Sheraro$NL + 13.6
Trait.BLUPs.Sheraro$NPP <- Trait.BLUPs.Sheraro$NPP + 0.063
Trait.BLUPs.Sheraro$PH <- Trait.BLUPs.Sheraro$PH + 279.6
Trait.BLUPs.Sheraro$GYP <- Trait.BLUPs.Sheraro$GYP + 38.204
Trait.BLUPs.Sheraro$GYP.adj <- Trait.BLUPs.Sheraro$GYP.adj + 38.4083
Trait.BLUPs.Sheraro$TSW <- Trait.BLUPs.Sheraro$TSW + 19.67
Trait.BLUPs.Sheraro$PPP <- Trait.BLUPs.Sheraro$PPP + 15.672

pairs.panels(Trait.BLUPs.Sheraro, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             cex.cor = 1.5)

#write.csv(Trait.BLUPs.Kobo, "1Trait BLUPs in Kobo.csv")
#write.csv(Trait.BLUPs.Meiso, "1Trait BLUPs in Meiso.csv")
#write.csv(Trait.BLUPs.Sheraro, "1Trait BLUPs in Sheraro.csv")


library(ggplot2)
library(plyr)
library(grid)
library(ggpubr)
library(reshape2)
setwd("~/Dropbox/PGML_Projects/Raw Pheno")

# Import blups from three environments
blup.Kobo <- read.csv("Trait BLUPs in Kobo.csv", header = TRUE)
blup.Meiso <- read.csv("Trait BLUPs in Meiso.csv", header = TRUE)
blup.Sheraro <- read.csv("Trait BLUPs in Sheraro.csv", header = TRUE)
blup.all <- read.csv("Trait BLUPs of three environments.csv", header = TRUE)
Retained.1178 <- read.csv("1178 BCNAM lines.csv", header = TRUE)

# Create Grain filling period
blup.Kobo$GFP <- blup.Kobo$DM - blup.Kobo$DF
blup.Meiso$GFP <- blup.Meiso$DM - blup.Meiso$DF
blup.Sheraro$GFP <- blup.Sheraro$DM - blup.Sheraro$DF
blup.all$GFP <- blup.all$DM - blup.all$DF

# Subset 1178 lines
blup.Kobo <- blup.Kobo[blup.Kobo$Taxa%in%Retained.1178$Taxa,]
blup.Meiso <- blup.Meiso[blup.Meiso$Taxa%in%Retained.1178$Taxa,]
blup.Sheraro <- blup.Sheraro[blup.Sheraro$Taxa%in%Retained.1178$Taxa,]
blup.all <- blup.all[blup.all$Taxa%in%Retained.1178$Taxa,]

colnames(blup.Kobo)[c(11,12,14,15)] <- c("SGY","SGY.adj","CPP","PGY")
colnames(blup.Meiso)[c(11,12,14,15)] <- c("SGY","SGY.adj","CPP","PGY")
colnames(blup.Sheraro)[c(8,9,11,12)] <- c("SGY","SGY.adj","CPP","PGY")
#############################################################################################################################
# get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

# Heatmap in ggplot2
# correlation heatmap in Kobo
mydata.kobo <- blup.Kobo[,-c(1,2,8,12,16)]
cormat <- round(cor(mydata.kobo),2)

# Melt the correlation matrix with reshape2 package
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
cor.kobo <- ggplot(data = melted_cormat, aes(Var2, Var1, fill=value)) + 
  geom_tile(col="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 0.9999), space = "Lab", name = "Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.title.x = element_text(vjust = 1, size = 12, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label=value), col="black") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1,0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal") + 
  guides(fill=guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

# correlation heatmap in Meiso
mydata.Meiso <- blup.Meiso[,-c(1,2,8,12,16)]
cormat <- round(cor(mydata.Meiso),2)

# Melt the correlation matrix with reshape2 package
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap in ggplot2
cor.Meiso <- ggplot(data = melted_cormat, aes(Var2, Var1, fill=value)) + 
  geom_tile(col="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 0.9999), space = "Lab", name = "Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.title.x = element_text(vjust = 1, size = 12, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label=value), col="black") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1,0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal") + 
  guides(fill=guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

# correlation heatmap in Sheraro
mydata.Sheraro <- blup.Sheraro[,-c(1,2,6,9,13)]
cormat <- round(cor(mydata.Sheraro),2)

# Melt the correlation matrix with reshape2 package
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap in ggplot2
cor.Sheraro <- ggplot(data = melted_cormat, aes(Var2, Var1, fill=value)) + 
  geom_tile(col="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 0.9999), space = "Lab", name = "Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.title.x = element_text(vjust = 1, size = 12, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label=value), col="black") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1,0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal") + 
  guides(fill=guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

tiff("Figure 2. Trait correlations.tiff", units = "in", width = 10, height = 10, res = 300)
ggarrange(cor.kobo, cor.Meiso, cor.Sheraro, nrow=2, ncol=2, labels = c("A","B","C"), common.legend = TRUE, legend = "top")
dev.off()


#############################################################################################################################
# boxplot
box.df <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=DF, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 77.41), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 83.88), color="#00BA38", linetype="dashed") +
  ylab("Days to 50% flowering") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.title.x=element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.dm <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=DM, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 125.58), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 120.95), color="#00BA38", linetype="dashed") +
  ylab("Days to maturity") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.title.x=element_blank(),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.nl <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=NL, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 8.78), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 9.90), color="#00BA38", linetype="dashed") +
  ylab("Number of tillers") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.he <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234" & !blup.all$Environment=="Sheraro",], aes(x=Population, y=HE, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#F8766D","#00BA38"))+
  scale_colour_manual(values=c("#F8766D","#00BA38"))+
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 1.12), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 2.23), color="#00BA38", linetype="dashed") +
  ylab("Head exsertion") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.ls <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234" & !blup.all$Environment=="Sheraro",], aes(x=Population, y=LS, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#F8766D","#00BA38"))+
  scale_colour_manual(values=c("#F8766D","#00BA38"))+
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 2.67), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 2.86), color="#00BA38", linetype="dashed") +
  ylab("Leaf senescence") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.nt <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234" & !blup.all$Environment=="Sheraro",], aes(x=Population, y=NT, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("#F8766D","#00BA38"))+
  scale_colour_manual(values=c("#F8766D","#00BA38"))+
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 0.66), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 0.20), color="#00BA38", linetype="dashed") +
  ylab("Number of tillers") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.ph <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=PH, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 179.14), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 181.48), color="#00BA38", linetype="dashed") +
  ylab("Plant height") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.gyp <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=GYP.adj, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 69.16), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 54.95), color="#00BA38", linetype="dashed") +
  ylab("Single plant-based grain yield") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.pgy <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=Plot.Yield, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 401), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 539), color="#00BA38", linetype="dashed") +
  ylab("Plot-based grain yield") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.cpp <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=PPP, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 8.8), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 9.3), color="#00BA38", linetype="dashed") +
  ylab("Count of plants per plot") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

box.tsw <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Population, y=TSW, col=Environment, fill=Environment)) +
  geom_boxplot() + 
  theme_bw(base_size = 12) + 
  geom_hline(aes(yintercept = 29.27), color="#F8766D", linetype="dashed") +
  geom_hline(aes(yintercept = 19.63), color="#00BA38", linetype="dashed") +
  ylab("1000-seed weight") +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1), 
    axis.text.y = element_text(color = "black", size = 11))

# Eight traits: NL, HX, LS, NT, SGY, TSW, CPP, PGY
tiff("Figure_8 Trait boxplots by pop at three envs.tiff", units = "in", width = 8, height = 10, res = 300)
ggarrange(box.nl, box.he, box.ls, box.nt, box.gyp, box.tsw, box.cpp, box.pgy, nrow=3, ncol=3, common.legend = TRUE)
dev.off()

# The other three traits: DF, DM, PH
tiff("Figure_3 Trait boxplots by pop at three envs.tiff", units = "in", width = 8, height = 6, res = 300)
ggarrange(box.df, box.dm, box.ph, nrow=3, ncol=1, common.legend = TRUE)
dev.off()

################################################################################################################################
# Density plot for each trait
# Consider the whole BCNAM population together, across three envs 

# DF: days to 50% flowering
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(DF))
head(mu)

p.df <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=DF, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("Days to 50% flowering")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# DM: days to maturity
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(DM))
head(mu)

p.dm <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=DM, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("Days to maturity")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# HE: head exsertion
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(HE))
head(mu)

p.he <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234" & !blup.all$Environment=="Sheraro",], aes(x=HE, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  scale_fill_manual(values=c("#F8766D","#00BA38"))+
  scale_colour_manual(values=c("#F8766D","#00BA38"))+
  geom_vline(data=mu[-3,], aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("Head exsertion (cm)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# LS: leaf senescence
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(LS))
head(mu)

p.ls <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234" & !blup.all$Environment=="Sheraro",], aes(x=LS, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  scale_fill_manual(values=c("#F8766D","#00BA38"))+
  scale_colour_manual(values=c("#F8766D","#00BA38"))+
  geom_vline(data=mu[-3,], aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("Leaf senescence (score)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))


# NL: number of leaves
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(NL))
head(mu)

p.nl <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=NL, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("No. of leaves (count)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# NT: number of tillers
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(NT))
head(mu)

p.nt <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234" & !blup.all$Environment=="Sheraro",], aes(x=NT, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  scale_fill_manual(values=c("#F8766D","#00BA38"))+
  scale_colour_manual(values=c("#F8766D","#00BA38"))+
  geom_vline(data=mu[-3,], aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("No. of tillers (count)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))


# PH: plant height
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(PH))
head(mu)

p.ph <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=PH, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("Plant height (cm)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# SGY: single plant-based grain yield
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(GYP))
head(mu)

p.gyp <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=GYP, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("Singe plant-based yield (g)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# PGY: Plot-based grain yield
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(Plot.Yield))
head(mu)

p.plotyield <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=Plot.Yield, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("Plot-based yield (g)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# TSW: 1000-seed weight
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(TSW))
head(mu)

p.tsw <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=TSW, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("1000-seed weight (g)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))

# CPP: count of plants per plot
mu <- ddply(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], "Environment", summarise, grp.mean=mean(PPP))
head(mu)

p.ppp <- ggplot(blup.all[!blup.all$Population=="Parent" & !blup.all$Population=="IS32234",], aes(x=PPP, fill=Environment)) +
  geom_density(aes(y=..scaled..), alpha=0.6) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Environment), linetype="dashed", lwd=1) +
  ylab("Density") +
  xlab("No. of plants per plot (count)")+
  theme_bw(base_size = 12) +
  theme(
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 11))


#ggarrange(p.df, p.dm, p.nl, p.he, p.ls, p.nt, p.ph, p.gyp, p.plotyield, p.ppp, p.tsw, nrow=4, ncol=3, common.legend = TRUE)
#dev.off()

tiff("Figure_Density plot of each trait.tiff", units = "in", width = 8, height = 6, res = 300)
ggarrange(p.nl, p.he, p.ls, p.nt, p.gyp, p.tsw, p.ppp, p.plotyield, nrow=3, ncol=3, common.legend = TRUE)
dev.off()

# summary statistics by population
library(psych)
Kobo.summary <- describeBy(blup.all[blup.all$Environment=="Kobo",], group = blup.all[blup.all$Environment=="Kobo",]$Population, mat = TRUE, digits = 2)
Meiso.summary <- describeBy(blup.all[blup.all$Environment=="Meiso",], group = blup.all[blup.all$Environment=="Meiso",]$Population, mat = TRUE, digits = 2)
Sheraro.summary <- describeBy(blup.all[blup.all$Environment=="Sheraro",], group = blup.all[blup.all$Environment=="Sheraro",]$Population, mat = TRUE, digits = 2)

write.csv(Kobo.summary, "Pheno summary by population in Kobo.csv")
write.csv(Meiso.summary, "Pheno summary by population in Meiso.csv")
write.csv(Sheraro.summary, "Pheno summary by population in Sheraro.csv")



setwd("~/Dropbox/PGML_Projects/Raw Pheno")
# Only progeny, no parents included
Pheno.raw <- read.csv("Sorghum BCNAM Raw Data across three environments.csv", header = TRUE)
Retained.1178 <- read.csv("1178 BCNAM lines.csv", header = TRUE)
Pheno <- Pheno.raw[!is.na(Pheno.raw$ID) & Pheno.raw$ID%in%Retained.1178$Taxa,]

Kobo <- Pheno[Pheno$Location=="Kobo",]
Meiso <- Pheno[Pheno$Location=="Meiso",]
Sheraro <- Pheno[Pheno$Location=="Sheraro",]

length(unique(Kobo$ID)) # 1031 genotypes
length(unique(Meiso$ID)) # 1155 genotypes
length(unique(Sheraro$ID)) # 1046 genotypes
length(unique(Pheno$ID)) # 1242

library(lme4)
library(MASS)
library(rcompanion)
library(car)

# ANOVA
m <- lm(LS ~ Population+ID+Location+Location:Population+Location:ID, 
        data = Pheno)
anova(m)

m1 <- lmer(NT ~ 1+ (1|Population) + (1|ID), 
           data = Meiso)
summary(m1)
