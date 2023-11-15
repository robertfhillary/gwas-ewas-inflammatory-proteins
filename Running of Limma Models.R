## Install Prerequisite Packages

if(!require(limma)){
  install.packages("limma")
}

if(!require(dplyr)){
  install.packages("dplyr")
}


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

## Load Packages 

library(limma)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

## Preparation of Files 

pData <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Limma Inflammatory/Input/pData.csv") 

meth <- load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")

a = which(colnames(dat) %in% pData$Basename)
dat1 = dat[,a]
ids = pData$Basename
dat1 = dat1[,match(ids, colnames(dat1))]

rownames(pData) <- pData$Basename 


all.equal(rownames(pData), colnames(dat1))


saveRDS(dat1, "/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Limma Inflammatory/Input/Methylation_Inflam.rds") 
saveRDS(pData, "/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Limma Inflammatory/Input/Phenotypic_Inflam.rds") 


## Limma Models

IL8 <- model.matrix(~IL8 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

VEGFA <- model.matrix(~VEGFA + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

MCP3 <- model.matrix(~MCP.3 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

GDNF <- model.matrix(~GDNF + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

CDCP1 <- model.matrix(~CDCP1 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

CD244 <- model.matrix(~CD244 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

IL7 <- model.matrix(~IL7 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

OPG <- model.matrix(~OPG + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

LAPTGFBeta1 <- model.matrix(~LAP.TGF.beta.1 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

uPA <- model.matrix(~uPA + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

IL6 <- model.matrix(~IL6 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

IL17C <- model.matrix(~IL.17C + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

MCP1 <- model.matrix(~MCP.1 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

IL17A <- model.matrix(~IL.17A + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

CXCL11 <- model.matrix(~CXCL11 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)



fits <- list()
models <- c("IL8", "VEGFA", "MCP3", "GDNF", 
           "CDCP1", "CD244", "IL7", "OPG", "LAPTGFBeta1", 
           "uPA", "IL6", "IL17C", "MCP1", "IL17A", "CXCL11")

for(i in models) { 
print(i)
timestamp()
fits[[i]] <- lmFit(dat1[,rownames(get(i))], get(i))
}

saveRDS(fits, file="/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Limma Inflammatory/Output/lmFit_First15_Proteins.rds")

fits2 <- list()
for(i in models) {
print(i)
timestamp()
fits2[[i]] <- eBayes(fits[[i]])
}
saveRDS(fits2, file="/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Limma Inflammatory/Output/eBayes_First15_Proteins.rds")

library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
genes <- data.frame(ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island)
rownames(genes) <- genes$ID

tt <- list()
for(i in models){
tt[[i]] <- topTable(fits2[[i]], coef=2, adjust.method="none", number=nrow(fits2[[i]]), 
                  genelist=genes[rownames(fits2[[i]]),c("geneSymbol", "CHR", "MAPINFO", "FEATURE", "CpGISLAND")])
				  }
saveRDS(tt, file="/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Limma Inflammatory/Output/TopTables_First15_Proteins.rds")

