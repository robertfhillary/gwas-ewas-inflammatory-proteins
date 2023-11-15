## Processing of Summary Files 

setwd("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/")
x <- read.csv("Normalised_data.csv")
x = x[-c(1048:1050),]
x$QC.Warning <- NULL
x$BDNF <- NULL 

seventy <- read.csv("Seventy_proteins.csv")
x = x[,c(2,94,which(colnames(x) %in% seventy$Protein))]


pc <- read.csv("LBC1936_PCA_from_Gail.csv")
phen <- read.table("example_phen.phen", header = T)
target <- read.csv("target_QC_age_sex_date.csv")

target = target[target$WAVE == 1,]
target = target[target$cohort %in% "LBC36",]
target = target[which(target$ID %in% phen$IID),]
target <- target[,c("age", "sex", "ID")]
x=x[which(x$lbc36no %in% target$ID),]
names(x)[1] <- "ID"

file = merge(target, x, by = "ID")

names(pc)[1] <- "ID"
file1 = merge(file, pc, by = "ID")

## Age 


list <- list()
for(i in colnames(file1)[5:74]){ 
  list[[i]] <- summary(lm(scale(file1[,i]) ~ scale(file1$age)))
  }

age <- do.call(rbind,lapply(list,function(x)coefficients(x)[2,c(1,2,4)]))



## Sex 

list <- list()
for(i in colnames(file1)[5:74]){ 
  list[[i]] <- summary(lm(scale(file1[,i]) ~ file1$sex))
}

sex <- do.call(rbind,lapply(list,function(x)coefficients(x)[2,c(1,2,4)]))




## PC1 

list <- list()
for(i in colnames(file1)[5:74]){ 
  list[[i]] <- summary(lm(scale(file1[,i]) ~ scale(file1$C1)))
}

PC1 <- do.call(rbind,lapply(list,function(x)coefficients(x)[2,c(1,2,4)]))



## PC2

list <- list()
for(i in colnames(file1)[5:74]){ 
  list[[i]] <- summary(lm(scale(file1[,i]) ~ scale(file1$C2)))
}

PC2 <- do.call(rbind,lapply(list,function(x)coefficients(x)[2,c(1,2,4)]))



## PC3

list <- list()
for(i in colnames(file1)[5:74]){ 
  list[[i]] <- summary(lm(scale(file1[,i]) ~ scale(file1$C3)))
}

PC3 <- do.call(rbind,lapply(list,function(x)coefficients(x)[2,c(1,2,4)]))



## PC4

list <- list()
for(i in colnames(file1)[5:74]){ 
  list[[i]] <- summary(lm(scale(file1[,i]) ~ scale(file1$C4)))
}

PC4 <- do.call(rbind,lapply(list,function(x)coefficients(x)[2,c(1,2,4)]))


## Plate 

anovalist <- list()
for(i in colnames(file1)[5:74]){
anovalist[[i]] = summary(aov(scale(file1[,i]) ~ as.factor(file1$Plate.ID)))
} 

anovaout <- do.call(rbind,lapply(anovalist,function(x)x[[1]][1,c(4,5)]))
anovaout <- as.data.frame(anovaout)

## Combine All Regressions for 70 Proteins 

total <- cbind(age, cbind(sex, cbind(PC1, cbind(PC2, cbind(PC3, cbind(PC4, cbind(anovaout)))))))

