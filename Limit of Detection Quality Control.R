
## Quality Control Proteins to Those Where Samples Surpassing Limit of Detection >= 40% 

raw = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Raw_data.csv")

lod = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Limits_of_detection.csv") 

raw$Plate.ID <- NULL 
raw = raw[-which(raw$QC.Warning %in% "Warning"),]
raw$QC.Warning <- NULL
raw <- raw[-c(1018:1020),] 
row.names(raw) <- raw$lbc36no
raw$lbc36no <- NULL


beyond_limit <- vector()
for(i in colnames(raw)){ 
  beyond_limit[i] <- length(which(raw[,i] <= lod[,i]))
}

beyond_limit<-as.data.frame(beyond_limit)
beyond_limit$Prop <- beyond_limit$beyond_limit/1017
beyond_limit$Proteins <- row.names(beyond_limit)

list <- beyond_limit[which(beyond_limit$Prop >= 0.4),"Proteins"] 
raw1=raw[,-which(colnames(raw) %in% list)]
raw1$BDNF=NULL

