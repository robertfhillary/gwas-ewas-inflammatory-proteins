## Install Prerequisite Packages 

if(!require(ggplot2)){
  install.packages("ggplot2")
}


if(!require(coloc)){
  install.packages("coloc")
}

if(!require(corrplot)){
  install.packages("corrplot")
}


if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
}

if(!require(wesanderson)){
  install.packages("wesanderson")
}

## Load Library Packages 

library(ggplot2)
library(coloc)
library(corrplot)
library(RColorBrewer)
library(wesanderson)

install.packages("devtools")
library(devtools)
install_github("phenoscanner/phenoscanner")
library(phenoscanner)

### Genetic Figures ###

    ## Genetic Atlas 

gen = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Genetic_data_prep1.csv")
p = ggplot(gen, aes(gen$Relative_Hit_Position, gen$Relative_Gene_Position))
q = p + geom_jitter(aes(colour = as.factor(gen$Type))) + xlab("pQTL Position") + ylab("Protein Position")
q = q + scale_x_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + scale_y_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + theme_bw()
q = q + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic")) + theme(panel.background = element_rect(colour = "black")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.1, color = "darkslategray4")
q = q + annotate(x = 3.85, y = 16.7, label = expression(paste("rs12075 (", italic("DARC"), ")")), geom = "text", size = 2.95)
q + annotate(x = 3.85, y = 15.9, label = paste("for MCP4 levels"), geom = "text", size = 2.95)

  ## Distance v P Values 

cis = gen[gen$Type %in% "Cis",]
cis$diff = as.numeric(as.character(cis$diff))
cis$diff = cis$diff/1e6
q <- ggplot(data = cis, aes(cis$diff, -log10(cis$P.Value))) 
q1 = q + geom_point() + xlab("Distance of Cis Variant to TSS (Mb)") + labs(y = expression(-log[10](P))) + ylim(8.55,100) + xlim(-0.3, 0.3) 
q1 + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic")) + theme(legend.position=c(0.88, 0.84)) + theme(legend.key.size = unit(0.6, "cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


    ## MAF v Beta 

beta = ggplot(gen, aes(gen$X, abs(gen$Beta)))
beta1 = beta + geom_point(aes(colour = factor(gen$Type))) + xlab("Minor Allele Frequency") + ylab("Effect Size") +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5, 2.0, 2.5)) 
beta2 = beta1  + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic")) + theme(legend.position=c(0.88, 0.84)) + theme(legend.key.size = unit(0.6, "cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
beta2 + theme(axis.title.x = element_text(size=13.5), axis.title.y = element_text(size=13.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))


    ## Variance explained by SNPs
y <- matrix(ncol = 2, nrow = 4)
colnames(y) <- c("Group", "Freq")
y <- data.frame(y)   
y$Group <- c("1","2","3","4")
y$Freq <- c(7, 5, 0, 1)
colour = c("#f9a65a", "#79c36a", "#9e66ab", "#d77fb3")
bar = ggplot(data = y, aes(y$Group, y$Freq)) + geom_bar(stat="identity", fill = colour) + xlab("Variance explained by SNPs") + ylab("Number of SNPs") + scale_x_discrete(labels = c("0-10%", "11-20%", "21-30%", ">30%"))
bar + theme(axis.title.x = element_text(size=13.5), axis.title.y = element_text(size=13.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7))

    ## Proportion of Variants and their Function 
a = as.data.frame(table(gen$Predicted.Function))
names(a)[1] <- "Function"
colours = c("#E78AC3", "#8DA0CB", "#66C2A5", "#FC8D62")
a$Proportion_of_Variants <- a$Freq/13
bar = ggplot(data = a, aes(a$Function, a$Proportion_of_Variants)) +
  geom_bar(stat="identity", fill = colours) + ylab("Proportion of Variants")
bar = bar + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
bar + theme(axis.title.x=element_blank()) + theme(axis.text.x = element_text(size = 11, colour = "grey20")) + theme(axis.title.y = element_text(size = 14)) + theme(axis.text.y = element_text(size = 11))

  
## Correlation Plots 

prot = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Normalised_data.csv")
prot <- prot[,c("ADA", "CCL25", "CD6", "CST5", "CXCL5", "CXCL6", "FGF.5", "IL.10RB", "IL.12B", "IL.18R1", "MCP.2", "MCP.4", "TNFB")]
corp = cor(prot, use="complete.obs")
corp = round(corp,2)
corrplot(corp, type="upper", order = "hclust", tl.col = "black",addCoef.col = "black")

### Epigenetic Figures ###

  ## Epigenetic Atlas

epi = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Epigenetic_data_prep1.csv")
epi$Present <- gsub("BOTH", "Linear and BayesR+", epi$Present)
epi$Present <- gsub("LINEAR", "Linear Only", epi$Present)
epi$Present <- gsub("BAYES", "BayesR+ Only", epi$Present)
epi$Present <- gsub("ALL", "All Models", epi$Present)
p = ggplot(epi, aes(epi$Relative_Hit_Position, epi$Relative_Gene_Position))
q = p + geom_jitter(aes(colour = as.factor(epi$Type), shape = as.factor(epi$Present)),width =0.21, size = 2.2) + xlab("CpG Position") + ylab("Protein Position")
q = q + scale_x_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + scale_y_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + theme_bw()
q = q + theme(legend.text = element_text(face = "italic")) + theme(panel.background = element_rect(colour = "black")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.1, color = "darkslategray4")
q = q + annotate(x = 5, y = 16.5, label = expression(paste("cg05575921 (", italic("AHRR"), ")")), geom = "text", size = 2.75)
q = q + annotate(x = 5, y = 16.0, label = expression(paste("for CCL11 levels")), geom = "text", size = 2.75)
q = q + annotate(x = 2.4, y = 5, label = expression(paste("cg03938978 (", italic("IL18RAP"), ")")), geom = "text", size = 2.6)
q = q + annotate(x = 2.4, y = 4.5, label = expression(paste("for IL18R1 levels")), geom = "text", size = 2.6)
q = q + annotate(x = 16.9, y = 3.6, label = expression(paste("cg07839457 (", italic("NLRC5"), ")")), geom = "text", size = 2.75)
q = q + annotate(x = 16.9, y = 3.1, label = expression(paste("for CXCL9 levels")), geom = "text", size = 2.75)
q = q  + labs(colour = "Type") + labs(shape = "Concordance across Models") + theme(legend.title = element_text(hjust = 0.5))
q + scale_shape_manual(values=c(18, 3, 0, 17))

  ## Correlation between existing betas and those in Hillary et al 

betas = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Replication_review3.csv")
cor.test(betas$beta, betas$beta_hillary, use = "complete.obs")


beta = ggplot(betas, aes(betas$beta, betas$beta_hillary))
beta1 = beta + geom_point(aes(color = factor(betas$study)))

beta2 = beta1 + 
  geom_smooth(method = "lm", se = FALSE) + 
  xlab("Betas from Literature") + 
  ylab(expression(paste("Betas from Hillary ", italic("et al.")))) +
  ggtitle("Correlation of Betas for pQTLs") +
  theme(plot.title= element_text(hjust =0.5)) + 
  labs(colour="Protein", shape="Study")
beta3 = beta2 + geom_abline(intercept = 0, slope = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
beta3 + annotate(x = 1.65, y = 0.52, label = expression(paste(italic("r"), " = 0.70")), size = 5.2, geom = "text") + theme(legend.title = element_text(hjust = 0.5))


## BayesR+ - Plotting of Variance Explained According To Mixture Variances

## Genetics

herit = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Heritability_1.csv")
list = c("ADA", "CCL25", "CD6", "CST5", "CXCL5", "CXCL6", "FGF.5", "IL.10RB", "IL.12B", "IL.18R1", "MCP.2", "MCP.4", "TNFB")
herit=herit[which(herit$Biomarker %in% list),]
protein <- rep(herit$Biomarker, 2)
protein <- protein[order(protein)]
level <- rep(c("Small" , "Large") , 13)
herit[,5] <- gsub("%", "", herit[,5])
herit[,6] <- gsub("%", "", herit[,6])
herit[,2] <- gsub("%", "", herit[,2])
herit[,3] <- gsub("%", "", herit[,3])
herit[,2] <- as.numeric(as.character(herit[,2]))
herit[,3] <- as.numeric(as.character(herit[,3]))
herit[,5] <- herit[,2]*(as.numeric(herit[,5])/100) 
herit[,6] <- herit[,2]*(as.numeric(herit[,6])/100) 
value <- rbind(herit[,5], herit[,6])
herit[,2] <- gsub("%", "", herit[,2])
herit[,4] <- gsub("%", "", herit[,4])
herit[,2] <- as.numeric(as.character(herit[,2]))
herit[,4] <- as.numeric(as.character(herit[,4]))
value <- as.numeric(value)
dat<-data.frame(protein,level,value)
dat$total <- rep(herit[,2],each=2)
dat$min <- rep(herit[,3], each=2)
dat$max = rep(herit[,4], each=2)
dat <- dat[order(dat$total),]
dat$level = gsub("Large", "Large Effects", dat$level)
dat$level = gsub("Small", "Medium Effects", dat$level)
dat$protein <- factor(dat$protein, levels = unique(as.character(dat$protein))[order(dat$total)])
bar = ggplot(data = dat, aes(dat$protein, dat$value, fill = dat$level)) 
bar = bar + geom_bar(position="stack", stat="identity") + labs(fill = "Proprortion of Variance Attributed to") + theme(legend.title.align=0.5) + geom_errorbar(ymin = dat$min, ymax=dat$max, width = 0.25)
bar = bar  + xlab("Protein Biomarker") + ylab("% Variance Explained in Levels by Common Variants") + scale_y_continuous(limits=c(0,80), breaks = c(0,10,20,30,40,50,60,70,80))
bar + theme(axis.text.y = element_text(size=13.5))+ scale_fill_brewer(palette="Accent") + theme(legend.position="top") 

## Epigenetics

herit = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Heritability_Inflam.csv")
list = c("CXCL9", "CCL11", "IL.18R1")
herit=herit[which(herit$Biomarker %in% list),]
protein <- rep(herit$Biomarker, each=3)
protein <- protein[order(protein)]
level <- rep(c("Small Effects", "Medium Effects", "Large Effects") , 3)
herit[,5] <- gsub("%", "", herit[,5])
herit[,6] <- gsub("%", "", herit[,6])
herit[,7] <- gsub("%", "", herit[,7])
herit[,3] <- gsub("%", "", herit[,3])
herit[,2] <- gsub("%", "", herit[,2])
herit[,4] <- gsub("%", "", herit[,4])
herit[,2] <- as.numeric(as.character(herit[,2]))
herit[,4] <- as.numeric(as.character(herit[,4]))
herit[,3] <- as.numeric(as.character(herit[,3]))
herit[,5] <- herit[,2]*(as.numeric(herit[,5])/100) 
herit[,6] <- herit[,2]*(as.numeric(herit[,6])/100) 
herit[,7] <- herit[,2]*(as.numeric(herit[,7])/100) 
value <- rbind(herit[,5], herit[,6])
value <- rbind(value, herit[,7])
value <- as.numeric(value)
dat<-data.frame(protein,level,value)
dat$total <- rep(herit[,2],each=3)
dat$LCI <- rep(herit[,3],each =3)
dat$HCI <- rep(herit[,4],each=3)
dat$protein <- factor(dat$protein, levels = c("CXCL9", "CCL11", "IL.18R1"))
bar = ggplot(data = dat, aes(factor(dat$protein), dat$value, fill = dat$level)) 
bar = bar + geom_bar(position="stack", stat="identity") + labs(fill = "Proprortion of Variance Attributed to") + theme(legend.title.align=0.5) + geom_errorbar(ymin = dat$LCI, ymax = dat$HCI, width = 0.25) 
bar = bar  + xlab("Protein Biomarker") + ylab("% Variance Explained in Levels by DNA Methylation") + scale_y_continuous(limits = c(0,70), breaks = c(0,10,20,30,40,50,60,70,80))
bar + theme(axis.text.y = element_text(size=13.5)) + scale_fill_manual(values = wes_palette(n=3,"FantasticFox1")) + theme(legend.position="top") 


## Correlation between Heritability in Hillary et al versus Ahsan et al 

## How to Calculate Correlation to Account for Variance 

cov1 = cov(herit[,2], herit[,3])
l = herit[,2]
sd_x = sd(l)/sqrt(29)
n = herit[,3]
sd_y = sd(n)/sqrt(29)
varx = var(herit[,2])
vary = var(herit[,3])
cov1/sqrt((vary-sd_y^2)*(varx-sd_x^2))

## Plot 

herit = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Heritability_Correlation.csv")
herit = herit[-30,]
plot = ggplot(data = herit, aes(herit$Estimate.from.Ahsan.et.al., herit$Estimate.from.Hillary.et.al.))  
plot = plot + geom_point() + xlim(0,0.6) + ylim(0,0.6) + xlab(expression(paste("Heritability Estimate from Ahsan",italic(" et al.")), sep= "")) + ylab(expression(paste("Heritability Estimate from Hillary",italic(" et al.")), sep= ""))  
plot = plot + annotate("text", x = 0.50, y = 0.265, label = expression(paste(italic("r"), "= 0.71"))) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_smooth(method = "lm", se = F)
plot + annotate("text", x = 0.57, y = 0.3, label = "MMP-1", size = 3.2)


## Correlation between Heritability - BayesR+ vs OSCA 

herit = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Comparison_Meth_Variance.csv")
herit[,2] <- as.numeric(gsub("%", "", herit[,2]))
herit[,3] <- as.numeric(gsub("%", "", herit[,3]))
herit[,4] <- as.numeric(gsub("%", "", herit[,4]))
herit[,5] <- as.numeric(gsub("%", "", herit[,5]))
herit[,6] <- as.numeric(gsub("%", "", herit[,6]))

herit[,2] = herit[,2]/100
herit[,3] = herit[,3]/100
herit[,4] = herit[,4]/100
herit[,5] = herit[,5]/100
herit[,6] = herit[,6]/100


## Calculation of Correlation Coefficient as Above 

cov = cov(herit$Mean.Variance.Explained, herit$Variance.Explained)
cov/sqrt((var(herit$Variance.Explained)-(sd(n)/sqrt(70))^2)*(var(herit$Mean.Variance.Explained)-(sd(l)/sqrt(70))^2))

## Plot 

plot = ggplot(data = herit, aes(herit$Mean.Variance.Explained, herit$Variance.Explained))  
plot = plot + geom_point() + xlim(0,0.6) + ylim(0,0.5) + xlab(expression(paste("H"^2, " BayesR+", sep = "  "))) + ylab(expression(paste("H"^2, " OSCA", sep= "  ")))  
plot = plot + annotate("text", x = 0.35, y = 0.475, label = "SCF", size = 4, colour = "purple") 
plot = plot + annotate("text", x = 0.51, y = 0.26, label = expression(paste(italic("r"), " = 0.73", sep = " "))) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_smooth(method = "lm", se = F)
plot + ggtitle("Estimate of Variance Explained in Protein Levels", subtitle = "by Genome-Wide DNA Methylation") + theme(plot.title = element_text(hjust=0.5))+ theme(plot.subtitle = element_text(hjust = 0.5))


## Combined BayesR Variance Output

var = read.csv("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/BayesR_Combined_Variance_Output_Infissue.csv")
protein <- rep(var$Biomarker, each=3)
condition <- rep(c("Epigenetics" , "Genetics" , "Combined"), 70)
data <- data.frame(protein,condition)
data$value <- 0
var[,2] <- as.numeric(gsub("%", "", var[,2]))
var[,3] <- as.numeric(gsub("%", "", var[,3]))
var[,4] <- as.numeric(gsub("%", "", var[,4]))
data$value[data$condition %in% "Epigenetics"] <- var$Variance.Estimate.from.Methylation
data$value[data$condition %in% "Genetics"] <- var$Variance.Estimate.from.Genetics
data$value[data$condition %in% "Combined"] <- var$Combined.Genetics.and.Epigenetic.Estimate
data$condition <- factor(data$condition, levels = c("Epigenetics", "Genetics", "Combined"))
# Grouped

## Set 1 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins1-7.pdf", width = 9, height = 6.5)
data1 <- data[1:21,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 2
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins8-14.pdf", width = 9, height = 6.5)
data1 <- data[22:42,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 3 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins15-21.pdf", width = 9, height = 6.5)
data1 <- data[43:63,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 4
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins22-28.pdf", width = 9, height = 6.5)
data1 <- data[64:84,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 5 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins29-35.pdf", width = 9, height = 6.5)
data1 <- data[85:105,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 6 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins36-42.pdf", width = 9, height = 6.5)
data1 <- data[106:126,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 7 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins43-49.pdf", width = 9, height = 6.5)
data1 <- data[127:147,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 8 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins50-56.pdf", width = 9, height = 6.5)
data1 <- data[148:168,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 9 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins57-63.pdf", width = 9, height = 6.5)
data1 <- data[169:189,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()

## Set 10 
pdf("U:/Datastore/IGMM/marioni-lab/Rob/Inflammatory pQTL/Combined_Variance_Plots/CombinedVariance_Proteins64-70.pdf", width = 9, height = 6.5)
data1 <- data[190:210,]
ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Protein Biomarker") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
dev.off()




## Bayesian Colocalisation Analysis Example - colocalisation of expression QTL and protein QTL data

  ## Example Protein: SIGLEC9

  ## Take 400kb region around sentinel cis variant 

  ## Extract known eQTLs in this region from phenoscanner 

res_eqtl <- phenoscanner(regionquery="chr19:51422766-51827766", proxies = "EUR", catalogue = "eQTL") 

  ## Fix Phenoscanner Results so no missing MAF, BETA present and eqtl results are in correct format

snps_from_phenoscanner <- res_eqtl$regions
results_from_phenoscanner <- res_eqtl$results

## Missing MAF 
if(length(which(results_from_phenoscanner$eur %in% "-")) >= 1) { 
  results_from_phenoscanner <- results_from_phenoscanner[-which(results_from_phenoscanner$eur %in% "-"),]
} 

## Missing BETA 
if(length(which(results_from_phenoscanner$beta %in% "NA")) >= 1) { 
  results_from_phenoscanner <- results_from_phenoscanner[-which(results_from_phenoscanner$beta %in% "NA"),]
} 

## Subset to Gene of Interest 
results_from_phenoscanner1 <- results_from_phenoscanner[which(results_from_phenoscanner$exp_gene %in% "SIGLEC9"),]
eqtls1 <- results_from_phenoscanner1

eqtls1$eur <- as.numeric(eqtls1$eur)
eqtls1$MAF <- ifelse(eqtls1$eur > 0.5, 1-eqtls1$eur, eqtls1$eur)
eqtls1$beta <- as.numeric(eqtls1$beta)

  ## Import GWAS Summary Statistics subset to 400 kb (+/- 200 kb) Region of Interest 

gwas_hits1 <- read.csv("siglec9_gwas.csv")

  
  ## Coloc analysis - calculation of abf (approximate Bayes Factor) and posterior probabilities 

dataset1 = list(snp = as.character(eqtls1$rsid), N = 31684, pvalues = as.numeric(eqtls1$p), MAF = as.numeric(eqtls1$MAF), type = "quant")
dataset2 = list(snp = as.character(gwas_hits1$marker), N = 750, pvalues = gwas_hits1$p, MAF = gwas_hits1$freq, type = "quant")
coloc = coloc.abf(dataset2, dataset1)


## Mendelian Randomisation analysis 
   
   ## swap exposure and outcome for bidirectional analysis 


## Open Exposure - replace with relevant datasets and assign column names appropriately 

exposure_dat <- read_exposure_data(
  filename = '/Proteomics/MR/ad_snps.txt',
  sep = ',',
  snp_col = 'SNP',
  beta_col = 'BETA',
  se_col = 'SE',
  effect_allele_col = 'A1',
  phenotype_col = '',
  units_col = '',
  other_allele_col = 'A2',
  eaf_col = '',
  samplesize_col = '',
  ncase_col = '',
  ncontrol_col = '',
  gene_col = '',
  pval_col = 'P'
)

## Clump Exposures

exposure_dat <- clump_data(exposure_dat)

## Open Outcome 

outcome_dat <- read_outcome_data(
  filename = '/Proteomics/MR/pvr1.txt',
  sep = ',',
  snp_col = 'marker',
  beta_col = 'beta',
  se_col = 'stderr',
  effect_allele_col = 'effect_allele',
  phenotype_col = '',
  units_col = '',
  other_allele_col = 'other_allele',
  eaf_col = 'freq',
  samplesize_col = 'n',
  ncase_col = '',
  ncontrol_col = '',
  gene_col = '',
  pval_col = 'p'
)

## Harmonisation of Alleles

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

## Run MR

mr_results <- mr(dat)

