## Mixed Model Method 

 ## Preparation for OSCA run in Terminal 
 
 ## (run in R)
 osca_dat <- data.frame(IID=colnames(dat), FID=colnames(dat))
 osca_dat <- cbind(osca_dat, t(dat))
 write.table(osca_dat, file="myprofile.txt", row.names=F, sep=' ')


## #Youll have to manually make the .opi annotation file:
## Manual Annotation 
 opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
 opi$chr <- gsub("chr", "", opi$chr)
 opi$chr <- gsub("X", "23", opi$chr)
 opi$chr <- gsub("Y", "24", opi$chr)
 opi$chr <- as.numeric(opi$chr)
 opi$UCSC_RefGene_Name	 <- as.factor(opi$UCSC_RefGene_Name	)
 opi$strand <- as.factor(opi$strand)
 opi[which(opi$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA
 write.table(opi, file="output_filenam.opi", 
 	             col.names=F, 
 	             row.names=F, 
 	             quote=F, sep='\t')


## Make phenotype file. 
 ## 3 columns: IID, FID, Phenotype
 pheno <- data.frame(FID = phenos$Basename,
 	                IID = phenos$Basename,
 	                simd = scale(phenos$ pheno))

 ## Make quantitative covariates file (e.g. age, pack years, cellcounts etc) 
 ## First two columns are IID and FID as before

 quant_cov <- data.frame(FID = phenos$Barcode, 
 	                    IID = phenos$Barcode,
 	                    CD8T = phenos$CD8T, 
 	                    CD4T = phenos$CD4T,
 	                    Mono = phenos$Mono,
 	                    Bcell = phenos$Bcell,
 	                    Gran = phenos$Gran,
 	                    Age = phenos$Age,
 	                    pack_years = phenos$pack_years)
 write.table(quant_cov, file="quant.cov", row.names=F, sep=' ')


# Make factor covariates file (e.g. smoking status, sex, batch etc) 
# First two columns are IID and FID as before
 fact_cov <- data.frame(FID =phenos$Barcode,
 	                   IID = phenos$Barcode,
 	                   ever_smoke = phenos$ever_smoke,
 	                   batch = phenos$plate_processing_batch)
 write.table(fact_cov, file="factors.cov", row.names=F, sep=' ')


## Run OSCA in terminal ##
# Create binary methylation data --methylation-beta if beta values, --methylation-m if M-values

system("osca_Linux --efile myprofile.txt --methylation-beta --make-bod --out output_filename")

osca_Linux \
	   --moment \
 	   --befile ../output_filename \
 	   --pheno ../simd.phen \
 	   --qcovar ../quant.cov \
 	   --covar ../factors.cov \
 	   --fast-linear \
 	   --out osca_simd_ewas \
 	   --methylation-m
 	   
 	   
 	   
### Extraction of Top Results ### in R 

setwd("/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/OSCA Inflammatory/") 

L <- list.files(".", ".mlma")

pip_files <- lapply(L, read.csv, header = T) 
names <- as.character(L) 
names <- gsub("_.*", "", names) 
names(pip_files) <- names 
pip_files = Map(cbind, pip_files, "Biomarker"=names) 

osca_top <- do.call(rbind, lapply(osca_files, function(x)x[x$p < 1e-7,]))

