## Install Prerequisite Package and Load Library

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")

library("biomaRt")


## obtain transcription start sites for 15 Olink Inflam EWAS Proteins 

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl, page = "structure")
attributes[grep("transcript", attributes$description, ignore.case = TRUE), ]

## list of proteins 

list = c("CSF1", "TNFSF12", "TNFRSF9", "CXCL10", "FLT3LG", "HGF", "IL18R1", "FGF21", "FGF5", "CXCL8", "TNFRSF11B", "CXCL11", "CXCL9", "OSM", "KITLG", "TGFA")



tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values = c(list),
             mart = ensembl)

## set up code to create biomart dataframe with gene, start and end columns
biomart_df <- tss

names(biomart_df)[8] <- "gene" 
names(biomart_df)[3] <- "start" 
names(biomart_df)[4] <- "end" 
names(biomart_df)[2] <- "chromosome_ensembl" 


## separate data into genes on positive strand and genes on negative strand 

biomart_df_pos <- biomart_df[biomart_df$strand == 1,] 
biomart_df_neg <- biomart_df[biomart_df$strand == -1,]

# sanity check to see if overlap is present/correct number of genes are present 

a = unique(biomart_df_pos$gene) 
b = unique(biomart_df_neg$gene) 
which(a %in% b) 
length(a) 
length(b) 

## create output data frames for genes on positive strand and those on negative strand to be r-bound later 
    
    
    ##out_df1 == genes on positive strand 
    
out_df1 <- matrix(nrow=length(unique(biomart_df_pos$gene)), ncol=2)
colnames(out_df1) <- c("start", "end")
rownames(out_df1) <- unique(biomart_df_pos$gene) 


    ##out_df2 == genes on negative strand 
    
out_df2 <- matrix(nrow=length(unique(biomart_df_neg$gene)), ncol=2)
colnames(out_df2) <- c("start", "end")
rownames(out_df2) <- unique(biomart_df_neg$gene) 

## loop through gene names to fill in min 5'end (start) and max 3'end (end) as on positive strand you want furthest away 5'(tss) and biggest distance to 3' end 
## fill out out_df1 based on the above

for(gene in unique(biomart_df_pos$gene)) {
tmp <- biomart_df_pos[which(biomart_df_pos$gene==gene), ]
out_df1[gene,"start"] <-  min(tmp$start)
out_df1[gene,"end"] <-  max(tmp$end)
}

## loop through gene names to fill in min 3'end (end) and max 5'end (start) as on negative strand you want furthest away 3'(tss (relative to positive strand)) and biggest distance to 5' end 
## fill out out_df2 based on the above

for(gene in unique(biomart_df_neg$gene)) {
tmp <- biomart_df_neg[which(biomart_df_neg$gene==gene), ]
out_df2[gene,"start"] <-  max(tmp$end)
out_df2[gene,"end"] <-  min(tmp$start)
}

## row bind the positive and negative strand dataframes to give one out_df 
out_df <- rbind(out_df1, out_df2)


## Tidy the files before saving for new row with protein names (out_df), then one for length of transcripts in biomart_df for future reference 

    ## set as data frames
    out_df <- as.data.frame(out_df) 
    biomart_df <- as.data.frame(biomart_df) 
    
    
    
out_df$Protein <- row.names(out_df) 
biomart_df$transcript_length <- abs(biomart_df$start - biomart_df$end) 

## save files 

write.csv(biomart_df, "/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Biomart_Output_Olink_Inflam.csv", row.names = F) 

write.csv(out_df, "/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/TSS_Olink_Inflamm_EWAS.csv", row.names = F) 
