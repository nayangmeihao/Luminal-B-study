# download TCGA_BRCA
rm(list = ls())
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(EDASeq)
library(SummarizedExperiment)

clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
write.csv(clinical,file = paste("01_01_TCGA_BRCA", "clinical.csv",sep = "-"))

query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

GDCdownload(query = query, 
            files.per.chunk=6, 
            method ="api", 
            directory ="BRCA_cancer")

dataPrep <- GDCprepare(query = query, save = TRUE, 
                       save.filename = "BRCA_cases.rda", 
                       directory ="BRCA_cancer")

count_matrix <- assay(dataPrep)
write.csv(count_matrix,file ="TCGA_BRCA_all_Counts.csv")

dataPrep1 <- TCGAanalyze_Preprocessing(object = dataPrep, 
                                       cor.cut = 0.6,
                                       datatype = "HTSeq - Counts",
                                       filename = "AAIC.png")


## ID transformation of RNA-seq data using genecode 
library(tidyr)
library(dplyr)
library(rtracklayer)
data <- read.csv("04_01_TCGA_LIHC_Filt_Counts.csv", head=FALSE) #42144*421
colnames(data) <- data[1,]
data <- data[-1,]

names(data)[names(data) == ""] <- "Ensembl_ID"

gtf <- rtracklayer::import("Genecode/gencode.v38.annotation.gtf.gz")
gtf <- as.data.frame(gtf)
dim(gtf)

a <- gtf[which(gtf$gene_type == "protein_coding" ), ]
b <- dplyr::select(a,c(gene_name,gene_id,gene_type))
library(stringr)
b$Ensembl_ID <- unlist(str_split(b$gene_id,"[.]",simplify=T))[,1]
join.data <- dplyr::inner_join(b, data, by ="Ensembl_ID")

join.dat <- select(join.data, -gene_id, -gene_type, -Ensembl_ID)
mRNAdata <- distinct(join.dat, gene_name,.keep_all = T)
rownames(mRNAdata)<- mRNAdata[,1]
mRNAdata <- mRNAdata[,-1] 

write.csv(mRNAdata,"Genecode/04_01_TCGA_BRCA_filt_mRNA_data.csv",quote = F,row.names = T)
save(mRNAdata,file = "Genecode/04_01_TCGA_BRCA_filt_mRNA_data.Rda")



library(TCGAbiolinks)
library(plyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)
library(DT)
library(data.table)
library(openxlsx)
sample.data <- read.xlsx("type1_type2_name_data.xlsx")

type1.sample <- sample.data[which(sample.data$label == "type1"), "sample"]
type2.sample <- sample.data[which(sample.data$label == "type2"), "sample"]


query <- GDCquery(project = "TCGA-BRCA", # Cancer type
                  data.category = "DNA Methylation",
                  data.type = "Methylation Beta Value",
                  platform = "Illumina Human Methylation 450",
                  barcode = c(type1.sample, type2.sample)
)

GDCdownload(query,
            method = "api",
            directory = "BRCAdata",
            files.per.chunk = 6)

# dataPrep1 <- GDCprepare(query = query, save = TRUE, save.filename = "BRCA_case.rda")
dataPrep1 <- GDCprepare(query = query)
dataPrep <- assay(dataPrep1)
save(dataPrep, file = "dataPrep1.Rdata")
write.csv(dataPrep, "dataPrep_data.csv")

dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1,
                                       cor.cut = 0.6,
                                       datatype ="Methylation Beta Value")

write.csv(dataPrep2, file = "BRCA_dataPrep.csv", quote = FALSE)



