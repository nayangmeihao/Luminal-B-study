


library(TCGAbiolinks)
library(dplyr)
library(limma)
library(biomaRt)
library(SummarizedExperiment)
library(data.table)
library(openxlsx)

sample.data <- read.xlsx("DNA_Methylation/type1_type2_name_data.xlsx")
type1.sample <- sample.data[which(sample.data$label == "type1"), "sample"]
type2.sample <- sample.data[which(sample.data$label == "type2"), "sample"]
group.data <- sample.data[, c("names", "label")]

dataPrep2 <- read.csv("DNA_Methylation/raw_data/dataPrep_data.csv")
dataPrep2[1:5, 1:5]
DNA.data <- dataPrep2
col.names <- read.xlsx("DNA_Methylation/raw_data/colnames.xlsx")
colnames(DNA.data) <- union("id", col.names$sample)
row.names(DNA.data) <- DNA.data[, 1]
DNA.data <- DNA.data[-1]
DNA.data[1:5, 1:5]
write.csv(DNA.data, "DNA_Methylation/diff_analysis/DNA_recolnames_data.csv")

# 匹配col.names和group.data
new.group.data <- left_join(col.names, group.data, by = c("sample" = "names"))
row.names(new.group.data) <- new.group.data[, 1]
new.group.data <- new.group.data[-1]
write.csv(new.group.data, "DNA_Methylation/diff_analysis/new_group_data.csv")

# remove rows with NA
fun <- function(df) {
  drop.rna <- c()
  df.copy <- df
  for (i in 1:nrow(df)) {
    # print(sum(is.na(df[i, ])))
    if (sum(is.na(df[i, ])) > 100) {
      drop.rna <- union(drop.rna, i)
    }
  }
  # print(drop.rna)
  df.copy <- df.copy[-drop.rna, ]
  return(df.copy)
}

dim(DNA.data) # 485577    129
DNA.df <- fun(DNA.data)
dim(DNA.df) # 15112  1183
table(is.na(DNA.df))  

library("impute")
beta <- as.matrix(DNA.df)
beta <- impute.knn(beta) # filled Na by knn
betaData <- beta$data
betaData <- betaData + 0.00001 
write.csv(betaData, "DNA_Methylation/diff_analysis/DNA_recolnames_RemoveNa_data.csv")
write.csv(colnames(betaData), "DNA_Methylation/diff_analysis/DNA_recolnames_RemoveNa_data_colnames.csv")

# It must be ensured that there is a one-to-one correspondence between the methylation signal value matrix and the phenotype information
identical(colnames(betaData), rownames(new.group.data))

library(ChAMP)
library(minfi)
library(methylationArrayAnalysis)
library(wateRmelon)

# beta Signal value matrix cannot have NA values
# Make champ objects
myLoad <- champ.filter(beta = betaData, pd = new.group.data)
myLoad
save(myLoad,file = "DNA_Methylation/diff_analysis/step1-output.Rdata")

######  difference analysis
## ChAMP 
load(file =  "DNA_Methylation/diff_analysis/step1-output.Rdata")
myLoad  # Stores the methylation signal matrix and phenotype information.
myNorm <- champ.norm(beta = myLoad$beta, arraytype = "450K", cores = 5)
dim(myNorm)
group.list <- myLoad$pd
table(group.list)
myDMP <- champ.DMP(beta = myNorm, pheno = group.list)
head(myDMP[[1]])
DMP.GUI()
DMP.GUI(DMP = myDMP[[1]], beta = myNorm, pheno = group.list)  # Automatically detect whether covariates are numeric or sub-types.
save(myDMP, file = "DNA_Methylation/diff_analysis/step2-output-myDMP.Rdata")
write.csv(myDMP[[1]], "DNA_Methylation/diff_analysis/DNA_recolnames_RemoveNa_data_ChAMP_result.csv")

#  minfi 
load(file =  "DNA_Methylation/diff_analysis/step1-output.Rdata")
beta.m <- myLoad$beta
# Follow-up difference analysis for beta.m, such as minfi package
grset <- makeGenomicRatioSetFromMatrix(beta.m, what = "Beta")
M <- getM(grset)
# Because the methylation chip is 450K or 850K, 
# with hundreds of thousands of lines of methylation sites, statistical testing is usually very slow.
dmp <- dmpFinder(M, pheno = group.list, type = "categorical")
dmpDiff <- dmp[(dmp$qval < 0.05) & (is.na(dmp$qval) == F), ]
dim(dmpDiff)

load(file = "DNA_Methylation/diff_analysis/step2-output-myDMP.Rdata")
champDiff <- myDMP[[1]]

dim(dmpDiff)
dim(champDiff)
length(intersect(rownames(dmpDiff), rownames(champDiff)))


### Select the first 1000 differential methylation sites to draw heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)
setwd("D:/study_data/RNA")
# colnames(log2_counts_df) are sample,rownames are genes
log2.counts.df <- read.csv("DNA_Methylation/diff_analysis/DNA_recolnames_RemoveNa_data.csv") 
rownames(log2.counts.df) <- log2.counts.df[, 1]
log2.counts.df <- log2.counts.df[-1]
log2.counts.names.df <- read.csv("DNA_Methylation/diff_analysis/DNA_recolnames_RemoveNa_data_colnames.csv") 
colnames(log2.counts.df) <- log2.counts.names.df$x

gene.df <- read.csv("DNA_Methylation/diff_analysis/DNA_recolnames_RemoveNa_data_ChAMP_result.csv")
gene.df <- gene.df[order(gene.df[, "adj.P.Val"], decreasing = T), ] 
select.genes <- gene.df
heatmap.data <- log2.counts.df[select.genes$X, ]
diff.data <- log2.counts.df[gene.df$X, ]
dim(diff.data)
head(diff.data)
# gene.df$X

all(rownames(diff.data) == gene.df$X)
diff.data$symbol <- gene.df$gene
dim(diff.data)
head(diff.data)



## DNA methylation data visualization
#### my data analysis #####
library(data.table)
library(openxlsx)
library(genoset)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(methyAnalysis)

data1 <- read.csv("DNA_ChAMP_data_join_keygenes_sorted.csv") # 301*132
data1 <- subset(data1, select = -reSymbol)
data1.symbol.data <- data1[, c("X", "symbol")]
data1.data <- subset(data1, select = -(symbol))
rownames(data1.data) <- data1.data[, 1]
data1.data <- data1.data[-1]
# scale data
data1.data.scale <- scale(t(data1.data))
data1.data.scale.new <- ifelse(data1.data.scale > 1, 1, ifelse(data1.data.scale < -1, -1, data1.data.scale))
# data.scale.new <- na.omit(data.scale.new)
data1.data.new <- as.data.frame(t(data1.data.scale.new))
all(colnames(data1.data) == colnames(data1.data.new))
all(rownames(data1.data) == rownames(data1.data.new))

data1.data.new$X <- rownames(data1.data.new) 
data1.new <- left_join(data1.data.new, data1.symbol.data, by = "X")
# rownames(data1.new) <- data1.new$X
# data1.new <- subset(data1.new, select = -X)

sub.key.gene.data <- read.csv("21key_genes.csv") # 21
sub.key.gene <- sub.key.gene.data$Gene
sub.data <- data1.new[which(data1.new$symbol %in% sub.key.gene), ] # 101*130

group.text <- read.csv("new_group_data.csv")
rownames(group.text) <- group.text[, 1]
group.text <- group.text[-1]

# # # test.1
# sub.data.copy <- subset(sub.data, select=-symbol)
# rownames(sub.data.copy) <- sub.data.copy$X
# sub.data.copy <- subset(sub.data.copy, select = -X)
# all(rownames(group.text) == colnames(sub.data.copy)) # TRUE

# Chromosome
Chromosome.data <- read.csv("one_sample_gdc_hg38.csv") # 485577*11
sub.Chromosome.data <- Chromosome.data[, c("ID", "Chromosome", "Start", "End")]

fun1 <- function (rowRanges, exprs = NULL, methylated = NULL, unmethylated = NULL, 
                  detection = NULL, pData = NULL, annotation = "", universe = NULL, 
                  assays = NULL, ...) 
{
  if (is(rowRanges, "RangedData")) 
    rowRanges <- as(rowRanges, "GRanges")
  if (!is.null(universe)) 
    genome(rowRanges) <- universe
  if (!is.null(assays)) {
    if (!all(c("exprs", "methylated", "unmethylated") %in% names(assays))) 
      stop("'exprs', 'methylated' and 'unmethylated are required in assayData!")
    object <- GenoSet(rowRanges = rowRanges, colData = pData, 
                      assays = assays, ...)
  }
  else {
    assays <- list(exprs = exprs, methylated = methylated, unmethylated = unmethylated, detecton = detection)
    # assays <- list(exprs = exprs)
    object <- GenoSet(rowRanges = rowRanges, colData = pData, assays = assays)
  }
  object <- new("MethyGenoSet", object)
  object@annotation <- annotation
  return(object)
}


#### Construct a GenoSet object ####
# array
funAnalysis <- function(gene, sub.data, sub.Chromosome.data, group.text) {
  plot.array <- sub.data[which(sub.data$symbol == gene), ]
  rownames(plot.array) <- plot.array$X
  plot.array <- subset(plot.array, select = -c(X, symbol))
  
  # plot.array <- list(as.matrix(plot.array, dimnames = list(2, 129)))
  dim(plot.array)[1]
  
  plot.array <- as.matrix(plot.array, dimnames = list(dim(plot.array)[1], 129))
  all(rownames(group.text) == colnames(plot.array)) # TRUE
  
  plot.Chr.data <- sub.Chromosome.data[which(sub.Chromosome.data$ID %in% rownames(plot.array)), ]
  plot.rowRanges <- GRanges(ranges = IRanges(start = plot.Chr.data$Start, end = plot.Chr.data$End, names = plot.Chr.data$ID), 
                            seqnames = plot.Chr.data$Chromosome)
  # colData
  plot.colData <- group.text
  # gs <- GenoSet(plot.rowRanges, plot.array, plot.colData)
  gc <- fun1(rowRanges = plot.rowRanges, exprs = plot.array, pData = plot.colData, methylated = plot.array, unmethylated = plot.array, detection = plot.array) 
  return(gc)
}

res.obj <- funAnalysis("BRCA2", sub.data, sub.Chromosome.data, group.text)
res.obj
rowRanges(res.obj)
methylated(res.obj)

pdf("BRCA2_DNA_methylation_heatmap_by_gene_of_selected_Ranges_scaled_hg38.pdf", width = 20, height = 20)
## plot the DNA methylation heatmap by gene of selected GRanges
plotMethylationHeatmapByGene('675', methyGenoSet = res.obj, 
                             phenoData = colData(res.obj), 
                             includeGeneBody = TRUE,
                             gradient = c("navy", "black", "yellow"),
                             phenoColor = list(c("red", "blue")),
                             genomicFeature = 'TxDb.Hsapiens.UCSC.hg38.knownGene')
dev.off()


pdf("FN1_DNA_methylation_heatmap_by_chromosome_scaled_hg38.pdf")
## plot the DNA methylation heatmap by chromosome location   	gene:a Entrez Gene ID
heatmapByChromosome(res.obj, gene = '675', genomicFeature = "TxDb.Hsapiens.UCSC.hg38.knownGene", includeGeneBody = TRUE)
dev.off()


