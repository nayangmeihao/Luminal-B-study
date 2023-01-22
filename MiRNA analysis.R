
# miRNA analysis
library(data.table)
library(openxlsx)

# Input data 
# DEseq2 analysis between tumor and normal
# Gene expression data，rownames(raw.data) are genes ,colnames are sample
raw.data <- read.csv("micRNA_data_drop_No_same_rename.csv")
rownames(raw.data) <- raw.data$miRNA_ID
raw.data <- subset(raw.data, select = -miRNA_ID)
raw.data.t <- as.data.frame(t(raw.data))
head(raw.data.t)

tumor.names.data <- read.xlsx("study_dat/tumor_names.xlsx")
normal.names.data <- read.xlsx("study_dat/normal_names.xlsx")

length(tumor.names.data$sample)
length(normal.names.data$sample)

# get tumor data
raw.data.t$sample <- rownames(raw.data.t)
tumor.data <- raw.data.t[which(raw.data.t$sample %in% tumor.names.data$sample ), ]
tumor.data <- subset(tumor.data, select = -sample)
tumor.data$label <- "tumor"

#get normal data
normal.data <- raw.data.t[which(raw.data.t$sample %in% normal.names.data$sample ), ]
normal.data <- subset(normal.data, select = -sample)
normal.data$label <- "normal"

all(colnames(normal.data)==colnames(tumor.data)) # TRUE

write.table(tumor.data, "study_dat2/tumor_data_data.csv", row.names = TRUE, col.names = TRUE, sep = ",") 
write.table(normal.data, "study_dat2/normal_data_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")


## difference analysis by DEseq2
library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(DESeq2)
library(data.table)
library(openxlsx)

tumor.data <- read.csv("study_dat2/tumor_data_data.csv")
normal.data <- read.csv("study_dat2/normal_data_data.csv")

all(colnames(normal.data)==colnames(tumor.data))

DEseq2.data <- rbind(tumor.data, normal.data)
DEseq2.data$sample <- rownames(DEseq2.data)

express.dat <- subset(DEseq2.data, select = -c(sample, label))
express.rec <- as.data.frame(t(express.dat)) # colnames(express.rec) are sample, rownames(express.rec) are genes

colData <- DEseq2.data[, c("sample", "label")]
colData <- subset(colData, select = -sample)# colnames(colData) are label, rownames(colData) are samples

label <- colData$label
all(rownames(colData) == colnames(express.rec))

dds <- DESeqDataSetFromMatrix(countData = express.rec, colData = colData, design = ~label, tidy = F)

# Retain the genes that can be found in 30% of samples, 348=1163*30%
keep <- rowSums(counts(dds) == 0) < 348
dds <- dds[keep, ]
counts(dds)[1:10, 1:3]
dim(dds) 

# Difference analysis
dds <- DESeq(dds)
suppressMessages(dds)
resultsNames(dds)

# 'label_tumor_vs_normal'
res <- results(dds, contrast = c("label", "tumor", "normal"))
summary(res)
plotMA(res, ylim = c(-2, 2))

# # shrink log2 fold change by lfcShrink
# library(apeglm)
# res <- lfcShrink(dds, coef = "label_tumor_vs_normal", type = "apeglm")
# summary(res)
# plotMA(res, ylim = c(-2,2))

# save result
resOrdered <- as.data.frame(res)
# In resOrdered, select genes by padj<0.05, |log2FoldChange|>1
resOrdered2 <- resOrdered[(abs(resOrdered[, "log2FoldChange"]) > 1)&(resOrdered[, "padj"] < 0.05), ]
dim(resOrdered2)
resOrdered2
# write.table(resOrdered, "study_dat2/tumor_normal_deseq2_All_result.csv", row.names = TRUE, col.names = TRUE, sep = ",")
# write.table(resOrdered2, "study_dat2/tumor_normal_deseq2_After_screening_result.csv", row.names = TRUE, col.names = TRUE, sep = ",") 


# Data Transformation
# this gives log2(n + 1)
ntd <- normTransform(dds)
log2.df <- as.data.frame(assay(ntd)) 
log2.dat <- log2.df[rownames(resOrdered2), ]
tumor.deseq2.dat <- log2.dat[, rownames(tumor.data)]
normal.deseq2.dat <- log2.dat[, rownames(normal.data)]
write.table(tumor.deseq2.dat, "study_dat2/tumor_deseq2_log2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")
write.table(normal.deseq2.dat, "study_dat2/normal_deseq2_log2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")


## difference analysis between type1 and type2
library(data.table)
library(openxlsx)
library(dplyr)

setwd("D:/study_data/RNA/micRNA")
tumor.data <- read.csv("study_dat2/tumor_data_data.csv")
type1.type2.name.data <- read.xlsx("study_dat/type1_type2_name_data.xlsx")

sub.tumor.data <- tumor.data[, union("label", rownames(resOrdered2))]
sub.tumor.data$sample <- rownames(sub.tumor.data)
dim(sub.tumor.data) # 1073*157
head(sub.tumor.data)

deseq2.data <- inner_join(type1.type2.name.data, sub.tumor.data, by = "sample")
rownames(deseq2.data) <- deseq2.data$sample
dim(deseq2.data) # 181 161
deseq2.data <- subset(deseq2.data, select = -c(sample_names, raw_sample_names, sample1, label.y))
names(deseq2.data)[names(deseq2.data) == "label.x"] <- "label"
# dim(deseq2.data) # 181*157
# head(deseq2.data)
write.table(deseq2.data, "study_dat2/B_counts_data.csv", row.names = TRUE, col.names = TRUE, sep = ",") 

table(deseq2.data$label)
# head(deseq2.data)

### diifference analysis between type1 and type2 by DEseq2 
library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(DESeq2)
library(data.table)
library(openxlsx)

deseq2.data <- read.csv("study_dat2/B_counts_data.csv")
express.dat <- subset(deseq2.data, select = -c(sample, label))
express.rec <- as.data.frame(t(express.dat))
# dim(express.rec) # 155*181
# head(express.rec)

colData <- deseq2.data[, c("sample", "label")]
colData <- subset(colData, select = -sample)
# dim(colData) # 181*1
label <- colData$label
all(rownames(colData) == colnames(express.rec))

dds <- DESeqDataSetFromMatrix(countData = express.rec, colData = colData, design = ~label, tidy = F)
# Retain the genes that can be found in 30% of samples, 54=181*30%
keep <- rowSums(counts(dds) == 0) < 54
dds <- dds[keep, ]
counts(dds)[1:10, 1:3]
dim(dds) 

# Difference analysis
dds <- DESeq(dds)
suppressMessages(dds)
resultsNames(dds)

# "label_type2_vs_type1"
res <- results(dds, contrast = c("label", "type2", "type1"))
summary(res)
plotMA(res, ylim = c(-2, 2))

# save result
resOrdered <- as.data.frame(res)

# In resOrdered, select genes by padj<0.05, |log2FoldChange|>1
resOrdered2 <- resOrdered[(abs(resOrdered[, "log2FoldChange"]) > 1)&(resOrdered[, "padj"] < 0.05), ]
dim(resOrdered2)

# Data Transformation
# this gives log2(n + 1)
ntd <- normTransform(dds)
log2.df <- as.data.frame(assay(ntd)) 
log2.dat <- log2.df[rownames(resOrdered2), ]

coldata <- colData
coldata$label <- as.character(coldata$label)
coldata$sample <- rownames(coldata)
type1.data <- coldata[which(coldata$label == "type1"), ]
type2.data <- coldata[which(coldata$label == "type2"), ]

type1.deseq2.dat <- log2.dat[, rownames(type1.data)]
type2.deseq2.dat <- log2.dat[, rownames(type2.data)]
write.table(type1.deseq2.dat, "study_dat2/type1_deseq2_log2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")
write.table(type2.deseq2.dat, "study_dat2/type2_deseq2_log2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")

# Data Transformation
# this gives log2(n + 1)
ntd <- normTransform(dds)
log2.df <- as.data.frame(assay(ntd))  
head(log2.df)
write.table(log2.df, "study_dat2/B_log2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",") 


library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)
library(dplyr)

log2.counts.df <- read.csv("B_log2_data.csv")  # 696*188
dim(log2.counts.df) # 147*182
head(log2.counts.df)

miRNA.name.data <- read.csv("type1_type2_deseq2_After_screening_result.csv")
miRNA.data <- log2.counts.df[which(log2.counts.df$micRNA %in% miRNA.name.data$micRNA), ]
dim(miRNA.data) #  9*182
rownames(miRNA.data) <- miRNA.data$micRNA
miRNA.data <- subset(miRNA.data, select = -micRNA)
write.table(miRNA.data, "miRNA.data_log2_data.csv",row.names = TRUE, col.names = TRUE, sep = ",") 

# Add label
miRNA.t.data <- as.data.frame(t(miRNA.data))
miRNA.t.data$sample <- rownames(miRNA.t.data)
type1.type2.name.data <- read.xlsx("type1_type2_name_data.xlsx")

miRNA.heatmap.data <- left_join(type1.type2.name.data, miRNA.t.data, by = c("sample" = "sample"))
rownames(miRNA.heatmap.data) <- miRNA.heatmap.data$sample
miRNA.heatmap.data <- na.omit(miRNA.heatmap.data)
miRNA.heatmap.data <- subset(miRNA.heatmap.data, select = -c(sample_names, sample1, raw_sample_names))
head(miRNA.heatmap.data)

heatmap.data <- miRNA.heatmap.data
# heatmap.data <- lncRNA.heatmap.data

# label
group.text <- heatmap.data[, c("sample","label")]
rownames(group.text) <- group.text$sample
group.text <- subset(group.text, select = -sample)
# head(group.text)

# expre data
expre.data <- subset(heatmap.data, select = -c(sample, label))
all(rownames(group.text) == rownames(expre.data))

# data scale
for (i in colnames(expre.data)) {
  expre.data[, i] <- as.numeric(as.character(expre.data[, i]))
}
all(rownames(group.text) == rownames(expre.data))
data.scale <- scale(expre.data)
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
# data.scale.new <- na.omit(data.scale.new)
dat <- t(data.scale.new)
all(rownames(group.text) == colnames(dat))

# plot Heatmap
mat <- as.matrix(dat)
group.mat <-  as.matrix(group.text)
all(rownames(group.mat) == colnames(dat))

pdf("miRNA_heatmap.pdf", width = 9, height = 9)
column.ha <- HeatmapAnnotation(df = group.text, 
                               col = list(label = c("type1" = "red", "type2" = "blue")))
ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#84C1FF", "black", "#decb00")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter",
                                          legend_height = unit(6, "cm"),
                                          legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = TRUE,
              #right_annotation = row_anno,
              row_title_gp = gpar(fontsize = 20),       
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = T,
              column_split = group.mat,
              column_title = NULL)
ht
dev.off()



## pathway analysis
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)

# micRNA enrichment analysis 
Enrichment.barplot2 <- function(data, gene_sets = NULL, title = NULL){
  enrich <- data
  enrich$FDR_q <- -log10(enrich[, "FDR.q.value"])
  # Sort the enrichment results in ascending order according to qvalue to ensure that the most significant pathway is on the top
  enrich1 <- enrich[order(enrich$FDR_q, decreasing = T), ]   
  print(max(enrich1$FDR_q))
  print(length(enrich1$FDR_q))
  enrich3 <- data.frame(enrich1[, "Gene.Set.Name"], enrich1[, "FDR_q"])
  colnames(enrich3) <- c("ID", "FDR_q")
  p <- ggplot(data = enrich3, aes(x = ID, y = FDR_q))   # fill = qvalue  
  p1 <- p + geom_bar(stat = "identity", width = 0.5) + coord_flip()            #coord_flip()  Reverse axis
  p2 <- p1 + theme(panel.background = element_rect(fill = "transparent", color = "gray"),
                   axis.text.y = element_text(color = "black", size = 12))
  # ylim(0,200) Change the range of abscissa. The coordinate axis here is reversed. 
  # Although it looks like the x-axis, it is actually the y-axis.
  p3 <- p2 + ylim(0, 4) + scale_fill_gradient(low = "red", high = "blue")   
  p4 <- p3 + scale_x_discrete(limits = rev(enrich3[,1])) +labs(x = "", y = "-lg(FDR_q_value)", 
                                                               title = "Enrichment Analysis of hsa-miR-6510")
  p5 <- p4+ theme(plot.margin = unit(rep(1,4), "lines"))
  p5
}


hsa.miR.144.data <- read.csv("hsa_mir_144_hallmark_enrichment_analysis.csv")  # 1*7
dim(hsa.miR.144.data)
head(hsa.miR.144.data)
Enrichment.barplot2(hsa.miR.144.data)  # ylim(0, 6)




## correlation analysis between micRNA and mRNA
# Organize micRNA target genes data
library(data.table)
library(openxlsx)
log2.mRNA.data <- read.csv("tumor_log2_without_label_data.csv")
all.genes <- as.character(log2.mRNA.data[, 1])
length(all.genes) # 16875
rownames(log2.mRNA.data) <- log2.mRNA.data[, 1]
log2.mRNA.data <- log2.mRNA.data[-1]
log2.mRNA.data.t <- as.data.frame(t(log2.mRNA.data)) # colnames(log2.mRNA.data.t) are genes, rownames(log2.mRNA.data.t) are genes

log2.micRNA.data <- read.csv("B_log2_data.csv")
log2.micRNA.data.t <- as.data.frame(t(log2.micRNA.data)) # colnames(log2.micRNA.data.t) are genes, rownames(log2.micRNA.data.t) are samplea
inter.sample <- intersect(rownames(log2.micRNA.data.t), rownames(log2.mRNA.data.t))
length(inter.sample) # 181

type.gene.data <- read.csv("type1_type2_deseq2_result.csv")
type.genes <- type.gene.data$genes
length(type.genes) #537

###  hsa.mir.144  ###
hsa.mir.144.gene.data.3p <- read.xlsx("hsa-mir-144/hsa-miR-144-3p/miRDB.xlsx")
hsa.mir.144.genes.3p <- hsa.mir.144.gene.data.3p[which(hsa.mir.144.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.144.genes.3p)

hsa.mir.144.gene.data.5p <- read.xlsx("hsa-mir-144/hsa-miR-144-5p/miRDB.xlsx")
hsa.mir.144.genes.5p <- hsa.mir.144.gene.data.5p[which(hsa.mir.144.gene.data.5p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.144.genes.5p)

hsa.mir.144.genes <- union(hsa.mir.144.genes.3p, hsa.mir.144.genes.5p) # 531
write.csv(hsa.mir.144.genes, "hsa_mir_144_target_mRNA_select_result.csv")
length(hsa.mir.144.genes) # 531

inter.all.mir144 <- intersect(hsa.mir.144.genes, all.genes) # 298
length(inter.all.mir144)

inter.type.mir144 <- intersect(hsa.mir.144.genes, type.genes) # 10
length(inter.type.mir144)
write.csv(inter.type.mir144, "inter_type_mir144_result.csv")


###  hsa.mir.147b  ###
hsa.mir.147b.gene.data.3p <- read.xlsx("hsa-mir-147b/hsa-miR-147b-3p/miRDB.xlsx")
hsa.mir.147b.genes.3p <- hsa.mir.147b.gene.data.3p[which(hsa.mir.147b.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.147b.genes.3p) 
hsa.mir.147b.genes <- hsa.mir.147b.genes.3p
write.csv(hsa.mir.147b.genes, "hsa_mir_147b_target_mRNA_select_result.csv")

inter.all.mir147b <- intersect(hsa.mir.147b.genes, all.genes)
length(inter.all.mir147b) # 5

inter.type.mir147b <- intersect(hsa.mir.147b.genes, type.genes)
length(inter.type.mir147b) # 0
write.csv(inter.type.mir147b, "inter_type_mir147b_result.csv")


### hsa.mir.184 ###
hsa.mir.184.gene.data <- read.xlsx("hsa-mir-184/miRDB.xlsx")
hsa.mir.184.genes <- hsa.mir.184.gene.data[which(hsa.mir.184.gene.data$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.184.genes) # 8

inter.all.mir184 <- intersect(hsa.mir.184.genes, all.genes)
length(inter.all.mir184) # 7
inter.type.mir184 <- intersect(hsa.mir.184.genes, type.genes)
length(inter.type.mir184) # 0


### hsa.mir.224 ###
hsa.mir.224.gene.data.3p <- read.xlsx("hsa-mir-224/hsa-mir-224-3p/miRDB.xlsx")
hsa.mir.224.genes.3p <- hsa.mir.224.gene.data.3p[which(hsa.mir.224.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.224.genes.3p)

hsa.mir.224.gene.data.5p <- read.xlsx("hsa-mir-224/hsa-mir-224-5p/miRDB.xlsx")
hsa.mir.224.genes.5p <- hsa.mir.224.gene.data.5p[which(hsa.mir.224.gene.data.5p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.224.genes.5p)

hsa.mir.224.genes <- union(hsa.mir.224.genes.3p, hsa.mir.224.genes.5p)

inter.all.mir224 <- intersect(hsa.mir.224.genes, all.genes)
length(inter.all.mir224) # 336
inter.type.mir224 <- intersect(hsa.mir.224.genes, type.genes)
length(inter.type.mir224) # 10


### hsa.mir.486.1 ###
hsa.mir.486.1.gene.data.3p <- read.xlsx("hsa-mir-486-1/hsa-miR-486-3p/miRDB.xlsx")
hsa.mir.486.1.genes.3p <- hsa.mir.486.1.gene.data.3p[which(hsa.mir.486.1.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.486.1.genes.3p)

hsa.mir.486.1.gene.data.5p <- read.xlsx("hsa-mir-486-1/hsa-miR-486-5p/miRDB.xlsx")
hsa.mir.486.1.genes.5p <- hsa.mir.486.1.gene.data.5p[which(hsa.mir.486.1.gene.data.5p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.486.1.genes.5p)

hsa.mir.486.1.genes <- union(hsa.mir.486.1.genes.3p, hsa.mir.486.1.genes.5p)

inter.all.mir486.1 <- intersect(hsa.mir.486.1.genes, all.genes)
length(inter.all.mir486.1) # 239
inter.type.mir486.1 <- intersect(hsa.mir.486.1.genes, type.genes)
length(inter.type.mir486.1) # 8


### hsa.mir.486.2 ###
hsa.mir.486.2.gene.data.3p <- read.xlsx("hsa-mir-486-2/hsa-miR-486-3p/miRDB.xlsx")
hsa.mir.486.2.genes.3p <- hsa.mir.486.2.gene.data.3p[which(hsa.mir.486.2.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.486.2.genes.3p)

hsa.mir.486.2.gene.data.5p <- read.xlsx("hsa-mir-486-2/hsa-miR-486-5p/miRDB.xlsx")
hsa.mir.486.2.genes.5p <- hsa.mir.486.2.gene.data.5p[which(hsa.mir.486.2.gene.data.5p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.486.2.genes.5p)

hsa.mir.486.2.genes <- union(hsa.mir.486.2.genes.3p, hsa.mir.486.2.genes.5p)

inter.all.mir486.2 <- intersect(hsa.mir.486.2.genes, all.genes)
length(inter.all.mir486.2) # 239
inter.type.mir486.2 <- intersect(hsa.mir.486.2.genes, type.genes)
length(inter.type.mir486.2) # 8


### hsa.mir.3065 ###
hsa.mir.3065.gene.data.3p <- read.xlsx("hsa-mir-3065/hsa-miR-3065-3p/miRDB.xlsx")
hsa.mir.3065.genes.3p <- hsa.mir.3065.gene.data.3p[which(hsa.mir.3065.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.3065.genes.3p)

hsa.mir.3065.gene.data.5p <- read.xlsx("hsa-mir-3065/hsa-miR-3065-5p/miRDB.xlsx")
hsa.mir.3065.genes.5p <- hsa.mir.3065.gene.data.5p[which(hsa.mir.3065.gene.data.5p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.3065.genes.5p)

hsa.mir.3065.genes <- union(hsa.mir.3065.genes.3p, hsa.mir.3065.genes.5p)

inter.all.mir3065 <- intersect(hsa.mir.3065.genes, all.genes)
length(inter.all.mir3065) # 656
inter.type.mir3065 <- intersect(hsa.mir.3065.genes, type.genes)
length(inter.type.mir3065) # 17


### hsa.mir.4728 ###
hsa.mir.4728.gene.data.3p <- read.xlsx("hsa-mir-4728/hsa-miR-4728-3p/miRDB.xlsx")
hsa.mir.4728.genes.3p <- hsa.mir.4728.gene.data.3p[which(hsa.mir.4728.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.4728.genes.3p)

hsa.mir.4728.gene.data.5p <- read.xlsx("hsa-mir-4728/hsa-miR-4728-5p/miRDB.xlsx")
hsa.mir.4728.genes.5p <- hsa.mir.4728.gene.data.5p[which(hsa.mir.4728.gene.data.5p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.4728.genes.5p)

hsa.mir.4728.genes <- union(hsa.mir.4728.genes.3p, hsa.mir.4728.genes.5p)

inter.all.mir4728 <- intersect(hsa.mir.4728.genes, all.genes)
length(inter.all.mir4728) # 638
inter.type.mir4728 <- intersect(hsa.mir.4728.genes, type.genes)
length(inter.type.mir4728) # 17


### hsa.mir.6510 ###
hsa.mir.6510.gene.data.3p <- read.xlsx("hsa-mir-6510/hsa-miR-6510-3p/miRDB.xlsx")
hsa.mir.6510.genes.3p <- hsa.mir.6510.gene.data.3p[which(hsa.mir.6510.gene.data.3p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.6510.genes.3p)

hsa.mir.6510.gene.data.5p <- read.xlsx("hsa-mir-6510/hsa-miR-6510-5p/miRDB.xlsx")
hsa.mir.6510.genes.5p <- hsa.mir.6510.gene.data.5p[which(hsa.mir.6510.gene.data.5p$Target.Score > 80), ]$Gene.Symbol
length(hsa.mir.6510.genes.5p)

hsa.mir.6510.genes <- union(hsa.mir.6510.genes.3p, hsa.mir.6510.genes.5p)

inter.all.mir6510 <- intersect(hsa.mir.6510.genes, all.genes)
length(inter.all.mir6510) # 82
inter.type.mir6510 <- intersect(hsa.mir.6510.genes, type.genes)
length(inter.type.mir6510) # 1


hsa.mir.union.genes <- Reduce(union, list(hsa.mir.144.genes, hsa.mir.147b.genes, 
                                          hsa.mir.184.genes, hsa.mir.224.genes, 
                                          hsa.mir.486.1.genes, hsa.mir.486.2.genes, 
                                          hsa.mir.3065.genes, hsa.mir.4728.genes, 
                                          hsa.mir.6510.genes))
length(hsa.mir.union.genes) # 2318


union.genes <- Reduce(union, list(inter.type.mir144, inter.type.mir147b, 
                                  inter.type.mir184, inter.type.mir224, 
                                  inter.type.mir486.1, inter.type.mir486.2, 
                                  inter.type.mir3065, inter.type.mir4728, 
                                  inter.type.mir6510))
length(union.genes) # 62


# correlation analysis
re <- c()
res <- c()
fun <- function(micRNA, genes) {
  dat1 <- log2.micRNA.data.t[inter.sample, which(colnames(log2.micRNA.data.t) == micRNA)]
  for (j in genes) {
    dat2 <- log2.mRNA.data.t[inter.sample, which(colnames(log2.mRNA.data.t) == j)]
    correlation <- cor.test(dat1, dat2, alternative = "two.side", method = "spearman", conf.level = 0.95, exact = FALSE)
    re <- cbind(micRNA, j, correlation$estimate, correlation$p.value)
    res <- rbind(res, re)
  }
  res <- as.data.frame(res)
  colnames(res) <- c("micRNA", "mRNA", "correlation", "p.value")
  res$p.value <- as.numeric(as.character(res$p.value))
  res$adj.pval <- p.adjust(res$p.value, method = 'BH')
  res <- res[order(res$adj.pval), ]
  return(res)
}

hsa.mir.144.res <- fun("hsa-mir-144", inter.all.mir144)

hsa.mir.147b.res <- fun("hsa-mir-147b", inter.all.mir147b)

hsa.mir.184.res <- fun("hsa-mir-184", inter.all.mir184)

hsa.mir.224.res <- fun("hsa-mir-224", inter.all.mir224)

hsa.mir.486.1.res <- fun("hsa-mir-486-1", inter.all.mir486.1)

hsa.mir.486.2.res <- fun("hsa-mir-486-2", inter.all.mir486.2)

hsa.mir.3065.res <- fun("hsa-mir-3065", inter.all.mir3065)

hsa.mir.4728.res <- fun("hsa-mir-4728", inter.all.mir4728)

hsa.mir.6510.res <- fun("hsa-mir-6510", inter.all.mir6510)

max(as.numeric(as.character(hsa.mir.486.1.res$correlation)))
min(as.numeric(as.character(hsa.mir.486.1.res$correlation)))


# correlation analysis
fun2 <- function(micRNA, inter.gene) {
  res <- c()
  dat1 <- log2.micRNA.data.t[inter.sample, micRNA]
  for (j in inter.gene) {
    re <- c()
    dat2 <- log2.mRNA.data.t[inter.sample, j]
    correlation <- cor.test(dat1, dat2, alternative = "two.side", method = "spearman", conf.level = 0.95, exact = FALSE)
    re <- cbind(micRNA, j, correlation$estimate, correlation$p.value, p.adjust(correlation$p.value, method = "BH"))
    res <- rbind(res, re)
  }
  res <- as.data.frame(res)
  colnames(res) <- c("micRNA", "mRNA", "cor", "pvalue", "padj")
  for (t in 3:ncol(res)) {
    res[, t] <- as.numeric(as.character(res[, t]))
  }
  return(res)
}

hsa.mir.144.dat <- fun2("hsa-mir-144", inter.type.mir144)
head(hsa.mir.144.dat)

# hsa.mir.147b.dat <- fun2("hsa-mir-184", inter.type.mir147b) # NULL

# hsa.mir.184.dat <- fun2("hsa-mir-184", inter.type.mir184) # NULL

hsa.mir.224.dat <- fun2("hsa-mir-224", inter.type.mir224)
head(hsa.mir.224.dat)

hsa.mir.486.1.dat <- fun2("hsa-mir-486-1", inter.type.mir486.1)
head(hsa.mir.486.1.dat)

hsa.mir.486.2.dat <- fun2("hsa-mir-486-2", inter.type.mir486.2)
head(hsa.mir.486.2.dat)

hsa.mir.3065.dat <- fun2("hsa-mir-3065", inter.type.mir3065)
head(hsa.mir.3065.dat)

hsa.mir.4728.dat <- fun2("hsa-mir-4728", inter.type.mir4728)
head(hsa.mir.4728.dat)

hsa.mir.6510.dat <- fun2("hsa-mir-6510", inter.type.mir6510)  # NULL
head(hsa.mir.6510.dat)  

union.cor <- Reduce(rbind, list(hsa.mir.144.dat, hsa.mir.224.dat, 
                                hsa.mir.486.1.dat, hsa.mir.486.1.dat, 
                                hsa.mir.486.2.dat, hsa.mir.3065.dat,
                                hsa.mir.4728.dat, hsa.mir.6510.dat))
write.csv(union.cor, "union_diff_micRNA_cor_diff_genes_result.csv")

union.cor <- read.csv("union_diff_micRNA_cor_diff_genes_result.csv")


for (i in 1:nrow(union.cor)) {
  if (union.cor[i, "padj"] > 0.05) {
    union.cor[i, "padj"] <- NA
  }
}

union.cor <- union.cor[complete.cases(union.cor[,5]),]  # Delete rows with empty values in the fifth column
union.cor[, "-lg(Pvalue)"] <- -log(10, union.cor[, "pvalue"])
require(ggplot2)
pdf("cor_plot.pdf", width = 12, height = 5)
p1 <-  ggplot(union.cor, aes(y = micRNA,  x = mRNA)) + scale_fill_continuous(low = "LightSteelBlue", high = "LightSlateGray")
P2 <- p1 + geom_point(aes(colour = union.cor[, 4], size = union.cor[, 5])) + 
  scale_color_gradient(low = "yellow", high = "DarkOrange") + 
  scale_size(range = c(1, 10)) + theme_bw()
p3 <- P2 + theme(axis.text.x = element_text(size = 8, vjust = 0.5, hjust = 0.5, angle = 30), legend.position = "right")
p3 + labs( colour = "Correlation",size = "-lg(Pvalue)") + guides(size = guide_legend(override.aes = list(shape = 1)))
dev.off()




