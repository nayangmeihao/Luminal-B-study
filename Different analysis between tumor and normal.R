library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(DESeq2)
library(data.table)
library(openxlsx)

# Input data 
# DEseq2 analysis between tumor and normal
# Gene expression data，rownames(express) are genes ,colnames are sample
express <- read.xlsx('D:/study_data/tumor_normal_data/tumor_normal_data.xlsx')

# Label datas,rownames(group_text) are sample ,colnames are label
group.text <- read.xlsx('D:/study_data/tumor_normal_data/tumor_normal_label.xlsx')

# data processing
express.rec <- as.data.frame(express)
group.text <- as.data.frame(group.text)

rownames(express.rec) <- express.rec$sample
express.rec <- subset(express.rec, select = -sample)

rownames(group.text) <- group.text$sample
group.text <- subset(group.text, select = -sample)

# Whether the column names of the expression matrix and the row names of the grouping matrix are consistent
all(rownames(group.text) == colnames(express.rec))

label <- group.text$label
table(label)
colData <- group.text

express <- read.xlsx('D:/study_data/tumor_normal_data/tumor_normal_data.xlsx')
dim(express)
head(express)

# DEseq2 analysis
dds <- DESeqDataSetFromMatrix(countData = express.rec, colData = colData, design = ~label, tidy = F)
# Filter out low-expression data (can be used when it is necessary to remove)
table(rowSums(counts(dds) == 0))
# Retain genes that can be found in 70% of samples,351=1176*30%
keep <- rowSums(counts(dds) == 0) < 351
dds <- dds[keep, ]
dds1 <- DESeq(dds)
suppressMessages(dds1)
resultsNames(dds1)

# label_tumor_vs_normal
res <- results(dds1, contrast = c("label", "tumor", "normal"))
# View summary information of the results
summary(res)
plotMA(res, ylim = c(-2,2))

# Correct p
library(apeglm)
res <- lfcShrink(dds1, coef = "label_tumor_vs_normal", type = 'apeglm')
summary(res)
plotMA(res, ylim = c(-2,2))

# Save result
# Sort res by padj (ascending order)
resOrdered <- res[order(res$padj), ]

resOrdered <- as.data.frame(resOrdered)
head(resOrdered)

# # According to padj<0.05, |log2FoldChange|>1, select genes in resOrdered
# resOrdered2 <- resOrdered[resOrdered[, "padj"] < 0.05, ]
# resOrdered3 <- resOrdered2[abs(resOrdered[, "log2FoldChange"]) > 1, ]
# head(resOrdered3)
# write.table(resOrdered3, "D:/study_data/tumor_normal_data/tumor_normal_deseq2_result.csv", row.names = TRUE, col.names = TRUE, sep = ",") 

# According to pvalue<0.05, |log2FoldChange|>1, select genes in resOrdered
resOrdered2 <- resOrdered[resOrdered[, "pvalue"] < 0.05, ]
resOrdered3 <- resOrdered2[abs(resOrdered[, "log2FoldChange"]) > 1, ]
head(resOrdered3)
write.table(resOrdered3, "D:/study_data/tumor_normal_data/tumor_normal_deseq2_pvalue_result.csv", row.names = TRUE, col.names = TRUE, sep = ",") 


# Data Transformation
# this gives log2(n + 1)
ntd <- normTransform(dds1)
log2.counts.df <- as.data.frame(assay(ntd))
write.table(log2.counts.df, "D:/study_data/tumor_normal_data/tumor_normal_log2_counts.csv", row.names = TRUE, col.names = TRUE, sep = ",") 


library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)
data <- read.csv("D:/study_data/tumor_normal_data/tumor_normal_log2_counts.csv")  # rownames(data) are samples,colnames are genes
key.genes.data <- read.csv("D:/study_data/tumor_normal_data/tumor_normal_deseq2_result.csv") 
key.genes <- key.genes.data$gene
y <- "label"
key.genes.union <- union(y, key.genes)

dat <- data[which(data$Gene %in% key.genes.union), ]
rownames(dat) <- dat$Gene
dat <- subset(dat, select = -Gene)
dat.t <- as.data.frame(t(dat))

for(i in 2:ncol(dat.t)) {
  dat.t[, i] <- as.numeric(as.character(dat.t[, i]))
}

dat.t.copy <- dat.t
dat.t.copy$sample <- rownames(dat.t.copy)
head(dat.t.copy)

# tumor
# rownames(data1) are samples
data1 <- dat.t[which(dat.t$label == "tumor"), ]

# cluster
data1 <- subset(data1, select = -1)
hc <- hclust(dist(data1, method = "euclidean"), method = "ward.D2")
new.label <- cutree(hc, k = 3)
new.df.label <- data.frame(new.label)

new.df.label["sample"] <- rownames(new.df.label)
la.t1.df <- new.df.label[order(-new.df.label$new.label), ]

# normal
data2 <- dat.t[which(dat.t$label == "normal"), ]
data1 <- subset(data2, select = -1)
hc <- hclust(dist(data1, method = "euclidean"), method = "ward.D2")
new.label <- cutree(hc, k = 3)
new.df.label <- data.frame(new.label)
new.df.label["sample"] <- rownames(new.df.label)
la.t2.df <- new.df.label[order(new.df.label$new.label), ]

new.label.data <- rbind(la.t1.df, la.t2.df)
m.data <- left_join(new.label.data, dat.t.copy, by = "sample")

# combin la.t1.df and la.t2.df
new.label.data <- rbind(la.t1.df, la.t2.df)
m.data <- left_join(new.label.data, dat.t.copy, by = "sample")
m.data2 <- subset(m.data, select = -new.label)
label.data <- m.data2[, c("sample", "label")]
rownames(label.data) <- label.data$sample
label.data <- subset(label.data, select = -sample)

plot.data <- m.data2
rownames(plot.data) <- plot.data$sample
plot.data <- subset(plot.data, select = -c(sample, label))

all(rownames(label.data) == rownames(plot.data))

log2.data <- plot.data
# scale()
data.scale <- scale(log2.data)
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
# data.scale.new <- na.omit(data.scale.new)
dat <- t(data.scale.new)
all(rownames(label.data) == colnames(dat))
table(is.na(dat))  

#  plot Heatmap
pdf("tumor_normal_heatmap.pdf", width = 9, height = 10)
mat <- as.matrix(dat)
group.mat <-  as.matrix(label.data)
column.ha <- HeatmapAnnotation(df = label.data, 
                               col = list(label = c("tumor" = "#aa2116", "normal" = "#b2d235")))

ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("navy", "black", "yellow")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter",
                                          legend_height = unit(6, "cm"),
                                          legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = FALSE,
              #               right_annotation = row_anno,
              row_title_gp = gpar(fontsize = 20),       
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = FALSE,
              column_split = group.mat,
              column_title = NULL)
ht
dev.off()








