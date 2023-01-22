
# LncRNA (tumor ,normal) Difference analysis by DEseq2
library(limma)
library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(DESeq2)
library(data.table)
library(openxlsx)
library(dplyr)

# Gene expression data，rownames(express) are genes ,colnames are sample
express <- read.csv("raw_data/LncRNA_count_data.csv")
sample.label.data <- read.xlsx("tumor_normal_deseq2/tumor_normal_names.xlsx")
LncRNA.data <- express[, union("Ensembl_ID", sample.label.data$sample)]

gene.data <- express[, c("Ensembl_ID", "Symbol")]
gene.data$Ensembl_ID <- as.character(gene.data$Ensembl_ID)
gene.data$Symbol <- as.character(gene.data$Symbol)

# LncRNA.data <- LncRNA.data[!duplicated(LncRNA.data$Ensembl_ID), ]
rownames(LncRNA.data) <- LncRNA.data[, 1]
LncRNA.data <- LncRNA.data[-1]

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

dim(LncRNA.data) # 17948  1183
LncRNA.df <- fun(LncRNA.data)
dim(LncRNA.df) # 15112  1183
table(is.na(LncRNA.df)) 
write.csv(LncRNA.data, "tumor_normal_deseq2/LncRNA_tumor_normal_No_na_data.csv")

# Label datas,rownames(group_text) are sample ,colnames are label
group.text <- sample.label.data[, c("sample", "label")]
rownames(group.text) <- group.text[, 1]
group.text <- group.text[-1]

# Whether the column names of the expression matrix and the row names of the grouping matrix are consistent
all(rownames(group.text) == colnames(LncRNA.df))
label <- group.text$label
table(label)
colData <- group.text
express.rec <- LncRNA.df
table(is.na(express.rec))  

# DEseq2 analysis
dds <- DESeqDataSetFromMatrix(countData = express.rec, colData = colData, design = ~label, tidy = F)
# Filter out low-expression data (can be used when it is necessary to remove)
table(rowSums(counts(dds) == 0))
# Retain genes that can be found in 70% of samples, 56=1183*30%
keep <- rowSums(counts(dds) == 0) < 354
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
res <- lfcShrink(dds1, coef = "label_tumor_vs_normal", type = "apeglm")
summary(res)
plotMA(res, ylim = c(-2, 2))

# Save result
# Sort res by padj (ascending order)
resOrdered <- as.data.frame(res[order(res$padj), ])
resOrdered$Ensembl_ID <- row.names(resOrdered)
res.data <- left_join(resOrdered, gene.data, by = "Ensembl_ID")
res.data <- res.data[, c("Ensembl_ID", "Symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
dim(res.data) # 6222*7
write.table(res.data, "tumor_normal_deseq2/LncRNA_tumor_vs_normal_deseq2_result.csv", 
            row.names = TRUE, col.names = TRUE, sep = ",") # 6135*7


# According to padj<0.05, |log2FoldChange|>1, select genes in resOrdered
resOrdered2 <- resOrdered[resOrdered[, "padj"] < 0.05, ]
resOrdered3 <- resOrdered2[abs(resOrdered[, "log2FoldChange"]) > 1, ]
res.data2 <- left_join(resOrdered3, gene.data, by = "Ensembl_ID")
res.data2 <- res.data2[, c("Ensembl_ID", "Symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
res.data2 <- na.omit(res.data2)
dim(res.data2) # 1521*7
write.table(res.data2, "tumor_normal_deseq2/LncRNA_tumor_vs_normal_deseq2_select_result.csv", 
            row.names = TRUE, col.names = TRUE, sep = ",")  # 393*7


# Data Transformation
# this gives log2(n + 1)
ntd <- normTransform(dds)
log2.df <- as.data.frame(assay(ntd)) 
log2.dat <- log2.df[res.data2$Ensembl_ID, ]

coldata <- colData
coldata$sample <- rownames(coldata)
tumor.data <- coldata[which(coldata$label == "tumor"), ]
normal.data <- coldata[which(coldata$label == "normal"), ]

tumor.deseq2.dat <- log2.dat[, rownames(tumor.data)]
normal.deseq2.dat <- log2.dat[, rownames(normal.data)]
write.table(tumor.deseq2.dat, "tumor_normal_deseq2//tumor_deseq2_log2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")
write.table(normal.deseq2.dat, "tumor_normal_deseq2//normal_deseq2_log2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")



## LncRNA(type1,type2) difference analysis  by DEseq2
library(limma)
library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(DESeq2)
library(data.table)
library(openxlsx)
library(dplyr)
# Gene expression data，rownames(express) are genes ,colnames are sample
express <- read.csv("raw_data/LncRNA_count_data.csv")
gene.data <- express[, c("Ensembl_ID", "Symbol")]
gene.data$Ensembl_ID <- as.character(gene.data$Ensembl_ID)
gene.data$Symbol <- as.character(gene.data$Symbol)
rownames(express) <- express[, 1]
express <- express[-1]

sample.names.data <- read.xlsx("type1_type2_deseq2/join_name_data.xlsx")
select.gene.data <- read.csv("tumor_normal_deseq2/LncRNA_tumor_vs_normal_deseq2_select_result.csv")
table(is.na(select.gene.data$Ensembl_ID))

LncRNA.data <- express[select.gene.data$Ensembl_ID, sample.names.data$raw_sample_names]
dim(LncRNA.data) # 1521*187
write.csv(LncRNA.data, "type1_type2_deseq2/LncRNA_type1_type2_No_na.csv")

# Label datas,rownames(group_text) are sample ,colnames are label
group.text <- sample.names.data[, c("raw_sample_names", "label")]
rownames(group.text) <- group.text[, 1]
group.text <- group.text[-1]

# Whether the column names of the expression matrix and the row names of the grouping matrix are consistent
all(rownames(group.text) == colnames(LncRNA.data))
label <- group.text$label
table(label)
colData <- group.text
express.rec <- LncRNA.data
table(is.na(express.rec))  

# DEseq2 analysis
dds <- DESeqDataSetFromMatrix(countData = express.rec, colData = colData, design = ~label, tidy = F)
# Filter out low-expression data (can be used when it is necessary to remove)
table(rowSums(counts(dds) == 0))
# Retain genes that can be found in 70% of samples, 56=187*30%
keep <- rowSums(counts(dds) == 0) < 56
dds <- dds[keep, ]

dds1 <- DESeq(dds)
suppressMessages(dds1)
resultsNames(dds1)

# label_type2_vs_type1
res <- results(dds1, contrast = c("label", "type2", "type1"))
# View summary information of the results
summary(res)
plotMA(res, ylim = c(-2,2))

# Correct p
library(apeglm)
res <- lfcShrink(dds1, coef = "label_type2_vs_type1", type = "apeglm")
summary(res)
plotMA(res, ylim = c(-2,2))

# Save result
# Sort res by padj (ascending order)
resOrdered <- as.data.frame(res[order(res$padj), ])
resOrdered$Ensembl_ID <- row.names(resOrdered)
res.data <- left_join(resOrdered, gene.data, by = "Ensembl_ID")
res.data <- res.data[, c("Ensembl_ID", "Symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
dim(res.data) # 1425*7
write.table(res.data, "type1_type2_deseq2/LncRNA_type2_vs_type1_deseq2_result.csv", 
            row.names = TRUE, col.names = TRUE, sep = ",") # 6135*7


# According to padj<0.05, |log2FoldChange|>1, select genes in resOrdered
resOrdered2 <- resOrdered[resOrdered[, "padj"] < 0.05, ]
resOrdered3 <- resOrdered2[abs(resOrdered[, "log2FoldChange"]) > 1, ]
res.data2 <- left_join(resOrdered3, gene.data, by = "Ensembl_ID")
res.data2 <- res.data2[, c("Ensembl_ID", "Symbol", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
res.data2 <- na.omit(res.data2)
dim(res.data2) # 153*7
write.table(res.data2, "type1_type2_deseq2/LncRNA_type2_vs_type1_deseq2_select_result.csv", 
            row.names = TRUE, col.names = TRUE, sep = ",") # 393*7


# Data Transformation
# this gives log2(n + 1)
ntd <- normTransform(dds1)
log2.counts.df <- as.data.frame(assay(ntd))
# log2.counts.df$Ensembl_ID <- row.names(log2.counts.df)
# log2.counts.df <- left_join(log2.counts.df, gene.data, by = "Ensembl_ID")
# log2.counts.df <- log2.counts.df[-18]
write.table(log2.counts.df, "type1_type2_deseq2/type1_type2_log2_counts.csv", row.names = TRUE, col.names = TRUE, sep = ",") 






## Volcano plot
library(ggplot2)
library(ggthemes)
library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(data.table)
library(openxlsx)

res <- read.csv("type1_type2_deseq2/LncRNA_type2_vs_type1_deseq2_result.csv")  
vol.res <- as.data.frame(res)
vol.res$Symbol <- as.character(vol.res$Symbol)

# If p_value<0.05 and LF>1, it is UP; if p_value<0.05 and LF<1, it is down;
# p_value<0.05 and |LF|<1, it is not (not significant)
vol.res$threshold <- as.factor(ifelse(vol.res$padj < 0.05 & abs(vol.res$log2FoldChange) >= log2(2),
                                      ifelse(vol.res$log2FoldChange > log2(2) , "type2_Up", "type1_Up"), "Not"))
# Sort the p-values of the differential genes in ascending 
vol.res <- vol.res[order(vol.res$padj), ]

# Select the 10 genes with the largest p value from the highly expressed genes
# vol.res$Symbol <- rownames(vol.res)
up.genes <- head(vol.res$Symbol[which(vol.res$threshold == "type2_Up")], 10)

# Among the low-expressed genes, select the 10 genes with the smallest p-value
down.genes <- head(vol.res$Symbol[which(vol.res$threshold == "type1_Up")], 10)

# Merge up genes and down genes 
deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
vol.res$label <- ifelse(vol.res$Symbol %in% c(deg.top10.genes), vol.res$Symbol, "")

#  Process the data for drawing the volcano map to make the volcano map symmetrical
# The processing method is: log2FoldChange>3 is assigned a value of 3; log2FoldChange<-3 is assigned a value of -3
vol.res$log2FoldChange <- ifelse(vol.res$log2FoldChange > 3, 3, ifelse(vol.res$log2FoldChange < -3, -3, vol.res$log2FoldChange))

# plot
plot2 <- ggplot(data = vol.res, aes(x = log2FoldChange, y = -log10(padj), colour = threshold, fill = threshold)) +
  scale_color_manual(values = c("grey", "red", "blue")) +
  geom_point(alpha = 0.4, size = 1.2) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "grey", lwd = 0.6) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "grey", lwd = 0.6) +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", color = "black", family = "Times", size = 8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold", color = "black", size = 12),
        axis.text.y = element_text(face = "bold",  color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        axis.title.y = element_text(face = "bold", color = "black", size = 12)) +
  labs(x= "log2 (Fold Change)", y= "-log10 (p-value)")
plot2


library(ggrepel)
# Label gene, add layer
gg <- plot2 + geom_text_repel(data = vol.res, aes(x = log2FoldChange, y = -log10(padj), label = label),
                              size = 3, box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.8, "lines"), 
                              segment.color = "black", 
                              show.legend = FALSE) 
ggsave("type1_type2_deseq2/lncRNA_volcano_plot.pdf", plot = print(gg), width = 8, height = 8, units = "in")




#Heatmap
### heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)

# colnames(log2_counts_df) are sample,rownames are genes
log2.counts.df <- read.csv("type1_type2_deseq2/type1_type2_log2_counts.csv") 
gene.df <- read.csv("type1_type2_deseq2/LncRNA_type2_vs_type1_deseq2_select_result.csv") 
heatmap.data <- log2.counts.df[gene.df$Ensembl_ID, ]

sample.names.data <- read.xlsx("type1_type2_deseq2/join_name_data.xlsx")
sample.names <- sample.names.data$raw_sample_names
# Label datas,rownames(group_text) are sample ,colnames are label
group.text <- sample.names.data[, c("raw_sample_names", "label")]
# group.text2 <- group.text[!duplicated(group.text$raw_sample_names), ]

rownames(group.text) <- group.text[, 1]
group.text <- group.text[-1]
log2.data <- heatmap.data

# scale()
data.scale <- scale(t(log2.data))
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
# data.scale.new <- na.omit(data.scale.new)
dat <- t(data.scale.new)
all(rownames(group.text) == colnames(dat))

table(is.na(dat))  

#  plot Heatmap
pdf("type1_type2_deseq2/lncRNA_heatmap.pdf", width = 9, height = 10)
mat <- as.matrix(dat)
group.mat <-  as.matrix(group.text)
column.ha <- HeatmapAnnotation(df = group.text, 
                               col = list(label = c("type1" = "red", "type2" = "blue")))

ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("navy", "black", "yellow")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter",
                                          legend_height = unit(6, "cm"),
                                          legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = TRUE,
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


####  lncRNA enrichment analysis ####
Enrichment.barplot2 <- function(data, gene_sets = NULL, title = NULL){
  enrich <- data
  enrich$FDR_q <- -log10(enrich[, "FDR.q.value"])
  enrich1 <- enrich[order(enrich$FDR_q,decreasing = T), ]  
  print(max(enrich1$FDR_q))
  print(length(enrich1$FDR_q))
  enrich3 <- data.frame(enrich1[, "Gene.Set.Name"],enrich1[, "FDR_q"])
  colnames(enrich3) <- c("ID", "FDR_q")
  p <- ggplot(data = enrich3,aes(x = ID, y = FDR_q))  
  p1 <- p + geom_bar(stat = "identity", width = 0.5) + coord_flip()           
  p2 <- p1 + theme(panel.background=element_rect(fill = "transparent", color = "gray"),
                   axis.text.y = element_text(color = "black", size = 12))
  
  p3 <- p2 + ylim(0, 7) + scale_fill_gradient(low = "red", high = "blue")  
  p4 <- p3 + scale_x_discrete(limits = rev(enrich3[,1])) + labs(x = "",y = "-lg(FDR_q_value)", 
                                                                title = "Enrichment Analysis of lncRNA(type2_KEGG)")
  p5 <- p4 + theme(plot.margin = unit(rep(1, 4), "lines"))
  p5
}

## type1 ##
micRNA.type1.data <- read.csv("lncRNA_type1_pos0.5_select_pathway.csv")  # 13*7
dim(micRNA.type1.data)
micRNA.type1.data
Enrichment.barplot2(micRNA.type1.data)  # ylim(0, 3)

## type2  hallmark ##
micRNA.type2.hallmark.data <- read.csv("lncRNA_type2_pos0.7_hallmark_select_pathway.csv")  # 11*7
dim(micRNA.type2.hallmark.data)
micRNA.type2.hallmark.data
Enrichment.barplot2(micRNA.type2.hallmark.data)  # ylim(0, 62)

## type2  KEGG ##
micRNA.type2.KEGG.data <- read.csv("lncRNA_type2_pos0.7_KEGG_select_pathway.csv")  # 15*7
dim(micRNA.type2.KEGG.data)
micRNA.type2.KEGG.data
Enrichment.barplot2(micRNA.type2.KEGG.data)  # ylim(0, 7)


## correlation analysis between lncRNA and micRNA
library(data.table)
library(openxlsx)
library(gdata)

# Batch process micRNA target lncRNA data
files <- list.files(pattern = "*.xls")
micRNA.lncRNA.myfiles <- do.call(rbind, lapply(files, function(x) read.table(x, sep = "\t", header = T)))
write.csv(micRNA.lncRNA.myfiles, "starbase_micRNA_lncRNA.csv")
micRNA.lncRNA.geneID <- unique(micRNA.lncRNA.myfiles$geneID)
micRNA.lncRNA.Symbol <- unique(micRNA.lncRNA.myfiles$geneName)

lncRNA.data <- read.csv("LncRNA_type2_vs_type1_deseq2_select_result.csv")
inter.Ensembl.ID <- intersect(lncRNA.data$Ensembl_ID,micRNA.lncRNA.geneID)
inter.Symbol <- intersect(lncRNA.data$Symbol,micRNA.lncRNA.Symbol)


####### Correlation analysis between lncRNA and micRNA, mRNA
### micRNA data (9) ###
#log2 data
diff.micRNA.data <- read.csv("type1_type2_deseq2_After_screening_result.csv")
diff.micRNA <- diff.micRNA.data$micRNA
micRNA.log2.data <- read.csv("B_log2_data.csv")
rownames(micRNA.log2.data) <- micRNA.log2.data$micRNA
micRNA.log2.data <- subset(micRNA.log2.data, select = -micRNA)
diff.micRNA.log2.data <- micRNA.log2.data[diff.micRNA, ] # 9*181
write.csv(diff.micRNA.log2.data, "diff_micRNA_log2_data.csv")

### lncRNA data（103） ###
# log2 data
diff.lncRNA.data <- read.csv("LncRNA_type2_vs_type1_deseq2_select_result.csv")
diff.lncRNA <- diff.lncRNA.data$Ensembl_ID  # 153
lncRNA.log2.data <- read.csv("type1_type2_log2_counts.csv")
# Extract characters from fistr to 16th position.
log2.new.colnames <- as.character(lapply(colnames(lncRNA.log2.data), function(x) substring(x, 1, 16)))
colnames(lncRNA.log2.data) <- log2.new.colnames
diff.lncRNA.log2.data <- lncRNA.log2.data[diff.lncRNA, ] # 153*187
write.csv(diff.lncRNA.log2.data, "diff_lncRNA_log2_data.csv")

### mRNA data (1073) ###
# tumor_normal_deseq2_result_log2_data
mRNA.log2.data <- read.csv("tumor_log2_without_label_data.csv")
rownames(mRNA.log2.data) <- mRNA.log2.data$Gene
mRNA.log2.data <- subset(mRNA.log2.data, select = -Gene)

# lncRNA(cor_mRNA)_data
mRNA.log2.data.lncrna <- mRNA.log2.data[, colnames(diff.lncRNA.log2.data)] # 10645*187
all(colnames(diff.lncRNA.log2.data) == colnames(mRNA.log2.data.lncrna)) # true
write.csv(diff.lncRNA.log2.data, "lncRNA(cor_mRNA)_data.csv")

# lncRNA(cor_micRNA)_data
diff.lncRNA.log2.data.micrna <- diff.lncRNA.log2.data[, colnames(diff.micRNA.log2.data)]  # 153*181
write.csv(diff.lncRNA.log2.data.micrna, "lncRNA(cor_micRNA)_data.csv")

# micRNA(cor_mRNA)
mRNA.log2.data.micrna <- mRNA.log2.data[, colnames(diff.micRNA.log2.data)]  # 10645*181
write.csv(mRNA.log2.data.micrna, "mRNA(cor_micRNA)data.csv")

# correction analysis
fun <- function(data1, data2) {
  result <- c()
  result.all <- c()
  for (i in rownames(data1)) {
    res <- c()
    print(i)
    dat1 <- as.numeric(data1[i, ])
    for (j in rownames(data2)) {
      re <- c()
      dat2 <- as.numeric(data2[j, ])
      correlation <- cor.test(dat1, dat2, alternative = "two.side", method = "spearman", conf.level = 0.95, exact = FALSE)
      if(p.adjust(correlation$p.value, method = "BH") < 0.05) {
        re <- cbind(i, j, correlation$estimate, correlation$p.value, p.adjust(correlation$p.value, method = "BH"))
        re <- as.data.frame(re)
        colnames(re) <- c("lncRNA", "target.rna", "cor", "pvalue", "padj")
      }
      else {
        next
      }
      res <- rbind(res, re)
      res <- as.data.frame(res)
    }
    if(is.null(res)) {
      next
    }
    else {
      res$cor <- as.numeric(as.character(res$cor))
      max.index <- which(res$cor == max(res$cor), arr.ind = TRUE)
      print(max.index)
      result <- rbind(result, res[max.index[1], ])
      min.index <- which(res$cor == min(res$cor))
      print(min.index)
      result <- rbind(result, res[min.index[1], ])
      result.all <- rbind(result.all, res)
    }
  }
  out <- list(one = result, two = result.all)
  return(out)
}


# lncRNA_cor_micRNA
all(colnames(diff.lncRNA.log2.data.micrna) == colnames(diff.micRNA.log2.data)) # TRUE
lncRNA.micRNA.out <- fun(diff.lncRNA.log2.data.micrna, diff.micRNA.log2.data)
lncRNA.micRNA.cor.res <- lncRNA.micRNA.out$one
write.csv(lncRNA.micRNA.cor.res, "lncRNA_micRNA_log2_cor_res.csv")
lncRNA.micRNA.cor.res.2 <- lncRNA.micRNA.out$two
write.csv(lncRNA.micRNA.cor.res.2, "lncRNA_micRNA_log2_cor_res_all.csv")


#### lncRNA_cor_mRNA(537), micRNA_cor_mRNA(537) ####
# lncRNA_cor_mRNA
all(colnames(diff.lncRNA.log2.data) == colnames(mRNA.log2.data)) # TRUE
lncRNA.mRNA.out <- fun(diff.lncRNA.log2.data, mRNA.log2.data)
lncRNA.mRNA.cor.res <- lncRNA.mRNA.out$one
lncRNA.mRNA.cor.res.2 <- lncRNA.mRNA.out$two


# micRNA_cor_mRNA
mirRNA.mRNA.out <- fun(diff.micRNA.log2.data, mRNA.log2.data.2)
mirRNA.mRNA.cor.res <- mirRNA.mRNA.out$one
mirRNA.mRNA.cor.res.2 <- mirRNA.mRNA.out$two


#### lncRNA_cor_mRNA(10645), micRNA_cor_mRNA(10645) ####
# lncRNA_cor_mRNA
all(colnames(diff.lncRNA.log2.data) == colnames(mRNA.log2.data.lncrna)) #TRUE
lncRNA.mRNA.out <- fun(diff.lncRNA.log2.data, mRNA.log2.data.lncrna)
lncRNA.mRNA.cor.res <- lncRNA.mRNA.out$one
lncRNA.mRNA.cor.res.2 <- lncRNA.mRNA.out$two


# micRNA_cor_mRNA
all(colnames(diff.micRNA.log2.data) == colnames(mRNA.log2.data.micrna)) #TRUE
mirRNA.mRNA.out <- fun(diff.micRNA.log2.data, mRNA.log2.data.micrna)
mirRNA.mRNA.cor.res <- mirRNA.mRNA.out$one
mirRNA.mRNA.cor.res.2 <- mirRNA.mRNA.out$two


## lncRNA_micRNA scatter plot##
lncRNA.micRNA.cor.res.2 <- read.csv("lncRNA_micRNA_log2_cor_res_all.csv")
names(lncRNA.micRNA.cor.res.2)[names(lncRNA.micRNA.cor.res.2) == 'target.rna'] <- 'micRNA'
lncRNA.micRNA.cor.res.2 <- read.csv("lncRNA_micRNA_log2_cor_res_all.csv")
lncRNA.names <- read.csv("lncRNA(12)_names.csv")

require(ggplot2)
data <- lncRNA.micRNA.cor.res.2[which(lncRNA.micRNA.cor.res.2$lncRNA %in% lncRNA.names$lncRNA), ]
data[, "-lg(Pvalue)"] <- -log(10, data[, "pvalue"])
pdf("lncRNA_cor_micRNA_plot.pdf", width = 12, height = 10)
pl1 <-  ggplot(data, aes(y = lncRNA,  x = micRNA)) +  
  scale_fill_continuous(low = "LightSteelBlue", high = "LightSlateGray")
P2 <- pl1 + geom_point(aes(colour = data[, 4],  size = data[, 5])) +
  scale_color_gradient(low = "yellow", high = "DarkOrange") + scale_size(range = c(1, 10)) + theme_bw()
p3 <- P2 + theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 30), legend.position = "right", 
                 legend.text = element_text(size = 7.5), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
p3 + labs(colour = "Correlation", size = "-lg(Pvalue)") + guides(size = guide_legend(override.aes = list(shape = 1)))
dev.off()


## lncRNA_mRNA scatter plot ##
lncRNA.diffmRNA.cor <- read.csv("lncRNA_cor_diff_mRNA.csv")
lncRNA.names <- read.csv("lncRNA(12)_names.csv")
key.mRNA.names <- read.csv("import_genes_withoutsame.csv")

require(ggplot2)
data.1 <- lncRNA.diffmRNA.cor[which(lncRNA.diffmRNA.cor$lncRNA %in% lncRNA.names$lncRNA), ]
data <- data.1[which(data.1$mRNA %in% key.mRNA.names$gene), ]
data[, "-lg(Pvalue)"] <- -log(10, data[, "pvalue"])
pdf("lncRNA_cor_mRNA_plot.pdf", width = 12, height = 10)
pl1 <-  ggplot(data, aes(y = lncRNA,  x = mRNA)) +  
  scale_fill_continuous(low = "LightSteelBlue", high = "LightSlateGray")
P2 <- pl1 + geom_point(aes(colour = data[, 3],  size = data[, 4])) +
  scale_color_gradient(low = "yellow", high = "DarkOrange") + scale_size(range = c(1, 10)) + theme_bw()
p3 <- P2 + theme(axis.text.x = element_text(size = 8, vjust = 0.5, hjust = 0.5, angle = 30), legend.position = "right", 
                 legend.text = element_text(size = 7.5), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
p3 + labs(colour = "Correlation",size = "-lg(Pvalue)") + guides(size = guide_legend(override.aes = list(shape = 1)))
dev.off()







