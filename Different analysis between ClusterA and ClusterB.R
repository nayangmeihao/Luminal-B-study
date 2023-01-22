library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(DESeq2)
library(data.table)
library(openxlsx)

# colnames(express) are samples,rownames(express) are genes
express <- read.xlsx("cluster2_B_count_data.xlsx", sheet = 1)  # 4248*188
group.text <- read.xlsx("s_cluster_2_res.xlsx", sheet = 1)  # 188*2
express.rec <- as.data.frame(express)
group.text <- as.data.frame(group.text)

rownames(express.rec) <- express.rec$Gene
express.rec <- subset(express.rec, select = -Gene)
rownames(group.text) <- group.text$sample
group.text<- subset(group.text, select = -sample)

all(rownames(group.text) == colnames(express.rec))
label <- group.text$label

colData <- group.text
dds <- DESeqDataSetFromMatrix(countData = express.rec, colData = colData, design = ~label, tidy = F)

# Retain the genes that can be found in 30% of samples, 58=187*30%
keep <- rowSums(counts(dds) == 0)< 56
dds <- dds[keep, ]
counts(dds)[1:10, 1:3]
dim(dds) 

# Difference analysis
dds <- DESeq(dds)
suppressMessages(dds)
resultsNames(dds)

res <- results(dds, contrast = c("label", "ClusterB", "ClusterA"))
summary(res)
plotMA(res, ylim = c(-2, 2))

# Correct q
library(apeglm)
res <- lfcShrink(dds, coef = "label_ClusterB_vs_ClusterA", type = 'apeglm')
summary(res)
plotMA(res, ylim = c(-2,2))

# save result
resOrdered <- as.data.frame(res)
# In resOrdered, select genes by padj<0.05, |log2FoldChange|>1
resOrdered2 <- resOrdered[(abs(resOrdered[, "log2FoldChange"]) > 1)&(resOrdered[, "padj"] < 0.05), ]
write.table(resOrdered2, "cluster2_deseq2_result.csv", row.names=TRUE, col.names = TRUE, sep = ",") 

# Save unfiltered difference analysis resultsb
write.table(resOrdered, "cluster2_deseq2_result.csv", row.names = TRUE, col.names = TRUE, sep = ",") 


# Data Transformation
# this gives log2(n + 1)
ntd <- normTransform(dds)
log2.df <- as.data.frame(assay(ntd))  
head(log2.df)
write.table(log2.df, "B_log2_data.csv",row.names=TRUE,col.names=TRUE,sep=",") 

# Take the genes with pvalue<0.05 and |log2FoldChange|>1 in resOrdered for WGCAN analysis
re2 <- resOrdered[(abs(resOrdered[, "log2FoldChange"]) > 1) & (resOrdered[, "pvalue"] < 0.05), ]
genes <- rownames(re2)
log2.data <- log2.df[genes, ]
write.table(log2.data, "B_DEseq2_log2_pvalue_data.csv", row.names = TRUE, col.names = TRUE, sep = ",")




#Heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)

# 0605
# colnames(log2_counts_df) are sample,rownames are genes
log2.counts.df <- read.xlsx("B_DEseq2_log2_data.xlsx")  # 537*188
rownames(log2.counts.df) <- log2.counts.df$gene
log2.counts.df <- subset(log2.counts.df, select = -gene)

clinical.data <- read.xlsx("final.data.xlsx")
rownames(clinical.data) <- clinical.data$sample
clinical.data <- subset(clinical.data, select = -sample)


group.text <- read.xlsx("cluster_2_res.xlsx", sheet = 1)
group.text <- as.data.frame(group.text)
rownames(group.text) <- group.text$sample
group.text <- subset(group.text, select = -sample)

log2.data <- log2.counts.df

# scale()
data.scale <- scale(t(log2.data))
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
# data.scale.new <- na.omit(data.scale.new)
dat <- t(data.scale.new)
all(rownames(group.text) == colnames(dat))

# plot Heatmap
mat <- as.matrix(dat)
mark.gene <- c("APOB", "CCND1", "CDH11", "COL1A1", "COL2A1", "COL5A1", "COX6C", "DDR2", "ETV1", "FCGR3A",
               "FGFR3", "HIST1H3B", "INHBA", "KLK2", "KRT15", "KRT31", "KRT75", "KRT81", "NRG1", "PDGFRA", 
               "RET", "ROS1", "RSPO3", "S100A7", "S100A8", "S100A9", "SOX11", "TNFRSF18", "TNFSF4", "ZNF703")

at <- c()
for (i in mark.gene) {
  at <- c(at, c(which(rownames(mat) == i)))  
}
row.anno <-  rowAnnotation(foo = anno_mark(at = at,
                                           labels = mark.gene,
                                           link_width = unit(3, "mm")))


column.ha <- HeatmapAnnotation(df = clinical.data,
                               col = list( Age = c("<=60" = "#feeeed", ">60" = "#f47920"),
                                           Race = c("White" =  "#feeeed", "Nonwhite" = "#f47920"),
                                           Tumor.stage = c("<=II" =  "#feeeed", ">II" = "#f47920"),
                                           Neoplasm.status = c("Tumor Free" =  "#f47920", "With Tumor" = "#feeeed"),
                                           Metastasis.status = c("YES" =  "#feeeed", "NO" = "#f47920"),
                                           Surgical.margin.status = c("Positive" = "#feeeed", "Negative" = "#f47920"),
                                           Menopause.status = c("<=12 mo since LMP" =  "#feeeed", ">12 mo since LMP" = "#f47920"),
                                           Lymph.node.stage = c("1" =  "#feeeed", "0" = "#f47920"),
                                           Tumor.location = c("Right" =  "#feeeed", "Left" = "#f47920"),
                                           HER2 = c("Positive" =  "#feeeed", "Negative" = "#f47920"),
                                           ER = c("Positive" =  "#feeeed", "Negative" = "#f47920"),
                                           PR = c("Positive" =  "#feeeed", "Negative" = "#f47920"),
                                           CNA = c("<=0.251" =  "#feeeed", ">0.251" = "#f47920"),
                                           Mutation.count = c("<=30" =  "#feeeed", ">30" = "#f47920"),
                                           Clinical.stage = c("<=II" =  "#feeeed", ">II" = "#f47920"),
                                           label = c("type1" =  "red", "type2" = "blue")), na_col = "gray")



ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#7fb80e", "white", "red")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter", 
                                          legend_height = unit(6, "cm"),
                                          legend_direction = "vertical"),
              top_annotation = column.ha,
              right_annotation = row.anno,
              show_row_dend = FALSE,
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = FALSE,
              column_split = group.text,
              column_title = NULL)
# ht

pdf("Cluster_heatmap.pdf", width = 9, height = 12)
lgd <- Legend(labels =  "NA",  legend_gp = gpar(fill = 8))
draw(ht, annotation_legend_list = list(lgd))
dev.off()


# volcano
# visualization
library(ggplot2)
library(ggthemes)
library(parallel)
library(BiocGenerics)
library(stats4)
library(S4Vectors)
library(data.table)
library(openxlsx)

vol.res <- read.csv("cluster2_deseq2_result.csv")

# If p_value<0.05 and LF>1, it is UP; if p_value<0.05 and LF<1, it is down;
# p_value<0.05 and |LF|<1, it is not (not significant)
vol.res$threshold <- as.factor(ifelse(vol.res$padj < 0.05 & abs(vol.res$log2FoldChange) >= log2(2),
                                      ifelse(vol.res$log2FoldChange > log2(2), "type2_Up", "type1_Up"), "Not"))

vol.res <- vol.res[order(vol.res$padj), ]
vol.res$Symbol <- vol.res$gene
# up.genes <- head(vol_res$Symbol[which(vol_res$threshold == 'Up')],10)
# down.genes <- head(vol_res$Symbol[which(vol_res$threshold == 'Down')],10)

# deg.top10.genes <- c(as.character(up.genes),as.character(down.genes))
key.genes <- c("APOB", "CCND1", "CDH11", "COL1A1", "COL2A1", "COL5A1", "COX6C", "DDR2",
               "ETV1", "FCGR3A", "FGFR3", "HIST1H3B", "INHBA", "KLK2", "KRT15", "KRT31", 
               "KRT75", "KRT81", "NRG1", "PDGFRA", "RET", "ROS1", "RSPO3", "S100A7", "S100A8",
               "S100A9", "SOX11", "TNFRSF18", "TNFSF4", "ZNF703")

vol.res$label <- ifelse(vol.res$Symbol %in% c(key.genes), vol.res$Symbol, "")
#  Process the data for drawing the volcano map to make the volcano map symmetrical
# The processing method is: log2FoldChange>3 is assigned a value of 3; log2FoldChange<-3 is assigned a value of -3
new.vol.res <- vol.res
new.vol.res$log2FoldChange <- ifelse(new.vol.res$log2FoldChange > 3, 3, 
                                     ifelse(new.vol.res$log2FoldChange < -3, -3, new.vol.res$log2FoldChange))

pdf("cluster2_volcano_plot.pdf", width = 8, height = 8)

plot2 <- ggplot(data = new.vol.res, aes(x = log2FoldChange, y = -log10(padj), colour = threshold, fill = threshold)) + 
  scale_color_manual(values = c("grey", "red", "blue"))+
  geom_point(alpha = 0.4, size = 1.2) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "grey", lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "grey", lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face = "bold", color = "black", family = "Times", size = 8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold", color = "black", size = 12),
        axis.text.y = element_text(face = "bold",  color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        axis.title.y = element_text(face = "bold",color = "black", size = 12)) +
  labs(x = "log2(Fold Change)", y = "-log10(pvalue)")
plot2

library(ggrepel)
gg <- plot2 + geom_text_repel(data = new.vol.res, aes(x = log2FoldChange, y = -log10(padj), label = label),
                              size = 3, box.padding = unit(0.5, "lines"),
                              point.padding = unit(0.8, "lines"), 
                              segment.color = "black", 
                              show.legend = FALSE) 
gg
dev.off()



