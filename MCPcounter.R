# Need to study whether MCPcount needs read counts or can use normalization values
# But Xcell is very clear. The value of read counts or normalization can be used, because it only does ranking
# When using the TCGA database to verify the original text of MCPcounter Methodology, normalized results are used

library(MCPcounter)
# The data passed in by MCPcounter is the data after log2 standardization
# rownames(all.rna) are genes, colnames are patients
estimates <- MCPcounter.estimate(all.rna, featuresType = "HUGO_symbols")


#Heatmap
library(data.table)
library(openxlsx)

# clinical data
clinical.data <- read.xlsx("final.data.xlsx")
clinical.data <- as.data.frame(clinical.data)

clinical.dat.hot <- clinical.data[, c("sample", "label", "HER2", "CNA")]
sub.clinical.data <- clinical.dat.hot[, c("sample", "label")]


# MCPcounter data
mcp.data <- read.xlsx("mcp_result.xlsx")
mcp.data <- data.frame(mcp.data)
raw.mcpdata <- mcp.data
rownames(mcp.data) <- mcp.data$sample
mcp.data <- subset(mcp.data, select = -sample)


# Combine clinical data with mcp data for t test and heatmap drawing
dat <- merge(sub.clinical.data, raw.mcpdata, by = "sample", sort = FALSE)
write.table(dat,"clinical_mcp_data.csv",row.names = TRUE,col.names = TRUE, sep = ",")

# Scale mcp.data
data.scale <- scale(mcp.data)
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
data.scale.new <- na.omit(data.scale.new)

data.scale.new <- data.frame(data.scale.new)
data.scale.new["sample"] <- rownames(data.scale.new)

# Add a label to the scaled data
mer.data <- merge(sub.clinical.data, data.scale.new, by = "sample", sort = FALSE)
mer.data.copy <- mer.data
rownames(mer.data) <- mer.data$sample
mer.data <- subset(mer.data, select = -sample)


# type1
# rownames(data1) are samples
data1 <- mer.data[which(mer.data$label == "type1"), ]

# cluster
data1 <- subset(data1, select = -1)
hc <- hclust(dist(data1, method = "euclidean"), method = "complete")
new.label <- cutree(hc, k = 3)
new.df.label <- data.frame(new.label)

new.df.label["sample"] <- rownames(new.df.label)
la.t1.df <- new.df.label[order(-new.df.label$new.label), ]


# type2
data2 <- mer.data[which(mer.data$label == "type2"), ]
data1 <- subset(data2, select = -1)
hc <- hclust(dist(data1, method = "euclidean"), method = "complete")
new.label <- cutree(hc, k = 3)
new.df.label <- data.frame(new.label)
new.df.label["sample"] <- rownames(new.df.label)
la.t2.df <- new.df.label[order(new.df.label$new.label), ]


# combin la.t1.df and la.t2.df
new.label.data <- rbind(la.t1.df, la.t2.df)
m.data <- merge(new.label.data, mer.data.copy, by = "sample", sort = FALSE)

m.data2 <- subset(m.data, select = c(-2, -3))
new.data <- merge(m.data2, clinical.dat.hot, by = "sample", sort = FALSE)

df <- new.data
rownames(df) <- df[, 1]
df <- df[, -1]

da1 <- df[, 9:11]
annotation.col <- da1[, c("CNA", "HER2", "label")]
an.col <- annotation.col
an.col$sample <- rownames(an.col)
group.text <- an.col[, c("sample", "label")]
rownames(group.text) <- group.text$sample
group.text<- subset(group.text, select = -sample)

mcp.data.hot <- df[, 1:8]
mcp.data.hot <- t(mcp.data.hot)
mat <- as.matrix(mcp.data.hot)


library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)

column.ha <- HeatmapAnnotation(df = annotation.col,
                               col = list(label = c("type1" =  "red", "type2" = "blue"),
                                          HER2 = c("Positive" =  "#feeeed", "Negative" = "#f47920"),
                                          CNA = c(">0.251" =  "#feeeed", "<=0.251" = "#f47920")),
                               na_col = "gray" )

ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#7fb80e", "white", "red")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter",
                                          legend_height = unit(6, "cm"), legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = FALSE,
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              column_split = group.text,
              column_title = NULL
)
ht


pdf("mcpcounter.pdf", width = 9, height = 9)
lgd <- Legend(labels =  "NA", legend_gp = gpar(fill = 8))
draw(ht, annotation_legend_list = list(lgd))
dev.off()





