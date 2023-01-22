
# search std_cut mean_cut 
library(data.table)
library(openxlsx)

# input data
B.data <- read.xlsx("D:\\study_data\\B_cluster_data\\0601\\B_log2_data.xlsx", sheet = 1)
rownames(B.data) <- B.data$Gene
B.data <- subset(B.data, select = -Gene )
dim(B.data)
head(B.data)

st <- seq(0, 3, by = 0.001)
mean <- seq(0, 11, by = 0.001)
svp <- data.frame(log2Exp = apply(B.data, 1, mean), sd = apply(B.data, 1, sd))
rownames(svp) <- rownames(B.data)
svp$selected <- rep('black', nrow(svp))

for (i in st) {
  for(j in mean) {
    sd.cut <- i
    mean.cut <- j
    R1 <- sum(svp$log2Exp > mean.cut & svp$sd > sd.cut)
    if (R1==1000) { 
      print(i)
      print(j)
    }
  }
}


library(data.table)
library(openxlsx)

# input data
B.data <- read.xlsx("D:\\study_data\\B_cluster_data\\0601\\B_log2_data.xlsx")
rownames(B.data) <- B.data$Gene
B.data <- subset(B.data, select = -Gene)

# Calculate the sd of each gene
svp <- data.frame(log2Exp = apply(B.data, 1, mean), sd = apply(B.data, 1, sd))
rownames(svp) <- rownames(B.data)
svp$selected <- rep("black", nrow(svp))
sd.cut <- 1.13
mean.cut <- 6.735
svp$selected[svp$log2Exp > mean.cut & svp$sd > sd.cut] <- "red"

pdf("Extract_genes.pdf")
plot(x = svp$log2Exp,
     y = svp$sd,
     col = svp$selected,
     pch = ".",
     xlab = "log2 expression",
     ylab = "Standard deviation")
abline(h = sd.cut)
abline(v = mean.cut)
R1 <- sum(svp$selected == "red")
R2 <- nrow(filter(svp, svp$log2Exp > mean.cut & svp$sd < sd.cut ))
L1 <- nrow(filter(svp, svp$log2Exp < mean.cut & svp$sd > sd.cut ))
L2 <- nrow(filter(svp, svp$log2Exp < mean.cut & svp$sd < sd.cut ))
text(x = 12, y = 4, label = paste0(R1, "genes"), col = "red")
text(x = 2, y = 0.5, label = paste0(L2, "genes"), col = "blue")
text(x = 2, y = 4, label = paste0(L1, "genes"), col = "blue")
text(x = 12, y = 0.5, label = paste0(R2, "genes"), col = "blue")
dev.off()

# character genes
svp$gene <- rownames(B.data)
write.xlsx(svp, "D:\\study_data\\B_cluster_data\\0603\\s1_extract_feature_res.xlsx")
gene.df <- svp[which(svp$selected == "red"), ]
da <- B.data[gene.df$gene,]  # rownames(da) are genes,colnames are sample
da$gene <- rownames(da)
write.xlsx(da, "D:\\study_data\\B_cluster_data\\0603\\cluster_2_res_data.xlsx")


## plot
library(data.table)
library(openxlsx)

# input data
B.data <- read.xlsx("D:\\study_data\\B_cluster_data\\0601\\B_log2_data.xlsx")
rownames(B.data) <- B.data$Gene
B.data <- subset(B.data, select = -Gene)

# Calculate the sd of each gene
svp <- data.frame(log2Exp = apply(B.data, 1, mean), sd = apply(B.data, 1, sd))
rownames(svp) <- rownames(B.data)
svp$selected <- rep("black", nrow(svp))
sd.cut <- 1.13
mean.cut <- 6.735
svp$selected[svp$log2Exp > mean.cut & svp$sd > sd.cut] <- "red"

pdf("Extract_genes_0929_2020.pdf", useDingbats = FALSE)
plot(x = svp$log2Exp,
     y = svp$sd,
     col = svp$selected,
     pch = ".",
     xlab = "log2 expression",
     ylab = "Standard deviation")
abline(h = sd.cut)
abline(v = mean.cut)
R1 <- sum(svp$selected == "red")

R2 <- nrow(svp[which(svp$log2Exp > mean.cut & svp$sd < sd.cut), ])
L1 <- nrow(svp[which(svp$log2Exp < mean.cut & svp$sd > sd.cut), ])
L2 <- nrow(svp[which(svp$log2Exp < mean.cut & svp$sd < sd.cut), ])
# R2 <- nrow(filter(svp, svp$log2Exp > mean.cut & svp$sd < sd.cut ))
# L1 <- nrow(filter(svp, svp$log2Exp < mean.cut & svp$sd > sd.cut ))
# L2 <- nrow(filter(svp, svp$log2Exp < mean.cut & svp$sd < sd.cut ))
text(x = 12, y = 4, label = paste0(R1, "genes"), col = "red")
text(x = 2, y = 0.5, label = paste0(L2, "genes"), col = "blue")
text(x = 2, y = 4, label = paste0(L1, "genes"), col = "blue")
text(x = 12, y = 0.5, label = paste0(R2, "genes"), col = "blue")
dev.off()


## clustering
library(data.table)
library(openxlsx)
# input data
B.data <- read.xlsx("D:\\study_data\\B_cluster_data\\0602\\s1_B_log2_data.xlsx", sheet = 1)

# data processing
rownames(B.data) <- B.data$Gene
B.data <- subset(B.data, select = -Gene )

data <- t(B.data)
# data was normalized ,normalized gene
df <- scale(data) 

# PCA analysis
# rownames(df) are patients, colnames are genes
# pdf("Scree_plot.pdf")
ttce.pca <- prcomp(df, center = TRUE, scale. = TRUE)
summary(ttce.pca)  # patients cluster
plot(ttce.pca, type = "l")  # Scree plot
dev.off()

# Cluster with ggfortify
library(ggfortify)
autoplot(ttce.pca, colour = 'label', data = df.pca)

# ConsensusClusterPlus
library(ConsensusClusterPlus)

# Data preprocessing
# Consensus Clustering with ConsensusClusterPlus
ccw <- t(df)  # rownames(ccw) are genes, colnames are patients

# Pre-process the data before clustering
ccw[is.na(ccw)] <- 0
ccw <- as.matrix(ccw)
# sweep function minus median to normalize
ccw <- sweep(ccw, 1, apply(ccw, 1, median, na.rm = TRUE))

# Graphic Output Description
# pdf("Census clustering with k-means_Surface.pdf", width = 6, height = 6)
title <- "Consensus Clustering for Breast cancer patients"
dt <- as.dist(1 - cor(ccw, method = "pearson"))
myr <- ConsensusClusterPlus(dt,
                            maxK = 3,
                            reps = 5,
                            pItem = 0.8,
                            pFeature = 1,
                            title = title,
                            distance = "euclidean",
                            clusterAlg = "km",
                            seed = 10100)
dev.off()


# Display 2 clusters
c2 <- data.frame(myr[[2]][["consensusClass"]], stringsAsFactors = FALSE)
c2 <- tibble::rownames_to_column(c2)
colnames(c2) <- c("Patients", "group")
c2$Patients <- paste0(gsub("\\_TotalPep", "", c2$Patients))
label.data  <- c2[sort(c2$group,index.return=TRUE)$ix, ]
write.xlsx(label.data, "D:\\study_data\\B_cluster_data\\0603\\cluster_2_res.xlsx")


## plot the result of clustering
library(data.table)
library(openxlsx)

# Input data
# colnames(df) are genes,rownames are samples
data <- read.xlsx("D:/study_data/PCA_data/cluster_2_res_data_order.xlsx")
rownames(data) <- data$sample
data <- subset(data, select = -sample)

df <- data[-1]
group <- factor(data[, 1])

# rownames are sample,colnames(df) are genes
pca <- prcomp(df, center = TRUE, retx = T) 
df.pca <- data.frame(pca$x, label = group)  
head(df.pca, 3)  

df.pca2 <- df.pca
df.pca2$label <- ifelse(df.pca2$label == "type1", -1, 1)
head(df.pca2)

# alpha = 1/2, levles = 0.95
library(ggplot2)
p1<- ggplot(data = df.pca, aes(x = df.pca$PC1, y = df.pca$PC2, color = df.pca$label))
p1.1 <- p1 + stat_ellipse(aes(fill = df.pca$label), type = "norm", geom = "polygon", alpha = 0.2, levles = 0.1)
p1.2 <- p1.1 + geom_point() + labs(x = xlab, y = ylab, color = "") + guides(fill = F)
p1.2



