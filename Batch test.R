# install.packages('devtools')
library(devtools)
install_github('fawda123/ggord')
library(ggord)
library(data.table)
library(openxlsx)
library(grid)
library(checkmate)

getwd()  # Check Path
setwd("D:/Experiment/PCA_data")  # Change the working directory to PCA_data
# input data
pca.data <- read.csv("readscount.csv")

# counts data,rownames(df) are sample,colnames are genes
df <- pca.data[-1]

# batch data
pca.group <- factor(pca.data[, 1])

pca1 <- prcomp(df, center = TRUE, retx = T)
df.pcs <- data.frame(pca1$x, batch_num = pca.group)  
head(df.pcs)  # View the results of principal component analysis

# plot with ggord
p1 <- ggord(pca1, grp_in = pca.group, arrow = 0, vec_ext = 0, txt = NULL)
p1


# plot with ggord
library(ggplot2)

# Add confidence ellipse and add percentage of PC1 and PC2
percentage <- round(pca1$sdev / sum(pca1$sdev) * 100, 2)
percentage <- paste(colnames(df.pcs), "(", paste(as.character(percentage), "%", ")", sep = ""))
ggplot(df_pcs, aes(x = PC1, y = PC2, color = batch_num))+
  geom_point()+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+ geom_point()+stat_ellipse(level = 0.95, show.legend = F) 


# Batch processing of data
# devtools::install_github("zhangyuqing/sva-devel")
library(sva)

df2 <- as.matrix(df)
df3 <- t(df2)
batch2 <- as.numeric(pca.group)

# rownames (df3) are genes,colnames are sample. batch is character data
adjusted <- ComBat_seq(df3, batch = batch2, group = NULL)
adj <- t(adjusted)
pca2 <- prcomp(adj, center = TRUE, retx = T)
df.pcs2 <- data.frame(pca2$x, batch.num = pca.group)  
head(df.pcs2)

# plot
percentage <- round(pca2$sdev / sum(pca2$sdev) * 100, 2)
percentage <- paste(colnames(df.pcs2), "(", paste(as.character(percentage), "%", ")", sep = ""))
ggplot(df_pcs2, aes(x = PC1, y = PC2, color = batch.num))+
  geom_point()+ 
  xlab(percentage[1]) +
  ylab(percentage[2]) + geom_point() + stat_ellipse(level = 0.95, show.legend = F) 

save.image()