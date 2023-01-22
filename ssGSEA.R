
rm(list = ls())
library(data.table)
library(openxlsx)
library(dplyr)
library(tibble)

log2.counts.df <- read.csv("TP53_log2_data.csv",header = FALSE) 
colnames(log2.counts.df) <- union("gene",log2.counts.df[1,][1:112])
log2.counts.df <- log2.counts.df[-1,]
rownames(log2.counts.df) <- log2.counts.df$gene
log2.counts.df <- log2.counts.df[,-1]
for (t in 1:ncol(log2.counts.df)) {log2.counts.df[, t] <- as.numeric(log2.counts.df[, t])}

label.data <- read.csv("Same_sample_data.csv",head = TRUE) #112*2

expe.data <- log2.counts.df[,label.data$sample] 
#ssGSEA analysis
library(GSEABase)
library(GSVA)

geneSet <- read.csv("/Users/nana/mycode/breast_cancer/add_exp/data/CellReports.txt",header = F,sep = "\t",) # 用EXCEL打开删除NA列
class(geneSet)
geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
save(l,file = "gene_set.Rdata")


dat <- expe.data
dat <- as.matrix(dat)
library(GSVA)
ssgsea<- gsva(dat, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

write.csv(ssgsea.1, "ssGSEA_result_plot.csv")


## Heatmap
label.data <- read.csv("Same_sample_data.csv",head = TRUE) #112*2
raw.oncodata <- read.xlsx("oncodata_data_processed.xlsx")
rownames(raw.oncodata) <- raw.oncodata$sample
label.df <- raw.oncodata[label.data$sample,c("sample","TP53_label")]
colnames(label.df) <- c("sample","label")
label.df <- label.df[order(label.df$label,decreasing = F),]
label.df <- subset(label.df,select=-sample)
head(label.df)
table(label.df$label) 


gsva_matrix <- ssgsea[,rownames(label.df)]
pheno <- label.df
library(pheatmap)
gsva_matrix1<- t(scale(t(gsva_matrix)))
gsva_matrix1[gsva_matrix1< -2] <- -2
gsva_matrix1[gsva_matrix1>2] <- 2
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 
                'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 
                'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
anti<- gsub('^ ','',rownames(gsva_matrix1))%in%anti_tumor
pro<- gsub('^ ','',rownames(gsva_matrix1))%in%pro_tumor
non <- !(anti|pro)
gsva_matrix1<- rbind(gsva_matrix1[anti,],gsva_matrix1[pro,],gsva_matrix1[non,])
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1)
annotation_col = pheno
rownames(annotation_col)<-colnames(nor_gsva_matrix1)
bk = unique(c(seq(0,1, length=100)))
pheatmap(nor_gsva_matrix1,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = annotation_col,
         breaks=bk,
         legend = T,
         cutree_rows = 3,
         gaps_col =sum(label.df == "mutated"),
         cutree_cols = 2,
         # cellwidth=5,
         # cellheight=5,
         # fontsize=5,
         # gaps_row = c(12,20),
         filename = 'PIK3CA_ssGSEA_heatmep.pdf',
         width = 8,
         height=6)



dat <- as.data.frame(t(gsva_matrix))
dat2 <- dat
dat2$sample <- rownames(dat2)
label.text <- label.df
label.text$sample <- rownames(label.text)
dat.df <- inner_join(label.text, dat2,by= "sample") #22283*162


res <- c()
for(i in 3:ncol(dat.df)){
  print(colnames(dat.df)[i])
  p <- t.test(dat.df[,i]~dat.df$label) 
  re <- cbind(colnames(dat.df)[i], p$p.value)
  res <- rbind(res, re)
}
write.csv(res, "PIK3CA_ssGSEA_Ttest_result.csv")
res[which(res[,2]<0.05),]


table(dat.df$label)
names(dat.df)[names(dat.df) == 'Activated CD8 T cell'] <- 'Activated_CD8_T_cell'
names(dat.df)[names(dat.df) == 'Central memory CD8 T cell'] <- 'Central_memory_CD8_T_cell'
names(dat.df)[names(dat.df) == 'Activated CD8 T cell'] <- 'Activated_CD8_T_cell'


# Shapiro-Wilk normality test for Men's weights
with(dat.df, shapiro.test(Activated_CD8_T_cell[label == "Cluster A"]))# p = 0.1901
# Shapiro-Wilk normality test for Women's weights
with(dat.df, shapiro.test(Activated_CD8_T_cell[label == "Cluster B"])) # p = 0.05518

res2 <- c()
for(i in 3:ncol(dat.df)){
  p <- wilcox.test(dat.df[,i] ~ dat.df$label, data = dat.df, var.equal = TRUE)
  re <- cbind(colnames(dat.df)[i], p$p.value)
  res2 <- rbind(res, re)
}
res2[which(res2[,2]<0.05),]




