
# Protein analysis
library(data.table)
library(openxlsx)
library(dplyr)

names.data <- read.xlsx("join_name_data.xlsx")  # 187*3
# level4 data
L4.data <- read.csv("level4/TCGA-BRCA-L4.csv")  # 901*224
L4.data$Sample_ID <- as.character(L4.data$Sample_ID)
rownames(L4.data) <- L4.data$Sample_ID
L4.data <- subset(L4.data, select = -Sample_ID)

name.data <- read.xlsx("level4/protein_names.xlsx")  
col.names <- name.data$gene
L4.dat <- L4.data[, col.names]

summary(as.vector(as.matrix(L4.dat)))
# survival data
sur.df <- read.xlsx("cluster_2_sur_dat.xlsx") # 187*4


library(ggpubr)
boxplotFun <- function(s.dat, group) { 
  df <- data.frame(gene = s.dat, stage = group)  # Combine group and data together
  p <- ggboxplot(df, x = "stage", y = "gene",  # x axis
                 color = "stage", palette = "jco",add = "jitter")
  p + stat_compare_means(method = "t.test")  # Add p-value
}


# intersect pationts
length(intersect(sur.df$sample, rownames(L4.dat)))
inter.id <- intersect(sur.df$sample, rownames(L4.dat))

m.df <- L4.dat[inter.id, ]
m.df$sample <- rownames(m.df)

# join survival data and protein data by sample
protein.df <- left_join(m.df, sur.df, by = "sample")
dim(protein.df)
head(protein.df, n = 2)
write.table(protein.df, "protein.df.csv", row.names = TRUE, col.names = TRUE, sep = ",")



# Evaluate the overall distribution of the data
protein.df <- read.csv("level4/protein.df.csv")  # 187*179
protein.df.copy <- protein.df
study.data <- protein.df[, 1:176]
rownames(study.data) <- study.data$sample
# rownames(protein.dat) are samples, colnames(protein.dat) are protein
study.data <- subset(study.data, select = -sample)  # 160*179

# Process the missing value of the data according to the missing rate of each protein
studyData1 <- study.data
summary(as.vector(as.matrix(studyData1)))

# input: data;    rownames(data) are samples, colnames(data) are genes
#  First, count the missing rate in each column, 
#  If the missing rate is less than 0.7, replace the missing value with the mean of valid data in each column
#  If the missing rate of a column is greater than 0.7, delete the column
# return  Processed data

MissValueProcess <- function(data) {
  res <- c()
  for (i in colnames(data)) {
    number.na <- sum(is.na(data[, i]))
    percentage <- number.na/nrow(data)
    res <- rbind(res, cbind(i, number.na, percentage))
  }
  colnames(res) <- c("gene", "number.na", "percentage")
  res <- as.data.frame(res)
  res$number.na <- as.numeric(as.character(res$number.na))
  res$percentage <- as.numeric(as.character(res$percentage))
  rownames(res) <- res$gene
  res <- subset(res, select = -gene)
  data.copy <- data
  for (j in rownames(res)) {
    if(res[j, "percentage"] == 0) {
      next
    } 
    if(res[j, "percentage"] < 0.3) {
      data.copy[, j][is.na(data.copy[, j])] <- 0
      mean <- sum(data.copy[, j]) / (nrow(data.copy) - res[j, "number.na"])
      data[, j][is.na(data[, j])] <- mean
    }
    else {
      data <- data[, -which(colnames(data) == j)]
    }
  }
  out <- list(one = res, two = data)
  return(out)
}

result <- MissValueProcess(studyData1) 
proteinWithoutNa <- result$two
write.table(proteinWithoutNa, "level4/proteinWithoutNa.csv", row.names = TRUE, col.names = TRUE, sep = ",")

# protein.df
summary(as.vector(as.matrix(proteinWithoutNa)))


# View the distribution of all proteins in a sample
studyData2 <- as.data.frame(t(proteinWithoutNa))
# input: data;    rownames(data) are genes, colnames(data) are sample
#  Count the minimum, mean, median, maximum, and variance of all gene expressions of each patient 
# return  Processed data

ExpressionDistribution <- function(data) {
  res <- c()
  for (i in colnames(data)) {
    series <- data[, i]
    min.value <- min(series)
    mean.value <- mean(series)
    median.value <- median(series)
    max.value <- max(series)
    var.value <- var(series)
    res <- rbind(res, cbind(i, min.value, mean.value, median.value, max.value, var.value))
  }
  res <- as.data.frame(res)
  colnames(res) <- c("sample", "min.value", "mean.value", "median.value", "max.value", "var.value")
  return(res)
}

res <- ExpressionDistribution(studyData2)
min(as.numeric(as.character(res$var.value)))


## Difference analysis
protein.df <- read.csv("level4/protein.df.csv")  
proteinWithoutNa <- read.csv("level4/proteinWithoutNa.csv")  
studyData3 <- proteinWithoutNa
studyData3$label <- protein.df$label  # 160*176

library(limma)
# Note that the limma package is an analysis method for the gene chip expression matrix. 
# It is not possible to analyze the reversed RNAseq expression matrix (because of the different data characteristics). 
# RNAseq requires DEseq2.

# Expression matrix, colnames(exprSet) are samples, rownames(exprSet) are genes
exprSet <- t(studyData3[, 1:175])

# Construct experimental design matrix
group <- studyData3$label
design <- model.matrix(~0 + factor(group))
colnames(design) <- c("type1", "type2")
rownames(design) <- colnames(exprSet)

# Construct a comparative model to compare the expression data under two experimental conditions
cont.matrix <- makeContrasts("type2-type1", levels = design)
cont.matrix  # This matrix statement, we compare the difference between type1 and type2

# difference analysis
# step1: Linear model fitting
fit <- lmFit(exprSet, design)

# step2: Calculate the difference according to the comparison model
fit2 <- contrasts.fit(fit, cont.matrix)

# step3: Bayesian test
fit2 <- eBayes(fit2)

# step4: Generate test results report for all genes
tempOutput <- topTable(fit2, coef = 1, n = Inf, adjust = "BH")
write.table(tempOutput, "level4/limma.result.csv", row.names = TRUE, col.names = TRUE, sep = ",")

# step5: Screen with adj.P.Val to get all differentially expressed genes
limma.dif <- tempOutput[tempOutput[, "adj.P.Val"] < 0.05, ]
dim(limma.dif)
head(limma.dif, n = 2)
diff.protein <- rownames(limma.dif)
length(diff.protein)

# culculate foldchange
data <- studyData3
fc <- c()
j <- levels(factor(data$label))
for (i in colnames(data)) {
  if(i %in% c("label", "sample")){
    next
  }
  else {
    log2FC <- mean(data[data$label == j[2], i]) - mean(data[data$label == j[1], i])  #type2 - type1
    fc <- rbind(fc, cbind(i, log2FC))
  }   
}
colnames(fc) <- c("gene", "log2FC")
fc <- as.data.frame(fc)
fc$log2FC <- as.numeric(as.character(fc$log2FC))
rownames(fc) <- fc$gene
fc <- subset(fc, select = -gene)




## Enrichment analysis
# join data 
nameData <- read.xlsx("level4/protein_join_gene_name.xlsx") 
limma.res <- read.csv("level4/limma.result.csv") 
limma.res$protein <- rownames(limma.res)

enrich.data <- left_join(limma.res, nameData, by = c("protein" = "RPPA_level4_raw_protein"))
names(enrich.data)[names(enrich.data) == "RPPA_level4_raw_protein_name"] <- "raw_protein_name"
names(enrich.data)[names(enrich.data) == "RPPA_level4_gene_names"] <- "gene_names"
write.xlsx(enrich.data, "level4/Enrichment/enrichData.xlsx")



## Heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)

# final.genes <- final.genes.t.u
final.genes <- diff.protein
y <- c("sample", "label")
final.genes.1 <- union(y, final.genes)
studyData3$sample <- rownames(studyData3)
p.dat <- studyData3[, final.genes.1]
p.dat$sample <- p.dat$sample
dim(p.dat)
head(p.dat, n = 2)

## plot 2 ##
# final.genes <- final.genes.t.u
final.genes <- c("DJ1", "YAP_pS127", "MAPK_pT202Y204", "FIBRONECTIN", "ERALPHA", 
                 "ARAF_pS299", "YB1", "PKCALPHA_pS657", "CLAUDIN7", "P27", "SRC_pY416",
                 "MEK1_pS217S221", "MYOSINIIA", "AKT_pS473", "SHP2_pY542", "YB1_pS102", 
                 "BAP1C4", "RAB11", "S6_pS235S236")
y <- c("sample", "label")
final.genes.1 <- union(y, final.genes)
studyData3$sample <- rownames(studyData3)
p.dat <- studyData3[, final.genes.1]
p.dat$sample <- p.dat$sample

# label
group.text <- p.dat[, c("sample", "label")]
rownames(group.text) <- group.text$sample
group.text <- subset(group.text, select = -sample)

# heatmap data
heatmap.data <- subset(p.dat, select = -c(label))
heatmap.data <- subset(heatmap.data, select = -sample)

all(rownames(group.text) == rownames(heatmap.data))
# scale
data.scale <- scale(heatmap.data)
# data.scale <- as.matrix(heatmap.data)
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))

# not do scale
# data.scale.new <- as.matrix(heatmap.data)

# remove rows with NAs
# data.scale.new <- na.omit(data.scale.new)
dat <- t(data.scale.new)

all(rownames(group.text)==colnames(dat))

pdf("protein_heatmap.pdf")
# plot Heatmap
mat <- as.matrix(dat)

column.ha <- HeatmapAnnotation(df = group.text,
                               col = list(label = c("type1" = "red", "type2" = "blue"))
)

ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#33a3dc", "white", "#f36c21")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter",
                                          legend_height = unit(6,"cm"), 
                                          legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = TRUE,
              #right_annotation = row_anno,
              row_title_gp = gpar(fontsize = 20),       
              cluster_rows = TRYE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              #column_split = group.mat,
              column_title = NULL
)
ht
dev.off()


## Randomforest
final.genes.1 <- c("DJ1", "YAP_pS127", "MAPK_pT202Y204", "FIBRONECTIN", "ERALPHA", "ARAF_pS299", 
                   "YB1", "PKCALPHA_pS657", "CLAUDIN7", "P27", "SRC_pY416", "MEK1_pS217S221", 
                   "MYOSINIIA", "AKT_pS473", "SHP2_pY542", "YB1_pS102", "BAP1C4", "RAB11", "S6_pS235S236")
y <- c("label")
final.genes.1 <- union(y, final.genes.1)
ran.dat <- studyData3[, final.genes.1]

library(party)
library(grid)
library(mvtnorm)
library(modeltools)
library(stats4)
library(strucchange)
library(zoo)
library(randomForest)
library(data.table)
library(openxlsx)

# Choose the optimal mtry parameter value
set.seed(123456)
errRate <- c(1)
for (i in 1:16) {  
  m <- randomForest(label~., data = ran.dat, mtry = i, proximity = TRUE)  
  err <- mean(m$err.rate)  
  errRate[i] <- err  
}  
print(errRate)

# Choose m with the smallest average error  
m  <-  which.min(errRate)  
print(m)

# Select the value of ntree, the default is 500, select the size of ntree when the error is basically stable, ntree=2000
ntree_fit <- randomForest(label~ ., data = ran.dat, mtry = m, ntree = 5000)
plot(ntree_fit)

set.seed(123456)
df.rf <- randomForest(as.factor(label)~., data = ran.dat, mtry = 15, importance = TRUE, ntree = 4000, proximity = TRUE)
impor <- importance(df.rf)  # The larger the value, the greater the importance of the variable
impor <- as.data.frame(impor)

write.table(impor, "RF_result.csv", row.names = TRUE, col.names = TRUE, sep = ",")

pdf("protein_random.pdf", width = 6, height = 6)
varImpPlot(df.rf, n.var = min(34, nrow(impor)), main = "Variable importance") 
dev.off()


## boxplot
library(data.table)
library(openxlsx)
library(dplyr)
group.text <- read.xlsx("cluster_2_res.xlsx", sheet = 1)
# rownames(group.text) <- group.text$sample
# group.text <- subset(group.text, select = -sample)

mRNA.df <- read.xlsx("B_DEseq2_log2_data.xlsx", sheet = 1)  # 538*275
# rownames(mRNA.df) <- mRNA.df$sample
# mRNA.df <- subset(mRNA.df, select = -sample)
mRNA.data <- mRNA.df[, c("sample", "FN1", "NRG1", "CDH2", "CCND1")]
mRNA.data$variable <- "mRNA"
mRNA.dat <- left_join(mRNA.data, group.text, by = "sample")

protein.df <- read.csv("level4/proteinWithoutNa.csv")  
protein.df$sample <- rownames(protein.df)
protein.data <- protein.df[, c("sample", "FIBRONECTIN", "HEREGULIN", "NCADHERIN", "CYCLIND1")]
protein.data$variable <- "protein"
protein.dat <- left_join(protein.data, group.text, by = "sample")

pdf("boxplot_CCND1.pdf", width = 9, height = 9)
p1 <- ggplot(mRNA.dat, aes(x = Label, y = CCND1, color = Label)) + geom_boxplot() + 
  scale_color_manual(values = c('brown', 'steelblue'))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.8, size = 1)+
  labs(title = 'Gene expression of type1 and type2 on CCND1', x = 'Label', y = 'vaIue') + 
  annotate("text", x= 1.5, y = 19.2, size= 3, parse = TRUE, label = "'P = 6.94E-11'")+theme(legend.position = "none") 

# "FIBRONECTIN","HEREGULIN","NCADHERIN","CYCLIND1" 
p2 <- ggplot(protein.dat, aes(x = Label, y = CYCLIND1, color = Label)) + geom_boxplot() + 
  scale_color_manual(values = c('brown', 'steelblue'))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.8, size = 1)+
  labs(title = 'Protein expression of type1 and type2 on CCND1', x = 'Label', y = 'vaIue') + 
  annotate("text", x= 1.5, y = 1.25, size= 3, parse = TRUE, label = "'P = 0.182'")

library(ggpubr)
ggarrange(p1,p2,ncol=2, nrow = 1,labels=c("A","B"))
dev.off()


## Sankey plot
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggalluvial)
plot.data <- read.xlsx("result.xlsx", sheet = "all_protein") 
names(plot.data)[names(plot.data) == "raw_protein"] <- "protein_names"

gene.data <- subset(plot.data, select = -RPPA_raw_protein_name)
# gene.data$label <- "withDiff"
gene.data$label <- ifelse(gene.data$gene_names %in% c("FN1", "NRG1", "CDH2", "CCND1"), "withDiff", "withoutDiff")
gene.data$num <- ifelse(gene.data$label == "withDiff", 4, 171)
gene.data$group <- "gene"
gene.data$id <- seq(1, 175)
head(gene.data)

key.protein <- c("DJ1", "MYOSINIIA_pS1943", "YAP_pS127", "AXL", "MAPK_pT202Y204", "FIBRONECTIN",
                 "ERALPHA", "BRAF_pS445", "ANNEXIN1", "ARAF_pS299", "YB1", "P16INK4A", "PKCALPHA_pS657", 
                 "CLAUDIN7", "P27", "TRANSGLUTAMINASE", "PREX1", "P53", "BRCA2", "SRC_pY416", "MEK1_pS217S221", 
                 "MYOSINIIA", "AKT_pS473", "PKCALPHA", "SYK", "GCN5L2", "SHP2_pY542", "YB1_pS102", "BAP1C4", 
                 "RAB11", "S6_pS235S236")
protein.data <- subset(plot.data, select = -RPPA_raw_protein_name)
protein.data$label <- ifelse(protein.data$protein_names %in% key.protein, "withDiff", "withoutDiff")
protein.data$num <- ifelse(protein.data$label == "withDiff", 31, 144)
protein.data$group <- "protein"
protein.data$id <- seq(1, 175)
head(protein.data)

pData <- rbind(gene.data, protein.data)
dim(pData)
head(pData)
table(pData$group)

levels(pData$label) <- rev(levels(pData$label))
head(pData)

# plot
pdf("sankeyplot.pdf", width = 10, height = 10)
ggplot(pData,
       aes(x = group, stratum = label, alluvium = id, y = num,
           fill = label, label = label)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  ggtitle("Differences in gene and protein expression in type1 and type2")
dev.off()


##  Survival analysis
protein.df <- read.csv("level4/protein.df.csv")  
sub.protein.df <- protein.df[, c("sample", "label", "OS_months", "OS_status")]
sub.protein.df$sample <- as.character(sub.protein.df$sample)
head(sub.protein.df)

proteinWithoutNa <- read.csv("D:/study_data/protein/level4/proteinWithoutNa.csv")  
proteinWithoutNa$sample <- rownames(proteinWithoutNa)
coxData <- left_join(proteinWithoutNa, sub.protein.df, by = "sample")
final.genes <- c("DJ1", "YAP_pS127", "MAPK_pT202Y204", "FIBRONECTIN", "ERALPHA", "ARAF_pS299", 
                 "YB1", "PKCALPHA_pS657", "CLAUDIN7", "P27", "SRC_pY416", "MEK1_pS217S221", 
                 "MYOSINIIA", "AKT_pS473", "SHP2_pY542", "YB1_pS102", "BAP1C4", "RAB11", "S6_pS235S236")
y <- c("sample", "label", "OS_months", "OS_status")
col.names <- union(y, final.genes)
expre.dat <- coxData[, col.names]
head(expre.dat, n = 2)


# Transform continuous data to discrete data
cox.med.df <- expre.dat

lab <- c("sample", "OS_months", "OS_status", "label")
for (i in colnames(cox.med.df)) {
  if (i %in%lab) {
    print(i)
  }
  else {
    cox.med.df[, i] <- ifelse(cox.med.df[,i] > median(cox.med.df[, i]), 1, 0)
  }
}
head(cox.med.df, n = 2)

res <- data.frame()
for (i in 1:length(final.genes)) {
  group <- cox.med.df[, final.genes[i]]
  surv <- as.formula(paste("Surv(cox.med.df$OS_months, cox.med.df$OS_status)~", "group"))
  x  <-  survdiff(surv, data = cox.med.df)
  pValue <- 1 - pchisq(x$chisq, df = 1) 
  res[i, 1]  <-  final.genes[i]
  res[i, 2]  <-  pValue
}
res

library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(survMisc)

# supremum (Renyi) test; two-sided; two covariate groups
g1 <- ten(Surv(cox.med.df$OS_months, cox.med.df$OS_status) ~ cox.med.df$FIBRONECTIN, data = cox.med.df)
comp(g1)



























