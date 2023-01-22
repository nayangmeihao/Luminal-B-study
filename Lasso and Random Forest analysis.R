

#lasso anaysis
library(glmnet)
library(Matrix)
library(data.table)
library(openxlsx)

df <- read.xlsx("B_DEseq2_log2_data.xlsx", sheet = 1)  # 538*275
df <- as.data.frame(df)
rownames(df) <- df$sample
df <- subset(df, select = -sample)

data.outcome <- df[, 1]
data.x <- as.matrix(df[, 2:538])
set.seed(123456)
# Use cross-checking to observe model errors with different lambda values
cv.fit <- cv.glmnet(data.x, data.outcome, alpha = 1, family = "binomial", type.measure = "auc", nfolds = 10)
plot(cv.fit)
cv.fit$lambda.min
cv.fit$lambda.1se

# Add names of variables to each curve
library(plotmo)
pdf("lasso.pdf", width = 6, height = 6)
fit <- glmnet(data.x, data.outcome, alpha = 1, family = "binomial")
plot_glmnet(fit, var = "rlambda", label = 37, s = cv.fit$lambda.1se)
dev.off()

# lambda.1se
coefficients <- coef(cv.fit,s = cv.fit$lambda.1se)
Active.Index <- coefficients[which(coefficients != 0), ]
Active.Index <- round(Active.Index, 4)
length(names(Active.Index))
lasso.genes <- names(Active.Index)[2:38]
lasso.genes



# Random Forest
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

df <- read.xlsx("B_DEseq2_log2_data.xlsx", sheet = 1)  # 538*275
df <- as.data.frame(df)
rownames(df) <- df$sample
df <- subset(df, select = -sample)

lasso.genes <- c("MMP1", "FN1", "LRP1", "RSPO3", "LINGO1", "EEF1A2", "RASGRF2", "VGLL3", 
                 "RNF157", "CFH", "FCRLB", "SERPING1", "RAB3A", "COL5A2", "MORN3", "SLC16A6",
                 "CCDC24", "ACKR4", "EMX2", "EGFL6", "PRRT2", "CORIN", "FCGR3A", "NT5E", 
                 "GPR37L1", "DLL3", "NYAP1", "SLC39A6", "FAM234B", "COL3A1", "HIST4H4", "PAH", 
                 "FNDC1", "CACNA2D2", "KIF12", "SLC2A12", "UPK3B")
y <- c("label")
lasso.genes <- union(y, lasso.genes)
length(lasso.genes)
ran.dat <- df[ , lasso.genes]
ran.dat$label <- as.factor(ran.dat$label)

# Choose the optimal mtry parameter value
set.seed(123456)
errRate <- c(1)
for (i in 1:30) {  
  m <- randomForest(label~., data = ran.dat, mtry = i, proximity = TRUE)  
  err <- mean(m$err.rate)  
  errRate[i] <- err  
}  
print(errRate)
m  <-  which.min(errRate)  
print(m)

set.seed(123456)
ntree_fit <- randomForest(label~ ., data = ran.dat, mtry = 3, ntree = 5000)
plot(ntree_fit)

impor <- importance(df.rf)  
impor <- as.data.frame(impor)
# save result
pdf("randomforest_1.pdf", width = 8, height = 8, useDingbats = FALSE)
varImpPlot(df.rf, n.var = min(20, nrow(impor)), main = "Top 20 - variable importance") 
dev.off()

# Descending order of data according to MeanDecreaseGini
impor <- impor[order(impor[ , 4], decreasing = T), ] 
best.gini.gene <- impor[1:20, ]
gini.x <- rownames(best.gini.gene)
y <- c("label")
gini.genes <- union(y, gini.x)
length(gini.genes)


# Descending order of data according to MeanDecreaseAccuracy
Acc.impor <- impor[order(impor[, 3], decreasing = T), ] 
best.acc.gene <- Acc.impor[1:20, ]
acc.x <- rownames(best.acc.gene)
y <- c("label")
acc.genes <- union(y, acc.x)
length(acc.genes)

genes <- union(gini.x, acc.x)
y <- c("label")
union.genes.1 <- union(y, genes)
length(union.genes.1)


#plot Heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)

heatmap.data <- ran.dat[genes]

# label data, rownames(group.text) are samples, colnames(group.text) is label
group.text <- read.xlsx("cluster_2_res.xlsx", sheet = 1)
group.text <- as.data.frame(group.text)
rownames(group.text) <- group.text$sample
group.text <- subset(group.text, select = -sample)

all(rownames(group.text) == rownames(heatmap.data))

# scale
data.scale <- scale(heatmap.data)  # normalize genes
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))

data.scale.new <- na.omit(data.scale.new)  # remove the rows with na 
dat <- t(data.scale.new)

all(rownames(group.text)==colnames(dat))

# plot Heatmap
mat <- as.matrix(dat)
group.mat <-  as.matrix(group.text)

all(rownames(group.mat) == colnames(dat))
pdf("21gene_heatmap.pdf")
column.ha = HeatmapAnnotation(df = group.text,
                              col = list(Label = c("type1" =  "red", "type2" = "blue"))                            )
ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#7fb80e", "white", "red")),
              heatmap_legend_param = list(
                title = "legend", 
                title_position = "topcenter",
                legend_height = unit(6,"cm"), legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = T,
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

