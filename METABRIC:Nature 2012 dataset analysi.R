library(data.table)
library(openxlsx)
library(dplyr)

## Data pre-processing
# KIF12 is not included in this dataset
raw.data <- read.xlsx("mRNA expression (microarray).xlsx", sheet = 1)
rownames(raw.data) <- raw.data$SAMPLE_ID
raw.data <- subset(raw.data,select = -SAMPLE_ID)
ra.data <- as.data.frame(t(raw.data))
summary(as.vector(as.matrix(ra.data)))

# Outlier detection
r.data <- ra.data
set.seed(123)
summary(r.data$CCDC24)
plot(density(r.data$CCDC24)) # Print out the probability density function of data
boxplot(r.data$CCDC24) # print out boxplot
Outliers <- boxplot.stats(r.data$CCDC24)$out  # print out outlier

r.data$CCDC24 <- ifelse(r.data$CCDC24 %in% Outliers, NA, r.data$CCDC24)  # Replace outliers with empty
sub <- which(is.na(r.data$CCDC24))  # Identify the number of rows with missing values
inputfile1 <- r.data[-sub, ]  # Split the data set into two parts, complete data and missing data
inputfile2 <- r.data[sub, ]

# The mean replacement method handles the missing, and the results are transferred
# Idea: split into two, assign the missing value to the average, and then put it back together
avg.sales <- mean(inputfile1$CCDC24)  # Find the mean of the missing part of the variable
inputfile2$CCDC24 <- rep(avg.sales, nrow(inputfile2))  # Replace missing with mean
result2 <- rbind(inputfile1, inputfile2)  # Incorporate the completed data

raw.dat <- result2
# key genes
key.genes <- colnames(raw.dat)
raw.dat$PATIENT_ID <- rownames(raw.dat)
dim(raw.dat)  # 1904*22
head(raw.dat, n = 3)

# clinical data, rownames(clinical.data) are samples ,colnames(clinical.data) are genes
clinical.data <- read.xlsx("raw_clinical_data.xlsx", sheet = 1)
clinical.names <- colnames(clinical.data)
# head(clinical.data)

# join express data and clinical data by sample
merge.data <- left_join(raw.dat, clinical.data, by = "PATIENT_ID")
dim(merge.data)
head(merge.data, n = 2)

class.line <- median(B.data$CCDC24)
class.line
B.data$label <- ifelse(B.data$CCDC24 > class.line, "type1", "type2")
B.data <- B.data[order(B.data$label), ]
table(B.data$label)
head(B.data, n = 2)


## Survival analysis
library(data.table)
library(openxlsx)
library(dplyr)

sur.df <- read.xlsx("set1_set2_survival_data.xlsx", sheet = "set2")
head(sur.df)


library(openxlsx)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(data.table)

# sur.df <- B.data[, c("PATIENT_ID", "label", "OS_STATUS", "OS_MONTHS")]
sur.df$OS_STATUS <- ifelse(sur.df$os_event == "1:DECEASED", 1, 0)
data2 <- sur.df
kmfit <- survfit(Surv(data2$os_months, data2$os_event)~ data2$label, data = data2)

# plot survival curve
pdf("set1_KM_analysis.pdf", width = 6, height = 6)
gg <- ggsurvplot(kmfit, palette = c('red', 'blue'),
                 legend.title = "label", 
                 legend = c(0.85, 0.9),
                 pval = TRUE,
                 #pval = "P value: 0.001",
                 pval.coord = c(0.13, 0.03),
                 conf.int = TRUE,
                 #ylim = c(0.6, 1),
                 xlab = "Time in months",
                 risk.table = TRUE,
                 risk.table.height = 0.2,
                 ggtheme = theme( legend.title = element_blank(),
                                  panel.background = element_rect(fill = "white", color = "black"),
                                  legend.key = element_rect(fill = "transparent"),
                                  axis.title = element_text(color = "black", vjust = 2),
                                  axis.ticks = element_line(color = "black")))
gg
dev.off()


## Heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)

y <- c("PATIENT_ID", "label")
col.names <- union(y, key.genes)
B.dat <- B.df[, col.names]

# label
group.text <- B.dat[, c("PATIENT_ID", "label", "CCDC24")]
rownames(group.text) <- group.text$PATIENT_ID
group.text <- subset(group.text, select = -PATIENT_ID)
group.text$CCDC24 <- ifelse(group.text$CCDC24 > class.line, "high", "low")


# heatmap.data
c2 <- c('FCGR3A', 'MMP1', 'SERPING1', 'EGFL6', 'SLC2A12', 'CORIN', 'EMX2',
        'ACKR4', 'COL3A1', 'COL5A2', 'FN1', 'FNDC1', 'VGLL3', 'CFH', 'NT5E',
        'RASGRF2', 'LRP1', 'CACNA2D2', 'CCDC24', 'RAB3A')
heatmap.data <- subset(B.dat, select = -label)
rownames(heatmap.data) <- heatmap.data$PATIENT_ID
heatmap.data <- subset(heatmap.data, select = -PATIENT_ID)

heatmap.data <- heatmap.data[, c2]

all(rownames(group.text) == rownames(heatmap.data))

# scale
data.scale <- scale(heatmap.data)
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
# data.scale.new <- na.omit(data.scale.new)
dat <- t(data.scale.new)

all(rownames(group.text) == colnames(dat))

# plot heatmap
mat <- as.matrix(dat)
group.mat <-  as.matrix(group.text)
all(rownames(group.mat) == colnames(dat))

pdf("set1_heatmap.pdf", width = 6, height = 6)
column.ha <- HeatmapAnnotation(df = group.text,
                               col = list(label = c("type1" = "red", "type2" = "blue"),
                                          CCDC24 = c("low"=  "#ffce7b", "high"= "#f15a22")))

ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#7fb80e", "white", "red")),
              heatmap_legend_param = list(
                title = "legend", 
                title_position = "topcenter",
                legend_height = unit(6,"cm"), 
                legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = T,
              row_title_gp = gpar(fontsize = 20),       
              cluster_rows = F,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = T,
              column_split = group.mat,
              column_title = NULL)
ht
dev.off()


##Didderent analysis between ClusterA and ClusterB
study.data <- B.dat
rownames(study.data) <- study.data$PATIENT_ID
study.data <- subset(study.data, select = -PATIENT_ID)

data <- study.data
col_names <- colnames(data)
re <- c()
normalRes <- c()
for (i in col_names) {
  if(i == "label") {
    next
  }
  else {
    high <- data[which(data$label == "type2"), i]
    low <- data[which(data$label == "type1"), i]
    #  it is in a normal distribution. The larger the W, the more the data conforms to the normal distribution
    nor.high <- shapiro.test(high)  # If p>0.05
    nor.low <- shapiro.test(low)  
    x <- c(high, low)
    group <- c(rep("high", length(high)), rep("low",length(low)))  # Homogeneity test of variance
    bar <- bartlett.test(x~group)  # Variance homogeneity test, close to 1 indicates homogeneous variance
    re <- cbind(i, nor.high$p.value, nor.low$p.value, bar$p.value)
    normalRes <- rbind(normalRes, re)    
  }
}
colnames(normalRes) <- c("genes", "type2.p", "type1.p", "bar.p")
normalRes <- as.data.frame(normalRes)

normalRes$genes <- as.character(normalRes$genes)
normalRes$type2.p <- as.character(normalRes$type2.p)
normalRes$type1.p <- as.character(normalRes$type1.p)
normalRes$bar.p <- as.character(normalRes$bar.p)
dim(normalRes)
head(normalRes)


# User a t-test on data that conforms to a normal distribution
res1 <- normalRes[normalRes[, "type2.p"] > 0.05, ]
res2 <- res1[res1[, "type1.p", ] > 0.05, ]
tRes <- res2[res2[, "bar.p", ] > 0.7, ]
dim(tRes)

# t_test, type1-type2
data <- study.data
dif <- c()
j <- levels(factor(data$label))
for (i in tRes$genes) {
  ttest <- t.test(data[, i] ~ data$label)  # data$group = data[, 1]
  log2FC <- mean(data[data$label == j[1], i]) - mean(data[data$label == j[2], i])
  dif <- rbind(dif, cbind(i, ttest$p.value, log2FC))
}
colnames(dif) <- c("gene", "pvalue", "log2FC")
dif <- as.data.frame(dif)
dif$pvalue <- as.numeric(as.character(dif$pvalue))
dif$adj.pval <- p.adjust(dif$pvalue, method = "BH")
dif <- dif[order(dif$adj.pval), ]
t.dif <- dif[which(dif$adj.pval <= 0.05), ]
dim(t.dif)
t.dif


# U test for data that does not conform to the normal distribution
uRes <- setdiff(normalRes$genes, tRes$genes)
length(uRes)  # 196
y <- c("label")
uRes.1 <- union(y, uRes)
remain.data <- study.data[, uRes.1]
dim(remain.data) 
# head(remain.data, n=2)

# U_test
data <- remain.data
dif <- c()
col_names <- colnames(data)
j <- levels(factor(data$label))
for (i in 2:ncol(data)) {
  high <- data[which(data$label == "type2"), i]
  low <- data[which(data$label == "type1"), i]
  Utest <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE, correct = FALSE, conf.level = 0.95)
  log2FC <- mean(data[data$label == j[1], i]) - mean(data[data$label == j[2], i])
  dif <- rbind(dif, cbind(col_names[i], Utest$p.value, log2FC))
}

colnames(dif) <- c("gene", "pvalue", "log2FC")
dif <- as.data.frame(dif)
dif$pvalue <- as.numeric(as.character(dif$pvalue))
dif$log2FC <- as.numeric(as.character(dif$log2FC))
dif$adj.pval <- p.adjust(dif$pvalue, method = "BH")
dif <- dif[order(dif$adj.pval), ]

u.dif <- dif[which(dif$adj.pval <= 0.05), ]
dim(u.dif)
u.dif


final.genes.t.u <- union(t.dif$gene, u.dif$gene)  #  45
length(final.genes.t.u)
final.genes.t.u


### Fisher test of Clinical fasts
library(MASS)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(data.table)

clinical.data <- read.xlsx("type1_type2_clinical_data.xlsx")
names(clinical.data)[names(clinical.data) == "LYMPH_NODES_EXAMINED_POSITIVE"] <- 'LYMPH_NODES'
names(clinical.data)[names(clinical.data) == "AGE_AT_DIAGNOSIS"] <- 'AGE'

# Data processing
clinical.data$NPI <- ifelse(clinical.data$NPI> median(clinical.data$NPI), 1, 0)
clinical.data$AGE <- ifelse(clinical.data$AGE> 60, 1, 0)

car.df2 <- table(clinical.data$label, clinical.data$LYMPH_NODES) 
car.df2
fisher.test(car.df2)


### Clinical factor cox
library(openxlsx)
library(data.table)

clinical.df <- read.xlsx("type1_type2_clinical_data.xlsx")
names <- c("PATIENT_ID", "OS_MONTHS", "OS_STATUS", "HER2_SNP6", "INTCLUST", "THREEGENE")
clinical.df$HER2_SNP6 <- ifelse(clinical.df$HER2_SNP6  == "null", NA, clinical.df$HER2_SNP6)
clinical.df$INTCLUST <- ifelse(clinical.df$INTCLUST  == "null", NA, clinical.df$INTCLUST)
clinical.df$THREEGENE <- ifelse(clinical.df$THREEGENE  == "null", NA, clinical.df$THREEGENE)
clinical.cox.dat <- clinical.df[, names]

library(survival)
library(survminer)
library(plyr)
library(survivalROC)
library(grid)
library(MASS)
library(leaps)

dd <- clinical.cox.dat
UniCox <- function(x) {
  FML <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS)~", x))
  rex.cox <- coxph(FML, data = dd)
  sum <- summary(rex.cox)
  y <- cbind(sum$coefficients, sum$conf.int)
  return(y)
}
VarNames <- colnames(dd)[c(4:6)]  
UniVar <- lapply(VarNames, UniCox)  
UniVar <- ldply(UniVar, data.frame)  
factor <- VarNames
results <- cbind(factor, UniVar)
results <- data.frame(factor = results[, 1], HR = results[, 3], 
                      z = results[, 5], pvalue = results[, 6],
                      lower95_CI = results[, 9], upper95_CI = results[, 10])
results1 <- arrange(results, results[, 4])  
dim(results1)
results1


# Univariate cox  of gene 
cox.data <- read.xlsx("cox_data.xlsx", sheet = 1)

cox.med.df <- cox.data
for (i in colnames(cox.med.df)){
  if (i %in% c("sample", "label", "OS_months", "OS_status")) {
    print(i)
  }
  else{
    cox.med.df[, i] <- ifelse(cox.med.df[, i] > median(cox.med.df[, i]), 1, 0)
  }
}
# head(cox.med.df, n = 2)
# write.table(cox.med.df,"cox_data2.csv",row.names=TRUE,col.names=TRUE,sep=",")

library(survival)
library(survminer)

dd <- cox.med.df
UniCox <- function(x) {
  FML <- as.formula(paste0("Surv(OS_months, OS_status)~", x))
  rex.cox <-coxph(FML, data = dd)
  sum <- summary(rex.cox)
  y <- cbind(sum$coefficients, sum$conf.int)
  return(y)
}

VarNames <- colnames(dd)[c(5:24)] 
UniVar <- lapply(VarNames,UniCox)  
UniVar <- ldply(UniVar,data.frame)  
factor <- VarNames
results <- cbind(factor, UniVar)
results <- data.frame(factor = results[, 1], HR = results[, 3], 
                      z = results[, 5], pvalue = results[, 6],
                      lower95_CI = results[, 9], upper95_CI = results[, 10])

results1 <- arrange(results, results[, 4]) 
dim(results1)
results1
# write.table(results, "Genes_single_cox_result.csv", row.names = TRUE, col.names = TRUE, sep = ",")


### Multivariate  cox of gene
# use four key genes to conxstruct  multivariable cox model
cox.gene <- c("CCDC24", "NT5E", "COL5A2", "CORIN", "FN1")
fml <- as.formula(paste0("Surv(OS_months, OS_status)~", paste0(cox.gene, collapse = '+')))
mode <- coxph(fml, data = cox.med.df)
Multisum <- summary(mode)
Multisum


### Forestplot
# plot1
library(forestplot)
rs.forest <- read.csv("forest_data1.csv", header = FALSE)

pdf("set1.1.forestplot.pdf", width = 6, height = 6)
forestplot(labeltext = as.matrix(rs.forest[, 1:4]),
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 2:5, col = "#000044")),
           mean = rs.forest$V5, 
           lower = rs.forest$V6,
           upper = rs.forest$V7, 
           is.summary = c(T, T, rep(F, 7)),
           zero = 1,
           boxsize = 0.4,  
           lineheight = unit(8, "mm"),  
           colgap = unit(2, "mm"),
           lwd.zero = 2, 
           lwd.ci = 2, 
           col = fpColors(box = "#458B00", summary = "#8B008B", lines = "black", zero = "#7AC5CD", hrz_lines = "#444444"),
           xlab = "The estimates", 
           lwd.xaxis = 2, 
           lty.ci = "solid",
           graph.pos = 4 
)
dev.off()


# plot2
library(forestplot)
rs.forest <- read.csv("forest_data2.csv", header = FALSE)

df("set1.2.forestplot.pdf", width = 6, height = 6)
forestplot(labeltext = as.matrix(rs.forest[, 1:3]),
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044")),
           mean = rs.forest$V4, 
           lower = rs.forest$V5,
           upper = rs.forest$V6, 
           is.summary = c(T, T, rep(F, 7)),
           zero = 1,
           boxsize = 0.4,  
           lineheight = unit(8, "mm"),  
           colgap = unit(2, "mm"),
           lwd.zero = 2, 
           lwd.ci = 2, 
           col = fpColors(box = "#458B00", summary = "#8B008B", lines = "black", zero = "#7AC5CD", hrz_lines = "#444444"),
           xlab = "The estimates", 
           lwd.xaxis = 2, 
           lty.ci = "solid",
           graph.pos = 3 )
dev.off()


### ki67 difference test
library(data.table)
library(openxlsx)
library(dplyr)
ki67.data <- read.csv("ki67_data.csv")

data <- ki67.data[, c("MKI67", "label")]
high <- data[which(data$label == "type2"), "MKI67"]
low <- data[which(data$label == "type1"), "MKI67"]

#  it is in a normal distribution. The larger the W, the more the data conforms to the normal distribution
nor.high <- shapiro.test(high)  # If p>0.05
nor.low <- shapiro.test(low)  
x <- c(high, low)
group <- c(rep("high", length(high)), rep("low",length(low)))  # Homogeneity test of variance
bar <- bartlett.test(x~group)  # Variance homogeneity test, close to 1 indicates homogeneous variance
re <- cbind(i, nor.high$p.value, nor.low$p.value, bar$p.value)
colnames(re) <- c("type2.p", "type1.p", "bar.p")
normalRes <- as.data.frame(re)
# normalRes$genes <- as.character(normalRes$genes)
normalRes$type2.p <- as.character(normalRes$type2.p)
normalRes$type1.p <- as.character(normalRes$type1.p)
normalRes$bar.p <- as.character(normalRes$bar.p)
dim(normalRes)
head(normalRes)

# U_test
data <- ki67.data[, c("MKI67", "label")]
col_names <- colnames(data)
j <- levels(factor(data$label))

high <- data[which(data$label == "type2"), "MKI67"]
low <- data[which(data$label == "type1"), "MKI67"]
Utest <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE, correct = FALSE, conf.level = 0.95)
Utest

log2FC <- mean(data[data$label == j[1], "MKI67"]) - mean(data[data$label == j[2], "MKI67"])
log2FC


## Heatmap Line and Scatter plot of genes
library(data.table)
library(openxlsx)
data <- read.csv('cox_mode_data.csv')
head(data, n = 2)

library(ComplexHeatmap)
library(circlize)
library(grid)

df <- data[, c("sample", "label", "risk.label", "CCDC24", "NT5E", "COL5A2", "CORIN", "FN1")]

group.text <- df[, c("sample", "label", "risk.label")]
rownames(group.text) <- group.text$sample
group.df <- subset(group.text, select = -sample)
group.df$risk.label <- as.factor(group.df$risk.label)
group.mat <- as.matrix(group.df)

da <- subset(df, select = -c(risk.label, label))
rownames(da) <- da$sample
da <- subset(da, select = -sample)

data.scale <- scale(da)
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
dat <- t(data.scale.new)

# plot Heatmap
mat <- as.matrix(dat)

column.ha <- HeatmapAnnotation(df = group.df,
                               col = list(risk.label = c("Low risk" = "#ffce7b", "High risk" = "#f15a22"),
                                          label = c( "type1" = "red","type2" = "blue")))
ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#7fb80e", "white", "red")),
              heatmap_legend_param = list(title = "legend", title_position = "topcenter",
                                          legend_height = unit(6, "cm"), legend_direction = "horizontal"),
              bottom_annotation = column.ha,
              show_row_dend = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              #column_split = group_mat,
              column_title = NULL)

draw(ht)
pdf("set1_3-1.pdf", width = 9, height = 3)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom",  annotation_legend_side = "bottom")
dev.off()

library(ggplot2)

p2.data <- data[, 1:6]
p2.data$OS_status <- factor(p2.data$OS_status)

pdf("set1_3-2.pdf", width = 9, height = 3, useDingbats = FALSE)

p2 <- ggplot(p2.data, aes(x = index, y = OS_months, color = OS_status)) + geom_point(shape = 16)
p2.1 <- p2 + scale_x_discrete("") + theme(axis.title.x = element_blank(),
                                          axis.title.y = element_blank(),
                                          axis.text.x = element_blank(),
                                          axis.ticks.x = element_blank(),
                                          legend.position = c(0.95, 0.12),
                                          legend.title = element_blank(),
                                          legend.background = element_rect(fill = "transparent", colour = NA),           
                                          panel.background = element_rect(fill = "white", color = "black"), 
                                          axis.title = element_text(color = "black", vjust = 2),
                                          axis.ticks = element_line(color = "black"))

p2.1
dev.off()


library(gcookbook)
library(ggplot2)

p3.data <- data[,c("index", "mode", "risk.label")]
pdf("set1_3-3.pdf", width = 9, height = 3)
p3 <- ggplot(data = p3.data, aes(x = index, y = mode)) + geom_line(aes(color = risk.label))
p3.1 <- p3 + scale_x_discrete("") + theme(axis.title.x = element_blank(), 
                                          axis.title.y = element_blank(),
                                          axis.text.x = element_blank(),
                                          axis.ticks.x = element_blank(),
                                          legend.title = element_blank(),
                                          panel.background = element_rect(fill = "white", color = "black"),
                                          legend.position = c(0.1, 0.8),
                                          axis.title = element_text(color = "black", vjust = 2),
                                          axis.ticks = element_line(color = "black")) + labs(y = "predictive scores")
p3
dev.off()


