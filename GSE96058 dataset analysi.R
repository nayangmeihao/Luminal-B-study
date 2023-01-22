
library(data.table)
library(openxlsx)
library(dplyr)

raw.data <- read.xlsx("set3.xlsx", sheet = 1)  # colnames(data) are sample, rownames are genes
rownames(raw.data) <- raw.data$sample
raw.data <- subset(raw.data, select = -sample)  # 22*3409
ra.data <- as.data.frame(t(raw.data))
key.genes <- colnames(ra.data)  # key genes
ra.data$title <- rownames(ra.data)


r.data <- ra.data
set.seed(123)
summary(r.data$CCDC24)

# Outlier recognition
par(mfrow = c(1,2))  # Divide the drawing window into one row and two columns, and display two pictures at the same time
dotchart(r.data$CCDC24)  # Draw univariate scatterplot, Dolan plot
pc <- boxplot(r.data$CCDC24, horizontal=F)  # Draw a box plot
Outliers <- boxplot.stats(r.data$CCDC24)$out  
length(Outliers)


# Outlier handling
r.data$CCDC24 <- ifelse(r.data$CCDC24 %in% Outliers, NA, r.data$CCDC24)
sub <- which(is.na(r.data$CCDC24))
inputfile1 <- r.data[-sub, ]
inputfile2 <- r.data[sub, ]
r.data$CCDC24 [is.na(r.data$CCDC24 )] <- mean(inputfile1$CCDC24)

ra.data <- r.data
key.genes <- colnames(ra.data)  # key genes
ra.data$title <- rownames(ra.data)
dim(ra.data)  # 3409*23



## Select luminal B data
library(data.table)
library(openxlsx)
library(dplyr)

clinical.data <- read.xlsx("D:/study_data/set3/clinical_data.xlsx")
clinical.data <- as.data.frame(clinical.data)
clinical.names <- colnames(clinical.data)
dim(clinical.data)  # 3273*17

b.data <- clinical.data[which(clinical.data$subtype == "LumB"), ]  # 729*17
dim(b.data)
B.data <- left_join(b.data, ra.data, by = "title")
dim(B.data)
head(B.data, n = 3)
# write.table(B.data, "B.data.csv", row.names = TRUE, col.names = TRUE, sep = ",")

med <- median(B.data$CCDC24)
med

B.data$label<- ifelse(B.data$CCDC24 > med , "type1" , "type2")
B.data <- B.data[order(B.data$label),  ]
table(B.data$label)
head(B.data, n = 3)


# Survival analysis
library(ggtext)
library(data.table)
library(openxlsx)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(magrittr)

my.data <- read.csv("my.data.csv")
data2 <- my.data[, c("label", "os_months", "os_event")]
y <- Surv(data2$os_months, data2$os_event)
kmfit <- survfit(y~ data2$label, data = data2)

pdf("set3_sur.pdf", width = 6, height = 6)
gg <- ggsurvplot(kmfit, palette = c("red", "blue"),
                 legend.title = "label", 
                 legend = c(0.85, 0.9),
                 #                  pval =TRUE,
                 pval = "P value: 0.017",
                 pval.coord = c(0.13, 0.62),
                 conf.int = TRUE,
                 ylim = c(0.6, 1),
                 xlab = 'Time in months',
                 risk.table = TRUE,
                 risk.table.height = 0.2,
                 ggtheme = theme(legend.title = element_blank(),
                                 panel.background = element_rect(fill = "white", color = "black"),
                                 legend.key = element_rect(fill = "transparent"),
                                 axis.title = element_text(color = "black", vjust = 2),
                                 axis.ticks = element_line(color = "black") ))
gg
dev.off()

# supremum (Renyi) test; two-sided; two covariate groups
g1 <- ten(Surv(data2$os_months, data2$os_event) ~ data2$label, data = data2)
comp(g1)



##Heatmap
library(ComplexHeatmap)
library(circlize)
library(grid)
library(data.table)
library(openxlsx)

my.data <- read.csv("my.data.csv")

y <- c("X","label")
col.names <- union(y, key.genes)
B.dat <- my.data[, col.names]

# label
group.text <- B.dat[, c("label","CCDC24")]
group.text$CCDC24 <- ifelse(group.text$CCDC24 > median(group.text$CCDC24), "high", "low")

# expre data
c2 <- c('FCGR3A', 'MMP1', 'SERPING1', 'EGFL6', 'SLC2A12', 'CORIN', 'EMX2', 'ACKR4',
        'COL3A1', 'COL5A2', 'FN1', 'FNDC1', 'VGLL3', 'CFH', 'NT5E', 'RASGRF2', 'LRP1',
        'CACNA2D2', 'KIF12', 'CCDC24', 'RAB3A')
heatmap.data <- B.dat[,c2]

all(rownames(group.text) == rownames(heatmap.data))
data.scale <- scale(heatmap.data)
data.scale.new <- ifelse(data.scale > 1, 1, ifelse(data.scale < -1, -1, data.scale))
# data.scale.new <- na.omit(data.scale.new)
dat <- t(data.scale.new)

all(rownames(group.text) == colnames(dat))

# plot Heatmap
mat <- as.matrix(dat)
group.mat <-  as.matrix(subset(group.text, select = -CCDC24))

all(rownames(group.mat) == colnames(dat))

pdf("set3_heatmap.pdf", width = 6, height = 6)
column.ha <- HeatmapAnnotation(df = group.text, 
                               col = list(label = c("type1" = "red", "type2" = "blue"), 
                                          CCDC24 = c("low" = "#ffce7b", "high" = "#f15a22")))

ht <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#7fb80e", "white", "red")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter",
                                          legend_height = unit(6, "cm"),
                                          legend_direction = "vertical"),
              top_annotation = column.ha,
              show_row_dend = TRUE,
              #right_annotation = row_anno,
              row_title_gp = gpar(fontsize = 20),       
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              column_split = group.mat,
              column_title = NULL)
ht
dev.off()


## Different analysis
library(data.table)
library(openxlsx)
library(dplyr)
my.data <- read.csv("my.data.csv")

c2 <- c("label", 'FCGR3A', 'MMP1', 'SERPING1', 'EGFL6', 'SLC2A12', 'CORIN', 'EMX2',
        'ACKR4', 'COL3A1', 'COL5A2', 'FN1', 'FNDC1', 'VGLL3', 'CFH', 'NT5E', 'RASGRF2',
        'LRP1', 'CACNA2D2', 'KIF12', 'CCDC24', 'RAB3A')
study.data <- my.data[, c2]
ncol(study.data)
table(study.data$label)

data <- study.data
col.names <- colnames(data)
re <- c()
normalRes <- c()
for (i in col.names) {
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
res2 <- res1[res1[, "type1.p"] > 0.05, ]
tRes <- res2[res2[, "bar.p"] > 0.7, ]
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

# U_test
data <- remain.data
dif <- c()
col.names <- colnames(data)
j <- levels(factor(data$label))
for (i in 2:ncol(data)) {
  high <- data[which(data$label == "type2"), i]
  low <- data[which(data$label == "type1"), i]
  Utest <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE, correct = FALSE, conf.level = 0.95)
  log2FC <- mean(data[data$label == j[1], i]) - mean(data[data$label == j[2], i])
  dif <- rbind(dif, cbind(col.names[i], Utest$p.value, log2FC))
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

final.genes.t.u <- union(as.character(t.dif$gene), as.character(u.dif$gene))  #  45
length(final.genes.t.u)
final.genes.t.u



## Fisher test of Clinical factors
library(data.table)
library(openxlsx)
my.data <- read.csv("my.data.csv")
clinical.data <- my.data
clinical.data$lymphnodegroup <- as.character(clinical.data$lymphnodegroup)
clinical.data$lymph.node.status <- as.character(clinical.data$lymph.node.status)

clinical.data$age <- ifelse(clinical.data$age > 60, 1, ifelse(clinical.data$age <= 60, 0, "Null"))
sub <- which(is.na(clinical.data$tumor.size))
val <- clinical.data[-sub, ]

med2 <- median(val$tumor.size)
med2

clinical.data$tumor.size <- ifelse(clinical.data$tumor.size > 20, 1, ifelse(clinical.data$tumor.size <= 20, 0, "Null"))
clinical.data[is.na(clinical.data)] <- "Null"
table(clinical.data$tumor.size)

car.df2 <- table(clinical.data$label, clinical.data$pr) 
car.df2
fisher.test(car.df2)


## Clinical factor cox analysis
library(survival)
library(survminer)
library(plyr)
library(grid)
library(MASS)
library(data.table)
library(openxlsx)
my.data <- read.csv("my.data.csv")
my.data$chemo.treated <- as.numeric(as.character(my.data$chemo.treated))
dd <- my.data
UniCox <- function(x) {
  FML <- as.formula(paste0("Surv(os_months, os_event)~", x))
  rex.cox <- coxph(FML, data = dd)
  sum <- summary(rex.cox)
  y <- cbind(sum$coefficients, sum$conf.int)
  return(y)
}
# VarNames <- c("her2")
# VarNames <- c("chemo.treated")
VarNames <- c("her2", "chemo.treated")
UniVar <- lapply(VarNames, UniCox)
UniVar <- ldply(UniVar, data.frame) 
factor <- VarNames
results <- cbind(factor, UniVar)
results <- data.frame(factor = results[, 1], HR = results[, 3], z = results[, 5], pvalue = results[, 6],
                      lower95_CI = results[, 9], upper95_CI = results[, 10])
results1 <- arrange(results, results[, 4])  
dim(results1)
results1


## Univariate cox  of gene 
library(data.table)
library(openxlsx)
my.data <- read.csv("my.data.csv")

# Convert continuous data to discrete data
cox.med.df <- my.data

lab <- c("os_days", "X", "title", "endocrine.treated", "chemo.treated", 
         "label", "os_months", "os_event","lymphnodegroup", "chemo.treated", 
         "age", "subtype", "title", 'tumor.size', "lymph.node.status", "er",
         "pr", "her2", "ki67", 'label', 'geo_accession')
for (i in colnames(cox.med.df)) {
  if (i %in%lab) {
    print(i)
  }
  else {
    cox.med.df[, i] <- ifelse(cox.med.df[, i] > median(cox.med.df[, i]), 1, 0)
  }
}

names <- c("os_months", "os_event", "chemo.treated", "SERPING1", "FN1", "CFH", "ACKR4", "COL5A2", "COL3A1",
           "NT5E", "FNDC1", "RASGRF2", "SLC2A12", "CORIN", "LRP1", "KIF12", "VGLL3", "EGFL6",
           "CCDC24", "CACNA2D2", "RAB3A", "MMP1", "EMX2", "FCGR3A")
sub.cox.med.df <- cox.med.df[, names]
# head(sub.cox.med.df, n = 2)

library(survival)
library(survminer)
library(plyr)

dd <- sub.cox.med.df
UniCox <- function(x) {
  FML <- as.formula(paste0("Surv(os_months, os_event)~", x))
  rex.cox <- coxph(FML, data = dd)
  sum <- summary(rex.cox)
  y <- cbind(sum$coefficients, sum$conf.int)
  return(y)
}

VarNames <- colnames(dd)[c(3:23)]  
UniVar <- lapply(VarNames, UniCox)
UniVar <- ldply(UniVar, data.frame) 
factor <- VarNames
results <- cbind(factor, UniVar)
results <- data.frame(factor = results[, 1], HR = results[, 3], z = results[, 5], pvalue = results[, 6],
                      lower95_CI = results[, 9], upper95_CI = results[, 10])

results1 <- arrange(results, results[, 4])  
dim(results1)
results1
# write.table(results, "Genes_single_cox_result.csv", row.names = TRUE, col.names = TRUE, sep = ",")


## Multivariate cox of gene
# use four key genes to conxstruct  multivariable cox model
cox.gene <- c("chemo.treated", "CCDC24", "NT5E", "COL5A2", "CORIN", "FN1")
sub.cox.med.df$chemo.treated <- as.numeric(as.character(sub.cox.med.df$chemo.treated))
fml <- as.formula(paste0("Surv(os_months, os_event)~", paste0(cox.gene, collapse = "+")))
mode <- coxph(fml, data = sub.cox.med.df)
Multisum <- summary(mode)
Multisum


### Forestplot
# plot1
library(forestplot)
rs.forest <- read.csv("forest_data1.csv", header = FALSE)

pdf("set3.1.forestplot.pdf", width = 6, height = 6)
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
pdf("set3.2.forestplot.pdf", width = 6, height = 6)
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


## analysis of KI67 
library(data.table)
library(openxlsx)
ki67.data <- read.csv("ki67_data.csv")
# ki67.data$index <- as.numeric(ki67.data$index)
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
library(ComplexHeatmap)
library(circlize)
library(grid)
data <- read.csv('cox_mode_data.csv')

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
# data.scale.new - a.omit(data.scale.new)
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
pdf("set3_3-1.pdf", width = 9, height = 3)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom",  annotation_legend_side = "bottom")
dev.off()


library(ggplot2)
p2.data <- data[,1:6]
p2.data$os_event <- factor(p2.data$os_event)

pdf("set3_3-2.pdf", width = 9, height = 3, useDingbats = FALSE)
p2 <- ggplot(p2.data, aes(x = index, y = os_months, color = os_event)) + geom_point(shape = 16)
p2.1 <- p2 + scale_x_discrete("") + theme(axis.title.x = element_blank(),
                                          axis.title.y = element_blank(),
                                          axis.text.x = element_blank(),
                                          axis.ticks.x = element_blank(),
                                          legend.position = c(0.05, 0.12),
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
pdf("set3_3-3.pdf", width = 9, height = 3)
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



