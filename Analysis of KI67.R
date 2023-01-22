
library(data.table)
library(openxlsx)
library(dplyr)

ki67.data <- read.xlsx("ki67_data.xlsx")
cox.med.df <- read.csv("med.data.csv")
label.data <- cox.med.df[, c("sample", "label")]

ki67.dat <- left_join(label.data, ki67.data, by = "sample")
ki67.dat$index <- rownames(ki67.dat)
ki67.dat$index <- as.numeric(ki67.dat$index)
dim(ki67.dat)
head(ki67.dat, n = 3)

# Normality test
data <- ki67.dat[, c("MKI67", "label")]
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

# 正态性检验（符合）
data <- ki67.dat[, c("MKI67", "label")]
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

