
## Sankeyplot
library(data.table)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggalluvial)

genes.data <- read.xlsx("B_DEseq2_log2_data.xlsx", sheet = 1)  # 187*538
genes.data <- as.data.frame(genes.data)

thr.df <- read.xlsx("ER_PR_her2_data.xlsx", sheet = 1)  #538*275
thr.df <- as.data.frame(thr.df)
her2.data <- thr.df[, c("sample_id", "HER2")]
merge.data <- left_join(her2.data, genes.data, by = c("sample_id" = "sample"))

# type1
data1 <- merge.data[which(merge.data$label == 1), ]

HER2.p.dat <- data1[which(data1$HER2 == "HER2+"), ]  # 10*540
HER2.n.dat <- data1[which(data1$HER2 == "HER2-"), ]  # 72*540

# Calculate the median of each gene of HER2.p.dat and HER2.n.dat
# HER2.p.dat
result <- c()
for (i in colnames(HER2.p.dat)) {
  re <- c()
  if (i %in% c("sample_id", "HER2", "label")) {
    print(i)
  }
  else {
    med <- median(HER2.p.dat[, i])
    re <- cbind(i, med)
    result <- rbind(result, re)
  }
}
result.p <- as.data.frame(result)
names(result.p)[names(result.p) == "i"] <- "genes"
result.p <- result.p[order(result.p$med), ]


# HER2.n.dat
result <- c()
for (i in colnames(HER2.n.dat)) {
  re <- c()
  if (i %in% c("sample_id", "HER2", "label")) {
    print(1)
  }
  else{
    med <- median(HER2.n.dat[, i])
    re <- cbind(i, med)
    result <- rbind(result, re)
  }
}
result.n <- as.data.frame(result)
names(result.n)[names(result.n) == "i"] <- "genes"
result.n <- result.n[order(result.n$med), ]


# HER2+\HER2-
# result.p
dat.p <- result.p
dat.p$med <- as.numeric(as.character(dat.p$med))

# Dividing numbers according to quartiles
q <- quantile(dat.p[, "med"])
dat.p[, "med"] <- ifelse(dat.p[, "med"] > q[4], "Q4",
                         ifelse(dat.p[, "med"] > q[3], "Q3",
                                ifelse(dat.p[, "med"] > q[2], "Q2", "Q1")))

# Dividing HER2.p.dat by quartile
her2.p.dat<- HER2.p.dat
result <- c()
for (i in 4:ncol(her2.p.dat)) {
  her2.p.dat[, i] <- ifelse(her2.p.dat[, i] > q[4], "Q4",
                            ifelse(her2.p.dat[, i] > q[3], "Q3",
                                   ifelse(her2.p.dat[, i] > q[2], "Q2", "Q1")))
}

# Count the number of patients in Q1\Q2\Q3\Q4 for each gene
result <- c()
col.name <- colnames(her2.p.dat)
for (i in 4:ncol(her2.p.dat)) {
  re <- c()
  gene <- col_name[i]
  Q <- dat.p[which(dat_p$gene == gene), "med"]
  num <- as.numeric(sum(unlist(her2.p.dat[, 4]) == Q))
  re <- rbind(gene, Q, num)
  re <- as.data.frame(t(re))
  result <- rbind(result, re)
}

result.p <- result
result.p$group <- c(rep("HER2+", 537))
result.p$id <- seq(1, 537)

result.p$Q <- ifelse(result.p$Q == "Q4", "Q4(134)",
                     ifelse(result.p$Q == "Q3", "Q3(134)",
                            ifelse(result.p$Q == "Q2", "Q2(134)", "Q1(135)")))


plot.data <- rbind(result.p, result.n)
write.table(plot.data, "type1_her2_data.csv", row.names = TRUE, col.names = TRUE, sep = ",") 

levels(plot.data$Q) <- rev(levels(plot.data$Q))

pdf("type1_sankeyplot.pdf")
ggplot(plot.data,
       aes(x = group, stratum = Q, 
           alluvium = id, y = num,
           fill = Q, label = Q)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  ggtitle("Gene expression of HER2+ and HER2- in type1")
dev.off()