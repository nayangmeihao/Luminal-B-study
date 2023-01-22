library(data.table)
library(openxlsx)

df <- read.xlsx("D:/study_data/plot_data/4248_genes_log2_data.xlsx", sheet = 1)
df <- as.data.frame(df)  # rownames(df) are samples,colnames are genes
rownames(df) <- df$sample
df <- subset(df, select = -sample)

f <- function(x) median(x)
median.df <- as.data.frame(apply(df, 2, f))
names(median.df)[names(median.df) == "apply(df, 2, f)"] <- "median"
median.df$gene <- rownames(median.df)

# Count the number of patients with expression of each gene
f <- function(x) sum(x != 0)
no.zero.df <- as.data.frame(apply(df, 2, f))
names(no.zero.df)[names(no.zero.df) == "apply(df, 2, f)"] <- "num_zero"
no.zero.df$gene <- rownames(no.zero.df)

# Combine the two dataframes obtained above. rownames(mer.data) are genes
mer.data <- merge(median.df, no.zero.df, by = "gene", sort = FALSE) # 4248*3
dim(mer.data)


# Divide the sequence composed of the median value corresponding to all genes 
# into 10 equal parts from small to large, and then determine the interval 
# to which each median value belongs.
max.d <- max(mer.data$median)  
min.d <- min(mer.data$median) 
s <- (max.d - min.d) / 10

re <- c()
for (i in mer.data$median) {
  if (i >= 0 & i < s){re <- rbind(re, 1) }
  if (i >= s & i < 2*s){re <- rbind(re, 2) }
  if (i >= 2*s & i < 3*s){re <- rbind(re, 3) }
  if (i >= 3*s & i < 4*s){re <- rbind(re, 4) }
  if (i >= 4*s & i < 5*s){re <- rbind(re, 5) }
  if (i >= 5*s & i < 6*s){re <- rbind(re, 6) }
  if (i >= 6*s & i < 7*s){re <- rbind(re, 7) }
  if (i >= 7*s & i < 8*s){re <- rbind(re, 8) }
  if (i >= 8*s & i < 9*s){re <- rbind(re, 9) }
  if (i >= 9*s ){re <- rbind(re, 10) }
}

length(re)
mer.data$Decile <- re
# mer.data$label <- as.numeric(mer.data$label)

# Descending
mer.data <- mer.data[order(-mer.data$num_zero), ]
# write.csv(mer.data, "mer.data.csv")

# Mark Key.genes
key.genes1 <- c("APOB", "CCND1", "CDH11", "COL1A1", "COL2A1", "COL5A1", "COX6C", "DDR2", 
                "ETV1", "FCGR3A", "FGFR3", "HIST1H3B", "INHBA", "KLK2", "KRT15", "KRT31",
                "KRT75", "KRT81", "NRG1", "PDGFRA", "RET", "ROS1", "RSPO3", "S100A7", 
                "S100A8", "S100A9", "SOX11", "TNFRSF18", "TNFSF4", "ZNF703")

mer.data$key.gene  <-  ifelse(mer.data$gene %in% c(key.genes1), mer.data$gene, "")   

library(ggplot2)
p.data <- mer.data
p.data$Decile <- factor(p.data$Decile)
# head(p.data)

# pdf("key_gene.pdf", width = 15, height = 10, useDingbats = FALSE)
p <- ggplot(p.data, aes(x = -num_zero, y = median)) + 
  geom_point(aes(colour = Decile), shape = 16) + scale_y_continuous("log2 Median intensity") + 
  theme_classic() + theme(legend.position = c(0.8, 0.75), 
                          panel.background = element_rect(fill = "white")) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(name = "Number of observations", 
                     breaks = c(-187, -175, -150, -125, -100), 
                     labels = c(187, 175, 150, 125, 100))
# p

library(ggrepel)
# Label gene, add layer
gg <- p + geom_text_repel(data = p.data, aes(x = -num_zero, y = median, label = key.gene), 
                          size = 3, box.padding = unit(0.5, "lines"),
                          point.padding = unit(0.8, "lines"), 
                          segment.color = "black", 
                          show.legend = FALSE,
                          check_overlap = TRUE) 
gg
# dev.off()