
library(data.table)
library(openxlsx)

data <- read.csv("cox_mode_data_2.csv")
data$index <- seq(1:187)
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
mat <- dat
column.ha <- HeatmapAnnotation(df = group.df,
                               col = list(risk.label = c( "Low risk" = "#ffce7b", "High risk" =  "#f15a22"),
                                          label = c( "type1" = "red", "type2" =  "blue")))
p1 <- Heatmap(mat,
              col = colorRamp2(c(-1, 0, 1), c("#7fb80e", "white", "red")),
              heatmap_legend_param = list(title = "legend", 
                                          title_position = "topcenter",
                                          legend_height = unit(6,"cm"), 
                                          legend_direction = "horizontal"),
              bottom_annotation = column.ha,
              show_row_dend = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              #column_split = group.mat,
              column_title = NULL)

draw(p1)
pdf("tcga_3-1.pdf", width = 9, height = 3)
draw(p1, merge_legend = TRUE, heatmap_legend_side = "bottom",  annotation_legend_side = "bottom")
dev.off()


library(ggplot2)
p2.data <- data[, c("index", "OS_months", "OS_status")]
p2.data$OS_status <- factor(p2.data$OS_status)
# head(p2.data)

pdf("tcga_3-2.pdf", width = 9, height = 3, useDingbats = FALSE)
p2 <- ggplot(p2.data, aes(x = index, y = OS_months, color = OS_status)) + geom_point(shape = 16)
p2.1 <- p2 + scale_x_discrete("") + theme(axis.title.x = element_blank(),
                                          axis.title.y = element_blank(),
                                          axis.text.x = element_blank(),
                                          axis.ticks.x = element_blank(),
                                          legend.position = c(0.93, 0.9278),
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
pdf("tcga_3-3.pdf", width = 9, height = 3)
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