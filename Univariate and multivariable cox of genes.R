
## Univariate cox of gene
library(survival)
library(survminer)

dd <- cox.med.df

UniCox <- function(x) {
  FML <- as.formula(paste0("Surv(OS_months, OS_status)~", x))
  rex.cox <- coxph(FML, data = dd)
  sum <- summary(rex.cox)
  y <- cbind(sum$coefficients, sum$conf.int)
  return(y)
}

VarNames <- colnames(dd)[c(5:25)] 
UniVar <- lapply(VarNames, UniCox) 
UniVar <- ldply(UniVar, data.frame) 
factor <- VarNames
results <- cbind(factor, UniVar)
results <- data.frame(factor = results[, 1], HR = results[, 3], z = results[, 5], pvalue = results[, 6],
                      lower95_CI = results[, 9], upper95_CI = results[, 10])

results1 <- arrange(results, results[, 4])  
dim(results1)
results1
write.table(results, "Genes_single_cox_result.csv",row.names=TRUE,col.names=TRUE,sep=",")


# use four key genes to conxstruct  multivariable cox model
cox.gene <- c("CCDC24", "NT5E", "COL5A2", "CORIN", "FN1")

fml <- as.formula(paste0("Surv(OS_months, OS_status)~", paste0(cox.gene, collapse = "+")))
mode <- coxph(fml, data = cox.med.df)
Multisum <- summary(mode)
Multisum

# mode = -1.3305 *CCDC24
cox.expre.df$mode <- (-1.3305 ) * cox.expre.df[, "CCDC24"]
write.table(cox.expre.df, "cox_mode_data2.csv", row.names = TRUE, col.names = TRUE, sep = ",")


# Forestplot
# plot1
library(forestplot)
rs.forest <- read.csv("forest_plot.csv", header = FALSE)
head(rs.forest)

pdf("tcga.1.forestplot.pdf", width = 6, height = 6)
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
# dev.off()

# plot2
library(forestplot)
rs.forest <- read.csv("forest_plot_2.csv", header = FALSE)

pdf("tcga.2.forestplot.pdf", width = 6, height = 6)
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



