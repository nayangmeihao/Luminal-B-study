library(openxlsx)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(data.table)
library(survMisc)

label.data <- read.xlsx("cluster_2_res.xlsx", sheet = 1)
table(label.data$label)
label.data$label[label.data$group == 1] <- "Cluster A" 
label.data$label[label.data$group == 2] <- "Cluster B" 

# rownames (dat) are sample,colnames are Patients,OS_months,OS_status
dat <- read.xlsx("sur_data2.xlsx", sheet = 1)
dat$OS_status <- ifelse(dat$OS_status == 'Alive', 0, 1)
sur.data <- dat[, c("Patients", "OS_months", "OS_status")]
sur.data$OS_months <- as.numeric(sur.data$OS_months)
df <- merge(label.data, sur.data, by.x = "Patients", sort = FALSE)

data2 <- df[, 1:5]
# data2 <- data2[- which(data2$OS_month >= 120), ]
y <- Surv(data2$OS_months, data2$OS_status)
kmfit <- survfit(y~ data2$label, data = data2)

#plot survival curve
# pdf("cluster2_KM_plot.pdf", width = 6, height = 6)
gg <- ggsurvplot(kmfit, palette = c("red", "blue"),
                 legend.title = "label",
                 legend = c(0.85, 0.9),
                 #pval =TRUE,
                 pval = "P value: 0.017",
                 pval.coord = c(0.13, 0.23),
                 conf.int = TRUE,
                 ylim = c(0.2, 1),
                 xlab = "Time in months",
                 risk.table = TRUE,
                 risk.table.height = 0.2,
                 ggtheme = theme(legend.title = element_blank(), 
                                 panel.background = element_rect(fill = "white", color = "black"),           
                                 legend.key = element_rect(fill = "transparent"),
                                 axis.title = element_text(color = "black", vjust = 2),
                                 axis.ticks = element_line(color = "black")) )
gg
# dev.off()

# Survival curves difference test
# logrank test
data2 <- data2[-which(data2$OS_month >= 120), ]
surv <- as.formula(paste("Surv(df$OS_months, df$OS_status)~", "label"))
x <- survdiff(surv, data = df)
# pvalue
pValue <- 1 - pchisq(x$chisq, df = 1) 
pValue

# supremum (Renyi) test; two-sided; two covariate groups
g1 <- ten(Surv(data2$OS_months, data2$OS_status) ~ group, data = data2)
comp(g1)

