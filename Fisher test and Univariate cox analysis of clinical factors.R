#Fisher Test of Clinical factors
library(data.table)
library(openxlsx)

clinical.data <- read.xlsx("final.data.xlsx", sheet = 1)
clinical.data <- as.data.frame(clinical.data)

# fisher test
car.df2  <- table(clinical.data$label, clinical.data$CNA) 
fisher.test(car.df2)


library(data.table)
library(openxlsx)
library(dplyr)

# clinical data 
clinical.data <- read.xlsx("final.data.xlsx", sheet = 1)
clinical.data <- as.data.frame(clinical.data)
clinical.dat <- clinical.data[, c("sample", "HER2", "CNA")]

# survival data
sur.df <- read.xlsx("cluster_2_sur_dat.xlsx") # 187*4

# join survival data and clinical data by sample
clinical.cox.dat <- left_join(sur.df, clinical.dat, by = "sample")
# head(clinical.cox.dat)


# cox analysis of clinical factor 
library(survival)
library(survminer)
library(plyr)
library(survivalROC)
library(grid)
library(MASS)
library(leaps)

dd <- clinical.cox.dat
# Define the function, input the gene expression data to be analyzed by cox single factor; return the analysis result
UniCox <- function(x) {
  FML <- as.formula(paste0("Surv(OS_months, OS_status)~", x))
  rex.cox <- coxph(FML, data = dd)
  sum <- summary(rex.cox)
  y <- cbind(sum$coefficients, sum$conf.int)
  return(y)
}

VarNames <- colnames(dd)[c(5:6)] 
UniVar <- lapply(VarNames, UniCox)
UniVar <- ldply(UniVar, data.frame) 
factor <- VarNames
results <- cbind(factor, UniVar)
results <- data.frame(factor = results[, 1], HR = results[, 3], z = results[, 5], pvalue = results[, 6],
                      lower95_CI = results[, 9], upper95_CI = results[, 10])

results1 <- arrange(results,results[, 4])  
dim(results1)
results1



