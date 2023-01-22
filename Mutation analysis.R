rm(list = ls())
library(data.table)
library(openxlsx)
library(SummarizedExperiment)
label.data <- read.csv("data/join_name_data.csv",head = TRUE)

library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

GDCdownload(query)
GDCprepare(query, save = T,save.filename = "TCGA-BRCA_SNP.Rdata")

library(maftools)
load(file = "TCGA-BRCA_SNP.Rdata")
maf.coad <- data
write.table(maf.coad, "raw_mutation_data.csv", row.names = TRUE, col.names = TRUE, sep = ",") 
intersect(label.data$raw_sample_names,maf.coad$Tumor_Sample_Barcode)

mutation.df <- read.xlsx("data/mutation_data/4keyGene_Mutation_data_from_TCGA_v2.xlsx")
cna.df <- read.xlsx("data/mutation_data/Sheet2.xlsx")

labelAB.data <- read.xlsx("cluster_2_sur_dat.xlsx")
length(intersect(labelAB.data$sample,cna.df$sample)) #187
length(intersect(labelAB.data$sample,mutation.df$sample)) #126

same.name <- intersect(labelAB.data$sample,cna.df$sample)
cna.df.2 <- cna.df[which(cna.df$sample %in% same.name),] #187*5

cna.df.3 <- cna.df.2
mutation.df.2 <- mutation.df
same.genes <- c("ERBB2","GATA3","PIK3CA","TP53")
for (i in 1:nrow(mutation.df.2)){
  ge <- mutation.df.2[i, "Gene"]
  sa <- mutation.df.2[i, "sample"]
  mu.val <- mutation.df.2[i, "Variant_Type"]
  mu.val.2 <- ifelse(mu.val%in%c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site-"), "Truncating",
                     ifelse(mu.val== "Missense_Mutation", "Missense", 
                            ifelse(mu.val=="Nonsense_Mutation", "", "Fusion")))
  
  can.val <- cna.df.3[which(cna.df.3$sample == sa), ge]
  if(mu.val==""){
    print("without change!")
    }
  else{
    new.val <- paste(mu.val.2, can.val, sep = ",")
    cna.df.3[which(cna.df.3$sample == sa), ge] <- new.val}
}
write.table(cna.df.3, "data/mutation_data/oncodata_data_raw.csv", row.names = FALSE, col.names = TRUE, sep = ",") 


label.data <- read.csv("mutation_data/Same_sample_data.csv",head = TRUE) #112*2
raw.oncodata <- read.xlsx("oncodata_data_processed.xlsx")
rownames(raw.oncodata) <- raw.oncodata$sample
info <- raw.oncodata[label.data$sample,c("sample","TP53_label")]
oncodata <- raw.oncodata[label.data$sample,c("sample","TP53")]

info <- raw.oncodata[label.data$sample,c("sample","ERBB2_label")]
oncodata <- raw.oncodata[label.data$sample,c("sample","ERBB2")]

info <- raw.oncodata[label.data$sample,c("sample","GATA3_label")]
oncodata <- raw.oncodata[label.data$sample,c("sample","GATA3")]

info <- raw.oncodata[label.data$sample,c("sample","PIK3CA_label")]
oncodata <- raw.oncodata[label.data$sample,c("sample","PIK3CA")]


colnames(info) <- c("sample","label")
info <- info[order(info$label,decreasing = F),]

oncodata.new <- oncodata[info$sample,]
colnames(oncodata.new) <- c("sample","gene_exp")
oncodata.new <- subset(oncodata.new,select=-c(sample))
print <- t(oncodata.new)

library(ComplexHeatmap)
col = c("Amplification" = "#99004b","LowcopyGain"="#FFD9EC", "ShallowDeletion"="#E0E0E0", 
        "Missense" = "#00E3E3", "Truncating"="#000000", "Fusion" = "#AAAAFF")

alter_fun = list( #horiz_margin = unit(0.2,"mm"), vertical_margin = unit(0.5,"mm"),
  background = alter_graphic("rect", horiz_margin = unit(0.5,"mm"),fill = "#eeeeee",width = 0.9), 
  Amplification = alter_graphic("rect", fill = col["Amplification"],width = 1.1),
  LowcopyGain = alter_graphic("rect", fill = col["LowcopyGain"],width = 1.1),
  ShallowDeletion = alter_graphic("rect", fill = col["ShallowDeletion"],width = 1.1),
  Missense = alter_graphic("rect", height = 0.33, fill = col["Missense"],width = 1.1),
  Truncating = alter_graphic("rect", height = 0.4, fill = col["Truncating"],width = 1.1),
  Fusion = alter_graphic("rect", height = 0.5, fill = col["Fusion"],width = 1.1))

heatmap_legend_param = list(title = "Alternations", at = c("Amplification","LowcopyGain","ShallowDeletion","Missense","Truncating","Fusion"), 
                            labels = c("Amplification","Low-copy Gain", "Shallow deletion", "Missense", "Truncating","Fusion"))

top_anno <- HeatmapAnnotation(Frequency = anno_oncoprint_barplot(), 
                              Group = info$label,
                              annotation_name_side = "left",
                              col = list(Group = c("mutated" = "red", "wild" = "blue")),
                              annotation_name_gp = gpar(fontsize=8))
p <-oncoPrint(print, 
             alter_fun = alter_fun, 
             col = col, 
             heatmap_legend_param = heatmap_legend_param,
             # row_order = colnames(oncodata),
             column_split = info$label,
             # remove_empty_columns = TRUE,
             row_names_side = "left", 
             pct_side = "right",
             top_annotation = top_anno,
             alter_fun_is_vectorized = FALSE,
             row_names_gp = gpar(fontsize=7))

pdf(file = "data/mutation_data/TP53_mutation_oncoPrint_plot.pdf", width = 9, height = 3,useDingbats = FALSE)
draw(p)
dev.off()





