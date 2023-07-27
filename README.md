
# Project Title
Integrative Analysis Identifies Two Molecular and Clinical Subsets in Luminal B Breast Cancer

# Project Description
Integrative analysis with multi-omics data identifies two molecular subtypes with distinct enriched signaling pathways in Luminal B breast cancer patients.

# How to use
Dependency
* R version
```
R version: 4.0.2
```

* Installing R dependency 
```
The required R package can be installed in the following ways：

install.packages(c("data.table","openxlsx","dplyr","Matrix"))


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RMassBank")


#load githubinstall at first
library(githubinstall)
install_github('dplyr') #install from github, e.g. dplyr
```

Experimental data

* The training set data was downloaded by：
```
RNA-seq within the TCGA database can be downloaded and ID transformed using the code Download dataset.R, in addition the code can also download DNA Methylation data. Therefore we recommend running the code sentence by sentence in R studio or r software.
Alternatively you can modify the code and run the file with the following command：
Rscript Download dataset.R
```
* Validation set data can be downloaded by.
```
GSE96058, GSE20685 and GSE54275 can be downloaded via the GEO(https://www.ncbi.nlm.nih.gov/geo/) website.
The E-MTAB-6703 can be downloaded via the ArrayExpress(https://www.ebi.ac.uk/arrayexpress/) website.
METABRIC/Nature 2012 be downloaded from the cBioPortal (https://www.cbioportal.org/).
```
* Run
```
As the data obtained needs to be pre-processed prior to the experiment, if you wish to run the code within this project please read the code 
comments in detail or please email Huina Wang AT nayangmeihao@outlook.com to request the pre-processed data.

# Clustering of the LuminalB expression profiles corresponds to Figures S1C-D in the paper.

    Rscript Clustering.R
 
# A batch effect test was performed on the data used in the paper, corresponding to Figure S1A in the paper.

    Rscript Batch test.R
    
# The following two codes are for the two subtypes of LuminalB and for the differential analysis of tumor and normal, respectively.

    Rscript Different analysis between ClusterA and ClusterB.R
    Different analysis between tumor and normal.R

```