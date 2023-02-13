
# Project Description
Exploring the relationship between the number of cattle cases and climate variables

# How to use
Dependency
* Installing R dependency 
```
R version: 4.2.0

install.packages(c("data.table","openxlsx","dplyr","Matrix"))
```
* Installing python dependency
```
python version: 3.8.13

pip freeze >requirements.txt

```
* Run
```
Firstly, the distribution of the data is viewed and the methods used are stored in main.py. The correlation between climate variables and counts is explored using Spearman and is stored in main.py. 
Filter_data.R is used to visualise, filter and discretise the data obtained. The Apriori.py file is used to explore the association rules between the strong correlation variables and count.

# Get data description and correlation between count and climate variables

    python main.py
 
# Filter out strongly correlated climate variables and discretize

    Rscript Filter_data.R
    
# Perform correlation rule analysis on the discretized data using Apriori algorithm

    python Apriori.py

```

the result of main.py