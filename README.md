# EDTox

www.edtox.fi

[![DOI](https://zenodo.org/badge/267520245.svg)](https://zenodo.org/badge/latestdoi/267520245)

https://academic.oup.com/bioinformatics/article/38/7/2066/6520808?login=true

## Description of EDTox
The EDTox is an R shiny application to predict endocrine disruption potential of a compound  using a machine learning approach based on toxicogenomics data. 
In this approach, a random walk with restart (RWR) is used on a gene co-expression network stating from the molecular initiating events (MIEs) to expand the list of genes perturbed by a certain chemical. The resulting gene set is then used for fast geneset enrichment analysis (FGSEA). The pathway activation scores obtained from FGSEA are used along with lists of MIEs from known EDCs and negative controls (non-EDCs) to train an elastic net GLM classifier to predict EDC probability of unknown compound from a set of primary perturbing known genes. Validation of the model is performed by k-fold cross-validation. The pipeline also allows an option to use pareto solution to select the optimal proportion of edges from the original network that will be used as input for RWR and the number of genes with highest probability after RWR to be taken into FGSEA.
This approach was tested to predict the EDC probability of about 12,000 compounds in the Comparative Toxicogenomics Database (CTD) using a large toxicogenomics data space composed of DrugMatrix, open TG-gates, LINCS and PPI. For this purpose, 21 gene co-expression networks were compiled from transcriptome data. The MIEs of compounds obtained from CTD were then used to obatain the EDC-class probability based on these 21 networks.
While users can train a new model based on their own network and training set composed of MIEs from known EDCs and negative controls, the application also provides a way to calculate EDC-class probability based on 21 pre-compiled classifiers mentioned above.The application can also be used to view the EDC-class probabilities of the 12,000 compounds in CTD predicted using  the EDTox pipeline and compare it with the ToxPi scores from ToxCast. 

## Install Dependencies
#### Install CRAN dependencies

```r
cran_pkgs <- c("doParallel", "foreach", "igraph", "plyr", "shiny","shinyjs", "shinyBS", "shinydashboard", "DT", "ggplot2","RColorBrewer", "dnet","tidyr","ggpubr","rPref","caret","glmnet","shape")

cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]

if(length(cran_pkgs.inst)>0){ 
  print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."));  
    install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE);print("Installed!!!") 
  }
}
```


#### Install Bioconductor dependencies for R versions > 4
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
bioc_pkgs <- c( "supraHex","hexbin",  "fgsea","Rgraphviz")
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))];
if(length(bioc_pkgs.inst)>0){
  print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"));  
  for(pkg in bioc_pkgs.inst){
    print(paste0("Installing Package:'", pkg, "'..."));  
    BiocManager::install(pkg)
    print("Installed!!!")
  }
}
```
## Launch from GitHub

```r
runGitHub("EDC_shiny", "vittoriofortino84", subdir="EDC_shiny")
## Using the archived file
runUrl("https://github.com/vittoriofortino84/EDC_shiny/archive/master.tar.gz", subdir="EDC-shiny")
runUrl("https://github.com/vittoriofortino84/EDC_shiny/archive/master.zip", subdir="EDC-shiny")
```
## Launch locally
```
git clone https://github.com/vittoriofortino84/EDC_shiny EDC_shiny_clone
```

## Launch from R

```r
## Run by using runApp()
setwd("~/EDC-shiny")
shiny::runApp(appDir = getwd(), port = getOption("shiny.port"),
              launch.browser = getOption("shiny.launch.browser", interactive()),
              host = getOption("shiny.host", "127.0.0.1"), workerId = "",
              quiet = FALSE, display.mode = c("auto", "normal", "showcase"), 
              test.mode = getOption("shiny.testmode", FALSE))
```

## EDTox Usage

### Summary
<img src="/example/figure/Summary_1.jpg" width="600">
  The summary tab presents an overview of the 21 toxicogenomics classifiers, the pathways used in this application as well as the distribution of the predicted toxicity scores across all 21 classifiers for 12k compounds in Comparative Toxicogenomics Database. It also includes a box plot of the F1 accuracy scores of the pre-compiled different classifiers in this application.The plot also allows for importing F1 scores of additional models to compare with the existing classifiers. 

### Toxicogenomics Pipeline
<img src="/example/figure/Toxicogenomics_pipeline_1.jpg" width="600">
  This is the primary tab of the application and used for generating classification models from toxicogenomics based co-expression network. The tab is divided into five sections/modules.The first or topmost section is divided into two halves. The left part displays the system time, allows the user to set the number of CPUs to be used during execution of the pipeline, the classifier name and the current execution status.  The right part of this section can be used to upload the gene-gene co-expression network to be used for building the classification model, and the lists of EDCs and negative controls along with their MIEs to be used as training set. All these three files should be uploaded in RDS format. The network file must be a three column dataframe with the first two columns containing the interacting nodes and the third column containing the edge parameter like edge weight or  or topological overlap generated from WGCNA or wTO package. The training set MIEs should be uploaded as R lists with names containing EDCs or negative control to which the MIEs belong. 

<img src="/example/figure/Toxicogenomics_pipeline_2.jpg" width="600">
  The second module of this tab allows for finding the optimal parameters to be used for RWR and FGSEA. This module calculates a Pareto solution for three parameters: (1) the proportion of edges from the network to be used for RWR (as %) (defined as k), (2) the number to nodes with highest probability from RWR to be used as input for FGSEA (defined as N), and (3) the silhouette scores to measure the EDC-cluster cohesion based on the set of N genes selected with the RWR. The final combination of K and N values giving the highest silhouette score is selected by to be used for downstream analysis.

<img src="/example/figure/Toxicogenomics_pipeline_3.jpg" width="600">
  The third module executes the training and validation of GLM based classifier using RWR-FGSEA-GLM. This section contains four inputs: the proportion of edges to be used for RWR (as %), the number of genes for FGSEA, the k-fold cross validations to be performed and the number of cross validation repeats. The first two fields are automatically updated based on Pareto solution but users also have the option to manipulate them if needed. Clicking the "Run RWR-FGSEA-GLM" starts a series of tasks including RWR simulation of the optimized network utilizing the input training set MIEs as seed genes, followed by FGSEA on the selected number of genes with highest probability from RWR and finally using the pathway activation scores from FGSEA into elastic-net generalized linear approach to train models for classification of EDCs and negative controls. A k-fold cross-validation step is also performed to estimate the accuracy of the generated machine learning model. Upon completion of the tasks in this module, it displays a histogram of the F1-scores for each fold validation steps. Note: Cross validation repeat > 1 is not recommended on common desktop or laptop computers.

<img src="/example/figure/Toxicogenomics_pipeline_4.jpg" width="600">
  The parameters for the classifier model generated in the previous step can be saved for future uses using the export module of this tab. This section also allows for exporting the pathway activation scores and F1 score of the model as RDS file. The F1 scores of the generated model can be compared with the F1-scores of the pre-compiled models in the Summary tab of the application as mentioned earlier.  


<img src="/example/figure/Toxicogenomics_pipeline_5.jpg" width="600">
  The final module of the toxicogenomics pipeline tab can be used to predict the EDC class probabilities for a given set of compounds based on their MIEs. It takes as input a RDS file containing lists of compounds and its MIEs, and model parameters (like the one generated in previous sections) based on which prediction is to be made. The output is a table listing the input compounds with EDC class probabilities and class labels. The table can be exported using the "Predicted scores" option in the export module of the tab. Note: If the MIEs of a compound is not included in the network used to build the classifier model, it may not be included in the resultant table.

### Pathway activation scores
<img src="/example/figure/Pathway_activation_scores_1.jpg" width="600">
  The putative pathways as the possible mode of action for the EDCs are displayed as a bubble plot in this tab. In the plot, the bubble size represents the GLM coefficients while the color represent the pathway activation (NES) scores. The threshold values for both these parameters to be considered for plotting the bubble plot can be altered using the slider bars present on the left panel of the tab. Apart from this, the pathway sources and the classifiers to be visualized in the plot can also be selected. The model generated in the toxicogenomics pipeline is automatically updated in the classifiers list. It is also possible to add new classifiers for visualization using the bubble plot. The new classifier should be entered as an RDS file containing a four columns: name of the data layer, average of pathway activation scores across all EDCs for that pathway, GLM coefficients for the pathway and name of the pathway. Note: The order of the column should be maintained as mentioned.

### EDC-class probability scores
<img src="/example/figure/EDC_class_probability_1.jpg" width="600">
  The underlying EDTox pipeline was used to evaluate the EDC class probabilities for about 12,000 compounds present in CTD. This tab can be used to visualize the class probability  scores of these compounds across all classifiers as histogram and also find the average and harmonic sum values across the classifiers. It also allows for a visual comparison of EDC class probabilities of upto five compounds at once. For compounds not among the 12K compounds in CTD, their MIEs can be used to calculate the class probabilities based on the 21 precompiled classifiers.

### Comparison with ToxPi Scores
<img src="/example/figure/ToxPi_1.jpg" width="600">
  This tab provides for an interactive scatter-plot that can be used to compare the EDC-class probability scores computed using the EDTox pipeline with the ToxPi scores. Each dot in the plot represents the 12K compounds from CTD plotted for average EDC-class probability along x-asis and ToxPi score along y-axis. The left panel allows the user to select the classifiers which will be used for calculation of the average EDC-class probabilities. The compounds or dots in the scatter plot can be highlighted based on a minimum ToxPi score and EDC-class probability by using the slider bars on the left panel. Users can also click on individual dots to see details like compound name, average EDC-class probability, and ToxPi score on a table below the plot. The table also provides a link to CompTox database which contains additional information regarding each compound.
