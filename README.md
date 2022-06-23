[![<praneet1988>](https://circleci.com/gh/praneet1988/ELeFHAnt.svg?style=shield)](https://github.com/praneet1988/ELeFHAnt)

# ELeFHAnt
Ensemble Learning for Harmonization and Annotation of Single Cells (ELeFHAnt) provides an easy to use R package for users to annotate clusters of single cells, harmonize labels across single cell datasets to generate a unified atlas and infer relationship among celltypes between two datasets. It provides users with the flexibility of choosing a single machine learning based classifier or letting ELeFHAnt automatically use the power of  randomForest and SVM (Support Vector Machines) to make predictions. It has three functions 1) CelltypeAnnotation 2) LabelHarmonization 3) DeduceRelationship 4) Benchmark ELeFHAnt.

## Version 1.1.3 is now available
* Improved speed (5X faster classification)
* Validate Predictions using GSEA msigdb [C8] (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C8) and CellMarkers Experimentally verified markers (http://bio-bigdata.hrbmu.edu.cn/CellMarker/)
* Users can select species (human or mouse) and tissues of interest for predicted celltype Marker validation. [Dotplots and Feature plots are generated for both reference and query]
* Users now can compare ELeFHAnt predictions against Seurat's Label Transfer and scPred by using BenchmarkELeFHAnt function
* A complete tutorial using PBMC datasets. Reference [https://www.nature.com/articles/ncomms14049] and query [https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz]

## Where to find previous versions
Users can access ELeFHAnt previous releases from Releases section of GitHub [https://github.com/praneet1988/ELeFHAnt/releases]

## Installation
```
library(devtools)
devtools::install_github('praneet1988/ELeFHAnt')
library(ELeFHAnt)
```
```
If you encounter any problems, try manually installing from the .zip or .tar.gz file with "R CMD INSTALL". 
```

## Tutorial using PBMC datasets
https://rpubs.com/praneet1988/917589

## Developers
```
Praneet Chaturvedi (MS Bioinformatics) : Lead Analyst Bioinformatics, Cincinnati Children's Hospital and Medical Center, Cincinnati, OH - USA
GitHub Username: praneet1988

Konrad Thorner (MS Bioinformatics) : Analyst Bioinformatics, Cincinnati Children's Hospital and Medical Center, OH - USA
GitHub Username: kthorner
```

## ELeFHAnt Model
![Graph](ELeFHAnt.png)

## Random Forest
Random Forests is a powerful tool used extensively across a multitude of fields. It is based on generating a large number of decision trees, each constructed using a different subset of your training set. These subsets are usually selected by sampling at random and with replacement from the original data set. The decision trees are then used to identify a classification consensus by selecting the most common output.

## SVM (Support Vector Machines)
SVM is a well-established supervised machine learning algorithm. It finds a hyperplane that separates data points into classes, where training occurs through “maximizing the margin”. Training time tends to be longer as data increases. SVM was ranked highly for cell annotation when benchmarked against other tools.

## Ensemble Learning
ELeFHAnt provides users the ability to use ensemble learning for classifying celltypes for un-annotated datasets. In this mode ELeFHAnt uses the classification accuracy of both Random forest and SVM. It does so by assigning weights (accuracy while learning) to the predictions from each classifier. Weighted confusion matrices from each classifier are normalized based on highest number of cells shared among celltypes and assigned clusters, which are then added together to make the final ensemble predictions.

## Celltype Annotation Function
Celltype annotation is a function used to annotate celltypes in single cell datasets. It requires a reference dataset (a processed Seurat Object with Celltypes column in metadata) and a query dataset (a processed seurat object with seurat_clusters column in metadata). One can choose from randomForest, SVM or Ensemble classifiction methods to learn celltypes from the reference dataset and then predict celltypes for query dataset.

## Label Harmonization Function
Label Harmonization is a function used to harmonize cell labels (celltypes) across single cell datasets. It requires a list of processed Seurat Objects with a Celltypes column in metadata or an integrated seurat object (seurat object with Celltypes and seurat_clusters columns in the metadata). One can choose from randomForest, SVM or Ensemble classifiction methods.

## Deduce Relationship Function
Deduce Relationship is a function used to infer the similarity between celltypes across single cell datasets. The output is a heatmap that shows which celltype in one reference best corresponds to a celltype in another reference. It requires two reference datasets (both processed Seurat Objects with Celltypes columns in the metadata). One can choose from randomForest, SVM or Ensemble classifiction methods.
  
## BenchmarkELeFHAnt Function
Benchmark ELeFHAnt is a function to compare ELeFHAnt celltype predictions against scPred and Seurat's Label Transfer.

## Celltype Annotation in detail
Celltype annotation is a function to annotate celltypes in single cell datasets.

### Requirements
It requires a reference dataset (a processed Seurat Object with Celltypes column in the metadata) and a query dataset (a processed Seurat object with seurat_clusters column in the metadata). One can choose from randomForest, SVM or Ensemble classifiction methods to learn celltypes from the reference dataset and then predict celltypes for the query dataset.

### How to use the function?
#### Load Library ELeFHAnt
library(ELeFHAnt)
#### Assing parameters in the function
out = CelltypeAnnotation(reference = reference.object, query = mydata.object, downsample = TRUE, downsample_to = 500, classification.method = "Ensemble", crossvalidationSVM = 5, validatePredictions = TRUE, selectvarfeatures = 2000, ntree = 500, classification.approach = "ClassifyCells")

## Label Harmonization in detail
Label Harmonization is a function to harmonize cell labels (celltypes) across single cell datasets.

### Requirements
It requires a list of processed Seurat objects with a Celltypes column in the metadata or an integrated Seurat object (integrated Seurat object with Celltypes and seurat_clusters columns in the metadata). One can choose from randomForest, SVM or Ensemble classifiction methods to harmonize celltypes. Please note: DefaultAssay of each object should be set to "RNA".

### How to use the function?
### Load Library ELeFHAnt
library(ELeFHAnt)
### Assing parameters in the function
out = LabelHarmonization(seurat.objects = c(seuratbject1, seuratbject2, seuratbject3, ..), perform_integration = TRUE, downsample = TRUE, downsample_to = 500, classification.method = "Ensemble", crossvalidationSVM = 5, validatePredictions = TRUE, integrated.atlas = NULL, npcs = 30, resolution = 0.8, selectanchorfeatures = 2000, ntree = 500)
  
Or Provide an integrated seurat object
  
out = LabelHarmonization(perform_integration = FALSE, downsample = TRUE, downsample_to = 500, classification.method = "Ensemble", crossvalidationSVM = 5, validatePredictions = TRUE, integrated.atlas = atlas, npcs = 30, resolution = 0.8, selectanchorfeatures = 2000, ntree = 500)

## Deduce Relationship function in detail
Deduce Relationship is a function used primarily to infer the similarity between celltypes across single cell datasets. As the name suggests, any kind of relationship between cell metadata (seurat_clusters, Celltypes, idents etc.) could also be determined.
  
### Requirements
It requires two reference datasets (processed Seurat Object with Celltypes column in the metadata). One can choose from randomForest, SVM or Ensemble classifiction methods to learn celltypes from the reference dataset and find the best corresponding celltypes in the other reference dataset. 

### How to use the function?
### Load Library ELeFHAnt
library(ELeFHAnt)
### Assing parameters in the function
out = DeduceRelationship(reference1 = NULL, reference2 = NULL, downsample = TRUE, downsample_to = 500, classification.method = "Ensemble", crossvalidationSVM = 5, selectvarfeatures = 2000)
  
## BenchmarkELeFHAnt function in detail
Benchmark ELeFHAnt is a function to compare ELeFHAnt celltype predictions against scPred and Seurat's Label Transfer.

### Requirements
It requires a reference dataset (a processed Seurat Object with Celltypes column in the metadata) and a query dataset (a processed Seurat object with seurat_clusters column in the metadata).

### How to use the function?
### Load Library ELeFHAnt
library(ELeFHAnt)
### Assing parameters in the function
out = BenchmarkELeFHAnt(reference = reference, query = query, downsample = TRUE, downsample_to = 500)

## ELeFHAnt Reference datasets as plugins
Download pre-processed reference datasets for Celltype Annotation, Label Harmonization or DeduceRelationship or Benchmark ELeFHAnt here: https://www.dropbox.com/sh/6hd2skriqqlokwp/AAAVol-_qPlCdA4DpERWjkeJa?dl=0

## Gut Development Datasets used to showcase ELeFHAnt functionalities
To demonstrate all the functions of ELeFHAnt we utilize three datasets of early gut development as shown below, as well as an integrated dataset of all three. "GSE158702" refers to a subset of terminal ileum (TI) data from an atlas for human fetal intestinal development called "STAR-FINDer" (https://www.sciencedirect.com/science/article/pii/S009286742031686X). "E-MTAB-8901" refers to a subset of duojejunum cell data from the Gut Cell Atlas, which also examines intestinal development from 6-10 weeks post-conception (https://www.sciencedirect.com/science/article/pii/S1534580720308868). Lastly, "E-MTAB-10187" refers a subset of fetal intestinal data from a multi-endodermal organ atlas (https://www.sciencedirect.com/science/article/pii/S0092867421005316). 

![Graph](Examples/gut_datasets.png)
  
## Celltype Annotation Example

To demostrate Celltype Annotation using ELeFHAnt we used E-MTAB-8901 (~21k cells) as the reference and E-MTAB-10187 (~77k cells) as the query. We used three scenarios 1) Reference was downsampled to 300 cells per celltype and query was downsampled to 200 cells per seurat_cluster. classiification.approach was set to "ClassifyCells" [Downsampled] 2) using all cells in reference and query. classiification.approach was set to "ClassifyCells" [All cells] and 3) using all cells in reference and query but setting using classification.approach = "ClassifyCells_usingApproximation" and downsample = 300 [Approximation]

### ELeFHAnt Celltype Anntation
![Graph](Examples/CelltypeAnnotation_GutAtlas_Figure3.png)


## Label Harmonization Example
To demonstrate LabelHarmonization we used three datasets: 1) E-MTAB-8901 (~21k cells) 2) E-MTAB-10187 (~77k cells) 3) GSE158702 (~20k cells). We integrated three datasets using Seurat's CCA based integration and then ran Label Harmonization (downsample = TRUE, downsample_to = 500) on the integrated object (~112k cells). Left panel shows all ~120 cell labels whereas right panel shows 33 granular cell labels assigned after harmonization.
  
### Harmonized Atlas ~112k cells
![Graph](Examples/LabelHarmonization.png)

## Deduce Relationship Example
To demonstrate Deduce Relationship we used two datasets that were also uses in the harmonization example: 1) E-MTAB-8901 (~21k cells) 2) E-MTAB-10187 (~77k cells). Parameters used: downsample = TRUE, downsample_to=300, classification.method="Ensemble"

### Relative Similarity among celltypes between two datasets
![Graph](Examples/DeduceRelationship_Example.png)

The output of Deduce Relationship is a representation of the confusion matrix as a heatmap, with each square signifying how many cells of a given celltype in one reference were classified as a celltype in the other. It is normalized in such a way that each celltype in reference 2 has a red square that shows the most closely related celltype in reference 1. We see that the related celltypes all make biological sense, such as all immune cells and their subtypes in reference2 being assigned to "Immune cells" reference1, and similarly mesenchyme subtypes being assigned to "Mesoderm 2".

## Citation
Please cite our preprint: https://www.biorxiv.org/content/10.1101/2021.09.07.459342v1 when using ELeFHAnt in your research.

## Bugs/Issues
Please report bugs, issues and improvements using the issues section of the GitHub page.

## Contribute to Code
Please open a pull request to contribute to existing code in terms of improvements / new features / bug removal.
