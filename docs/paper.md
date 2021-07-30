---
title: 'ELeFHAnt: Ensemble Learning for Harmonization and Annotation of Single Cell Data'
tags:
  - R
  - single cell
  - sequencing
  - annotation
  - classification
  - machine learning
authors:
  - name: Konrad Thorner
    affiliation: 1
  -name: Aaron Zorn
    affiliation: 1
  - name: Praneet Chaturvedi
    affiliation: 1 
affiliations:
 - name: Cincinnati Children’s Hospital Medical Center
   index: 1
date: 19 July 2021
bibliography: paper.bib
---

# Summary 

Annotation of single cells has become an important step in the single cell analysis framework. With advances in sequencing technology thousands to millions of cells are processed to understand the intricacies of the biological system in question. There are currently ~200 computational tools available to help researchers automatically annotate single cells using supervised / unsupervised machine learning, cell type markers, or tissue-based markers from bulk RNA-seq to annotate cells. Another popular method to annotate cells is through manual curation of markers to assign cell types based on a priori knowledge, but this is cumbersome given the exponential growth in data that can be generated. With the expansion of publicly available data there is also a need for a tool which can integrate multiple references into a unified cellular atlas, and compare two datasets to help understand the different annotations between them. Here we present ELeFHAnt: Ensemble learning for harmonization and annotation of single cells. ELeFHAnt is an easy to use R package with three main functions: 1) CelltypeAnnotation 2) LabelHarmonization 3) DeduceRelationship. CelltypeAnnotation is a function to annotate cells in a query Seurat object using a reference Seurat object with annotated cell types. LabelHarmonization can be utilized to integrate multiple cell atlases (references) into an unified cellular atlas with harmonized cell types. Finally, DeduceRelationship is a function that compares cell types between two scRNA-seq datasets.  ELeFHAnt can be accessed from GitHub at https://github.com/praneet1988/ELeFHAnt.

# Statement of Need

Single cell sequencing has become a powerful method for understanding biological systems at an increasingly granular level. For single cell RNA-seq data specifically, one of the primary questions is using the transcriptome to determine cell identity. It is common to visualize such data and find clusters of cells with similarity in gene expression, but assigning each cluster a cell type is a much more open-ended task. Taking advantage of publicly-available, annotated datasets in combination with supervised learning is a powerful method for addressing this question. Ensemble Learning for Harmonization and Annotation of Single Cells (ELeFHAnt) provides an easy to use R package for users to annotate clusters of single cells, harmonize labels across single cell datasets to generate a unified atlas, and infer relationships among cell types between two datasets. It provides users with the flexibility of choosing between random forest and SVM (Support Vector Machine) based classifiers or letting ELeFHAnt apply both in combination to make predictions. 

As an alternative to manual annotation, there are many automatic cell annotation tools currently available that employ either gene marker, correlation, or machine learning-based methods, each with varying levels of performance [@Pasquini:2021]. The label transfer functionality of Seurat is among the most well-known, but occurs on an individual cell level rather than a community level, often leading to over-annotation [@Stuart:2019]. There are also deep learning-based tools emerging such as scANVI, which utilizes generative models but requires significantly more computation time [@Xu:2021]. 

ELeFHAnt is a supervised machine learning-based tool that enables researchers to identify cell types in their scRNA-seq data while providing additional unique features. ELeFHAnt gives users the ability to use and compare not just one but multiple classification algorithms simultaneously through its ensemble method, weighting their predictions to produce the best consensus among them. In this ensemble, SVM and random forest are our two classifiers selected for their superior accuracy and computation time in a benchmarking study [@Abdelaal:2019]. Additionally, selecting the optimal reference is a challenge addressed by harmonization, that allows users to integrate multiple datasets together into an atlas. A standardized set of labels is generated across all of them, which can subsequently be used to annotate new datasets. Relationships between two datasets can also be deduced to better understand how each was annotated. This is provided in an easy to interpret heatmap format that compares all the cell types between them. Finally, a subsampling procedure is used to enable faster predictions while being shown not to influence reproducibility.

ELeFHAnt has been tested on public datasets involving multiple time points in fetal development, where it was still able to identify cell types. It has also been used internally within Cincinnati Children’s Hospital Medical Center across different projects, including annotation of snRNA (single nucleus RNA) sequencing with a harmonized gut cell atlas. 

# Overview and Examples

ELeFHAnt makes use of Seurat for the initial input data and pre-processing. It will then generate the training and test sets from the reference and query respectively, with optional subsampling. SVM and Random Forest are the two classifiers that can be used seperately or as an ensemble. Classification accuracy of both are used to assign weights to the predictions from each classifier. These weighted confusion matrices are normalized based on the largest number of cells shared among celltypes and assigned clusters. They are then added together for the final ensemble predictions.

#### Figure 1

![](https://raw.githubusercontent.com/praneet1988/ELeFHAnt/main/ELeFHAnt.png) 

#### Table 1

![](https://raw.githubusercontent.com/praneet1988/ELeFHAnt/main/Examples/gut_datasets.png){ width=500px }

The attributes of our three example datasets of early gut development are shown below, as well as those of the integrated dataset. "Fetal" refers to a subset of terminal ileum (TI) data from an atlas for human fetal intestinal development called "STAR-FINDer" [@Fawkner-Corbett:2021]. "Gut" refers to a subset of duojejunum cell data from the Gut Cell Atlas, which also examines intestinal development from 6-10 weeks post-conception [@Elmentaite:2020]. Lastly, "Spence" refers a subset of fetal intestinal data from a multi-endodermal organ atlas [@Yu:2021].

### Celltype Annotation

Cell type annotation in ELeFHAnt is performed by two sub functions: 1) ApproximationBasedCelltypeAssignment 2) ClassifyCells. ELeFHAnt runs diagnostics on the inputted reference and query Seurat objects to calculate the ratio between number of cells in the query versus reference. If the ratio is greater than 1.5, then ELeFHAnt automatically deploys ApproximationBasedCelltypeAssignment to approximate the best suitable cell type for each cluster of cells in the query. If the ratio is < 1.5, or in other words the number of cells in the reference are larger than query cells, a per cell based annotation is carried out (please see methods for detail). To demonstrate approximation based assignment, we used Gut as reference (~21k cells) and Spence as query (~77k cells) (please see Table 1). ELeFHAnt automatically calculated the ratio and performed the approximation based assignment. Using approximation based assignment ELeFHAnt labeled all clusters with the best approximated cell type (Figure 2).
Similarly, to exhibit the per cell annotation, we downsampled 300 cells per cell type from Gut to use as reference and 200 cells per seurat clusters from Spence to use as query (please see Table 1). As the ratio < 1.5, ELeFHAnt started per-cell annotation to annotate each cell in the query (Figure 3). Comparing the predicted annotations to known cell types in the query, we see 1:1 consensus across broad cell types including endothelial, epithelial, mesenchymal, and neuronal.


#### Figure 2

#### A

![](https://raw.githubusercontent.com/praneet1988/ELeFHAnt/main/Examples/GutCell_Reference.png){ width=450px }

#### B

![](https://raw.githubusercontent.com/praneet1988/ELeFHAnt/main/Examples/FetalIntestine_SpenceLab_Query.png){ width=450px }

#### C
4
![](https://raw.githubusercontent.com/praneet1988/ELeFHAnt/main/Examples/CelltypeAnnotation_Example2.png){ width=450px }

(A) represents celltypes in the reference dataset displayed on a UMAP. (B) represents Seurat clusters displayed on the query, and (C) represents the predicted celltypes as determined by ELeFHAnt's ensemble approach.

### Label Harmonization

Integrating single cell RNA datasets with different sets of cell type assignments can become very difficult to infer the cell type label that an integrated cluster should be associated with. To solve this problem, we designed a function (LabelHarmonization) to harmonize cell types from multiple datasets into an unified atlas with cell types assigned to each cluster of cells in the integrated dataset. To demonstrate ELeFHAnt’s LabelHarmonization we used three datasets with cells profiled from fetal gut development (Table 1). Briefly, we integrate three atlases using Seurat’s canonical correlation analysis based integration algorithm and then create training and testing sets to learn and harmonize cell types (please see methods for details). In Figure 4, the left panel shows integrated cells labeled with 145 cell types that each reference contributed; whereas on the right panel integrated cells labeled with harmonized cell types using ensemble learning in ELeFHAnt are shown. We can clearly assess the cell type assigned to each cell community instead of needing to manually annotate cell types. For example, on the left panel the all neuronal cell types can be seen which were present in each dataset whereas on the right ELeFHAnt harmonized and even found subtypes to annotate the integrated clusters.


#### Figure 4

![](https://raw.githubusercontent.com/praneet1988/ELeFHAnt/main/Examples/HarmonizationExample_ELeFHAnt.png)

The UMAP on the left represents the 144 total celltypes from the three datasets after integration. The UMAP on the right is the result of ELeFHAnt's harmonization with the ensemble method, showing the resulting labels for each cluster. 

### Deduce Relationship

With the growing number of single cell RNA-Seq (scRNA-Seq) datasets, hypothesis generation can be piloted by utilizing multiple scRNA datasets. Comparing datasets to infer similarity or comparing in-house datasets to publicly available datasets is important for facilitating experimental design and finding similarity across cell types. The DeduceRelationship function in ELeFHAnt compares the cell types between two datasets to calculate relative similarity. In addition to hypothesis generation this function can be helpful in comparing an in-house sequencing dataset to multiple public datasets for assessing robustness of sequencing, annotation, or clustering. To demonstrate, we compared Gut vs Spence and found the Mesoderm2 cell type in Gut is similar to the Mesenchyme subtypes in Spence. Similarly, all immune sub-cell types were similar to the Immune cell type in Gut (Figure 5).

#### Figure 5

![](https://raw.githubusercontent.com/praneet1988/ELeFHAnt/main/Examples/DeduceRelationship_Example.png)

The heatmap depicts the relationship between celltypes for the two references, with each red square showing which cell type in reference 2 best matches a particular celltype in reference 1. 

# Acknowledgements

We would like to thank Drs. Emily Miraldi, Nathan Salomonis and Anil Jegga at Cincinnati Children's Hospital Medical Center for their valuable feedback. 

# Author Contributions

K.T. and P.C. designed, implemented, and tested the algorithm. A.Z. provided feedback for testing and improving the algorithm from a research perspective.

# Competing Interests

Authors have no competing interests. 

# References
