###################### Load R packages ####################
library(ELeFHAnt)
library(scPred)
library(magrittr)
library(harmony)
library(SingleCellExperiment)
library(scater)
################### #####################################

################################### Preparing Datasets for Benchmarking ###################

setwd('/Volumes/Praneet_Backup/ELeFHAnt_Benchmarking/Benchmarking_Inter')

##### Download Gut_Atlas_duojejunum.rds from ELeFHAnt reference datasets using https://www.dropbox.com/sh/6hd2skriqqlokwp/AAAVol-_qPlCdA4DpERWjkeJa?dl=0

##### Download FetalTissues_human_JasonSpenceLab.rds from ELeFHAnt reference datasets using https://www.dropbox.com/sh/6hd2skriqqlokwp/AAAVol-_qPlCdA4DpERWjkeJa?dl=0

GutData_duojejunum <- readRDS('Gut_Atlas_duojejunum.rds')
fetal_jason_intestine = readRDS('FetalTissues_human_JasonSpenceLab.rds')

Idents(fetal_jason_intestine) <- fetal_jason_intestine$Corrected_organ_group
fetal_jason_intestine <- subset(fetal_jason_intestine, idents = "Intestine")

DefaultAssay(GutData_duojejunum) <- "RNA"
DefaultAssay(fetal_jason_intestine) <- "RNA"

GutData_duojejunum$Celltypes <- GutData_duojejunum$cell_name_detailed
fetal_jason_intestine$Celltypes <- fetal_jason_intestine$Cell_type

Idents(GutData_duojejunum) <- GutData_duojejunum$Celltypes
Idents(fetal_jason_intestine) <- fetal_jason_intestine$seurat_clusters

reference <- subset(GutData_duojejunum, downsample = 300)
query <- subset(fetal_jason_intestine, downsample = 200)

query = FindNeighbors(query)
query = FindClusters(query, resolution = 1)

p4 = DimPlot(reference, group.by = "Celltypes", label=T, repel = T) + NoLegend() + ggtitle("reference Celltypes")
p5 = DimPlot(query, label=T, repel = T, reduction = "umap") + NoLegend() + ggtitle("query seurat_clusters")
p6 = DimPlot(query, label=T, repel = T, reduction = "umap", group.by = "Celltypes") + NoLegend() + ggtitle("query known Celltypes")

p4+p5+p6
ggsave('Benchmark_scData_reference_query.InterDataset.png', width = 15, height = 15)
###########################################################################################

################################## Running ELeFHAnt Celltype Annotation ###################

start_time <- Sys.time()

out.ELeFHAnt <- CelltypeAnnotation(reference = reference, query = query, downsample = FALSE, classification.method = "Ensemble", validatePredictions=FALSE)

end_time <- Sys.time()
p1 = DimPlot(out.ELeFHAnt, group.by = "PredictedCelltype_UsingEnsemble", label=T, repel = T, reduction = "umap", pt.size=1) + NoLegend() + ggtitle("ELeFHAnt Predictions")
p1

Total_time_ELeFHAnt = end_time - start_time
saveRDS(out.ELeFHAnt, file = 'InterDataset.query.annotated.ELeFHAnt.rds')
###########################################################################################

################################## Running scPred Celltype Annotation ######################
start_time <- Sys.time()

reference.scPred <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

reference.scPred <- getFeatureSpace(reference.scPred, "Celltypes")
reference.scPred <- trainModel(reference.scPred)
get_probabilities(reference.scPred) %>% head()
get_scpred(reference.scPred)

query.scPred <- NormalizeData(query)
query.scPred <- scPredict(query.scPred, reference.scPred)

end_time <- Sys.time()
p2 = DimPlot(query.scPred, group.by = "scpred_prediction", repel = TRUE, reduction = "umap", label = T, pt.size=1) + NoLegend() + ggtitle("scPred Predictions")
p2

Total_time_scPred = end_time - start_time
saveRDS(query.scPred, file = 'InterDataset.query.annotated.scPred.rds')
###########################################################################################

################################## Running seurat Label transfer ######################
start_time <- Sys.time()

anchors <- FindTransferAnchors(reference = reference, query = query)
predictions <- TransferData(anchorset = anchors, refdata = reference$Celltypes)
query <- AddMetaData(object = query, metadata = predictions)

end_time <- Sys.time()
p3 = DimPlot(query, group.by = "predicted.id", repel = TRUE, reduction = "umap", label = T, pt.size=1) + NoLegend() + ggtitle("LabelTransfer Predicted Celltypes")
p3

Total_time_LT = end_time - start_time
saveRDS(query, file = 'InterDataset.query.annotated.seuratLabelTransfer.rds')
###########################################################################################

p1+p2+p3
ggsave('Benchmark_CelltypeAnnotation_ELeFHAnt_scPred_LT.InterDataset.png', width = 15, height = 15)

print ('ELeFHAnt Time elapsed')
Total_time_ELeFHAnt

print ('scPred Time elapsed')
Total_time_scPred

print ('seurat LabelTransfer Time elapsed')
Total_time_LT