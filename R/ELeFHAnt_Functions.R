#' ELeFHAnt Celltype Annotation
#'
#' Celltype annotation is a function to annotate celltypes in a single cell datasets.
#' It requires a reference dataset (a processed Seurat Object with Celltypes column in metadata) and a query dataset (a processed 
#' seurat object with seurat_clusters column in metadata). One can choose from randomForest, SVM or Ensemble classifiction method
#' to learn celltypes from reference dataset and then predict celltypes for query dataset.
#' 
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import kernlab
#' @import dbscan
#' @import randomForest
#' @import caTools
#' @import e1071
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference and query enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 100] in reference and query for Celltypes and seurat_clusters resspectively
#' @param classification.method choose classification method for learning and predicting celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses estimation robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis
#' @param selectvarfeatures number of variable features to select while training (default: 2000)
#' @return query seurat object with predictions added to meta.data of the object
#' @export
CelltypeAnnotation <- function(reference = NULL, query = NULL, downsample = FALSE, downsample_to = 100, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 10, validatePredictions = TRUE, selectvarfeatures = 2000) {
    if(downsample == TRUE)
    {
        message ("Setting Assay of reference and query to RNA")
        DefaultAssay(reference) <- "RNA"
        DefaultAssay(query) <- "RNA"

        message ("Donwsampling Reference and Query")
        reference$Dataset <- rep("reference", ncol(reference))
        query$Dataset <- rep("query", ncol(query))

        Idents(reference) <- reference$Celltypes
        Idents(query) <- query$seurat_clusters

        reference_use <- subset(reference, downsample = downsample_to)
        query_use <- subset(query, downsample = downsample_to)
    }

    if(downsample == FALSE)
    {
        message ("Setting Assay of reference and query to RNA")
        DefaultAssay(reference) <- "RNA"
        DefaultAssay(query) <- "RNA"

        reference$Dataset <- rep("reference", ncol(reference))
        query$Dataset <- rep("query", ncol(query))

        Idents(reference) <- reference$Celltypes
        Idents(query) <- query$seurat_clusters

        reference_use <- reference
        query_use <- query
    }

    message ("Merging reference and query")
    combined <- merge(x=reference_use, y=query_use)
    message ("Normalization, Variable Feature Selection and scaling")
    combined <- NormalizeData(combined)
    combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = selectvarfeatures)
    combined <- ScaleData(combined)
    combined_exp <- combined[['RNA']]@scale.data
    combined_exp <- t(combined_exp)
    combined_exp <- data.frame(combined_exp)
    num_features <- ncol(combined_exp)
    message (paste0("Number of Features selected:", num_features))
    combined_exp$Dataset <- combined$Dataset
    message ("Generating train and test sets")
    X <- split(combined_exp, combined_exp$Dataset)
    X$reference$Celltypes <- reference_use$Celltypes
    X$query$Clusters <- query_use$seurat_clusters

    if(classification.method == "randomForest")
    {
        message ("Setting up randomForest classifier learning")
        message ("Training randomForest classifier")
        rf_Celltypes.1 = randomForest_predictor(train = X$reference[,1:num_features], test = X$query[,1:num_features], train_label = X$reference$Celltypes, test_label = X$query$Clusters, ntree = 1000)
        message ("Predicting using trained randomForest classifier")
        rf_pred.1 = predict(rf_Celltypes.1, newdata=X$query[,1:num_features])
        rf_cm.1 = table(X$query$Clusters, rf_pred.1)

        message ("Calculating weights for each randomForest classifier")
        rf_acccuracy_estimate.1 <- (1-tail(rf_Celltypes.1$err.rate[,1], 1))*100
        message (paste0("Accuray estimate of randomForest classifier:", rf_acccuracy_estimate.1))

        message ("Assigning weights to randomForest predictions")
        rf_cm.1 <- as.matrix(rf_cm.1) * rf_acccuracy_estimate.1

        rf_cm <- rf_cm.1

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("PredictedCelltype_UsingRF", "seurat_clusters")
        PredictedCelltype_UsingRF <- as.character(rf_celltype_pred[match(query$seurat_clusters, rf_celltype_pred$seurat_clusters), "PredictedCelltype_UsingRF"])
        query[["PredictedCelltype_UsingRF"]] <- PredictedCelltype_UsingRF
        message ("Added Predicted celltypes using randomForest to query")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message("randomForest based learning and celltype annotation completed. Starting validation of celltype assignments using GSEA")
            reference.validation.use <- subset(reference_use, idents = rf_celltype_pred$PredictedCelltype_UsingRF)
            validation = ValidatePredictions(reference = reference.validation.use, query = query_use)
            message ("Validation completed. Please see summary of GSEA below")
            print (validation)
            write.table(validation, "Summary_GeneSetEnrichmentAnalysis.txt", quote=F, sep="\t")
            return(query)
        }
        if(validatePredictions == FALSE)
        {
            message("randomForest based learning and celltype annotation completed")
            return(query)
        }
    }

    if(classification.method == "SVM")
    {
        message ("Setting up SVM classifier learning")
        message ("Training SVM classifier")
        svm_Celltypes.1 = svm_predictor(train = X$reference[,1:num_features], test = X$query[,1:num_features], train_label = X$reference$Celltypes, test_label = X$query$Clusters, crossvalidationSVM = crossvalidationSVM, cachesize = 100, cost = 10, kernel = "linear")
        
        message ("Predicting using trained SVM classifier")
        svm_pred.1 = predict(svm_Celltypes.1, newdata=X$query[,1:num_features])
        svm_cm.1 = table(X$query$Clusters, svm_pred.1)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate.1 <- svm_Celltypes.1$tot.accuracy
        message (paste0("Accuray estimate of SVM classifier:", svm_accuracy_estimate.1))

        message ("Assigning weights to SVM predictions")
        svm_cm.1 <- as.matrix(svm_cm.1) * svm_accuracy_estimate.1

        svm_cm = svm_cm.1

        svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
        colnames(svm_celltype_pred) <- c("PredictedCelltype_UsingSVM", "seurat_clusters")
        PredictedCelltype_UsingSVM <- as.character(svm_celltype_pred[match(query$seurat_clusters, svm_celltype_pred$seurat_clusters), "PredictedCelltype_UsingSVM"])
        query[["PredictedCelltype_UsingSVM"]] <- PredictedCelltype_UsingSVM
        message ("Added Predicted celltypes using SVM to query")
        write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message("SVM based learning and celltype annotation completed. Starting validation of celltype assignments using GSEA")
            reference.validation.use <- subset(reference_use, idents = svm_celltype_pred$PredictedCelltype_UsingSVM)
            validation = ValidatePredictions(reference = reference.validation.use, query = query_use)
            message ("Validation completed. Please see summary of GSEA below")
            print (validation)
            write.table(validation, "Summary_GeneSetEnrichmentAnalysis.txt", quote=F, sep="\t")
            return(query)
        }
        if(validatePredictions == FALSE)
        {
            message ("SVM based learning and celltype annotation completed")
            return(query)
        }
    }

    if(classification.method == "Ensemble")
    {
        message ("Ensemble learning using classification accuracy of both Random Forest and SVM classifiers")
        message ("Setting up randomForest classifier learning")
        message ("Training randomForest classifier")
        rf_Celltypes.1 = randomForest_predictor(train = X$reference[,1:num_features], test = X$query[,1:num_features], train_label = X$reference$Celltypes, test_label = X$query$Clusters, ntree = 1000)
        message ("Predicting using trained randomForest classifier")
        rf_pred.1 = predict(rf_Celltypes.1, newdata=X$query[,1:num_features])
        rf_cm.1 = table(X$query$Clusters, rf_pred.1)

        message ("Calculating weights for each randomForest classifier")
        rf_acccuracy_estimate.1 <- (1-tail(rf_Celltypes.1$err.rate[,1], 1))*100
        message (paste0("Accuray estimate of randomForest classifier:", rf_acccuracy_estimate.1))

        message ("Assigning weights to randomForest predictions")
        rf_cm.1 <- as.matrix(rf_cm.1) * rf_acccuracy_estimate.1

        rf_cm <- rf_cm.1

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("PredictedCelltype_UsingRF", "seurat_clusters")
        PredictedCelltype_UsingRF <- as.character(rf_celltype_pred[match(query$seurat_clusters, rf_celltype_pred$seurat_clusters), "PredictedCelltype_UsingRF"])
        query[["PredictedCelltype_UsingRF"]] <- PredictedCelltype_UsingRF
        message ("Added Predicted celltypes using randomForest to query")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        
        message ("Setting up SVM classifier learning")
        message ("Training SVM classifier")
        svm_Celltypes.1 = svm_predictor(train = X$reference[,1:num_features], test = X$query[,1:num_features], train_label = X$reference$Celltypes, test_label = X$query$Clusters, crossvalidationSVM = crossvalidationSVM, cachesize = 100, cost = 10, kernel = "linear")
        
        message ("Predicting using trained SVM classifier")
        svm_pred.1 = predict(svm_Celltypes.1, newdata=X$query[,1:num_features])
        svm_cm.1 = table(X$query$Clusters, svm_pred.1)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate.1 <- svm_Celltypes.1$tot.accuracy
        message (paste0("Accuray estimate of SVM classifier:", svm_accuracy_estimate.1))

        message ("Assigning weights to SVM predictions")
        svm_cm.1 <- as.matrix(svm_cm.1) * svm_accuracy_estimate.1

        svm_cm = svm_cm.1

        svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
        colnames(svm_celltype_pred) <- c("PredictedCelltype_UsingSVM", "seurat_clusters")
        PredictedCelltype_UsingSVM <- as.character(svm_celltype_pred[match(query$seurat_clusters, svm_celltype_pred$seurat_clusters), "PredictedCelltype_UsingSVM"])
        query[["PredictedCelltype_UsingSVM"]] <- PredictedCelltype_UsingSVM
        message ("Added Predicted celltypes using SVM to query")
        write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")

        message ("randomForest and SVM based learning and predictions completed. Using predictions from RF and SVM to make Ensemble Predictions")

        consensus_cm = rf_cm/max(rf_cm) + svm_cm/max(svm_cm)
        consensus_celltype_pred <- data.frame(colnames(consensus_cm)[apply(consensus_cm,1,which.max)], rownames(consensus_cm), check.names=F)
        colnames(consensus_celltype_pred) <- c("PredictedCelltype_UsingEnsemble", "seurat_clusters")
        PredictedCelltype_UsingEnsemble <- as.character(consensus_celltype_pred[match(query$seurat_clusters, consensus_celltype_pred$seurat_clusters), "PredictedCelltype_UsingEnsemble"])
        query[["PredictedCelltype_UsingEnsemble"]] <- PredictedCelltype_UsingEnsemble
        message ("Added Predicted celltypes using Ensemble learning to query")
        write.table(consensus_cm, "ConfusionMatrix_EnsembleLearning.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message("Ensembl celltype annotation completed. Starting validation of celltype assignments using GSEA")
            reference.validation.use <- subset(reference_use, idents = consensus_celltype_pred$PredictedCelltype_UsingEnsemble)
            validation = ValidatePredictions(reference = reference.validation.use, query = query_use)
            message ("Validation completed. Please see summary of GSEA below")
            print (validation)
            write.table(validation, "Summary_GeneSetEnrichmentAnalysis.txt", quote=F, sep="\t")
            return(query)
        }
        if(validatePredictions == FALSE)
        {
            message("Ensembl celltype annotation completed.")
            return(query)
        }
    }
 }

#' ELeFHAnt Label Harmonization
#'
#' Label Harmonization is a function to harmonize cell labels (celltypes) across single cell datasets.
#' It requires a list of processed Seurat Objects with Celltypes column in metadata or a integrated seurat object
#' (seurat object with Celltypes and seurat_clusters column in metadata). One can choose from randomForest, SVM or Ensemble classifiction method
#' to harmonize celltypes. Please See: DefaultAssay of each object should be set to "RNA"
#' 
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import kernlab
#' @import dbscan
#' @import randomForest
#' @import caTools
#' @import e1071
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @param seurat.objects a list of processed seurat objects (please set Default Assay to "RNA") with Celltypes column in their respective meta.data to perform integration on
#' @param perform_integration logical Indicator (TRUE or FALSE) to perform integration using list of seurat.objects
#' @param integrated.atlas an integrated seurat object with CellTypes and seurat_clusters column in meta.data. Required if perform_integration = FALSE
#' @param downsample logical Indicator (TRUE or FALSE) to downsample seurat objects enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 100]
#' @param npcs number of principal components to compute after integration
#' @param resolution value of the resolution parameter, decides size of cell communities.
#' @param classification.method choose classification method for learning and harmonizing celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses learning robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis
#' @param selectanchorfeatures number of anchor features to use for integrating datasets (Default: 2000)
#' @return integrated seurat object with harmonized celltypes added to meta.data of the object
#' @export
LabelHarmonization <- function(seurat.objects = c(), perform_integration = TRUE, integrated.atlas = NULL, downsample = TRUE, downsample_to = 100, npcs = 30, resolution = 0.5, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 10, validatePredictions = TRUE, selectanchorfeatures = 2000) {
    if(perform_integration == TRUE)
    {
        if(downsample == TRUE)
        {
            message ("Downsampling seurat objects")
            seurat.objects <- lapply(X = seurat.objects, FUN = function(x) {
                DefaultAssay(x) <- "RNA"
                Idents(x) <- x$Celltypes
                x <- subset(x, downsample = downsample_to)
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = selectanchorfeatures)
            })
        }
        if(downsample == FALSE)
        {
            seurat.objects <- seurat.objects
        }
        message ("Starting integration using Seurat")
        integration.anchors <- FindIntegrationAnchors(object.list = seurat.objects, anchor.features = selectanchorfeatures)
        integrated.atlas <- IntegrateData(anchorset = integration.anchors)
        DefaultAssay(integrated.atlas) <- "integrated"
        message ("Integration Completed. Performing Scaling, Dimension reduction and clustering")
        integrated.atlas <- ScaleData(integrated.atlas, verbose = FALSE)
        integrated.atlas <- RunPCA(integrated.atlas, npcs = npcs, verbose = FALSE)
        integrated.atlas <- RunUMAP(integrated.atlas, reduction = "pca", dims = 1:npcs)
        integrated.atlas <- FindNeighbors(integrated.atlas, reduction = "pca", dims = 1:npcs)
        integrated.atlas <- FindClusters(integrated.atlas, resolution = resolution)
        integrated.use <- integrated.atlas
    }
    if(perform_integration == FALSE)
    {
        if(downsample == TRUE)
        {
            message ("Setting Assay of integrated.atlas to integrated and downsampling")
            integrated.atlas <- integrated.atlas
            DefaultAssay(integrated.atlas) <- "integrated"
            Idents(integrated.atlas) <- integrated.atlas$Celltypes
            integrated.use <- subset(integrated.atlas, downsample = downsample_to)
        }
        if(downsample == FALSE)
        {
            message ("Setting Assay of integrated.atlas to integrated")
            integrated.atlas <- integrated.atlas
            DefaultAssay(integrated.atlas) <- "integrated"
            integrated.use <- integrated.atlas
        }
    }
    message ("Generating train and test datasets using stratification -- 60% for training & 40% for testing")
    scaled_data <- integrated.use[['integrated']]@scale.data
    scaled_data <- t(scaled_data)
    scaled_data <- data.frame(scaled_data)
    num_features <- ncol(scaled_data)
    message (paste0("Number of Features selected:", num_features))
    scaled_data$Celltypes <- integrated.use$Celltypes
    scaled_data$Clusters <- integrated.use$seurat_clusters
    scaled_data_stratified <- stratified(scaled_data, group = "Clusters", size = 0.6, bothSets = TRUE)
    train <- scaled_data_stratified$SAMP1[,1:num_features]
    test <- scaled_data_stratified$SAMP2[,1:num_features]
    train_label <- scaled_data_stratified$SAMP1$Celltypes
    test_label <- scaled_data_stratified$SAMP2$Clusters

    if(classification.method == "randomForest")
    {
        message ("Setting up randomForest classifier learning")
        message ("Training randomForest classifier")
        rf_Celltypes.1 = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = 1000)

        message ("Predicting using trained randomForest classifier")
        rf_pred.1 = predict(rf_Celltypes.1, newdata=test)
        rf_cm.1 = table(test_label, rf_pred.1)
    

        message ("Calculating weights for randomForest classifier")
        rf_acccuracy_estimate.1 <- (1-tail(rf_Celltypes.1$err.rate[,1], 1))*100
        message (paste0("Accuray estimate of randomForest classifier:", rf_acccuracy_estimate.1))

        message ("Assigning weights to randomForest predictions")
        rf_cm.1 <- as.matrix(rf_cm.1) * rf_acccuracy_estimate.1

        rf_cm <- rf_cm.1

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("HarmonizedLabels_UsingRF", "seurat_clusters")
        HarmonizedLabels_UsingRF <- as.character(rf_celltype_pred[match(integrated.atlas$seurat_clusters, rf_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingRF"])
        integrated.atlas[["HarmonizedLabels_UsingRF"]] <- HarmonizedLabels_UsingRF
        message ("Added Harmonized Labels using randomForest to integrated object")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message("randomForest based learning and harmonization completed. Starting validation of celltype assignments using GSEA")
            Idents(integrated.atlas) <- integrated.atlas$Celltypes
            integrated.validation.use <- subset(integrated.atlas, idents = rf_celltype_pred$HarmonizedLabels_UsingRF)
            validation = ValidatePredictions(reference = integrated.validation.use, query = integrated.atlas)
            message ("Validation completed. Please see summary of GSEA below")
            print (validation)
            write.table(validation, "Summary_GeneSetEnrichmentAnalysis.txt", quote=F, sep="\t")
            return(integrated.atlas)
        }
        if(validatePredictions == FALSE)
        {
            message("randomForest based learning and harmonization completed")
            return(integrated.atlas)
        }
    }

    if(classification.method == "SVM")
    {
        message ("Setting up SVM classifier learning")
        message ("Training SVM classifier")
        svm_Celltypes.1 = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cachesize = 100, cost = 10, kernel = "linear")
        
        message ("Predicting using trained SVM classifier")
        svm_pred.1 = predict(svm_Celltypes.1, newdata=test)
        svm_cm.1 = table(test_label, svm_pred.1)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate.1 <- svm_Celltypes.1$tot.accuracy
        message (paste0("Accuray estimate of SVM classifier:", svm_accuracy_estimate.1))

        message ("Assigning weights to SVM predictions")
        svm_cm.1 <- as.matrix(svm_cm.1) * svm_accuracy_estimate.1

        svm_cm <- svm_cm.1

        svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
        colnames(svm_celltype_pred) <- c("HarmonizedLabels_UsingSVM", "seurat_clusters")
        HarmonizedLabels_UsingSVM <- as.character(svm_celltype_pred[match(integrated.atlas$seurat_clusters, svm_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingSVM"])
        integrated.atlas[["HarmonizedLabels_UsingSVM"]] <- HarmonizedLabels_UsingSVM
        message ("Added harmonized labels using SVM to integrated object")
        write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message ("SVM based learning and harmonization completed. Starting validation of celltype assignments using GSEA")
            Idents(integrated.atlas) <- integrated.atlas$Celltypes
            integrated.validation.use <- subset(integrated.atlas, idents = svm_celltype_pred$HarmonizedLabels_UsingSVM)
            validation = ValidatePredictions(reference = integrated.validation.use, query = integrated.atlas)
            message ("Validation completed. Please see summary of GSEA below")
            print (validation)
            write.table(validation, "Summary_GeneSetEnrichmentAnalysis.txt", quote=F, sep="\t")
            return(integrated.atlas)
        }
        if(validatePredictions == FALSE)
        {
            message ("SVM based learning and harmonization completed")
            return(integrated.atlas)
        }
    }

    if(classification.method == "Ensemble")
    {
        message ("Ensemble learning using classification accuracy of both Random Forest and SVM classifiers")
        message ("Setting up randomForest classifier learning")
        message ("Training randomForest classifier")
        rf_Celltypes.1 = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = 1000)

        message ("Predicting using trained randomForest classifier")
        rf_pred.1 = predict(rf_Celltypes.1, newdata=test)
        rf_cm.1 = table(test_label, rf_pred.1)
    
        message ("Calculating weights for randomForest classifier")
        rf_acccuracy_estimate.1 <- (1-tail(rf_Celltypes.1$err.rate[,1], 1))*100
        message (paste0("Accuray estimate of randomForest classifier:", rf_acccuracy_estimate.1))

        message ("Assigning weights to randomForest predictions")
        rf_cm.1 <- as.matrix(rf_cm.1) * rf_acccuracy_estimate.1

        rf_cm <- rf_cm.1

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("HarmonizedLabels_UsingRF", "seurat_clusters")
        HarmonizedLabels_UsingRF <- as.character(rf_celltype_pred[match(integrated.atlas$seurat_clusters, rf_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingRF"])
        integrated.atlas[["HarmonizedLabels_UsingRF"]] <- HarmonizedLabels_UsingRF
        message ("Added Harmonized Labels using randomForest to integrated object")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        
        message ("Setting up SVM classifier learning")
        message ("Training SVM classifier")
        svm_Celltypes.1 = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cachesize = 100, cost = 10, kernel = "linear")
        
        message ("Predicting using trained SVM classifier")
        svm_pred.1 = predict(svm_Celltypes.1, newdata=test)
        svm_cm.1 = table(test_label, svm_pred.1)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate.1 <- svm_Celltypes.1$tot.accuracy
        message (paste0("Accuray estimate of SVM classifier:", svm_accuracy_estimate.1))

        message ("Assigning weights to SVM predictions")
        svm_cm.1 <- as.matrix(svm_cm.1) * svm_accuracy_estimate.1

        svm_cm <- svm_cm.1

        svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
        colnames(svm_celltype_pred) <- c("HarmonizedLabels_UsingSVM", "seurat_clusters")
        HarmonizedLabels_UsingSVM <- as.character(svm_celltype_pred[match(integrated.atlas$seurat_clusters, svm_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingSVM"])
        integrated.atlas[["HarmonizedLabels_UsingSVM"]] <- HarmonizedLabels_UsingSVM
        message ("Added harmonized labels using SVM to integrated object")
        write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")

        message ("randomForest and SVM based learning and harmonization completed. Using predictions from all models for Ensemble harmonization")

        consensus_cm = rf_cm/max(rf_cm) + svm_cm/max(svm_cm)
        
        consensus_celltype_pred <- data.frame(colnames(consensus_cm)[apply(consensus_cm,1,which.max)], rownames(consensus_cm), check.names=F)
        colnames(consensus_celltype_pred) <- c("HarmonizedLabels_UsingEnsemble", "seurat_clusters")
        HarmonizedLabels_UsingEnsemble <- as.character(consensus_celltype_pred[match(integrated.atlas$seurat_clusters, consensus_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingEnsemble"])
        integrated.atlas[["HarmonizedLabels_UsingEnsemble"]] <- HarmonizedLabels_UsingEnsemble
        message ("Added Harmonized labels using Ensemble learning to query")
        write.table(consensus_cm, "ConfusionMatrix_EnsembleLearning.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message ("Ensembl harmonization completed. Starting validation of celltype assignments using GSEA")
            Idents(integrated.atlas) <- integrated.atlas$Celltypes
            integrated.validation.use <- subset(integrated.atlas, idents = consensus_celltype_pred$HarmonizedLabels_UsingEnsemble)
            validation = ValidatePredictions(reference = integrated.validation.use, query = integrated.atlas)
            message ("Validation completed. Please see summary of GSEA below")
            print (validation)
            write.table(validation, "Summary_GeneSetEnrichmentAnalysis.txt", quote=F, sep="\t")
            return(integrated.atlas)
        }
        if(validatePredictions == FALSE)
        {
            message ("Ensembl harmonization completed.")
            return(integrated.atlas)
        }
    }
 }

#' ELeFHAnt Deduce Relationship
#'
#' Deduce Relationship is a function that hels deduce relationships among celltypes between two datasets.
#' It requires two datasets (both processed Seurat Objects with Celltypes column in metadata). One can choose from randomForest, SVM or Ensemble classifiction method
#' to obtain relationships among celltypes. Function outputs a heatmap with celltypes from dataset1 as rows and celltypes from dataset2 as columns.
#' 
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import kernlab
#' @import dbscan
#' @import randomForest
#' @import caTools
#' @import e1071
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @param reference1 a processed Seurat Object with Celltypes column in metadata
#' @param reference2 a processed seurat object with Celltypes column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference and query enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 100] in reference and query for Celltypes and seurat_clusters resspectively
#' @param classification.method choose classification method for learning and predicting celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses estimation robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model
#' @param selectvarfeatures number of variable features to select while training (default: 2000)
#' @return ggplot2 heatmap object and heatmap is automatically saved
#' @export
DeduceRelationship <- function(reference1 = NULL, reference2 = NULL, downsample = FALSE, downsample_to = 100, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 5, selectvarfeatures = 2000) {
  if(downsample == TRUE)
  {
    message ("Setting Assay of reference1 and reference2 to RNA")
    DefaultAssay(reference1) <- "RNA"
    DefaultAssay(reference2) <- "RNA"
    
    message ("Downsampling reference1 and reference2")
    reference1$Dataset <- rep("reference1", ncol(reference1))
    reference2$Dataset <- rep("reference2", ncol(reference2))
    
    Idents(reference1) <- reference1$Celltypes
    Idents(reference2) <- reference2$Celltypes
    
    reference1_use <- subset(reference1, downsample = downsample_to)
    reference2_use <- subset(reference2, downsample = downsample_to)
  }
  
  if(downsample == FALSE)
  {
    message ("Setting Assay of reference1 and reference2 to RNA")
    DefaultAssay(reference1) <- "RNA"
    DefaultAssay(reference2) <- "RNA"
    
    reference1$Dataset <- rep("reference1", ncol(reference1))
    reference2$Dataset <- rep("reference2", ncol(reference2))
    
    Idents(reference1) <- reference1$Celltypes
    Idents(reference2) <- reference2$Celltypes
    
    reference1_use <- reference1
    reference2_use <- reference2
  }
  
  message ("Merging reference1 and reference2")
  combined <- merge(x=reference1_use, y=reference2_use)
  message ("Normalization, Variable Feature Selection and scaling")
  combined <- NormalizeData(combined)
  combined <- FindVariableFeatures(combined)
  combined <- ScaleData(combined)
  combined_exp <- combined[['RNA']]@scale.data
  combined_exp <- t(combined_exp)
  combined_exp <- data.frame(combined_exp)
  num_features <- ncol(combined_exp)
  message (paste0("Number of Features selected:", num_features))
  combined_exp$Dataset <- combined$Dataset
  message ("Generating train and test sets")
  X <- split(combined_exp, combined_exp$Dataset)
  X$reference1$Celltypes <- reference1_use$Celltypes
  X$reference2$Celltypes <- reference2_use$Celltypes
  
  if(classification.method == "randomForest")
  {
    message ("Setting up randomForest classifier learning.")
    message ("Training randomForest classifier")
    rf_Celltypes.1 = randomForest_predictor(train = X$reference1[,1:num_features], test = X$reference2[,1:num_features], train_label = X$reference1$Celltypes, test_label = X$reference2$Celltypes, ntree = 1000)
    
    message ("Predicting using trained randomForest classifier")
    rf_pred.1 = predict(rf_Celltypes.1, newdata=X$reference2[,1:num_features])
    rf_cm.1 = table(X$reference2$Celltypes, rf_pred.1)
    
    message ("Calculating weight for randomForest classifier")
    rf_acccuracy_estimate.1 <- (1-tail(rf_Celltypes.1$err.rate[,1], 1))*100
    message (paste0("Accuray estimate of randomForest classifier:", rf_acccuracy_estimate.1))
    
    message ("Assigning weights to randomForest predictions")
    rf_cm.1 <- as.matrix(rf_cm.1) * rf_acccuracy_estimate.1
    
    message ("Generating confusion matrix and heatmap")
    rf_cm <- rf_cm.1
    write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
    rf_cm_norm <- round(rf_cm/apply(rf_cm,1,max),3)
    rf_df <- as.data.frame(rf_cm_norm)
    colnames(rf_df) <- c("reference2","reference1","Cells")
    plot = ggplot(data = rf_df, aes(x=reference2, y=reference1, fill=Cells)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90))
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_RandomForest.png", width = 10, height = 10)
    message("randomForest based learning and relationship inference completed")
    return(plot)
  }
  
  if(classification.method == "SVM")
  {
    message ("Setting up SVM classifier learning.")
    message ("Training SVM classifier")
    svm_Celltypes.1 = svm_predictor(train = X$reference1[,1:num_features], test = X$reference2[,1:num_features], train_label = X$reference1$Celltypes, test_label = X$reference2$Celltypes, crossvalidationSVM = crossvalidationSVM, cachesize = 100, cost = 10)
    
    message ("Predicting using trained SVM classifier")
    svm_pred.1 = predict(svm_Celltypes.1, newdata=X$reference2[,1:num_features])
    svm_cm.1 = table(X$reference2$Celltypes, svm_pred.1)
    
    message ("Calculating weight for SVM classifier")
    svm_accuracy_estimate.1 <- svm_Celltypes.1$tot.accuracy
    message (paste0("Accuray estimate of SVM classifier:", svm_accuracy_estimate.1))
    
    message ("Assigning weights to SVM predictions")
    svm_cm.1 <- as.matrix(svm_cm.1) * svm_accuracy_estimate.1
    
    message ("Generating confusion matrix and heatmap")
    svm_cm <- svm_cm.1
    write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
    svm_cm_norm <- round(svm_cm/apply(svm_cm,1,max),3)
    svm_df <- as.data.frame(svm_cm_norm)
    colnames(svm_df) <- c("reference2","reference1","Cells")
    plot = ggplot(data = svm_df, aes(x=reference2, y=reference1, fill=Cells)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90))
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_SVM.png", width = 10, height = 10)
    message ("SVM based learning and relationship inference completed")
    return(plot)
  }
  
  if(classification.method == "Ensemble")
  {
    message ("Ensemble learning using classification accuracy of both Random Forest and SVM classifiers")
    message ("Setting up randomForest classifier learning.")
    message ("Training randomForest classifier")
    rf_Celltypes.1 = randomForest_predictor(train = X$reference1[,1:num_features], test = X$reference2[,1:num_features], train_label = X$reference1$Celltypes, test_label = X$reference2$Celltypes, ntree = 1000)
    
    message ("Predicting using trained randomForest classifier")
    rf_pred.1 = predict(rf_Celltypes.1, newdata=X$reference2[,1:num_features])
    rf_cm.1 = table(X$reference2$Celltypes, rf_pred.1)
    
    message ("Calculating weight for randomForest classifier")
    rf_acccuracy_estimate.1 <- (1-tail(rf_Celltypes.1$err.rate[,1], 1))*100
    message (paste0("Accuray estimate of randomForest classifier:", rf_acccuracy_estimate.1))
    
    message ("Assigning weights to randomForest predictions")
    rf_cm.1 <- as.matrix(rf_cm.1) * rf_acccuracy_estimate.1
    
    message ("Generating confusion matrix and heatmap")
    rf_cm <- rf_cm.1
    write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
    rf_cm_norm <- round(rf_cm/apply(rf_cm,1,max),3)
    rf_df <- as.data.frame(rf_cm_norm)
    colnames(rf_df) <- c("reference2","reference1","Cells")
    plot = ggplot(data = rf_df, aes(x=reference2, y=reference1, fill=Cells)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90))
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_RandomForest.png", width = 10, height = 10)
    
    message ("Setting up SVM classifier learning.")
    message ("Training SVM classifier")
    svm_Celltypes.1 = svm_predictor(train = X$reference1[,1:num_features], test = X$reference2[,1:num_features], train_label = X$reference1$Celltypes, test_label = X$reference2$Celltypes, crossvalidationSVM = crossvalidationSVM, cachesize = 100, cost = 10)
    
    message ("Predicting using trained SVM classifier")
    svm_pred.1 = predict(svm_Celltypes.1, newdata=X$reference2[,1:num_features])
    svm_cm.1 = table(X$reference2$Celltypes, svm_pred.1)
    
    message ("Calculating weight for SVM classifier")
    svm_accuracy_estimate.1 <- svm_Celltypes.1$tot.accuracy
    message (paste0("Accuray estimate of SVM classifier:", svm_accuracy_estimate.1))
    
    message ("Assigning weights to SVM predictions")
    svm_cm.1 <- as.matrix(svm_cm.1) * svm_accuracy_estimate.1
    
    message ("Generating confusion matrix and heatmap")
    svm_cm <- svm_cm.1
    write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
    svm_cm_norm <- round(svm_cm/apply(svm_cm,1,max),3)
    svm_df <- as.data.frame(svm_cm_norm)
    colnames(svm_df) <- c("reference2","reference1","Cells")
    plot = ggplot(data = svm_df, aes(x=reference2, y=reference1, fill=Cells)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90))
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_SVM.png", width = 10, height = 10)
    
    message ("randomForest and SVM based learning and relationship inference completed. Using predictions from all models to make Ensemble Predictions")
    
    message ("Generating confusion matrix and heatmap")
    consensus_cm = rf_cm/max(rf_cm) + svm_cm/max(svm_cm)
    write.table(consensus_cm, "ConfusionMatrix_EnsembleLearning.txt", quote=F, sep="\t")
    consensus_cm_norm <- round(consensus_cm/apply(consensus_cm,1,max),3)
    consensus_df <- as.data.frame(consensus_cm_norm)
    colnames(consensus_df) <- c("reference2","reference1","Cells")
    plot = ggplot(data = consensus_df, aes(x=reference2, y=reference1, fill=Cells)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90))
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_Ensemble.png", width = 10, height = 10)
    message("Ensemble based learning and relationship inference completed")
    return(plot)
  }
}


#' ELeFHAnt Validate Predictions
#'
#' Validate Predictions is a function to validate celltype assignments. It uses Gene Set Enrichment Analysis (GSEA) to asses the enrichment.
#' It requires a reference dataset (a processed Seurat Object with Celltypes column in metadata) and a query dataset (a processed 
#' seurat object with seurat_clusters column in metadata).
#' 
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import kernlab
#' @import dbscan
#' @import randomForest
#' @import caTools
#' @import e1071
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @return GSEA result
#' @export

ValidatePredictions <- function(reference = NULL, query = NULL) {

    Idents(reference) <- reference$Celltypes
    message ("Obtaining markers per celltype")
    celltype_markers <- FindAllMarkers(reference)
    top_celltype <- celltype_markers %>% group_by(cluster) %>% slice(1:100)
    write.table(top_celltype, "temp.R.txt", quote = F, sep = "\t")
    top_celltype <- read.table('temp.R.txt', sep = "\t", header = T)
    message ("Generating Gene sets using top 100 markers per celltype from reference")
    fgsea_sets <- top_celltype %>% split(x = .$gene, f = .$cluster)
    

    message ("Obtaining markers per annotated cluster")
    Idents(query) <- query$seurat_clusters
    cluster_markers <- FindAllMarkers(query)
    top_cluster <- cluster_markers %>% group_by(cluster) %>% slice(1:100)
    write.table(top_cluster, "temp.R.1.txt", quote = F, sep = "\t")
    top_cluster <- read.table('temp.R.1.txt', sep = "\t", header = T)

    message ("Performing Gene Set Enrichment Analysis (GSEA) using gene sets from reference and top 100 markers per cluster in query")
    iter = length(unique(top_cluster$cluster))-1
    gsea.res.return <- c()
    for(f in 0:iter)
    {   
        data <- c()
        cluster.data <- c()
        cluster.data.use <- c()
        ranks <- c()
        gsea.res <- c()
        cluster_info <- c()
        cluster_number <- c()
        cluster_number <- paste0("seurat_cluster:", f)
        data <- top_cluster
        data %>% dplyr::filter(cluster == f) %>% arrange(desc(avg_log2FC), desc(p_val_adj)) %>% head(n = 100)
        cluster.data <- data %>% dplyr::filter(cluster == f) %>% arrange(desc(p_val_adj)) %>% dplyr::select(gene, avg_log2FC)
        cluster.data.use <- data.frame(cluster.data$gene, cluster.data$avg_log2FC)
        ranks <- deframe(cluster.data.use)
        gsea.res <- fgsea(fgsea_sets, ranks, minSize=10, maxSize = 1000, nperm = 500)
        gsea.res <- gsea.res %>% arrange(desc(abs(NES))) %>% top_n(10, -padj)
        gsea.res <- data.frame(gsea.res)
        cluster_info <- rep(cluster_number, length(rownames(gsea.res)))
        gsea.res <- data.frame(cluster_info, gsea.res)
        colnames(gsea.res) <- c("Cluster", "Celltypes", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
        gsea.res.return <- rbind(gsea.res.return, gsea.res)
        
    }
    gsea.res.return <- data.frame(gsea.res.return)
    gsea.res.return <- apply(gsea.res.return,2,as.character)
    unlink("temp.R.txt")
    unlink("temp.R.1.txt")
    return(gsea.res.return)
}


randomForest_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL, ntree = NULL) {
    rf_Celltypes <- randomForest(factor(train_label) ~ ., data=train, ntree = ntree)
    return(rf_Celltypes)
}

svm_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL, crossvalidationSVM = NULL, cachesize = NULL, cost = NULL, kernel = "linear") {
    svm_Celltypes <- svm(factor(train_label) ~ ., data=train, scale = FALSE, cross = crossvalidationSVM, cachesize = cachesize, cost = cost, kernel = "linear")
    return(svm_Celltypes)
}