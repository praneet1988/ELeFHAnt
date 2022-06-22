control_parsnip(verbosity=2)
tidymodels_prefer()

#' ELeFHAnt mouse tisssues
#' mouse tissues
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'
mouse_tissues = c("Bone marrow", "Taste bud", "Adipose tissue", "Undefined", "Meniscus", "Peyer patch", "Thyroid", "Lung", "Kidney", "Heart", "Skin", "Epidermis", "Spleen", "Brain", "Blood", "Small intestine", "Serum", "Peripheral blood", "Retina", "Cerebellum", "Esophagus", "Prostate", "Dermis", "Embryoid body", "Pancreas", "Ovary", "Inner nuclear layer of retina", "Ganglion cell layer of retina", "Umbilical cord", "White adipose tissue", "Lymph node", "Colon epithelium", "Cochlea", "Basilar membrane", "Hippocampus", "Colon", "Spinal cord", "Intestinal crypt", "Peritoneal cavity", "Submandibular gland", "Muscle", "Yolk sac", "Embryo", "Epithelium", "Umbilical cord blood", "Testis", "Skeletal muscle", "Heart muscle", "Ileum", "Inner Ear", "Mammary gland", "Embryonic stem cell", "Bone", "Thymus", "Lacrimal gland", "Liver", "Mesonephros", "Bronchiole", "Gonad", "Neural tube", "Hair follicle", "Stomach", "Intestine", "Mesenteric lymph node", "Carotid artery", "Lymphoid tissue", "Pancreatic islet", "Blood vessel", "Uterus", "Gastrointestinal tract", "Colorectum", "Artery", "Aorta", "Eye", "Corneal epithelium", "Fetal liver", "Breast", "Mammary epithelium", "Bladder")

#' ELeFHAnt human tisssues
#' human tissues
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'
human_tissues = c("Kidney", "Liver", "Endometrium", "Germ", "Corneal epithelium", "Placenta", "Periosteum", "Amniotic membrane", "Primitive streak", "Adipose tissue", "Scalp", "Heart", "Chorionic villus", "Brain", "Amniotic fluid", "Cartilage", "Bone", "Retina", "Peripheral blood", "Blood", "Tongue", "Gingiva", "Undefined", "Uterus", "Dental pulp", "Breast", "Stomach", "Nucleus pulposus", "Urine", "Pancreas", "Colon", "Gonad", "Testis", "Bone marrow", "Lung", "Embryo", "Skin", "Thymus", "Blood vessel", "Synovial fluid", "Periodontal ligament", "Lymph node", "Colorectum", "Pancreatic islet", "Synovium", "Oral mucosa", "Skeletal muscle", "Eye", "Esophagus", "Umbilical cord blood", "Hair follicle", "Ovarian follicle", "Limbal epithelium", "Oral cavity", "Nasal concha", "Oviduct", "Cornea", "Sclerocorneal tissue", "Intervertebral disc", "Plasma", "Embryoid body", "Inferior colliculus", "Ovary", "Ovarian cortex", "Subcutaneous adipose tissue", "Prostate", "Myocardium", "Deciduous tooth", "Salivary gland", "Tonsil", "Gastric gland", "Vagina", "Muscle", "Abdominal adipose tissue", "Venous blood", "Parotid gland", "Antecubital vein", "Intestinal crypt", "Endometrium stroma", "Foreskin", "Submandibular gland", "Retinal pigment epithelium", "Lymph", "Gastrointestinal tract", "Lymphoid tissue", "Renal glomerulus", "Pyloric gland", "Bladder", "Sputum", "Decidua", "Umbilical cord", "Sympathetic ganglion", "Alveolus", "Seminal plasma", "Serum", "Duodenum", "Rectum", "Fetal pancreas", "Embryonic stem cell")


#' ELeFHAnt Celltype Annotation
#'
#' Celltype annotation is a function to annotate celltypes in a single cell datasets.
#' It requires a reference dataset (a processed Seurat Object with Celltypes column in metadata) and a query dataset (a processed 
#' seurat object with seurat_clusters column in metadata). One can choose from randomForest, SVM or Ensemble classifiction method
#' to learn celltypes from reference dataset and then predict celltypes for query dataset.
#' 
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference, enabling fast computation. if classification.approach is set to "ClassifyCells_usingApproximation" query will be downsampled along with reference.
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 500]
#' @param classification.method choose classification method for learning and predicting celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses estimation robustness of both randomForest and SVM to predict)(
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model (Default: 5)
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis and assessing markers from CellMarker Database. 
#' @param selectvarfeatures number of variable features to select while training (default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param cost cost of constraints violation (Default: 1)
#' @param classification.approach apprach to classify cells 1) ClassifyCells 2) ClassifyCells_usingApproximation. Default: ClassifyCells. We recommend using ClassifyCells_usingApproximation
#' when reference has significantly less number of cells compared to query
#' @param species human or mouse if validatePredictions = TRUE
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @return query seurat object with predictions added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

CelltypeAnnotation <- function(reference = NULL, query = NULL, downsample = FALSE, downsample_to = 500, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 5, validatePredictions = TRUE, selectvarfeatures = 2000, ntree = 500, cost = 1, classification.approach = "ClassifyCells", species = NULL, tissue = NULL) {
    if(downsample == TRUE)
    {
        message ("Setting Assay of reference and query to RNA")
        DefaultAssay(reference) <- "RNA"
        DefaultAssay(query) <- "RNA"

        
        reference$Dataset <- rep("reference", ncol(reference))
        query$Dataset <- rep("query", ncol(query))

        Idents(reference) <- reference$Celltypes
        Idents(query) <- query$seurat_clusters

        message("Running Diagonistis on reference and query")
        reference_cells <- nrow(reference@meta.data)
        query_cells <- nrow(query@meta.data)

        message (paste0("Number of cells in reference:", reference_cells))
        message (paste0("Number of cells in query:", query_cells))

        message ("Downsampling reference")
        
        reference_use <- subset(reference, downsample = downsample_to)
        reference_cells_afterdownsampling <- nrow(reference_use@meta.data)
        message (paste0("Number of cells in reference after downsampling per celltype:", reference_cells_afterdownsampling))

        message ("Calculating ratio of number of cells in downsampled reference vs query")
        Ratio = query_cells / reference_cells_afterdownsampling
        message (paste0("Ratio of number of cells in query vs downsampled reference:", Ratio))

        if(classification.approach == "ClassifyCells")
        {
            query = ClassifyCells(reference = reference_use, query = query, downsample = downsample, downsample_to = downsample_to, classification.method = classification.method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost, species = species, tissue = tissue)
            return(query)
        }
        if(classification.approach == "ClassifyCells_usingApproximation")
        {
            query = ApproximationBasedCelltypeAssignment(reference = reference_use, query = query, downsample = downsample, downsample_to = downsample_to, classification.method = classification.method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost, species = species, tissue = tissue)
            return(query)
        }
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

        message("Running Diagonistis on reference and query")
        reference_cells <- nrow(reference@meta.data)
        query_cells <- nrow(query@meta.data)

        message (paste0("Number of cells in reference:", reference_cells))
        message (paste0("Number of cells in query:", query_cells))

        message ("Calculating ratio of number of cells in reference vs query")
        Ratio = query_cells / reference_cells
        message (paste0("Ratio of number of cells in query vs reference:", Ratio))

        if(classification.approach == "ClassifyCells")
        {
            query = ClassifyCells(reference = reference, query = query, downsample = downsample, downsample_to = downsample_to, classification.method = classification.method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost, species = species, tissue = tissue)
            return(query)
        }
        if(classification.approach == "ClassifyCells_usingApproximation")
        {
            query = ApproximationBasedCelltypeAssignment(reference = reference, query = query, downsample = downsample, downsample_to = downsample_to, classification.method = classification.method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost, species = species, tissue = tissue)
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
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @param seurat.objects a list of processed seurat objects (please set Default Assay to "RNA") with Celltypes column in their respective meta.data to perform integration on
#' @param perform_integration logical Indicator (TRUE or FALSE) to perform integration using list of seurat.objects
#' @param integrated.atlas an integrated seurat object with CellTypes and seurat_clusters column in meta.data. Required if perform_integration = FALSE
#' @param downsample logical Indicator (TRUE or FALSE) to downsample Seurat objects or integrated seurat object, enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 500]
#' @param npcs number of principal components to compute after integration [Default: 30]
#' @param resolution value of the resolution parameter, decides size of cell communities. [Default: 0.8]
#' @param classification.method choose classification method for learning and harmonizing celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses learning robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model [Default: 5]
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis
#' @param selectanchorfeatures number of anchor features to use for integrating datasets (Default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param cost cost of constraints violation (Default: 1)
#' @param species human or mouse if validatePredictions = TRUE
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @param k.anchor  How many neighbors (k) to use when picking anchors [Default: 5]
#' @param k.filter  How many neighbors (k) to use when filtering anchors [Default: 200]
#' @param k.score  How many neighbors (k) to use when scoring anchors [Default: 30]
#' @param dims  Which dimensions to use from the CCA to specify the neighbor search space [Default: 1:30]
#' @return integrated seurat object with harmonized celltypes added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

LabelHarmonization <- function(seurat.objects = c(), perform_integration = TRUE, integrated.atlas = NULL, downsample = TRUE, downsample_to = 500, npcs = 30, resolution = 0.8, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 5, validatePredictions = TRUE, selectanchorfeatures = 2000, ntree = 500, cost = 1, k.anchor = 5, k.filter = 200, k.score = 30, dims = 1:30, species = NULL, tissue = NULL) {
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
        integration.anchors <- FindIntegrationAnchors(object.list = seurat.objects, anchor.features = selectanchorfeatures, k.anchor = k.anchor, k.filter = k.filter, k.score = k.score, dims = dims)
        integrated.atlas <- IntegrateData(anchorset = integration.anchors)
        DefaultAssay(integrated.atlas) <- "integrated"
        message ("Integration Completed. Performing Scaling, Dimension reduction and clustering")
        integrated.atlas <- ScaleData(integrated.atlas, verbose = FALSE)
        integrated.atlas <- RunPCA(integrated.atlas, npcs = npcs, verbose = FALSE)
        integrated.atlas <- RunUMAP(integrated.atlas, reduction = "pca", dims = 1:npcs)
        integrated.atlas <- FindNeighbors(integrated.atlas, reduction = "pca", dims = 1:npcs)
        integrated.atlas <- FindClusters(integrated.atlas, resolution = resolution)
        integrated.use <- integrated.atlas
        
        num_cells <- nrow(integrated.use@meta.data)
        message (paste0("Number of cells in integrated atlas:", num_cells))
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
            num_cells <- nrow(integrated.use@meta.data)
            message (paste0("Number of cells in integrated atlas:", num_cells))
        }
        if(downsample == FALSE)
        {
            message ("Setting Assay of integrated.atlas to integrated")
            integrated.atlas <- integrated.atlas
            DefaultAssay(integrated.atlas) <- "integrated"
            integrated.use <- integrated.atlas
            num_cells <- nrow(integrated.use@meta.data)
            message (paste0("Number of cells in integrated atlas:", num_cells))
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
        message ("Training & Classifying using randomForest classifier")
        rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
        message ("Predicting using trained randomForest classifier")
        rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
        rf_cm = table(test_label, rf_pred)
    

        message ("Calculating weights for randomForest classifier")
        rf_accuracy_estimate <- rf[[2]]*100
        message (paste0("Accuracy estimate of randomForest classifier:", rf_accuracy_estimate))

        message ("Assigning weights to randomForest predictions")
        rf_cm <- as.matrix(rf_cm) * rf_accuracy_estimate

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("HarmonizedLabels_UsingRF", "seurat_clusters")
        HarmonizedLabels_UsingRF <- as.character(rf_celltype_pred[match(integrated.atlas$seurat_clusters, rf_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingRF"])
        integrated.atlas[["HarmonizedLabels_UsingRF"]] <- HarmonizedLabels_UsingRF
        message ("Added Harmonized Labels using randomForest to integrated object")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message("randomForest based learning and harmonization completed. Starting validation of celltype assignments using GSEA")
            validation = ValidatePredictions(species = species, query = integrated.atlas, tissue = tissue)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
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
        message ("Training & Classification using SVM classifier")
        svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
        message ("Predicting using trained SVM classifier")
        svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
        svm_cm = table(test_label, svm_pred)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate <- svm[[2]]*100
        message (paste0("Accuracy estimate of SVM classifier:", svm_accuracy_estimate))

        message ("Assigning weights to SVM predictions")
        svm_cm <- as.matrix(svm_cm) * svm_accuracy_estimate

        svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
        colnames(svm_celltype_pred) <- c("HarmonizedLabels_UsingSVM", "seurat_clusters")
        HarmonizedLabels_UsingSVM <- as.character(svm_celltype_pred[match(integrated.atlas$seurat_clusters, svm_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingSVM"])
        integrated.atlas[["HarmonizedLabels_UsingSVM"]] <- HarmonizedLabels_UsingSVM
        message ("Added harmonized labels using SVM to integrated object")
        write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message ("SVM based learning and harmonization completed. Starting validation of celltype assignments using GSEA")
            validation = ValidatePredictions(species = species, query = integrated.atlas, tissue = tissue)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
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
        message ("Training & Classifying using randomForest classifier")
        rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
        message ("Predicting using trained randomForest classifier")
        rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
        rf_cm = table(test_label, rf_pred)
    

        message ("Calculating weights for randomForest classifier")
        rf_accuracy_estimate <- rf[[2]]*100
        message (paste0("Accuracy estimate of randomForest classifier:", rf_accuracy_estimate))

        message ("Assigning weights to randomForest predictions")
        rf_cm <- as.matrix(rf_cm) * rf_accuracy_estimate

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("HarmonizedLabels_UsingRF", "seurat_clusters")
        HarmonizedLabels_UsingRF <- as.character(rf_celltype_pred[match(integrated.atlas$seurat_clusters, rf_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingRF"])
        integrated.atlas[["HarmonizedLabels_UsingRF"]] <- HarmonizedLabels_UsingRF
        message ("Added Harmonized Labels using randomForest to integrated object")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        
        message ("Setting up SVM classifier learning")
        message ("Training & Classification using SVM classifier")
        svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
        message ("Predicting using trained SVM classifier")
        svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
        svm_cm = table(test_label, svm_pred)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate <- svm[[2]]*100
        message (paste0("Accuracy estimate of SVM classifier:", svm_accuracy_estimate))

        message ("Assigning weights to SVM predictions")
        svm_cm <- as.matrix(svm_cm) * svm_accuracy_estimate

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
            validation = ValidatePredictions(species = species, query = integrated.atlas, tissue = tissue)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
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
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @param reference1 a processed Seurat Object with Celltypes column in metadata
#' @param reference2 a processed seurat object with Celltypes column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference1 and reference2, enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 500]
#' @param classification.method choose classification method for learning and predicting celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses estimation robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model (Default: 5)
#' @param selectvarfeatures number of variable features to select while training (Default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param cost cost of constraints violation (Default: 1)
#' @return ggplot2 heatmap object and heatmap is automatically saved
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

DeduceRelationship <- function(reference1 = NULL, reference2 = NULL, downsample = FALSE, downsample_to = 500, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 5, selectvarfeatures = 2000, ntree = 500, cost = 1) {
  if(downsample == TRUE)
  {
    message ("Setting Assay of reference1 and reference2 to RNA")
    DefaultAssay(reference1) <- "RNA"
    DefaultAssay(reference2) <- "RNA"
    
    message ("Downsampling reference1 and reference2")
    reference1$Dataset <- rep("reference1", ncol(reference1))
    reference2$Dataset <- rep("reference2", ncol(reference2))

    num_cells_reference1 <- nrow(reference1@meta.data)
    message (paste0("Number of cells in reference1:", num_cells_reference1))

    num_cells_reference2 <- nrow(reference2@meta.data)
    message (paste0("Number of cells in reference2:", num_cells_reference2))
    
    Idents(reference1) <- reference1$Celltypes
    Idents(reference2) <- reference2$Celltypes
    
    reference1_use <- subset(reference1, downsample = downsample_to)
    reference2_use <- subset(reference2, downsample = downsample_to)

    num_cells_reference1_use <- nrow(reference1_use@meta.data)
    message (paste0("Number of cells in reference1 after downsampling:", num_cells_reference1_use))

    num_cells_reference2_use <- nrow(reference2_use@meta.data)
    message (paste0("Number of cells in reference2 after downsampling:", num_cells_reference2_use))
  }
  
  if(downsample == FALSE)
  {
    message ("Setting Assay of reference1 and reference2 to RNA")
    DefaultAssay(reference1) <- "RNA"
    DefaultAssay(reference2) <- "RNA"
    
    reference1$Dataset <- rep("reference1", ncol(reference1))
    reference2$Dataset <- rep("reference2", ncol(reference2))

    num_cells_reference1 <- nrow(reference1@meta.data)
    message (paste0("Number of cells in reference1:", num_cells_reference1))

    num_cells_reference2 <- nrow(reference2@meta.data)
    message (paste0("Number of cells in reference2:", num_cells_reference2))
    
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
  train = X$reference1[,1:num_features]
  test = X$reference2[,1:num_features]
  train_label = X$reference1$Celltypes
  test_label = X$reference2$Celltypes

  if(classification.method == "randomForest")
  {
    message ("Setting up randomForest classifier learning.")
    message ("Training & Classifying using randomForest classifier")
    rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
    message ("Predicting using trained randomForest classifier")
    rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
    rf_cm = table(test_label, rf_pred)
    
    message ("Calculating weight for randomForest classifier")
    rf_accuracy_estimate <- rf[[2]]*100
    message (paste0("Accuracy estimate of randomForest classifier:", rf_accuracy_estimate))
    
    message ("Assigning weights to randomForest predictions")
    rf_cm <- as.matrix(rf_cm) * rf_accuracy_estimate
    
    message ("Generating confusion matrix and heatmap")
    write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
    row_order <- hclust(dist(rf_cm))$order
    col_order <- hclust(dist(t(rf_cm)))$order
    rf_cm <- rf_cm[match(rownames(rf_cm)[row_order], rownames(rf_cm)),match(colnames(rf_cm)[col_order], colnames(rf_cm))]
    rf_cm_norm <- round(rf_cm/apply(rf_cm,1,max),3)
    rf_df <- as.data.frame(rf_cm_norm)
    colnames(rf_df) <- c("reference2","reference1","Relative_Similarity")
    plot = ggplot(data = rf_df, aes(x=reference2, y=reference1, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_RandomForest.png", width = 10, height = 10, dpi = 800)
    message("randomForest based learning and relationship inference completed")
    return(plot)
  }
  
  if(classification.method == "SVM")
  {
    message ("Setting up SVM classifier learning.")
    message ("Training & Classification using SVM classifier")
    svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
    message ("Predicting using trained SVM classifier")
    svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
    svm_cm = table(test_label, svm_pred)
    
    message ("Calculating weight for SVM classifier")
    svm_accuracy_estimate <- svm[[2]]*100
    message (paste0("Accuracy estimate of SVM classifier:", svm_accuracy_estimate))
    
    message ("Assigning weights to SVM predictions")
    svm_cm <- as.matrix(svm_cm) * svm_accuracy_estimate
    
    message ("Generating confusion matrix and heatmap")
    write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
    row_order <- hclust(dist(svm_cm))$order
    col_order <- hclust(dist(t(svm_cm)))$order
    svm_cm <- svm_cm[match(rownames(svm_cm)[row_order], rownames(svm_cm)),match(colnames(svm_cm)[col_order], colnames(svm_cm))]
    svm_cm_norm <- round(svm_cm/apply(svm_cm,1,max),3)
    svm_df <- as.data.frame(svm_cm_norm)
    colnames(svm_df) <- c("reference2","reference1","Relative_Similarity")
    plot = ggplot(data = svm_df, aes(x=reference2, y=reference1, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_SVM.png", width = 10, height = 10, dpi = 800)
    message ("SVM based learning and relationship inference completed")
    return(plot)
  }
  
  if(classification.method == "Ensemble")
  {
    message ("Ensemble learning using classification accuracy of both Random Forest and SVM classifiers")
    message ("Setting up randomForest classifier learning.")
    message ("Training & Classifying using randomForest classifier")
    rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
    message ("Predicting using trained randomForest classifier")
    rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
    rf_cm = table(test_label, rf_pred)
    
    message ("Calculating weight for randomForest classifier")
    rf_accuracy_estimate <- rf[[2]]*100
    message (paste0("Accuracy estimate of randomForest classifier:", rf_accuracy_estimate))
    
    message ("Assigning weights to randomForest predictions")
    rf_cm <- as.matrix(rf_cm) * rf_accuracy_estimate
    
    message ("Generating confusion matrix and heatmap")
    write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
    row_order <- hclust(dist(rf_cm))$order
    col_order <- hclust(dist(t(rf_cm)))$order
    rf_cm <- rf_cm[match(rownames(rf_cm)[row_order], rownames(rf_cm)),match(colnames(rf_cm)[col_order], colnames(rf_cm))]
    rf_cm_norm <- round(rf_cm/apply(rf_cm,1,max),3)
    rf_df <- as.data.frame(rf_cm_norm)
    colnames(rf_df) <- c("reference2","reference1","Relative_Similarity")
    plot = ggplot(data = rf_df, aes(x=reference2, y=reference1, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_RandomForest.png", width = 10, height = 10, dpi = 800)
    
    message ("Setting up SVM classifier learning.")
    message ("Training & Classification using SVM classifier")
    svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
    message ("Predicting using trained SVM classifier")
    svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
    svm_cm = table(test_label, svm_pred)
    
    message ("Calculating weight for SVM classifier")
    svm_accuracy_estimate <- svm[[2]]*100
    message (paste0("Accuracy estimate of SVM classifier:", svm_accuracy_estimate))
    
    message ("Assigning weights to SVM predictions")
    svm_cm <- as.matrix(svm_cm) * svm_accuracy_estimate
    
    message ("Generating confusion matrix and heatmap")
    write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
    row_order <- hclust(dist(svm_cm))$order
    col_order <- hclust(dist(t(svm_cm)))$order
    svm_cm <- svm_cm[match(rownames(svm_cm)[row_order], rownames(svm_cm)),match(colnames(svm_cm)[col_order], colnames(svm_cm))]
    svm_cm_norm <- round(svm_cm/apply(svm_cm,1,max),3)
    svm_df <- as.data.frame(svm_cm_norm)
    colnames(svm_df) <- c("reference2","reference1","Relative_Similarity")
    plot = ggplot(data = svm_df, aes(x=reference2, y=reference1, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_SVM.png", width = 10, height = 10, dpi = 800)
    
    message ("randomForest and SVM based learning and relationship inference completed. Using predictions from all models to make Ensemble Predictions")
    
    message ("Generating confusion matrix and heatmap")
    consensus_cm = rf_cm/max(rf_cm) + svm_cm/max(svm_cm)
    write.table(consensus_cm, "ConfusionMatrix_EnsembleLearning.txt", quote=F, sep="\t")
    row_order <- hclust(dist(consensus_cm))$order
    col_order <- hclust(dist(t(consensus_cm)))$order
    consensus_cm <- consensus_cm[match(rownames(consensus_cm)[row_order], rownames(consensus_cm)),match(colnames(consensus_cm)[col_order], colnames(consensus_cm))]
    consensus_cm_norm <- round(consensus_cm/apply(consensus_cm,1,max),3)
    consensus_df <- as.data.frame(consensus_cm_norm)
    colnames(consensus_df) <- c("reference2","reference1","Relative_Similarity")
    plot = ggplot(data = consensus_df, aes(x=reference2, y=reference1, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
    plot
    ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_Ensemble.png", width = 10, height = 10, dpi = 800)
    message("Ensemble based learning and relationship inference completed")
    return(plot)
  }
}


#' ELeFHAnt Approximation based Celltype Annotation
#'
#' Approximation based Celltype Annotation is a function to annotate cells in a single cell datasets using approximation. Approximation sets highest voted cell type for a query cluster.
#' It requires a reference dataset (a processed Seurat Object with Celltypes column in metadata) and a query dataset (a processed 
#' seurat object with seurat_clusters column in metadata). One can choose from randomForest, SVM or Ensemble classifiction method
#' to learn celltypes from reference dataset and then predict celltypes for query dataset.
#' 
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference and query enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 500]
#' @param classification.method choose classification method for learning and predicting celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses estimation robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model (Default: 5)
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis
#' @param selectvarfeatures number of variable features to select while training (Default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param cost cost of constraints violation (Default: 1)
#' @param species human or mouse if validatePredictions = TRUE
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @return query seurat object with predictions added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

ApproximationBasedCelltypeAssignment <- function(reference = reference, query = query, downsample = downsample, downsample_to = downsample_to, classification.method = classification.method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, species = NULL, tissue = NULL, cost = cost) {
    message ("Downsampling query")
    query_use = subset(query, downsample = downsample_to)

    reference_use = reference

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
    train = X$reference[,1:num_features]
    test = X$query[,1:num_features]
    train_label = X$reference$Celltypes
    test_label = X$query$Clusters

    if(classification.method == "randomForest")
    {
        message ("Setting up randomForest classifier learning")
        message ("Training & Classifying using randomForest classifier")
        rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
        rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
        rf_cm = table(test_label, rf_pred)

        message ("Calculating weights for each randomForest classifier")
        rf_accuracy_estimate <- rf[[2]]*100
        message (paste0("Accuracy estimate of randomForest classifier:", rf_accuracy_estimate))

        message ("Assigning weights to randomForest predictions")
        rf_cm <- as.matrix(rf_cm) * rf_accuracy_estimate

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("PredictedCelltype_UsingRF", "seurat_clusters")
        PredictedCelltype_UsingRF <- as.character(rf_celltype_pred[match(query$seurat_clusters, rf_celltype_pred$seurat_clusters), "PredictedCelltype_UsingRF"])
        query[["PredictedCelltype_UsingRF"]] <- PredictedCelltype_UsingRF
        message ("Added Predicted celltypes using randomForest to query")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message("randomForest based learning and celltype annotation completed. Starting validation of celltype assignments using GSEA")
            validation = ValidatePredictions(species = species, query = query_use, tissue = tissue, reference = reference_use)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
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
        message ("Training & Classifying using SVM classifier")
        svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
        svm_cm = table(test_label, svm_pred)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate <- svm[[2]]*100
        message (paste0("Accuracy estimate of SVM classifier:", svm_accuracy_estimate))

        message ("Assigning weights to SVM predictions")
        svm_cm <- as.matrix(svm_cm) * svm_accuracy_estimate

        svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
        colnames(svm_celltype_pred) <- c("PredictedCelltype_UsingSVM", "seurat_clusters")
        PredictedCelltype_UsingSVM <- as.character(svm_celltype_pred[match(query$seurat_clusters, svm_celltype_pred$seurat_clusters), "PredictedCelltype_UsingSVM"])
        query[["PredictedCelltype_UsingSVM"]] <- PredictedCelltype_UsingSVM
        message ("Added Predicted celltypes using SVM to query")
        write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
        if(validatePredictions == TRUE)
        {
            message("SVM based learning and celltype annotation completed. Starting validation of celltype assignments using GSEA")
            validation = ValidatePredictions(species = species, query = query_use, tissue = tissue, reference = reference_use)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
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
        message ("Training & Classifying using randomForest classifier")
        rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
        rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
        rf_cm = table(test_label, rf_pred)

        message ("Calculating weights for each randomForest classifier")
        rf_accuracy_estimate <- rf[[2]]*100
        message (paste0("Accuracy estimate of randomForest classifier:", rf_accuracy_estimate))

        message ("Assigning weights to randomForest predictions")
        rf_cm <- as.matrix(rf_cm) * rf_accuracy_estimate

        rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
        colnames(rf_celltype_pred) <- c("PredictedCelltype_UsingRF", "seurat_clusters")
        PredictedCelltype_UsingRF <- as.character(rf_celltype_pred[match(query$seurat_clusters, rf_celltype_pred$seurat_clusters), "PredictedCelltype_UsingRF"])
        query[["PredictedCelltype_UsingRF"]] <- PredictedCelltype_UsingRF
        message ("Added Predicted celltypes using randomForest to query")
        write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
        
        message ("Setting up SVM classifier learning")
        message ("Training & Classifying using SVM classifier")
        svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
        svm_cm = table(test_label, svm_pred)

        message ("Calculating weights for each SVM classifier")
        svm_accuracy_estimate <- svm[[2]]*100
        message (paste0("Accuracy estimate of SVM classifier:", svm_accuracy_estimate))

        message ("Assigning weights to SVM predictions")
        svm_cm <- as.matrix(svm_cm) * svm_accuracy_estimate

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
            validation = ValidatePredictions(species = species, query = query_use, tissue = tissue, reference = reference_use)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
            return(query)
        }
        if(validatePredictions == FALSE)
        {
            message("Ensembl celltype annotation completed.")
            return(query)
        }
    }
}


#' ELeFHAnt Cell type Annotation per cell in query
#'
#' Classify Cells is a function to annotate cells in a single cell datasets.
#' It requires a reference dataset (a processed Seurat Object with Celltypes column in metadata) and a query dataset (a processed 
#' seurat object with seurat_clusters column in metadata). One can choose from randomForest, SVM or Ensemble classifiction method
#' to learn celltypes from reference dataset and then predict celltypes for query dataset.
#' 
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference, enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 500]
#' @param classification.method choose classification method for learning and predicting celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses estimation robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model (Default: 5)
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis
#' @param selectvarfeatures number of variable features to select while training (Default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param cost cost of constraints violation (Default: 1)
#' @param species human or mouse if validatePredictions = TRUE
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @return query seurat object with predictions added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

ClassifyCells <- function(reference = reference, query = query, downsample = downsample, downsample_to = downsample_to, classification.method = classification.method, crossvalidationSVM = crossvalidationSVM, validatePredictions = validatePredictions, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost, species = NULL, tissue = NULL) {
    query_use = query

    reference_use = reference

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
    train = X$reference[,1:num_features]
    test = X$query[,1:num_features]
    train_label = X$reference$Celltypes
    test_label = X$query$Clusters

    if(classification.method == "randomForest")
    {
        message ("Setting up randomForest classifier learning")
        message ("Training & Classifying using randomForest classifier")
        rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
        rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
        query$PredictedCelltype_UsingRF <- rf_pred
        message ("Added Predicted celltypes using randomForest to query")

        if(validatePredictions == TRUE)
        {
            message("randomForest based learning and celltype annotation completed. Starting validation of celltype assignments using GSEA")
            validation = ValidatePredictions(species = species, query = query_use, tissue = tissue, reference = reference_use)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
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
        message ("Training & Classifying using SVM classifier")
        svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
        query$PredictedCelltype_UsingSVM <- svm_pred
        message ("Added Predicted celltypes using SVM to query")

        if(validatePredictions == TRUE)
        {
            message("SVM based learning and celltype annotation completed. Starting validation of celltype assignments using GSEA")
            validation = ValidatePredictions(species = species, query = query_use, tissue = tissue, reference = reference_use)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
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
        message ("Training & Classifying using randomForest classifier")
        rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
        rf_pred <- predict(rf[[1]], new_data = test)$.pred_class
        query$PredictedCelltype_UsingRF <- rf_pred
        message ("Added Predicted celltypes using randomForest to query\n")
        
        message ("Setting up SVM classifier learning")
        message ("Training & Classifying using SVM classifier")
        svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM, cost = cost)
        svm_pred = predict(svm[[1]], new_data = test, type="raw")$predictions
        query$PredictedCelltype_UsingSVM <- svm_pred
        message ("Added Predicted celltypes using SVM to query\n")

        message ("randomForest and SVM based learning and predictions completed. Using predictions from RF and SVM to make Ensemble Predictions")

        message ("Calculating weights for randomForest classifier")
        rf_accuracy_estimate <- rf[[2]]*100
        message (paste0("Accuracy estimate of randomForest classifier:", rf_accuracy_estimate))

        message ("Calculating weights for SVM classifier")
        svm_accuracy_estimate <- svm[[2]]*100
        message (paste0("Accuracy estimate of SVM classifier:", svm_accuracy_estimate))

        predictions = data.frame(rf_pred, svm_pred)
        predictions$rf_accuracy = rep(rf_accuracy_estimate, nrow(predictions))
        predictions$svm_accuracy = rep(svm_accuracy_estimate, nrow(predictions))
        predictions$Match = predictions$rf_pred == predictions$svm_pred
        predictions$final_pred <- ifelse(as.character(predictions$Match) == 'TRUE', as.character(predictions$rf_pred), ifelse(predictions$rf_accuracy > predictions$svm_accuracy, as.character(predictions$rf_pred), as.character(predictions$svm_pred)))
        colnames(predictions) <- c("PredictedCelltype_UsingRF", "PredictedCelltype_UsingSVM", "RF_Accuracy", "SVM_accuracy", "Match", "PredictedCelltype_UsingEnsemble")
        query$PredictedCelltype_UsingEnsemble <- predictions$PredictedCelltype_UsingEnsemble
        message ("Added Predicted celltypes using Ensemble learning to query")
        
        if(validatePredictions == TRUE)
        {
            message("Ensembl celltype annotation completed. Starting validation of celltype assignments using GSEA")
            validation = ValidatePredictions(species = species, query = query_use, tissue = tissue, reference = reference_use)
            message ("Validation completed. Please see ValidatePredictions Folder for results")
            return(query)
        }
        if(validatePredictions == FALSE)
        {
            message("Ensembl celltype annotation completed.")
            return(query)
        }
    }
}


#' ELeFHAnt Validate Predictions
#'
#' Validate Predictions is a function to validate celltype assignments. It uses Gene Set Enrichment Analysis (GSEA) and cell type markers from CellMarker Database to let user assess celltype annotation/harmonization.
#' It requires species, tissue and query dataset. For GSEA, ELeFHAnt uses C8 Hallmark Gene Sets and from CellMarker Database it uses markers curated based on experiments.
#' 
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @param species human or mouse to select the C8 hallmark cell type gene sets
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param reference a processed seurat object with Celltypes column in metadata
#' @return GSEA result
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

ValidatePredictions <- function(species = NULL, tissue = NULL, query = NULL, reference = NULL) {

    DefaultAssay(query) = "RNA"
    message("\nSetting up Directory to write ValidatePredictions Results\n")
    current_WD = getwd()
    dir_create.temp = paste0(current_WD, "/", "ValidatePredictions/")
    dir.create(dir_create.temp)
    dir_create_GSEA = paste0(dir_create.temp, "/", "GSEA_usingC8HallmarkGeneSets/")
    dir.create(dir_create_GSEA)

    dir_create_Cellmarkers.temp = paste0(dir_create.temp, "/", "CellMarkers/")
    dir.create(dir_create_Cellmarkers.temp)
    

    message("Downloading Experiment based curated Cell type Markers from CellMarker Database and C8 Hallmark gene sets GSEA")
    
    if(species == "human")
    {
        download.file('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt', destfile = "all_cell_markers_Human.txt")
        fgsea_sets = msigdbr(species = "human", category = "C8")
        cellmarkers = read.table('all_cell_markers_Human.txt', sep="\t", header=T, check.names = F)
        cellmarkers_human = subset(cellmarkers, speciesType == "Human")
        cellmarkers_experiment = subset(cellmarkers_human, markerResource == "Experiment")
    }
    if(species == "mouse")
    {
        download.file('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt', destfile = "all_cell_markers_Mouse.txt")
        fgsea_sets = msigdbr(species = "mouse", category = "C8")
        cellmarkers = read.table('all_cell_markers_Mouse.txt', sep="\t", header=T, check.names = F)
        cellmarkers_mouse = subset(cellmarkers, speciesType == "Mouse")
        cellmarkers_experiment = subset(cellmarkers_mouse, markerResource == "Experiment")
    }

    message('\nGSEA BASED VALIDATION\n')
    msigdbr_list = split(x = fgsea_sets$gene_symbol, f = fgsea_sets$gs_name)
    message ("Obtaining markers per annotated cluster")
    Idents(query) <- query$seurat_clusters
    cluster_markers <- FindAllMarkers(query)
    cluster_markers = data.frame(cluster_markers)
    cluster_markers = subset(cluster_markers, p_val_adj <= 0.05)
    top_cluster <- cluster_markers %>% group_by(cluster) %>% dplyr::slice(1:50)
    write.table(top_cluster, "temp.R.1.txt", quote = F, sep = "\t")
    top_cluster <- read.table('temp.R.1.txt', sep = "\t", header = T)

    message ("Performing Gene Set Enrichment Analysis (GSEA) using gene sets from C8 Hallmark MsigDB and top 50 markers per cluster in query")
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
        gsea.res <- fgsea(msigdbr_list, ranks, minSize=10, maxSize = 1000, scoreType = "pos")
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
    gsea.res.return = data.frame(gsea.res.return)
    gsea.res.return.filtered = subset(gsea.res.return, padj <= 0.05)
    gsea.res.return.filtered = gsea.res.return.filtered %>% group_by(Cluster) %>% dplyr::slice(1:3)
    plot = ggplot(gsea.res.return.filtered, aes(fill=Celltypes, y=padj, x="")) + geom_bar(position="dodge", stat="identity") + ggtitle("Top 3 C8 Hallmark GSEA Celltypes Enrichment at pdj <= 0.05") + facet_wrap(~Cluster) + theme_ipsum() + xlab("") + theme_classic()
    filename_GSEA_plot = paste0(dir_create_GSEA, "/", 'Summary_GeneSetEnrichmentAnalysis.png')
    ggsave(filename_GSEA_plot, plot, width = 20, height = 20, dpi = 800)
    filename_GSEA_res = paste0(dir_create_GSEA, "/", 'Summary_GeneSetEnrichmentAnalysis.txt')
    write.table(gsea.res.return, filename_GSEA_res, quote=F, sep="\t")
    
    message('\nGSEA VALIDATION COMPLETED\n')

    message('\nCellMarker DATABASE BASED VALIDATION')

    message('\nCellMarker DATABASE BASED VALIDATION FOR QUERY')
    dir_create_Cellmarkers.temp.query = paste0(dir_create_Cellmarkers.temp, "query", "/")
    dir.create(dir_create_Cellmarkers.temp.query)
    for(i in 1:length(tissue))
    {
        message(paste0("Tissue of interest:", tissue[i]))
        dir_create_Cellmarkers = paste0(dir_create_Cellmarkers.temp.query, tissue[i], "/")
        dir.create(dir_create_Cellmarkers)
        cellmarkers_tissue = subset(cellmarkers_experiment, tissueType == tissue[i])
        celltypes = unique(cellmarkers_tissue$cellName)
        for(f in 1:length(celltypes))
        {
            tryCatch({
                message(paste0('Generating DotPlot/FeaturePlot for experimental evidence based markers for:', celltypes[f], " in ", tissue[i]))
                tempdata = subset(cellmarkers_tissue, cellName == celltypes[f])
                splitlist = strsplit(noquote(tempdata$cellMarker), ", ")
                print (unique(splitlist[[1]]))
                genes_temp = unique(splitlist[[1]])
                genes = genes_temp[!is.na(genes_temp)]
                if(length(genes) == 0)
                {
                    messsage(paste0("No Markers found for: ", celltypes[f]))
                }
                if(length(genes) > 0)
                {
                    check_status = intersect(genes, rownames(query))
                    if(length(check_status) > 0)
                    {
                        DotPlot(query, features = genes, group.by = "seurat_clusters") + RotatedAxis() + ggtitle(celltypes[f]) + theme_classic()
                        filename = paste0(dir_create_Cellmarkers, celltypes[f], " MarkerGenes DotPlot.png")
                        ggsave(filename, width = 10, height = 10, dpi = 800)
                        FeaturePlot(query, features = genes, order = T)
                        filename = paste0(dir_create_Cellmarkers, celltypes[f], " MarkerGenes FeaturePlot.png")
                        ggsave(filename, width = 10, height = 10, dpi = 800)
                    }
            
                }
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }

    }

    if(is.null(reference))
    {
        message('\nReference not provided')
        message('\nCellMarker DATABASE BASED VALIDATION COMPLETED')
    }
    else {
        message('\nCellMarker DATABASE BASED VALIDATION FOR REFERENCE')
        dir_create_Cellmarkers.temp.reference = paste0(dir_create_Cellmarkers.temp, "reference", "/")
        dir.create(dir_create_Cellmarkers.temp.reference)
        for(i in 1:length(tissue))
        {
            message(paste0("Tissue of interest:", tissue[i]))
            dir_create_Cellmarkers = paste0(dir_create_Cellmarkers.temp.reference, tissue[i], "/")
            dir.create(dir_create_Cellmarkers)
            cellmarkers_tissue = subset(cellmarkers_experiment, tissueType == tissue[i])
            celltypes = unique(cellmarkers_tissue$cellName)
            for(f in 1:length(celltypes))
            {
                tryCatch({
                    message(paste0('Generating DotPlot/FeaturePlot for experimental evidence based markers for:', celltypes[f], " in ", tissue[i]))
                    tempdata = subset(cellmarkers_tissue, cellName == celltypes[f])
                    splitlist = strsplit(noquote(tempdata$cellMarker), ", ")
                    print (unique(splitlist[[1]]))
                    genes_temp = unique(splitlist[[1]])
                    genes = genes_temp[!is.na(genes_temp)]
                    if(length(genes) == 0)
                    {
                        messsage(paste0("No Markers found for: ", celltypes[f]))
                    }
                    if(length(genes) > 0)
                    {
                        check_status = intersect(genes, rownames(reference))
                        if(length(check_status) > 0)
                        {
                            DotPlot(reference, features = genes, group.by = "Celltypes") + RotatedAxis() + ggtitle(celltypes[f]) + theme_classic()
                            filename = paste0(dir_create_Cellmarkers, celltypes[f], " MarkerGenes DotPlot.png")
                            ggsave(filename, width = 10, height = 10, dpi = 800)
                            FeaturePlot(reference, features = genes, order = T)
                            filename = paste0(dir_create_Cellmarkers, celltypes[f], " MarkerGenes FeaturePlot.png")
                            ggsave(filename, width = 10, height = 10, dpi = 800)
                        }
            
                    }
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            }

        }
        message('\nCellMarker DATABASE BASED VALIDATION COMPLETED')
    }   
}



#' Benchmark ELeFHAnt against scPred and Label Transfer
#'
#' Celltype annotation is a function to annotate celltypes in a single cell datasets.
#' Benchmark ELeFHAnt is an easy to use function if the user wants to compare the predictions from ELeFHAnt to Label Transfer and scPred
#' 
#' @import tidymodels
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import class
#' @import splitstackshape
#' @import fgsea
#' @import tibble
#' @import msigdbr
#' @import scPred
#' @import magrittr
#' @import harmony
#' @import scater
#' @import parsnip
#' @import hrbrthemes
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference, enabling fast computation.
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 500]
#' @param selectvarfeatures number of variable features to select while training (Default: 2000)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model (Default: 5)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param cost cost of constraints violation (Default: 1)
#' @return query seurat object with predictions from ELeFHAnt, Label Transfer and scPred added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

BenchmarkELeFHAnt <- function(reference = reference, query = query, downsample = TRUE, downsample_to = 500, crossvalidationSVM = 5, selectvarfeatures = 2000, ntree = 500, cost = 1) {

    if(downsample == TRUE)
    {
        message('\nDownsampling reference cells to enable fast computation\n')

        Idents(reference) = reference$Celltypes
        reference_use = subset(reference, downsample = downsample_to)
        reference_use = NormalizeData(reference_use)
        reference_use = FindVariableFeatures(reference_use)
        reference_use = ScaleData(reference_use)
    }
    
    if(downsample == FALSE)
    {
        reference_use = reference
    }

    message('\nDeploying ELeFHAnt: classification.method == Ensemble | classification.approach == ClassifyCells\n')

    out.ELeFHAnt = CelltypeAnnotation(reference = reference_use, query = query, downsample = downsample, downsample_to = downsample_to, classification.method = "Ensemble", crossvalidationSVM = crossvalidationSVM, validatePredictions = FALSE, selectvarfeatures = selectvarfeatures, ntree = ntree, cost = cost, classification.approach = "ClassifyCells", species = NULL, tissue = NULL)
    p1 = DimPlot(out.ELeFHAnt, group.by = "PredictedCelltype_UsingEnsemble", label=T, repel = T, label.size = 6) + NoLegend() + ggtitle("ELeFHAnt Predictions")

    message('\nDeploying Seurat Label Transfer\n')

    anchors <- FindTransferAnchors(reference = reference_use, query = query)
    predictions <- TransferData(anchorset = anchors, refdata = reference_use$Celltypes)
    out.seurat <- AddMetaData(object = query, metadata = predictions)
    p2 = DimPlot(out.seurat, group.by = "predicted.id", label=T, repel = T, label.size = 6) + NoLegend() + ggtitle("LabelTransfer Predictions")

    message('\nDeploying scPred\n')
    reference.scPred <- reference_use %>%
        RunPCA() %>% 
        RunUMAP(dims = 1:30)

    reference.scPred <- getFeatureSpace(reference.scPred, "Celltypes")
    reference.scPred <- trainModel(reference.scPred)
    get_probabilities(reference.scPred) %>% head()
    get_scpred(reference.scPred)
    out.scPred <- NormalizeData(query)
    out.scPred <- scPredict(out.scPred, reference.scPred)
    p3 = DimPlot(out.scPred, group.by = "scpred_prediction", label=T, repel = T, label.size = 6) + NoLegend() + ggtitle("scPred Predictions")

    plot = p1+p2+p3
    ggsave('Benchmarking_ELeFHAnt_LabelTransfer_scPred_CelltypeAnnotation_Plot.png', plot, width = 20, height = 20, dpi = 800)
    query$PredictedCelltype_UsingEnsemble = out.ELeFHAnt$PredictedCelltype_UsingEnsemble
    query$predicted.id = out.seurat$predicted.id
    query$scpred_prediction = out.scPred$scpred_prediction
    return(query)
}

randomForest_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL, ntree = NULL) {
    rf_model <- rand_forest(trees = ntree) %>% set_mode("classification") %>% set_engine("ranger")
    rf_Celltypes <- rf_model %>% fit(factor(train_label) ~ ., data=train)
    rf_acc <- 1-rf_Celltypes$fit$prediction.error
    return(list(rf_Celltypes,rf_acc))
}

svm_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL, crossvalidationSVM = NULL, cost = NULL) {
    cost_est <- heuristicC(train)
    svm_model <- svm_linear(cost = cost) %>% set_mode("classification") %>% set_engine("LiblineaR")
    svm_Celltypes <- svm_model %>% fit(factor(train_label) ~ ., data=train)
    svm_acc <- LiblineaR(data=train, target=train_label, type=1, cross = crossvalidationSVM, cost = cost)
    return(list(svm_Celltypes,svm_acc))
}