randomForest_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL) {
	message ("Training randomForest classifier")
	rf_Celltypes <- randomForest(factor(train_label) ~ ., data=train)
	return(rf_Celltypes)
}

svm_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL, crossvalidationSVM = NULL) {
	message ("Training SVM classifier")
 	svm_Celltypes <- svm(factor(train_label) ~ ., data=train, scale = FALSE, cross = crossvalidationSVM)
	return(svm_Celltypes)
}

ValidatePredictions <- function(reference = NULL, query = NULL, celltype_assignments = NULL) {
	Idents(reference) <- reference$Celltypes
	message ("Obtaining markers per celltype")

	celltype_markers <- FindAllMarkers(reference)
	top_celltype <- celltype_markers %>% group_by(cluster) %>% slice(1:100)

	message ("Obtaining markers per annotated cluster")

	Idents(query) <- query$seurat_clusters
	cluster_markers <- FindAllMarkers(query)
	top_cluster <- cluster_markers %>% group_by(cluster) %>% slice(1:100)

	colnames(celltype_assignments) <- c("label", "cluster")

	message ("Generating ratio of number of markers shared. Please see: Using top 100 markers for comparison")
	matches <- c()
	for (x in 1:nrow(celltype_assignments)) {
  		cluster <- celltype_assignments$cluster[x]
  		cell_type <- celltype_assignments$label[x]
  		subset_cluster <- top_cluster[top_cluster$cluster == cluster,]
  		subset_celltype <- top_celltype[top_celltype$cluster == cell_type,]
  		matches <- c(matches,(sum(subset_cluster$gene %in% subset_celltype$gene)/100))
	}
	celltype_assignments$SharedMarkers_Ratio <- matches
	return(celltype_assignments)
}

#' ELeFHAnt Celltype Annotation
#'
#' Celltype annotation is a function to annotate celltypes in a single cell datasets.
#' It requires a reference dataset (a processed Seurat Object with Celltypes column in metadata) and a query dataset (a processed 
#' seurat object with seurat_clusters column in metadata). One can choose from randomForest, SVM or Ensemble classifiction method
#' to learn celltypes from reference dataset and then predict celltypes for query dataset.
#' 
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference and query enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 100] in reference and query for Celltypes and seurat_clusters resspectively
#' @param classification.method choose classification method for learning and predicting celltypes. 
#' Choices: randomForest (decision trees), SVM (Support Vector Machines) or Ensemble (uses estimation robustness of both randomForest and SVM to predict)
#' @param crossvalidationSVM if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by computing number of markers shared between assigned celltype and annotated cluster
#' @return query seurat object with predictions added to meta.data of the object
#' @export
CelltypeAnnotation <- function(reference = NULL, query = NULL, downsample = FALSE, downsample_to = 100, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 10, validatePredictions = TRUE) {
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
	combined <- FindVariableFeatures(combined)
	combined <- ScaleData(combined)
	combined_exp <- combined[['RNA']]@scale.data
	combined_exp <- t(combined_exp)
	combined_exp <- data.frame(combined_exp)
	combined_exp$Dataset <- combined$Dataset
	message ("Generating train and test sets")
	X <- split(combined_exp, combined_exp$Dataset)
	X$reference$Celltypes <- reference_use$Celltypes
	X$query$Clusters <- query_use$seurat_clusters

	if(classification.method == "randomForest")
	{
		rf_Celltypes = randomForest_predictor(train = X$reference[,1:2000], test = X$query[,1:2000], train_label = X$reference$Celltypes, test_label = X$query$Clusters)
		message ("Predicting using trained randomForest classifier")
		rf_pred = predict(rf_Celltypes, newdata=X$query[,1:2000])
		rf_cm = table(X$query$Clusters, rf_pred)
		rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
		colnames(rf_celltype_pred) <- c("PredictedCelltype_UsingRF", "seurat_clusters")
 		PredictedCelltype_UsingRF <- as.character(rf_celltype_pred[match(query$seurat_clusters, rf_celltype_pred$seurat_clusters), "PredictedCelltype_UsingRF"])
 		query[["PredictedCelltype_UsingRF"]] <- PredictedCelltype_UsingRF
 		message ("Added Predicted celltypes using randomForest to query")
 		write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
 		if(validatePredictions == TRUE)
 		{
 			message("randomForest based learning and celltype annotation completed. Starting validation of celltype assignments")
 			reference.validation.use <- subset(reference_use, idents = rf_celltype_pred$PredictedCelltype_UsingRF)
 			validation = ValidatePredictions(reference = reference.validation.use, query = query_use, celltype_assignments = rf_celltype_pred)
 			message ("Validation completed. Please see summary of number of shared markers below")
 			print (validation)
 			write.table(validation, "Summary_numberofsharedmarkers_between_AssignedCelltypesandQueryClusters.txt", quote=F, sep="\t")
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
 		svm_Celltypes = svm_predictor(train = X$reference[,1:2000], test = X$query[,1:2000], train_label = X$reference$Celltypes, test_label = X$query$Clusters, crossvalidationSVM = crossvalidationSVM)
 		message ("Predicting using trained SVM classifier")
		svm_pred = predict(svm_Celltypes, newdata=X$query[,1:2000])
		svm_cm = table(X$query$Clusters, svm_pred)
		svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
		colnames(svm_celltype_pred) <- c("PredictedCelltype_UsingSVM", "seurat_clusters")
 		PredictedCelltype_UsingSVM <- as.character(svm_celltype_pred[match(query$seurat_clusters, svm_celltype_pred$seurat_clusters), "PredictedCelltype_UsingSVM"])
 		query[["PredictedCelltype_UsingSVM"]] <- PredictedCelltype_UsingSVM
 		message ("Added Predicted celltypes using SVM to query")
 		write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
 		if(validatePredictions == TRUE)
 		{
 			message("SVM based learning and celltype annotation completed. Starting validation of celltype assignments")
 			reference.validation.use <- subset(reference_use, idents = svm_celltype_pred$PredictedCelltype_UsingSVM)
 			validation = ValidatePredictions(reference = reference.validation.use, query = query_use, celltype_assignments = svm_celltype_pred)
 			message ("Validation completed. Please see summary of number of shared markers below")
 			print (validation)
 			write.table(validation, "Summary_numberofsharedmarkers_between_AssignedCelltypesandQueryClusters.txt", quote=F, sep="\t")
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
 		message ("Ensemble learning using classification accuracy of both Random Forest and SVM")

 		rf_Celltypes = randomForest_predictor(train = X$reference[,1:2000], test = X$query[,1:2000], train_label = X$reference$Celltypes, test_label = X$query$Clusters)
		message ("Predicting using trained randomForest classifier")
		rf_pred = predict(rf_Celltypes, newdata=X$query[,1:2000])
		rf_cm = table(X$query$Clusters, rf_pred)
		rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
		colnames(rf_celltype_pred) <- c("PredictedCelltype_UsingRF", "seurat_clusters")
 		PredictedCelltype_UsingRF <- as.character(rf_celltype_pred[match(query$seurat_clusters, rf_celltype_pred$seurat_clusters), "PredictedCelltype_UsingRF"])
 		query[["PredictedCelltype_UsingRF"]] <- PredictedCelltype_UsingRF
 		message ("Added Predicted celltypes using randomForest to query")
 		write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
 		
 		svm_Celltypes = svm_predictor(train = X$reference[,1:2000], test = X$query[,1:2000], train_label = X$reference$Celltypes, test_label = X$query$Clusters, crossvalidationSVM = crossvalidationSVM)
 		message ("Predicting using trained SVM classifier")
		svm_pred = predict(svm_Celltypes, newdata=X$query[,1:2000])
		svm_cm = table(X$query$Clusters, svm_pred)
		svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
		colnames(svm_celltype_pred) <- c("PredictedCelltype_UsingSVM", "seurat_clusters")
 		PredictedCelltype_UsingSVM <- as.character(svm_celltype_pred[match(query$seurat_clusters, svm_celltype_pred$seurat_clusters), "PredictedCelltype_UsingSVM"])
 		query[["PredictedCelltype_UsingSVM"]] <- PredictedCelltype_UsingSVM
 		message ("Added Predicted celltypes using SVM to query")
 		write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")

		rf_cm <- as.matrix(rf_cm)
 		svm_cm <- as.matrix(svm_cm)
 	
		message("Assigning weights based on Random Forest and SVM error estimation for Ensemble learning")
 		
 		rf_acccuracy_estimate <- (1-tail(rf_Celltypes$err.rate[,1], 1))*100
 		svm_accuracy_estimate <- svm_Celltypes$tot.accuracy

 		message (paste0("Accuracy estimate of randomForest based on training:", rf_acccuracy_estimate))
 		message (paste0("Accuracy estimate of SVM based on training:", svm_accuracy_estimate))

 		if((rf_acccuracy_estimate == "NaN") | (svm_accuracy_estimate == "NaN"))
 		{
 			consensus_cm = rf_cm + svm_cm
 		}
 		if((rf_acccuracy_estimate > 0) | (svm_accuracy_estimate > 0))
 		{
			rf_cm = rf_cm*rf_acccuracy_estimate
 			svm_cm = svm_cm*svm_accuracy_estimate
 			consensus_cm = rf_cm + svm_cm
 		}
 		if((rf_acccuracy_estimate == 0) | (svm_accuracy_estimate == 0))
 		{
 			consensus_cm = rf_cm + svm_cm
 		}
 		consensus_celltype_pred <- data.frame(colnames(consensus_cm)[apply(consensus_cm,1,which.max)], rownames(consensus_cm), check.names=F)
		colnames(consensus_celltype_pred) <- c("PredictedCelltype_UsingEnsemble", "seurat_clusters")
 		PredictedCelltype_UsingEnsemble <- as.character(consensus_celltype_pred[match(query$seurat_clusters, consensus_celltype_pred$seurat_clusters), "PredictedCelltype_UsingEnsemble"])
 		query[["PredictedCelltype_UsingEnsemble"]] <- PredictedCelltype_UsingEnsemble
 		message ("Added Predicted celltypes using Ensemble learning to query")
 		write.table(consensus_cm, "ConfusionMatrix_EnsembleLearning.txt", quote=F, sep="\t")
 		if(validatePredictions == TRUE)
 		{
 			message("Ensembl learning and celltype annotation completed. Starting validation of celltype assignments")
 			reference.validation.use <- subset(reference_use, idents = consensus_celltype_pred$PredictedCelltype_UsingEnsemble)
 			validation = ValidatePredictions(reference = reference.validation.use, query = query_use, celltype_assignments = consensus_celltype_pred)
 			message ("Validation completed. Please see summary of number of shared markers below")
 			print (validation)
 			write.table(validation, "Summary_numberofsharedmarkers_between_AssignedCelltypesandQueryClusters.txt", quote=F, sep="\t")
 			return(query)
 		}
 		if(validatePredictions == FALSE)
 		{
 			message("Ensembl learning and celltype annotation completed.")
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
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by computing number of markers shared between assigned celltype and annotated cluster
#' @return integrated seurat object with harmonized celltypes added to meta.data of the object
#' @export
LabelHarmonization <- function(seurat.objects = c(), perform_integration = TRUE, integrated.atlas = NULL, downsample = TRUE, downsample_to = 100, npcs = 30, resolution = 0.5, classification.method = c("randomForest", "SVM", "Ensemble"), crossvalidationSVM = 10, validatePredictions = TRUE) {
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
    				x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
			})
		}
		if(downsample == FALSE)
		{
			seurat.objects <- seurat.objects
		}
		message ("Starting integration using Seurat")
		integration.anchors <- FindIntegrationAnchors(object.list = seurat.objects)
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
	scaled_data$Celltypes <- integrated.use$Celltypes
	scaled_data$Clusters <- integrated.use$seurat_clusters
	scaled_data_stratified <- stratified(scaled_data, group = "Clusters", size = 0.6, bothSets = TRUE)
	train <- scaled_data_stratified$SAMP1[,1:2000]
	test <- scaled_data_stratified$SAMP2[,1:2000]
	train_label <- scaled_data_stratified$SAMP1$Celltypes
	test_label <- scaled_data_stratified$SAMP2$Clusters

	if(classification.method == "randomForest")
	{
		rf_Celltypes = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label)
		message ("Predicting using trained randomForest classifier")
		rf_pred = predict(rf_Celltypes, newdata=test)
		rf_cm = table(test_label, rf_pred)
		rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
		colnames(rf_celltype_pred) <- c("HarmonizedLabels_UsingRF", "seurat_clusters")
 		HarmonizedLabels_UsingRF <- as.character(rf_celltype_pred[match(integrated.atlas$seurat_clusters, rf_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingRF"])
 		integrated.atlas[["HarmonizedLabels_UsingRF"]] <- HarmonizedLabels_UsingRF
 		message ("Added Harmonized Labels using randomForest to integrated object")
 		write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
 		if(validatePredictions == TRUE)
 		{
 			message("randomForest based harmonization completed. Starting validation of celltype assignments")
 			Idents(integrated.atlas) <- integrated.atlas$Celltypes
 			integration.validation.use <- subset(integrated.atlas, idents = rf_celltype_pred$HarmonizedLabels_UsingRF)
 			validation = ValidatePredictions(reference = integration.validation.use, query = integrated.atlas, celltype_assignments = rf_celltype_pred)
 			message ("Validation completed. Please see summary of number of shared markers below")
 			print (validation)
 			write.table(validation, "Summary_numberofsharedmarkers_between_AssignedCelltypesandIntegratedClusters.txt", quote=F, sep="\t")
 			return(integrated.atlas)
 		}
 		if(validatePredictions == FALSE)
 		{
 			message ("randomForest based harmonization completed")
 			return(integrated.atlas)
 		}
 	}

 	if(classification.method == "SVM")
 	{
 		svm_Celltypes = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM)
 		message ("Predicting using trained SVM classifier")
		svm_pred = predict(svm_Celltypes, newdata=test)
		svm_cm = table(test_label, svm_pred)
		svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
		colnames(svm_celltype_pred) <- c("HarmonizedLabels_UsingSVM", "seurat_clusters")
 		HarmonizedLabels_UsingSVM <- as.character(svm_celltype_pred[match(integrated.atlas$seurat_clusters, svm_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingSVM"])
 		integrated.atlas[["HarmonizedLabels_UsingSVM"]] <- HarmonizedLabels_UsingSVM
 		message ("Added harmonized labels using SVM to integrated object")
 		write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")
 		if(validatePredictions == TRUE)
 		{
 			message ("SVM based harmonization completed. Starting validation of celltype assignments")
 			Idents(integrated.atlas) <- integrated.atlas$Celltypes
 			integrated.validation.use <- subset(integrated.atlas, idents = svm_celltype_pred$PredictedCelltype_UsingSVM)
 			validation = ValidatePredictions(reference = integrated.validation.use, query = integrated.atlas, celltype_assignments = svm_celltype_pred)
 			message ("Validation completed. Please see summary of number of shared markers below")
 			print (validation)
 			write.table(validation, "Summary_numberofsharedmarkers_between_AssignedCelltypesandIntegratedClusters.txt", quote=F, sep="\t")
 			return(integrated.atlas)
 		}
 		if(validatePredictions == FALSE)
 		{
 			message ("SVM based harmonization completed")
 			return(integrated.atlas)
 		}
 	}

 	if(classification.method == "Ensemble")
 	{
 		message ("Ensemble learning using classification accuracy of both Random Forest and SVM")

 		rf_Celltypes = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label)
		message ("Predicting using trained randomForest classifier")
		rf_pred = predict(rf_Celltypes, newdata=test)
		rf_cm = table(test_label, rf_pred)
		rf_celltype_pred <- data.frame(colnames(rf_cm)[apply(rf_cm,1,which.max)], rownames(rf_cm), check.names=F)
		colnames(rf_celltype_pred) <- c("HarmonizedLabels_UsingRF", "seurat_clusters")
 		HarmonizedLabels_UsingRF <- as.character(rf_celltype_pred[match(integrated.atlas$seurat_clusters, rf_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingRF"])
 		integrated.atlas[["HarmonizedLabels_UsingRF"]] <- HarmonizedLabels_UsingRF
 		message ("Added Harmonized Labels using randomForest to integrated object")
 		write.table(rf_cm, "ConfusionMatrix_RandomForest.txt", quote=F, sep="\t")
 		
 		svm_Celltypes = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label, crossvalidationSVM = crossvalidationSVM)
 		message ("Predicting using trained SVM classifier")
		svm_pred = predict(svm_Celltypes, newdata=test)
		svm_cm = table(test_label, svm_pred)
		svm_celltype_pred <- data.frame(colnames(svm_cm)[apply(svm_cm,1,which.max)], rownames(svm_cm), check.names=F)
		colnames(svm_celltype_pred) <- c("HarmonizedLabels_UsingSVM", "seurat_clusters")
 		HarmonizedLabels_UsingSVM <- as.character(svm_celltype_pred[match(integrated.atlas$seurat_clusters, svm_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingSVM"])
 		integrated.atlas[["HarmonizedLabels_UsingSVM"]] <- HarmonizedLabels_UsingSVM
 		message ("Added harmonized labels using SVM to integrated object")
 		write.table(svm_cm, "ConfusionMatrix_SVM.txt", quote=F, sep="\t")

		rf_cm <- as.matrix(rf_cm)
 		svm_cm <- as.matrix(svm_cm)
 	
		message ("Assigning weights based on Random Forest and SVM error estimation for Ensemble learning")
 		
 		rf_acccuracy_estimate <- (1-tail(rf_Celltypes$err.rate[,1], 1))*100
 		svm_accuracy_estimate <- svm_Celltypes$tot.accuracy

 		message (paste0("Accuracy estimate of randomForest based on training:", rf_acccuracy_estimate))
 		message (paste0("Accuracy estimate of SVM based on training:", svm_accuracy_estimate))

 		if((rf_acccuracy_estimate == "NaN") | (svm_accuracy_estimate == "NaN"))
 		{
 			consensus_cm = rf_cm + svm_cm
 		}
 		if((rf_acccuracy_estimate > 0) | (svm_accuracy_estimate > 0))
 		{
			rf_cm = rf_cm*rf_acccuracy_estimate
 			svm_cm = svm_cm*svm_accuracy_estimate
 			consensus_cm = rf_cm + svm_cm
 		}
 		if((rf_acccuracy_estimate == 0) | (svm_accuracy_estimate == 0))
 		{
 			consensus_cm = rf_cm + svm_cm
 		}
 		consensus_celltype_pred <- data.frame(colnames(consensus_cm)[apply(consensus_cm,1,which.max)], rownames(consensus_cm), check.names=F)
		colnames(consensus_celltype_pred) <- c("HarmonizedLabels_UsingEnsemble", "seurat_clusters")
 		HarmonizedLabels_UsingEnsemble <- as.character(consensus_celltype_pred[match(integrated.atlas$seurat_clusters, consensus_celltype_pred$seurat_clusters), "HarmonizedLabels_UsingEnsemble"])
 		integrated.atlas[["HarmonizedLabels_UsingEnsemble"]] <- HarmonizedLabels_UsingEnsemble
 		message ("Added Harmonized labels using Ensemble learning to query")
 		write.table(consensus_cm, "ConfusionMatrix_EnsembleLearning.txt", quote=F, sep="\t")
 		if(validatePredictions == TRUE)
 		{
 			message ("Ensembl learning based harmonization completed. Starting validation of celltype assignments")
 			Idents(integrated.atlas) <- integrated.atlas$Celltypes
 			integrated.validation.use <- subset(integrated.atlas, idents = consensus_celltype_pred$HarmonizedLabels_UsingEnsemble)
 			validation = ValidatePredictions(reference = integrated.validation.use, query = integrated.atlas, celltype_assignments = consensus_celltype_pred)
 			message ("Validation completed. Please see summary of number of shared markers below")
 			print (validation)
 			write.table(validation, "Summary_numberofsharedmarkers_between_AssignedCelltypesandQueryClusters.txt", quote=F, sep="\t")
 			return(integrated.atlas)
 		}
 		if(validatePredictions == FALSE)
 		{
 			message ("Ensembl learning based harmonization completed.")
 			return(integrated.atlas)
 		}
 	}
 }
