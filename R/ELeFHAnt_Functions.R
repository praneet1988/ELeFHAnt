#' ELeFHAnt mouse tissues
#'
#' Select tissues for CellMarker v2.0 based validation
#'
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'
mouse_tissues = c("Adipose tissue", "Adrenal gland", "Adventitia", "Afferent artery", "Airway", "Alveolar capillary", "Alveolus", "Amygdala", "Annulus fibrosus", "Anorectal junction", "Anterior lobule", "Aorta", "Aortic root", "Aortic valve", "Arm", "Artery", "Arthrosis", "Ascending aorta", "Auditory cortex", "Basilar membrane", "Bed nucleus of the stria terminalis", "Belly", "Bile duct", "Biliary tract", "Bladder", "Bladder mucosa", "Blood", "Blood brain barrier", "Blood vessel", "Bone", "Bone marrow", "Brain", "Breast", "Bronchiole", "Bronchus", "Buccal mucosal epithelium", "Caecum", "Capillary", "Cardiac neural crest", "Cardiovascular system", "Carotid artery", "Cartilage", "Cauda epididymis", "Caudal ganglionic eminence", "Central nervous system", "Cerebellar nuclei", "Cerebellum", "Cerebral cortex", "Cerebral motor cortex", "Cerebral organoid", "Cerebrospinal fluid", "Choroid", "Choroid plexus", "Choroid plexus capillary", "Cochlea", "Cochlear duct", "Colon", "Colon epithelium", "Colonic crypt", "Colorectum", "Conjunctiva", "Connective tissue", "Cornea", "Corneal epithelium", "Cornu ammonis 1", "Cornu ammonis 2", "Coronary artery", "Corpus callosum", "Cortex", "Dermal microvasculature", "Dermal papilla", "Dermis", "Diaphragm", "Distal limb mesenchyme", "Distal lung endoderm", "Dorsal forebrain", "Dorsal root ganglia", "Dorsal root ganglion", "Dorsal skin", "Dorsolateral prefrontal cortex", "Dorsomedial hypothalamus", "Ear", "Ectoderm", "Efferent artery", "Embryo", "Embryoid body", "Embryonic brain", "Embryonic breast", "Embryonic ectoderm", "Embryonic endoderm", "Embryonic heart", "Embryonic Kidney", "Embryonic mesoderm", "Embryonic stem cell", "Embryos", "Endocardium", "Endoderm", "Endodontium", "Endometrium", "Endothelium", "Enteric neural crest", "Epiblast", "Epidermis", "Epithelium", "Esophagus", "External genitalia", "Extra-embryonic ectoderm", "Extra-embryonic endoderm", "Extra-embryonic mesoderm", "Extra-embryonic tissue", "Eye", "Fat pad", "Femur bone", "Fetal hypothalamus", "Fetal kidney", "Fetal liver", "Fetal ovary", "Fetal skin", "First heart field(FHF)", "Flesh", "Focculus", "Foregut endoderm", "Ganglion cell layer of retina", "Gastric corpus", "Gastric epithelium", "Gastric gland", "Gastric isthmus", "Gastrointestinal tract", "Germinal center", "Gingiva", "Glomerular capillary", "Glomerulus", "Gonad", "Gut", "Hair canal", "Hair follicle", "Head and Neck", "Heart", "Heart muscle", "Heart valve", "Hind limb", "Hippocampus", "Hypothalamic brain slice", "Hypothalamic nucleus", "Hypothalamus", "Hypothalamus-POA", "Ileum", "Incisor", "Inferior colliculus", "Inner cell mass", "Inner Ear", "Inner nuclear layer of retina", "Interfollicular epidermis", "Intestinal crypt", "Intestine", "Juxta-cardiac field (JCF)", "Kidney", "Kidney cortex", "Knee", "Lacrimal gland", "Large intestine", "Large peritoneal", "Lateral hypothalamus", "Left ventricle", "limb", "Limb bud", "Liver", "Lobule VI", "Lower dermis", "Lower hair follicle", "Lung", "Lymph", "Lymph node", "Lymphatic vessel", "Lymphoid tissue", "Macrovessel", "Main olfactory epithelia", "Mammary epithelium", "Mammary gland", "Mandibular alveolar bone", "Meninge", "Meniscus", "Mesenteric lymph node", "Mesoderm", "Mesodermal precursor", "Mesonephros", "Microvessel", "Midbrain", "Molar", "Motor cortex", "Muscle", "Myenteric plexus", "Myocardium", "Nasal cavity", "Neocortex", "Nerve", "Neural tube", "Nodose", "Nodulus", "Non-Vasculature", "Nucleus accumbens", "Olfactory neuroepithelium", "Omentum", "Oral cavity", "Outflow tract", "Ovarian follicle", "Ovary", "Pancreas", "Pancreatic duct", "Pancreatic islet", "PeriBiliary cell gland", "Peribiliary gland", "Perichondrium", "Periosteum", "Peripheral blood", "Peritoneal cavity", "Peritoneum", "Peyer patch", "Pharynx", "Pituitary", "Placenta", "Pluripotent stem cell", "Polyp", "Posterior lobule", "Posterior second heart field", "Prefrontal cortex", "Presomitic mesoderm", "Primary motor cortex", "Primary visual cortex", "Primitive endoderm", "Primordial germ", "Prostate", "Proximal lung endoderm", "Pulmonary aorta", "Pulmonary arteriy", "Pylorus", "Red pulp", "Renal glomerulus", "Retina", "Retina vessel", "Retinal pigment epithelium", "Salivary duct", "Salivary gland", "Sciatic nerve", "Sebaceous gland", "Seminal plasma", "Serum", "Sinoatrial node", "Skeletal muscle", "Skin", "Skin of back", "Small intestinal crypt", "Small intestine", "Smooth muscle", "Soft palate", "Soft tissue", "Somatosensory cortex", "Spinal cord", "Spleen", "Stomach", "Striatum", "Subcutaneous adipose tissue", "Subgranular zone", "Submandibular gland", "Subventricular zone", "Superior cervical ganglion", "Suture mesenchyme", "Synovium", "Taste bud", "Tendon", "Testis", "Thoracic aorta", "Thymus", "Thyroid", "Tongue", "Tonsil", "Trachea", "Trophectoderm", "Umbilical cord", "Umbilical cord blood", "Undefined", "Upper hair follicle", "Urethra", "Uterine cervix", "Uterus", "Vein", "Ventral posterior hypothalamus (VPH)", "Ventral tegmental area", "Ventromedial hypothalamus (VMHvl)", "Vessel", "White adipose tissue", "White matter", "Yolk sac")

#' ELeFHAnt human tissues
#'
#' Select tissues for CellMarker v2.0 based validation
#'
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'
human_tissues = c("Abdomen", "Abdominal adipose tissue", "Abdominal fat pad", "Acinus", "Adipose tissue", "Adrenal gland", "Adventitia", "Airway", "Airway epithelium", "Allocortex", "Alveolus", "Amniotic fluid", "Amniotic membrane", "Ampullary", "Antecubital vein", "Anterior cruciate ligament", "Anterior presomitic mesoderm", "Aorta", "Aortic valve", "Artery", "Arthrosis", "Articular Cartilage", "Ascites", "Ascitic fluid", "Auditory cortex", "Basal airway", "Basilar membrane", "Bile duct", "Biliary tract", "Bladder", "Blood", "Blood vessel", "Bone", "Bone marrow", "Brain", "Breast", "Bronchial vessel", "Bronchiole", "Bronchoalveolar lavage", "Bronchoalveolar system", "Bronchus", "Brown adipose tissue", "Calvaria", "Capillary", "Cardiovascular system", "Carotid artery", "Carotid plaque", "Cartilage", "Caudal cortex", "Caudal forebrain", "Caudal ganglionic eminence", "Cavernosum", "Central amygdala", "Central nervous system", "Cerebellum", "Cerebral organoid", "Cerebrospinal fluid", "Chorionic villi", "Chorionic villus", "Choroid", "Choroid plexus", "Colon", "Colon epithelium", "Colorectum", "Cornea", "Corneal endothelium", "Corneal epithelium", "Coronary artery", "Corpus callosum", "Corpus luteum", "Cortex", "Cortical layer", "Cortical thymus", "Decidua", "Deciduous tooth", "Dental pulp", "Dermis", "Diencephalon", "Dorsal forebrain", "Dorsal root ganglion", "Dorsolateral prefrontal cortex", "Ductal tissue", "Duodenum", "Ectocervix", "Ectoderm", "Embryo", "Embryoid body", "Embryonic brain", "Embryonic heart", "Embryonic Kidney", "Embryonic stem cell", "Endocardium", "Endocrine", "Endoderm", "Endometrium", "Endometrium stroma", "Entorhinal cortex", "Epidermis", "Epithelium", "Esophagus", "Eye", "Fetal brain", "Fetal heart", "Fetal ileums", "Fetal kidney", "Fetal Leydig", "Fetal liver", "Fetal lung", "Fetal pancreas", "Fetal thymus", "Fetal umbilical cord", "Fetus", "Foreskin", "Frontal cortex", "Fundic gland", "Gall bladder", "Gastric corpus", "Gastric epithelium", "Gastric gland", "Gastrointestinal tract", "Germ", "Germinal center", "Gingiva", "Gonad", "Gut", "Hair follicle", "Head and neck", "Heart", "Heart muscle", "Hippocampus", "Ileum", "Iliac crest", "Inferior colliculus", "Interfollicular epidermis", "Intervertebral disc", "Intestinal crypt", "Intestine", "Intrahepatic cholangio", "Jejunum", "Kidney", "Lacrimal gland", "Large intestine", "Laryngeal squamous epithelium", "Larynx", "Lateral ganglionic eminence", "Left lobe", "Limb bud", "Limbal epithelium", "Liver", "Lumbar vertebra", "Lung", "Lymph", "Lymph node", "Lymphatic vessel", "Lymphoid tissue", "Malignant pleural effusion", "Mammary epithelium", "Mammary gland", "Medial ganglionic eminence", "Medullary thymus", "Meniscus", "Mesenchyme", "Mesoblast", "Mesoderm", "Microvascular endothelium", "Microvessel", "Midbrain", "Middle temporal gyrus", "Milk", "Molar", "Muscle", "Myenteric plexus", "Myocardium", "Myometrium", "Nasal concha", "Nasal epithelium", "Nasal mucosa", "Nasal polyp", "Nasopharyngeal mucosa", "Nasopharynx", "Neocortex", "Nerve", "Nose", "Nucleus pulposus", "Olfactory neuroepithelium", "Omentum", "Optic nerve", "Oral cavity", "Oral mucosa", "Osteoarthritic cartilage", "Ovarian cortex", "Ovarian follicle", "Ovary", "Oviduct", "Palatine tonsil", "Pancreas", "Pancreatic acinar tissue", "Pancreatic duct", "Pancreatic islet", "Parotid gland", "Periodontal ligament", "Periodontium", "Periosteum", "Peripheral blood", "Peritoneal fluid", "Peritoneum", "Pituitary", "Pituitary gland", "Placenta", "Plasma", "Pleura", "Pluripotent stem cell", "Polyp", "Posterior fossa", "Posterior presomitic mesoderm", "Prefrontal cortex", "Premolar", "Presomitic mesoderm", "Primitive streak", "Prostate", "Pulmonary arteriy", "Pyloric gland", "Rectum", "Renal glomerulus", "Respiratory tract", "Retina", "Retinal organoid", "Retinal pigment epithelium", "Right ventricle", "Saliva", "Salivary gland", "Scalp", "Sclerocorneal tissue", "Seminal plasma", "Septum transversum", "Serum", "Serum exosome", "Sinonasal mucosa", "Sinus tissue", "Skeletal muscle", "Skin", "Small intestinal crypt", "Small intestine", "Soft tissue", "Sperm", "Spinal cord", "Spleen", "Splenic red pulp", "Sputum", "Stomach", "Subcutaneous adipose tissue", "Submandibular gland", "Subpallium", "Subplate", "Subventricular zone", "Superior frontal gyrus", "Sympathetic ganglion", "Synovial fluid", "Synovium", "Taste bud", "Tendon", "Testis", "Thalamus", "Thymus", "Thyroid", "Tongue", "Tonsil", "Tooth", "Trachea", "Tracheal airway epithelium", "Transformed artery", "Trophoblast", "Umbilical cord", "Umbilical cord blood", "Umbilical vein", "Undefined", "Urine", "Urothelium", "Uterine cervix", "Uterus", "Vagina", "Vein", "Venous blood", "Ventral thalamus", "Ventricular and atrial", "Ventricular zone", "Vessel", "Visceral adipose tissue", "Vocal cord", "Vocal fold", "White adipose tissue", "White matter", "Yolk sac")


#' ELeFHAnt Celltype Annotation
#'
#' Celltype annotation is a function to annotate celltypes in a single cell datasets.
#' It requires a reference dataset and a query dataset. Ensemble classifiction method is used
#' to learn celltypes from reference dataset and then predict celltypes for query dataset.
#'
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference, enabling fast computation.
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 1000]
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis and assessing markers from CellMarker Database. 
#' @param selectvarfeatures number of variable features to select while training (default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param species human or mouse if validatePredictions = TRUE
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @param annotationCol Name of column in metdata which has Celltypes information in reference [Default: Celltypes]
#' @return query seurat object with predictions added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

CelltypeAnnotation <- function(reference = NULL, query = NULL, downsample = TRUE, downsample_to = 1000, validatePredictions = TRUE, selectvarfeatures = 2000, ntree = 500, species = NULL, tissue = NULL, annotationCol = "Celltypes") {
    if(downsample == TRUE)
    {
        message ("Setting Assay of reference and query to RNA")
        DefaultAssay(reference) <- "RNA"
        DefaultAssay(query) <- "RNA"

        Idents(reference) <- annotationCol
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

        reference_use = NormalizeData(reference_use)
        reference_use = FindVariableFeatures(reference_use, selection.method = "vst", nfeatures = selectvarfeatures)
        reference_use = ScaleData(reference_use)

        query_use = query
        query_use = NormalizeData(query_use)
        query_use = FindVariableFeatures(query_use, selection.method = "vst", nfeatures = selectvarfeatures)
        query_use = ScaleData(query_use)

        message('Finding common variable features between reference and query')
        features_common = intersect(VariableFeatures(reference_use), VariableFeatures(query_use))

        if(length(features_common) == 0)
        {
            stop('No common variable features found between reference and query! Please increase the number of variable featres using selectvarfeatures')
        }

        message('Subsetting reference and query for common variable features')
        reference_use = subset(reference_use, features = features_common)
        reference_use = NormalizeData(reference_use)
        query_use = subset(query_use, features = features_common)
        query_use = NormalizeData(query_use)
    }

    if(downsample == FALSE)
    {
        message ("Setting Assay of reference and query to RNA")
        DefaultAssay(reference) <- "RNA"
        DefaultAssay(query) <- "RNA"

        Idents(reference) <- annotationCol
        Idents(query) <- query$seurat_clusters

        message("Running Diagonistis on reference and query")
        reference_cells <- nrow(reference@meta.data)
        query_cells <- nrow(query@meta.data)

        message (paste0("Number of cells in reference:", reference_cells))
        message (paste0("Number of cells in query:", query_cells))

        message ("Calculating ratio of number of cells in reference vs query")
        Ratio = query_cells / reference_cells
        message (paste0("Ratio of number of cells in query vs reference:", Ratio))

        reference_use = reference
        reference_use = NormalizeData(reference_use)
        reference_use = FindVariableFeatures(reference_use, selection.method = "vst", nfeatures = selectvarfeatures)
        reference_use = ScaleData(reference_use)

        query_use = query
        query_use = NormalizeData(query_use)
        query_use = FindVariableFeatures(query_use, selection.method = "vst", nfeatures = selectvarfeatures)
        query_use = ScaleData(query_use)

        message('Finding common variable features between reference and query')
        features_common = intersect(VariableFeatures(reference_use), VariableFeatures(query_use))

        message('Subsetting reference and query for common variable features')
        reference_use = subset(reference_use, features = features_common)
        reference_use = NormalizeData(reference_use)
        query_use = subset(query_use, features = features_common)
        query_use = NormalizeData(query_use)
    }
    message('Preparing train and test datasets from reference and query')
    reference_matrix = reference_use[['RNA']]@data
    reference_matrix = Matrix::t(reference_matrix)
    reference_matrix = data.frame(reference_matrix, stringsAsFactors = FALSE)
    num_features = ncol(reference_matrix)
    reference_matrix$Celltypes = reference_use@meta.data[,annotationCol]
    message('Scaling reference to obtain training set')
    scaled_reference = scale(reference_matrix[,1:num_features],center=TRUE,scale=TRUE)
    train = data.frame(scaled_reference)
    train_label = reference_matrix$Celltypes

    query_matrix = query_use[['RNA']]@data
    query_matrix = Matrix::t(query_matrix)
    query_matrix = data.frame(query_matrix, stringsAsFactors = FALSE)
    num_features = ncol(query_matrix)
    query_matrix$seurat_clusters = query_use$seurat_clusters
    message('Scaling query to obtain test set')
    scaled_query = scale(query_matrix[,1:num_features],attr(scaled_reference,"scaled:center"),attr(scaled_reference,"scaled:scale"))
    test = data.frame(scaled_query)
    test_label = query_matrix$seurat_clusters

    message('\nSetting up three classifiers: randomForest, SVM and LR')
    rf = randomForest_predictor(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
    svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label)
    lr = lr_predictor(train = train, test = test, train_label = train_label, test_label = test_label)

    message('\nClassifying cells in query using each classifier')
    rf_pred = predict(rf[[1]], test)
    rf_Probablity = data.frame(rf_pred$predictions, check.names = F)
    rownames(rf_Probablity) = colnames(query)
    rf_Probablity_Predictions <- data.frame(colnames(rf_Probablity)[apply(rf_Probablity,1,which.max)], rownames(rf_Probablity), check.names=F)
    colnames(rf_Probablity_Predictions) = c("ELeFHAnt_RF_CelltypePrediction", "Cells")
    query$ELeFHAnt_RF_CelltypePrediction = rf_Probablity_Predictions$ELeFHAnt_RF_CelltypePrediction
    colnames(rf_Probablity) = paste0("ELeFHAnt_RF_", colnames(rf_Probablity), " Probability")
    query = AddMetaData(query, metadata = rf_Probablity)

    svm_pred = predict(svm[[1]], test, decisionValues = TRUE)
    svm_Probablity = data.frame(svm_pred$decisionValues, check.names = F)
    rownames(svm_Probablity) = colnames(query)
    svm_Predictions <- svm_pred$predictions
    query$ELeFHAnt_SVM_CelltypePrediction = svm_Predictions
    colnames(svm_Probablity) = paste0("ELeFHAnt_SVM_", colnames(svm_Probablity), " Decision Values")
    query = AddMetaData(query, metadata = svm_Probablity)

    lr_pred = predict(lr[[1]], test, prob = TRUE, decisionValues = TRUE)
    lr_Probablity = data.frame(lr_pred$probabilities, check.names = F)
    rownames(lr_Probablity) = colnames(query)
    lr_Predictions <- lr_pred$predictions
    query$ELeFHAnt_LR_CelltypePrediction = lr_Predictions
    colnames(lr_Probablity) = paste0("ELeFHAnt_LR_", colnames(lr_Probablity), " Probability")
    query = AddMetaData(query, metadata = lr_Probablity)

    predictions = data.frame(query$ELeFHAnt_RF_CelltypePrediction, query$ELeFHAnt_SVM_CelltypePrediction, query$ELeFHAnt_LR_CelltypePrediction)
    colnames(predictions) = c("RF", "SVM", "LR")

    message('\nObtaing Ensemble Predictions using RF, SVM and LR')
    ELeFHAnt_Ensemble_CelltypePrediction = c()
    for(f in 1:nrow(predictions))
    {
        value_frame = data.frame(table(c(as.character(predictions$RF[f]), as.character(predictions$SVM[f]), as.character(predictions$LR[f]))))
        colnames(value_frame) = c("Celltypes", "Frequency")
        rownames(value_frame) = value_frame$Celltypes
        value_frame = value_frame %>% arrange(desc(Frequency))
        if(max(value_frame$Frequency) > 1)
        {
            ELeFHAnt_Ensemble_CelltypePrediction = c(ELeFHAnt_Ensemble_CelltypePrediction, rownames(value_frame)[1])
        }   
        if(max(value_frame$Frequency) == 1)
        {
            ELeFHAnt_Ensemble_CelltypePrediction = c(ELeFHAnt_Ensemble_CelltypePrediction, as.character(predictions$SVM[f]))   ############## we break ties using SVM predictions.
        }
    }
    query$ELeFHAnt_Ensemble_CelltypePrediction = ELeFHAnt_Ensemble_CelltypePrediction
    message('\nCelltype predictions are stored in query metadata. Please see: ELeFHAnt_RF_CelltypePrediction, ELeFHAnt_SVM_CelltypePrediction, ELeFHAnt_LR_CelltypePrediction, ELeFHAnt_Ensemble_CelltypePrediction')

    if(validatePredictions == TRUE)
    {
        message("Ensembl celltype annotation completed. Starting validation of celltype assignments using GSEA")
        validation = ValidatePredictions(species = species, query = query, tissue = tissue)
        message ("Validation completed. Please see ValidatePredictions Folder for results")
        return(query)
    }
    if(validatePredictions == FALSE)
    {
        message("Ensembl celltype annotation completed.")
        return(query)
    }
}


#' ELeFHAnt Label Harmonization
#'
#' Label Harmonization is a function to harmonize cell labels (celltypes) across single cell datasets.
#' It requires a list of processed Seurat Objects or a integrated seurat object.
#'
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @param seurat.objects a list of processed seurat objects (please set Default Assay to "RNA") with Celltypes column in their respective meta.data to perform integration on
#' @param perform_integration logical Indicator (TRUE or FALSE) to perform integration using list of seurat.objects
#' @param integrated.atlas an integrated seurat object with CellTypes and seurat_clusters column in meta.data. Required if perform_integration = FALSE
#' @param downsample logical Indicator (TRUE or FALSE) to downsample Seurat objects or integrated seurat object, enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 1000]
#' @param npcs number of principal components to compute after integration [Default: 30]
#' @param resolution value of the resolution parameter, decides size of cell communities. [Default: 0.8]
#' @param validatePredictions logical indicator (TRUE or FALSE) to asses predictions by deploying Gene set enrichment analysis
#' @param selectanchorfeatures number of anchor features to use for integrating datasets (Default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param species human or mouse if validatePredictions = TRUE
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @param k.anchor  How many neighbors (k) to use when picking anchors [Default: 5]
#' @param k.filter  How many neighbors (k) to use when filtering anchors [Default: 200]
#' @param k.score  How many neighbors (k) to use when scoring anchors [Default: 30]
#' @param dims  Which dimensions to use from the CCA to specify the neighbor search space [Default: 1:30]
#' @param annotationCol Name of column in metdata which has Celltypes information [Default: Celltypes]. It has to be same across atlases you want to integrate or have integrated
#' @return integrated seurat object with harmonized celltypes added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

LabelHarmonization <- function(seurat.objects = c(), perform_integration = FALSE, integrated.atlas = NULL, downsample = TRUE, downsample_to = 1000, npcs = 30, resolution = 0.8, validatePredictions = TRUE, selectanchorfeatures = 2000, ntree = 500, k.anchor = 5, k.filter = 200, k.score = 30, dims = 1:30, species = NULL, tissue = NULL, annotationCol = "Celltypes") {
    if(perform_integration == TRUE)
    {
        if(downsample == TRUE)
        {
            message ("Downsampling seurat objects")
            seurat.objects <- lapply(X = seurat.objects, FUN = function(x) {
                DefaultAssay(x) <- "RNA"
                Idents(x) <- annotationCol
                x <- subset(x, downsample = downsample_to)
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = selectanchorfeatures)
            })
        }
        if(downsample == FALSE)
        {
            seurat.objects <- seurat.objects
        }
        message ("Starting integration using Seurat Canonical Correlation Algorithm")
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
            Idents(integrated.atlas) <- integrated.atlas$seurat_clusters
            integrated.use <- subset(integrated.atlas, downsample = downsample_to)
            num_cells <- nrow(integrated.use@meta.data)
            message (paste0("Number of cells in integrated atlas:", num_cells))
        }
        if(downsample == FALSE)
        {
            message ("Setting Assay of integrated.atlas to integrated")
            integrated.atlas <- integrated.atlas
            DefaultAssay(integrated.atlas) <- "integrated"
            Idents(integrated.atlas) <- integrated.atlas$seurat_clusters
            integrated.use <- integrated.atlas
            num_cells <- nrow(integrated.use@meta.data)
            message (paste0("Number of cells in integrated atlas:", num_cells))
        }
    }
    message ("Generating train and test datasets using stratification -- 70% for training & 30% for testing")
    integrated_data <- integrated.use[['integrated']]@scale.data
    integrated_data <- Matrix::t(integrated_data)
    integrated_data <- data.frame(integrated_data)
    num_features <- ncol(integrated_data)
    message (paste0("Number of Anchor Features selected:", num_features))
    integrated_data$Celltypes <- integrated.use@meta.data[,annotationCol]
    integrated_data$Clusters <- integrated.use$seurat_clusters
    integrated_data_stratified <- stratified(integrated_data, group = "Clusters", size = 0.7, bothSets = TRUE)
    train <- integrated_data_stratified$SAMP1[,1:num_features]
    test <- integrated_data_stratified$SAMP2[,1:num_features]
    train_label <- integrated_data_stratified$SAMP1$Celltypes
    test_label <- integrated_data_stratified$SAMP2$Clusters

    message('\nSetting up three classifiers: randomForest, SVM and LR')
    rf = randomForest_predictor_withoutProb(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
    svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label)
    lr = lr_predictor(train = train, test = test, train_label = train_label, test_label = test_label)

    message('\nClassifying cells in query using each classifier and Generating scaled confusion matrix')
    rf_pred = predict(rf[[1]], test)
    rf_cm = as.matrix(table(test_label, rf_pred$predictions))

    svm_pred = predict(svm[[1]], test, decisionValues = TRUE)
    svm_cm = as.matrix(table(test_label, svm_pred$predictions))

    lr_pred = predict(lr[[1]], test, prob = TRUE, decisionValues = TRUE)
    lr_cm = as.matrix(table(test_label, lr_pred$predictions))

    rf_cm_norm = t(apply(rf_cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
    write.table(rf_cm_norm, "ConfusionMatrix_RandomForest_LabelHarmonization.txt", sep="\t", quote=F)
    rf_predictions = data.frame(colnames(rf_cm_norm)[apply(rf_cm_norm,1,which.max)], rownames(rf_cm_norm), check.names=F)
    colnames(rf_predictions) <- c("ELeFHAnt_RF_HarmonizedCelltype", "seurat_clusters")
    ELeFHAnt_RF_HarmonizedCelltype <- as.character(rf_predictions[match(integrated.atlas$seurat_clusters, rf_predictions$seurat_clusters), "ELeFHAnt_RF_HarmonizedCelltype"])
    integrated.atlas$ELeFHAnt_RF_HarmonizedCelltype = ELeFHAnt_RF_HarmonizedCelltype

    svm_cm_norm = t(apply(svm_cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
    write.table(svm_cm_norm, "ConfusionMatrix_SVM_LabelHarmonization.txt", sep="\t", quote=F)
    svm_predictions = data.frame(colnames(svm_cm_norm)[apply(svm_cm_norm,1,which.max)], rownames(svm_cm_norm), check.names=F)
    colnames(svm_predictions) <- c("ELeFHAnt_SVM_HarmonizedCelltype", "seurat_clusters")
    ELeFHAnt_SVM_HarmonizedCelltype <- as.character(svm_predictions[match(integrated.atlas$seurat_clusters, svm_predictions$seurat_clusters), "ELeFHAnt_SVM_HarmonizedCelltype"])
    integrated.atlas$ELeFHAnt_SVM_HarmonizedCelltype = ELeFHAnt_SVM_HarmonizedCelltype

    lr_cm_norm = t(apply(lr_cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
    write.table(lr_cm_norm, "ConfusionMatrix_LR_LabelHarmonization.txt", sep="\t", quote=F)
    lr_predictions = data.frame(colnames(lr_cm_norm)[apply(lr_cm_norm,1,which.max)], rownames(lr_cm_norm), check.names=F)
    colnames(lr_predictions) <- c("ELeFHAnt_LR_HarmonizedCelltype", "seurat_clusters")
    ELeFHAnt_LR_HarmonizedCelltype <- as.character(lr_predictions[match(integrated.atlas$seurat_clusters, lr_predictions$seurat_clusters), "ELeFHAnt_LR_HarmonizedCelltype"])
    integrated.atlas$ELeFHAnt_LR_HarmonizedCelltype = ELeFHAnt_LR_HarmonizedCelltype

    cm = rf_cm_norm + svm_cm_norm + lr_cm_norm
    cm_norm = t(apply(cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
    write.table(cm_norm, "ConfusionMatrix_Ensemble_LabelHarmonization.txt", sep="\t", quote=F)
    cm_predictions = data.frame(colnames(cm_norm)[apply(cm_norm,1,which.max)], rownames(cm_norm), check.names=F)
    colnames(cm_predictions) <- c("ELeFHAnt_Ensemble_HarmonizedCelltype", "seurat_clusters")
    ELeFHAnt_Ensemble_HarmonizedCelltype <- as.character(cm_predictions[match(integrated.atlas$seurat_clusters, cm_predictions$seurat_clusters), "ELeFHAnt_Ensemble_HarmonizedCelltype"])
    integrated.atlas$ELeFHAnt_Ensemble_HarmonizedCelltype = ELeFHAnt_Ensemble_HarmonizedCelltype

    message('\nHarmonized Celltype predictions are stored in integrated metadata. Please see: ELeFHAnt_RF_HarmonizedCelltype, ELeFHAnt_SVM_HarmonizedCelltype, ELeFHAnt_LR_HarmonizedCelltype, ELeFHAnt_Ensemble_HarmonizedCelltype')

    if(validatePredictions == TRUE)
    {
        message("Ensembl celltype harmonization completed. Starting validation of celltype assignments using GSEA & CellMarker")
        validation = ValidatePredictions(species = species, query = integrated.atlas, tissue = tissue)
        message ("Validation completed. Please see ValidatePredictions Folder for results")
        return(integrated.atlas)
    }
    if(validatePredictions == FALSE)
    {
        message("Ensembl celltype harmonization completed.")
        return(integrated.atlas)
    }
 }


 #' ELeFHAnt Deduce Relationship
#'
#' Deduce Relationship is a function that hels deduce relationships among celltypes/annotationCol between two datasets.
#' It requires two datasets (both processed Seurat Objects). Obtain relationships among metadata of interest. Function outputs a heatmap with annotation from dataset1 y-axis and annotation from dataset2 on x-axis.
#'
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @param reference1 a processed Seurat Object with Celltypes column in metadata
#' @param reference2 a processed seurat object with Celltypes column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference1 and reference2, enabling fast computation
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 1000]
#' @param selectvarfeatures number of variable features to select while training (Default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param annotationCol_ref1 Name of column in metdata which has Celltypes information [Default: Celltypes] in reference1
#' @param annotationCol_ref2 Name of column in metdata which has Celltypes information [Default: Celltypes] in reference2
#' @return ggplot2 heatmap object and heatmap is automatically saved
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

DeduceRelationship <- function(reference1 = NULL, reference2 = NULL, downsample = TRUE, downsample_to = 1000, selectvarfeatures = 2000, ntree = 500, annotationCol_ref1 = "Celltypes", annotationCol_ref2 = "Celltypes") {
  if(downsample == TRUE)
  {
    message ("Setting Assay of reference1 and reference2 to RNA")
    DefaultAssay(reference1) <- "RNA"
    DefaultAssay(reference2) <- "RNA"

    num_cells_reference1 <- nrow(reference1@meta.data)
    message (paste0("Number of cells in reference1:", num_cells_reference1))

    num_cells_reference2 <- nrow(reference2@meta.data)
    message (paste0("Number of cells in reference2:", num_cells_reference2))
    
    Idents(reference1) <- annotationCol_ref1
    Idents(reference2) <- annotationCol_ref2
    
    reference1_use <- subset(reference1, downsample = downsample_to)
    reference1_use = NormalizeData(reference1_use)
    reference1_use = FindVariableFeatures(reference1_use, selection.method = "vst", nfeatures = selectvarfeatures)
    reference1_use = ScaleData(reference1_use)

    reference2_use <- subset(reference2, downsample = downsample_to)
    reference2_use = NormalizeData(reference2_use)
    reference2_use = FindVariableFeatures(reference2_use, selection.method = "vst", nfeatures = selectvarfeatures)
    reference2_use = ScaleData(reference2_use)

    num_cells_reference1_use <- nrow(reference1_use@meta.data)
    message (paste0("Number of cells in reference1 after downsampling:", num_cells_reference1_use))

    num_cells_reference2_use <- nrow(reference2_use@meta.data)
    message (paste0("Number of cells in reference2 after downsampling:", num_cells_reference2_use))

    message('Finding common variable features between reference and query')
    features_common = intersect(VariableFeatures(reference1_use), VariableFeatures(reference2_use))

    if(length(features_common) == 0)
    {
        stop('No common variable features found between reference1 and reference2! Please increase the number of variable featres using selectvarfeatures')
    }

    message('Subsetting reference1 and reference2 for common variable features')
    reference1_use = subset(reference1_use, features = features_common)
    reference1_use = NormalizeData(reference1_use)
    reference2_use = subset(reference2_use, features = features_common)
    reference2_use = NormalizeData(reference2_use)
  }
  
  if(downsample == FALSE)
  {
    message ("Setting Assay of reference1 and reference2 to RNA")
    DefaultAssay(reference1) <- "RNA"
    DefaultAssay(reference2) <- "RNA"

    num_cells_reference1 <- nrow(reference1@meta.data)
    message (paste0("Number of cells in reference1:", num_cells_reference1))

    num_cells_reference2 <- nrow(reference2@meta.data)
    message (paste0("Number of cells in reference2:", num_cells_reference2))
    
    Idents(reference1) <- annotationCol_ref1
    Idents(reference2) <- annotationCol_ref2
    
    reference1_use <- reference1
    reference1_use = NormalizeData(reference1_use)
    reference1_use = FindVariableFeatures(reference1_use, selection.method = "vst", nfeatures = selectvarfeatures)
    reference1_use = ScaleData(reference1_use)

    reference2_use <- reference2
    reference2_use = NormalizeData(reference2_use)
    reference2_use = FindVariableFeatures(reference2_use, selection.method = "vst", nfeatures = selectvarfeatures)
    reference2_use = ScaleData(reference2_use)

    message('Finding common variable features between reference and query')
    features_common = intersect(VariableFeatures(reference1_use), VariableFeatures(reference2_use))

    if(length(features_common) == 0)
    {
        stop('No common variable features found between reference1 and reference2! Please increase the number of variable featres using selectvarfeatures')
    }

    message('Subsetting reference1 and reference2 for common variable features')
    reference1_use = subset(reference1_use, features = features_common)
    reference1_use = NormalizeData(reference1_use)
    reference2_use = subset(reference2_use, features = features_common)
    reference2_use = NormalizeData(reference2_use)
  }
  
  message('Preparing train and test datasets from reference1 and reference2')
  reference1_matrix = reference1_use[['RNA']]@data
  reference1_matrix = Matrix::t(reference1_matrix)
  reference1_matrix = data.frame(reference1_matrix, stringsAsFactors = FALSE)
  num_features = ncol(reference1_matrix)
  reference1_matrix$Annotation1 = reference1_use@meta.data[,annotationCol_ref1]
  message('Scaling reference1 to obtain training set')
  scaled_reference1 = scale(reference1_matrix[,1:num_features],center=TRUE,scale=TRUE)
  train = data.frame(scaled_reference1)
  train_label = reference1_matrix$Annotation1

  reference2_matrix = reference2_use[['RNA']]@data
  reference2_matrix = Matrix::t(reference2_matrix)
  reference2_matrix = data.frame(reference2_matrix, stringsAsFactors = FALSE)
  num_features = ncol(reference2_matrix)
  reference2_matrix$Annotation2 = reference2_use@meta.data[,annotationCol_ref2]
  message('Scaling reference2 to obtain test set')
  scaled_reference2 = scale(reference2_matrix[,1:num_features],attr(scaled_reference1,"scaled:center"),attr(scaled_reference1,"scaled:scale"))
  test = data.frame(scaled_reference2)
  test_label = reference2_matrix$Annotation2

  message('\nSetting up three classifiers: randomForest, SVM and LR')
  rf = randomForest_predictor_withoutProb(train = train, test = test, train_label = train_label, test_label = test_label, ntree = ntree)
  svm = svm_predictor(train = train, test = test, train_label = train_label, test_label = test_label)
  lr = lr_predictor(train = train, test = test, train_label = train_label, test_label = test_label)
  
  message('\nClassifying cells in query using each classifier and Generating scaled confusion matrix')
  rf_pred = predict(rf[[1]], test)
  rf_cm = as.matrix(table(test_label, rf_pred$predictions))

  svm_pred = predict(svm[[1]], test, decisionValues = TRUE)
  svm_cm = as.matrix(table(test_label, svm_pred$predictions))

  lr_pred = predict(lr[[1]], test, prob = TRUE, decisionValues = TRUE)
  lr_cm = as.matrix(table(test_label, lr_pred$predictions))
  
  rf_cm_norm = t(apply(rf_cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
  write.table(rf_cm_norm, "ConfusionMatrix_RandomForest_DeduceRelationship.txt", quote=F, sep="\t")
  svm_cm_norm = t(apply(svm_cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
  write.table(svm_cm_norm, "ConfusionMatrix_SVM_DeduceRelationship.txt", quote=F, sep="\t")
  lr_cm_norm = t(apply(lr_cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
  write.table(lr_cm_norm, "ConfusionMatrix_LR_DeduceRelationship.txt", quote=F, sep="\t")
  cm = rf_cm_norm + svm_cm_norm + lr_cm_norm
  cm_norm = t(apply(cm, 1, function(x)(x-min(x))/(max(x)-min(x))))
  write.table(cm_norm, "ConfusionMatrix_Ensemble_DeduceRelationship.txt", quote=F, sep="\t")

  message('\nUsing Relative Similarity from Normalized Confusion matrices to generate Reference1 vs Reference2 similarity')
  ############################### RF ##################################
  row_order <- hclust(dist(rf_cm_norm))$order
  col_order <- hclust(dist(t(rf_cm_norm)))$order
  rf_cm_norm <- rf_cm_norm[match(rownames(rf_cm_norm)[row_order], rownames(rf_cm_norm)),match(colnames(rf_cm_norm)[col_order], colnames(rf_cm_norm))]
  rf_df <- melt(rf_cm_norm)
  colnames(rf_df) <- c("Reference2_Annotation","Reference1_Annotation","Relative_Similarity")
  plot = ggplot(data = rf_df, aes(x=Reference2_Annotation, y=Reference1_Annotation, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
  plot
  ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_RandomForest.png", width = 10, height = 10, dpi = 800)
  ####################################################################

  ############################### SVM ##################################
  row_order <- hclust(dist(svm_cm_norm))$order
  col_order <- hclust(dist(t(svm_cm_norm)))$order
  svm_cm_norm <- svm_cm_norm[match(rownames(svm_cm_norm)[row_order], rownames(svm_cm_norm)),match(colnames(svm_cm_norm)[col_order], colnames(svm_cm_norm))]
  svm_df <- melt(svm_cm_norm)
  colnames(svm_df) <- c("Reference2_Annotation","Reference1_Annotation","Relative_Similarity")
  plot = ggplot(data = svm_df, aes(x=Reference2_Annotation, y=Reference1_Annotation, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
  plot
  ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_SVM.png", width = 10, height = 10, dpi = 800)
  ####################################################################
  
  ############################### LR ##################################
  row_order <- hclust(dist(lr_cm_norm))$order
  col_order <- hclust(dist(t(lr_cm_norm)))$order
  lr_cm_norm <- lr_cm_norm[match(rownames(lr_cm_norm)[row_order], rownames(lr_cm_norm)),match(colnames(lr_cm_norm)[col_order], colnames(lr_cm_norm))]
  lr_df <- melt(lr_cm_norm)
  colnames(lr_df) <- c("Reference2_Annotation","Reference1_Annotation","Relative_Similarity")
  plot = ggplot(data = lr_df, aes(x=Reference2_Annotation, y=Reference1_Annotation, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
  plot
  ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_LR.png", width = 10, height = 10, dpi = 800)
  ####################################################################
  
  ############################### Ensemble ##################################
  row_order <- hclust(dist(cm_norm))$order
  col_order <- hclust(dist(t(cm_norm)))$order
  cm_norm <- cm_norm[match(rownames(cm_norm)[row_order], rownames(cm_norm)),match(colnames(cm_norm)[col_order], colnames(cm_norm))]
  cvm_df <- melt(cm_norm)
  colnames(cvm_df) <- c("Reference2_Annotation","Reference1_Annotation","Relative_Similarity")
  plot = ggplot(data = cvm_df, aes(x=Reference2_Annotation, y=Reference1_Annotation, fill=Relative_Similarity)) + geom_tile() + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=14),axis.title=element_text(size=18),legend.title=element_text(size=14),legend.text=element_text(size=12))    
  plot
  ggsave("Reference1_vs_Reference2_RelationshipInference_Heatmap_Ensemble.png", width = 10, height = 10, dpi = 800)
  return(plot)
  ####################################################################

}

#' ELeFHAnt Validate Predictions
#'
#' Validate Predictions is a function to validate celltype assignments. It uses Gene Set Enrichment Analysis (GSEA) and cell type markers from CellMarker Database to let user assess celltype annotation/harmonization.
#' It requires species, tissue and query dataset. For GSEA, ELeFHAnt uses C8 Hallmark Gene Sets and from CellMarker Database it uses markers curated based on experiments.
#'
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @param species human or mouse to select the C8 hallmark cell type gene sets
#' @param tissue please check human_tissues or mouse_tissues if validatePredictions = TRUE
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param reference a processed seurat object with Celltypes column in metadata
#' @return GSEA result
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

ValidatePredictions <- function(species = NULL, tissue = NULL, query = NULL) {

    DefaultAssay(query) = "RNA"
    message("\nSetting up Directory to write ValidatePredictions Results\n")
    current_WD = getwd()
    dir_create.temp = paste0(current_WD, "/", "ValidatePredictions/")
    dir.create(dir_create.temp)
    dir_create_GSEA = paste0(dir_create.temp, "/", "GSEA_usingC8HallmarkGeneSets/")
    dir.create(dir_create_GSEA)

    dir_create_Cellmarkers.temp = paste0(dir_create.temp, "/", "CellMarkers/")
    dir.create(dir_create_Cellmarkers.temp)
    

    message("Extracting Cell type Markers from CellMarker Database v2.0 [Experiment Based] and C8 Hallmark gene sets GSEA")
    
    if(species == "human")
    {
        fgsea_sets = msigdbr(species = "human", category = "C8")
        if(file.exists("Cell_marker_Human.xlsx"))
        {
            cellmarkers = read_excel('Cell_marker_Human.xlsx')
            cellmarkers = data.frame(cellmarkers)
            cellmarkers_human = subset(cellmarkers, Species == "Human")
            cellmarkers_experiment = subset(cellmarkers_human, Marker.source == "Experiment")
        }
        else
        {
            download.file('https://github.com/praneet1988/ELeFHAnt/raw/main/Cell_marker_Human.xlsx', destfile = "Cell_marker_Human.xlsx")
            cellmarkers = read_excel('Cell_marker_Human.xlsx')
            cellmarkers = data.frame(cellmarkers)
            cellmarkers_human = subset(cellmarkers, Species == "Human")
            cellmarkers_experiment = subset(cellmarkers_human, Marker.source == "Experiment")
        }
        
    }
    if(species == "mouse")
    {
        fgsea_sets = msigdbr(species = "mouse", category = "C8")
        if(file.exists("Cell_marker_Mouse.xlsx"))
        {
            cellmarkers = read_excel('Cell_marker_Mouse.xlsx')
            cellmarkers = data.frame(cellmarkers)
            cellmarkers_mouse = subset(cellmarkers, Species == "Mouse")
            cellmarkers_experiment = subset(cellmarkers_mouse, Marker.source == "Experiment")
        }
        else
        {
            download.file('https://github.com/praneet1988/ELeFHAnt/raw/main/Cell_marker_Mouse.xlsx', destfile = "Cell_marker_Mouse.xlsx")
            cellmarkers = read_excel('Cell_marker_Mouse.xlsx')
            cellmarkers = data.frame(cellmarkers)
            cellmarkers_mouse = subset(cellmarkers, Species == "Mouse")
            cellmarkers_experiment = subset(cellmarkers_mouse, Marker.source == "Experiment")
        }
    }

    message('\nGSEA BASED VALIDATION\n')
    msigdbr_list = split(x = fgsea_sets$gene_symbol, f = fgsea_sets$gs_name)
    message ("Obtaining markers per annotated cluster")
    DefaultAssay(query) = "RNA"
    Idents(query) <- query$seurat_clusters
    cluster_markers <- FindAllMarkers(query, max.cells.per.ident = 500)
    cluster_markers = data.frame(cluster_markers)
    cluster_markers = subset(cluster_markers, p_val_adj <= 0.05)
    top_cluster <- cluster_markers
    top_cluster_onlyTop25 <- cluster_markers %>% group_by(cluster) %>% dplyr::slice(1:25)
    marker_genes = unique(top_cluster_onlyTop25$gene)
    message ("Performing Gene Set Enrichment Analysis (GSEA) using gene sets from C8 Hallmark MsigDB")
    iter = length(unique(top_cluster$cluster))-1
    gsea.res.return <- c()
    for(f in 0:iter)
    {   
        message(paste0('Obtaing GSEA statistics for cluster:', f))
        data <- c()
        cluster.data <- c()
        cluster.data.use <- c()
        ranks <- c()
        gsea.res <- c()
        cluster_info <- c()
        cluster_number <- c()
        fgseaResTidy <- c()
        cluster_number <- paste0("seurat_cluster:", f)
        data <- top_cluster
        data %>% dplyr::filter(cluster == f) %>% arrange(desc(avg_log2FC), p_val_adj) %>% head(n = 100)
        cluster.data <- data %>% dplyr::filter(cluster == f) %>% arrange(p_val_adj) %>% dplyr::select(gene, avg_log2FC)
        cluster.data.use <- data.frame(cluster.data$gene, cluster.data$avg_log2FC)
        ranks <- deframe(cluster.data.use)
        gsea.res <- fgseaMultilevel(msigdbr_list, stats = ranks, scoreType = "pos", eps = 0, nPermSimple = 1000)
        fgseaResTidy <- gsea.res %>% as_tibble() %>% arrange(desc(NES))
        fgseaResTidy = drop_na(fgseaResTidy)
        if(min(fgseaResTidy$pval) <= 0.05)
        {
            message(paste0('Generating a Barplot with Normalized Enrichment Score for cluster:', f))
            filename_GSEA_plot = paste0(dir_create_GSEA, "/", 'Query_GeneSetEnrichmentAnalysis_NES_BarPlot_Cluster', '-', f, '.png')
            plot = ggplot(fgseaResTidy %>% filter(pval < 0.05) %>% head(n=10), aes(reorder(pathway, NES), NES)) + geom_col(aes(fill= pval < 0.001)) + coord_flip() + labs(x="Celltypes", y="Normalized Enrichment Score", title="C8: cell type signature gene sets GSEA") + theme_minimal() + theme_classic() + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
            plot
            ggsave(filename_GSEA_plot, plot, width = 10, height = 10, dpi = 800)
        }
        if(min(fgseaResTidy$pval) > 0.05)
        {
            message(paste0('No Significant Enrichment reported from GSEA at pval <= 0.05 for Cluster:', f))
        }  
    }
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
        cellmarkers_tissue = subset(cellmarkers_experiment, Tissue.type == tissue[i])
        cellmarkers_tissue <- cellmarkers_tissue[cellmarkers_tissue$Marker %in% marker_genes, ]
        if(is.null(cellmarkers_tissue))
        {
            stop('CellMarker Database Validation: Stopped! None query cluster markers present in Cellmarker')
        }
        celltypes = unique(cellmarkers_tissue$Cell.name)
        for(f in 1:length(celltypes))
        {
            tryCatch({
                message(paste0('Generating DotPlot/FeaturePlot for experimental evidence based markers for:', celltypes[f], " in ", tissue[i]))
                tempdata = subset(cellmarkers_tissue, Cell.name == celltypes[f])
                genes = unique(tempdata$Marker)
                genes = intersect(genes, marker_genes)
                genes = genes[!is.na(genes)]
                if(length(genes) == 0)
                {
                    message(paste0("No Markers found for: ", celltypes[f]))
                }
                if(length(genes) > 0)
                {
                    dir_create_Celltypes = paste0(dir_create_Cellmarkers, celltypes[f], "/")
                    dir.create(dir_create_Celltypes)
                    check_status = intersect(genes, rownames(query))
                    len_genes = length(check_status)
                    if(len_genes > 10)
                    {
                        gene_sets_found = split(check_status, ceiling(seq_along(check_status) / 10))
                        for (gs in 1:length(gene_sets_found))
                        {
                            DotPlot(query, features = gene_sets_found[[gs]], group.by = "seurat_clusters") + RotatedAxis() + ggtitle(celltypes[f]) + theme_classic()
                            filename = paste0(dir_create_Celltypes, celltypes[f], " Set", gs, " MarkerGenes DotPlot.png")
                            ggsave(filename, width = 10, height = 10, dpi = 800)
                            FeaturePlot(query, features = gene_sets_found[[gs]], order = T, reduction = "umap")
                            filename = paste0(dir_create_Celltypes, celltypes[f], " Set", gs, " MarkerGenes FeaturePlot.png")
                            ggsave(filename, width = 10, height = 10, dpi = 800)
                        }
                    }
                    if(len_genes <= 10)
                     {
                        DotPlot(query, features = check_status, group.by = "seurat_clusters") + RotatedAxis() + ggtitle(celltypes[f]) + theme_classic()
                        filename = paste0(dir_create_Celltypes, celltypes[f], " MarkerGenes DotPlot.png")
                        ggsave(filename, width = 10, height = 10, dpi = 800)
                        FeaturePlot(query, features = check_status, order = T, reduction = "umap")
                        filename = paste0(dir_create_Celltypes, celltypes[f], " MarkerGenes FeaturePlot.png")
                        ggsave(filename, width = 10, height = 10, dpi = 800)
                    }
            
                }
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }

    }
    message('\nCellMarker DATABASE BASED VALIDATION COMPLETED')
}

#' ELeFHAnt CrossSpecies_Conversion
#'
#' CrossSpecies_Conversion is a function to convert gene names of a or more seurat objects from one species to other. For example: 
#' if your data comes from mouse and you want to compare to human then you can use this function to convert. Currently we support following species: human, mouse, rhesus, zebrafish, chicken and rat
#'
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @param seurat.objects one or more seurat objects
#' @param species_from current species: human, mouse, rhesus, zebrafish, chicken and rat
#' @param format_from which format genes are present in the object: ensembl or symbol
#' @param species_to convert to: human, mouse, rhesus, zebrafish, chicken and rat
#' @param format_to which format genes should be converted to: ensembl or symbol
#' @return Cross-species converted Seurat Object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'


CrossSpecies_Conversion <- function(seurat.objects = c(), species_from = NULL, species_to = NULL, format_from = NULL, format_to = NULL) {
  if (is.null(species_to)) {
    species_from = species_to
  }
  print(paste0("Converting genes from ", species_from, " ", format_from, " to ", species_to, " ", format_to))
  if (species_from == "human" | species_to == "human") {
        human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
    }
    if (species_from == "mouse" | species_to == "mouse") {
        mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
    }
    if (species_from == "rhesus" | species_to == "rhesus") {
        rhesus <- useMart('ensembl', dataset = 'mmulatta_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
    }
  if (species_from == "zebrafish" | species_to == "zebrafish") {
    zebrafish <- useMart('ensembl', dataset = 'drerio_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
  }
  if (species_from == "chicken" | species_to == "chicken") {
    chicken <- useMart('ensembl', dataset = 'ggallus_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
  }
  if (species_from == "rat" | species_to == "rat") {
    rat <- useMart('ensembl', dataset = 'rnorvegicus_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
  }
  if (format_from == "ensembl") {
    format_from = "Gene.stable.ID"
  } else if (format_from == "symbol") {
    format_from = "Gene.name"
  }
  if (format_to == "ensembl") {
    format_to = "Gene.stable.ID.1"
  } else if (format_to == "symbol") {
    format_to = "Gene.name.1"
  }
    gene_table <- getLDS(mart = get(species_from), attributes = c('ensembl_gene_id','external_gene_name'), martL = get(species_to), attributesL = c('ensembl_gene_id','external_gene_name'))
    gene_table <- gene_table[,c(format_from,format_to)]
    gene_table <- gene_table[gene_table[,format_to] != "",]
    gene_table <- gene_table[gene_table[,format_from] != "",]
    gene_table <- distinct_at(gene_table,format_to,.keep_all = T)
    gene_table <- distinct_at(gene_table,format_from,.keep_all = T)
    new.seurat.objects <- c()
    seurat.objects <- list(seurat.objects)
    for (seurat.object in seurat.objects) {
        DefaultAssay(seurat.object) <- "RNA"
        matches <- gene_table
        genes <- rownames(seurat.object)
        genes_start <- length(genes)
        print(paste0("Genes before: ", genes_start))
        seurat.object <- seurat.object[genes %in% matches[,format_from],]
        genes <- rownames(seurat.object)
        genes_end <- length(genes)
        print(paste0("Genes after: ", genes_end))
        matches <- matches[match(genes,matches[,format_from]),]
        new_genes <- matches[,format_to]
        seurat.object@assays$RNA@counts@Dimnames[[1]] <- new_genes
        seurat.object@assays$RNA@data@Dimnames[[1]] <- new_genes
        rownames(seurat.object@assays$RNA@meta.features) <- new_genes
        if (genes_start != genes_end) {
        seurat.object <- NormalizeData(seurat.object)
        seurat.object <- FindVariableFeatures(seurat.object)
        seurat.object <- ScaleData(seurat.object)
        seurat.object <- RunPCA(seurat.object)
        seurat.object <- FindNeighbors(seurat.object, dims = 1:20)
        seurat.object <- FindClusters(seurat.object, resolution = 0.8)
        seurat.object <- RunUMAP(seurat.object, dims = 1:20)
        }
        new.seurat.objects <- c(new.seurat.objects,seurat.object)
    }
    return(new.seurat.objects)
}

#' Benchmark ELeFHAnt against scPred and Label Transfer
#'
#' Celltype annotation is a function to annotate celltypes in a single cell datasets.
#' Benchmark ELeFHAnt is an easy to use function if the user wants to compare the predictions from ELeFHAnt to Label Transfer and scPred
#' 
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
#' @import ranger
#' @import LiblineaR
#' @import caTools
#' @import biomaRt
#' @import reshape2
#' @import readxl
#' @import Matrix
#' @param reference a processed Seurat Object with Celltypes column in metadata
#' @param query a processed seurat object with seurat_clusters column in metadata
#' @param downsample logical Indicator (TRUE or FALSE) to downsample reference, enabling fast computation.
#' @param downsample_to a numerical value > 1 to downsample cells [Default: 1000]
#' @param selectvarfeatures number of variable features to select while training (Default: 2000)
#' @param ntree number of trees randomForest classifier should build (Default: 500)
#' @param annotationCol Name of column in metdata which has Celltypes information in reference [Default: Celltypes]
#' @return query seurat object with predictions from ELeFHAnt, Label Transfer and scPred added to meta.data of the object
#' @export
#' @author Praneet Chaturvedi & Konrad Thorner
#'

BenchmarkELeFHAnt <- function(reference = reference, query = query, downsample = TRUE, downsample_to = 1000, selectvarfeatures = 2000, ntree = 500, annotationCol = "Celltypes") {

    if(downsample == TRUE)
    {
        message('\nDownsampling reference cells to enable fast computation\n')

        Idents(reference) = annotationCol
        reference_use = subset(reference, downsample = downsample_to)
        reference_use = NormalizeData(reference_use)
        reference_use = FindVariableFeatures(reference_use)
        reference_use = ScaleData(reference_use)
    }
    
    if(downsample == FALSE)
    {
        reference_use = reference
    }

    message('\nDeploying ELeFHAnt\n')

    out.ELeFHAnt = CelltypeAnnotation(reference = reference_use, query = query, downsample = downsample, downsample_to = downsample_to, validatePredictions = FALSE, selectvarfeatures = selectvarfeatures, ntree = ntree, annotationCol = annotationCol)
    p1 = DimPlot(out.ELeFHAnt, group.by = "ELeFHAnt_Ensemble_CelltypePrediction", label=T, repel = T, label.size = 6, reduction = "umap") + NoLegend() + ggtitle("ELeFHAnt Predictions")

    message('\nDeploying Seurat Label Transfer\n')

    anchors <- FindTransferAnchors(reference = reference_use, query = query)
    predictions <- TransferData(anchorset = anchors, refdata = reference_use@meta.data[,annotationCol])
    out.seurat <- AddMetaData(object = query, metadata = predictions)
    p2 = DimPlot(out.seurat, group.by = "predicted.id", label=T, repel = T, label.size = 6, reduction = "umap") + NoLegend() + ggtitle("LabelTransfer Predictions")

    message('\nDeploying scPred\n')
    reference.scPred <- reference_use %>%
        RunPCA() %>% 
        RunUMAP(dims = 1:30)

    reference.scPred <- getFeatureSpace(reference.scPred, annotationCol)
    reference.scPred <- trainModel(reference.scPred)
    get_probabilities(reference.scPred) %>% head()
    get_scpred(reference.scPred)
    out.scPred <- NormalizeData(query)
    out.scPred <- scPredict(out.scPred, reference.scPred)
    p3 = DimPlot(out.scPred, group.by = "scpred_prediction", label=T, repel = T, label.size = 6, reduction = "umap") + NoLegend() + ggtitle("scPred Predictions")

    plot = p1+p2+p3
    ggsave('Benchmarking_ELeFHAnt_LabelTransfer_scPred_CelltypeAnnotation_Plot.png', plot, width = 20, height = 20, dpi = 800)
    query$ELeFHAnt_Ensemble_CelltypePrediction = out.ELeFHAnt$ELeFHAnt_Ensemble_CelltypePrediction
    query$predicted.id = out.seurat$predicted.id
    query$scpred_prediction = out.scPred$scpred_prediction
    return(query)
}


randomForest_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL, ntree = NULL) {
    message('Initializing randomForest')
    rf_Celltypes <- ranger(factor(train_label) ~ ., data=train, num.trees = ntree, classification = TRUE, probability = TRUE)
    message('randomForest Complete')
    return(list(rf_Celltypes))
}

randomForest_predictor_withoutProb <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL, ntree = NULL) {
    message('Initializing randomForest')
    rf_Celltypes <- ranger(factor(train_label) ~ ., data=train, num.trees = ntree, classification = TRUE, probability = FALSE)
    message('randomForest Complete')
    return(list(rf_Celltypes))
}

svm_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL) {
    message('Initializing SVM')
    cost_est <- heuristicC(train)
    svm_Celltypes = LiblineaR(data=train,target=factor(train_label),type=5,cost=cost_est,verbose=TRUE)
    message('SVM Complete')
    return(list(svm_Celltypes))
}

lr_predictor <- function(train = NULL, test = NULL, train_label = NULL, test_label = NULL) {
    message('Initializing LR')
    cost_est <- heuristicC(train)
    lr_Celltypes = LiblineaR(data=train,target=factor(train_label),type=7,cost=cost_est,verbose=TRUE)
    message('LR Complete')
    return(list(lr_Celltypes))
}