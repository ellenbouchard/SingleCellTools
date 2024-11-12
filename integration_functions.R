# Functions for integration and clustering
# Author: Ellen Bouchard
# Date Created: 11/12/2024

# Function remove_scale_data: 
# Function to remove the scale.data layer from every object in a list
# Useful for prepping objects for integration and merging
# Input: a list of Seurat objects that have a scale.data layer in the RNA assay
# Output: a list of Seurat objects without scale.data layers
remove_scale_data <- function(object_list) {
    object_list <- lapply(object_list, function(obj) {
        obj[["RNA"]]$scale.data <- NULL
        return(obj)
        })
    return(object_list)
    }

# Function merge_object_list
# Merges all of the objects in a list into one Seurat object
# Input: a list of Seurat objects
# Output: one Seurat object that consists of all the individual objects, merged
merge_object_list <- function(object_list){
    merged_object <- merge(x = object_list[[1]], y = object_list[-1])
    return(merged_object)
    }

# Function norm_scale_pca
# Runs NormalizeData, FindVariableFeatures, ScaleData, and RunPCA 
# In preparation to integrate
# Arguments:
    # object : Seurat object to process
    # nfeatures: number of variable features to calculate, default 5000
    # vars.to.regress: variables to regress during scaling, default none
    # normalization.method: method to use while normalizing, default "LogNormalize"
    # selection.method : method used to calculate variable features, default "vst"
    # scale.factor: factor to normalize by, default 10,000
norm_scale_pca <- function(object,
                           nfeatures = 5000,
                           vars.to.regress = NULL,
                           normalization.method = "LogNormalize",
                           selection.method = "vst",
                           scale.factor = 10000) {
    object <- NormalizeData(object,
                            normalization.method = normalization.method, 
                            scale.factor = scale.factor,
                            verbose = FALSE)
    object <- FindVariableFeatures(object, nfeatures = nfeatures,
                                  verbose = FALSE)
    object <- ScaleData(object, vars.to.regress = vars.to.regress,
                       verbose = FALSE)
    object <- RunPCA(object, verbose = FALSE)
    return(object)
    }

# Function integrate_object : Integrate using specified method
# Default is Harmony integration
# Input: a Seurat object with separated layers
# Output: a Seurat object with an integrated dimensionality reduction 
# Arguments: 
    # object : Seurat object to integrate, must have separated layers
    # method : Integration method to use, default HarmonyIntegration
    # orig.reduction : dimensional reduction to use for integration, default "pca"
    # new.reduction : name of new integration reduction, default "harmony"
    # verbose: default FALSE
integrate_object <- function(object,
                             method = HarmonyIntegration,
                             orig.reduction = "pca",
                             new.reduction = "harmony") {
    integrated_obj <- IntegrateLayers(
                      object = object, method = method,
                      orig.reduction = orig.reduction, new.reduction = new.reduction,
                      verbose = FALSE)
    return(integrated_obj)
    }

# Function join_layers
# Used to join layers after integration
# Just a wrapper for JoinLayers
# Input: object : a Seurat object with separated layers
# Output: a Seurat object with joined layers
join_layers <- function(object) {
    object[["RNA"]] <- JoinLayers(object[["RNA"]])
    return(object)
    }

# Function cluster_umap
# Used to run FindNeighbors, FindClusters, and RunUMAP on a freshly integrated object
# Input: a Seurat object that was created by integrated multiple objects
# Oputput: a Seurat object that has been clustered and has a UMAP dimensionality reduction
# Arguments:
    # reduction: the integrated dimensionality reduction, default "harmony"
    # ndims: the number of dimensions to use when calculating UMAP, default 20
    # resolution: the resolution to cluster at, default = 0.5,
    # cluster.name: what to name resulting clusters, default "harmony_clusters"
    # reduction.name = name of newly added UMAP reduction, default "umap.harmony"
cluster_umap <- function(object,
                         reduction = "harmony",
                         ndims = 20,
                         resolution = 0.5,
                         cluster.name = "harmony_clusters",
                         reduction.name = "umap.harmony") {
    object <- FindNeighbors(object, reduction = reduction, dims = 1:ndims, verbose = FALSE)
    object <- FindClusters(object, resolution = resolution, cluster.name = cluster.name, verbose = FALSE)
    object <- RunUMAP(object, reduction = reduction, dims = 1:ndims, reduction.name = reduction.name, verbose = FALSE) 
    }

# Finally, Function integration_pipeline
# Uses the above functions to fully normalize, scale, merge, integrate, join layers, and cluster a list of Seurat objects
# Input: A list of Seurat objects (should be QC'd already)
# Output: A single Seurat object that is made from integrating the input list of Seurat objects
# Arguments:
    # object_list : list of Seurat objects to integrate
integration_pipeline <- function(object_list,
                                 nfeatures = 5000,
                                 vars.to.regress = NULL,
                                 normalization.method = "LogNormalize",
                                 selection.method = "vst",
                                 scale.factor = 10000,
                                 method = HarmonyIntegration,
                                 orig.reduction = "pca",
                                 new.reduction = "harmony",
                                 verbose = FALSE,
                                 reduction = "harmony",
                                 ndims = 20,
                                 resolution = 0.5,
                                 cluster.name = "harmony_clusters",
                                 reduction.name = "umap.harmony"){
    message("Removing scale.data layers")
    object_list <- remove_scale_data(object_list)
    message("Scale.data layers removed successfully")
    message("Merging objects")
    merged_object <- merge_object_list(object_list)
    message("Objects merged successfully")
    message("Processing merged object")
    processed_object <- norm_scale_pca(merged_object, 
                           nfeatures = nfeatures,
                           vars.to.regress = vars.to.regress,
                           normalization.method = normalization.method,
                           selection.method = selection.method,
                           scale.factor = scale.factor)
    message("Merged object processed successfully")
    message("Integrating object")
    integrated_object <- integrate_object(processed_object,
                             method = method,
                             orig.reduction = orig.reduction,
                             new.reduction = new.reduction,
                             verbose = verbose)
    message("Object integrated successfully")
    message("Joining layers")
    joined_object <- join_layers(integrated_object)
    message("Layers joined successfully")
    message("Clustering final object")
    final_object <- cluster_umap(joined_object,
                             reduction = reduction,
                             ndims = ndims,
                             resolution = resolution,
                             cluster.name = cluster.name,
                             reduction.name = reduction.name)
    message("Final object clustered successfully, pipeline finished")
    return(final_object)
    }
