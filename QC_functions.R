# Author: Ellen Bouchard
# Date Created: November 11 2024
# Functions for applying quality-control pipeline to a set of Seurat objects

# Functions included:
# process_objects : applies NormalizeData, ScaleData, FindVariableFeatures,
# 	RunPCA, RunUMAP, FindNeighbors, and FindClusters to a list of Seurat objects
# add_pct_mt_rb : adds columns in metadata for percent mitochondrial genes, percent rps genes, percent rpl genes, and overall ribosomal gene percentage
# filter_mt : filters objects to remove cells with high mitochondrial gene percentage
# filter_umi : filters objeccts to remove cells with low or high UMI count


# Function for normalizing, scaling, and clustering a list of seurat objects
# Input: obj_list, a list of Seurat objects
# Output: list of Seurat object that have been normalized, scaled, calculated for variable features,
    # calculated for PCA, calculated for UMAP, and clustered
# Allows for specificaiton of the following variables: 
    # nfeatures = number of variable features to calculate, default 5,000
    # ndims = number of dimensions to use for UMAP, default 10
    # resolution = clustering resolution, default 0.3
    # normalization.method = method used for data normalization, default "LogNormalize"
    # vars.to.regress = variables to regress out during scaling, default none
    # scale.factor = scale factor for cell-level normalization, default 10,000
    # selection.method = method for finding top variable features, default "vst"
process_objects <- function(obj_list, 
                            nfeatures = 5000, 
                            ndims = 10, 
                            resolution = 0.3, 
                            vars.to.regress = NULL,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000, 
                            selection.method = "vst") {
    obj_list <- lapply(obj_list, function(obj) { 
        obj <- NormalizeData(obj, 
                             normalization.method = normalization.method, 
                             scale.factor = scale.factor,
			     verbose = FALSE)
        obj <- ScaleData(obj, vars.to.regress = vars.to.regress,
			verbose = FALSE)
        obj <- FindVariableFeatures(obj, nfeatures = nfeatures,
			verbose = FALSE)
        obj <- RunPCA(obj, verbose = FALSE)
        obj <- RunUMAP(obj, dims = 1:ndims, verbose = FALSE)
        obj <- FindNeighbors(obj, verbose = FALSE)
        obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
        return(obj)
    })
    return(obj_list)
}

# Function for adding percent mitochondrial content and percent ribosomal content as metadata
# Input: a list of Seurat objects
# Output: a list of Seurat objects with columns in the metadata for percent mitochondrial gene percentage,
    # percent RPS gene percentage, percent RPL gene percentage, and overall ribosomal gene percentage
# Arguments: mt.name = name of metadata column to add that has mitochondrial gene percentage, default "percent.mt"
    # rps.name = name of metadata column to add that has RPS gene percentage, default "percent.rps"
    # rpl.name = name of metadata column to add that has RPL gene percentage, default "percent.rpl"
    # rb.name = name of metadata column to add that has ribosomal gene percentage, default "percent.rb"
    # species = either "human" or "mouse", default "human" 
# Function calculates percent of genes per cell that start with "MT-", "RPS-" or "RPL-" for human;
    # "mt-", "rps-", "rpl-" for mouse
# If any calculations result in all cells being "zero," will warn you that you may have chosen wrong species
add_pct_mt_rb <- function(obj_list,
                          species = "human",
                          mt.name = "percent.mt",
                          rps.name = "percent.rps",
                          rpl.name = "percent.rpl",
                          rb.name = "percent.rb") {
    if(!species %in% c("human","mouse")) {stop("Species must be either `human` or `mouse`")}
    pattern_list <- switch(
                        species,
                          "human" = c("^MT-", "^RPS", "^RPL"),
                          "mouse" = c("^mt-", "^Rps", "^Rpl"))
    obj_list <- lapply(obj_list, function(obj) {
        obj[[mt.name]] <- PercentageFeatureSet(obj, pattern = pattern_list[1])
        obj[[rps.name]] <- PercentageFeatureSet(obj, pattern = pattern_list[2])
        obj[[rpl.name]] <- PercentageFeatureSet(obj, pattern = pattern_list[3])
        obj[[rb.name]] <- obj[[rps.name]] + obj[[rpl.name]]
        if(all(obj$percent.mt == 0) == TRUE) {warning("All cells have a mitochondrial gene percentage of zero (0). 
                                                      Did you choose the correct species? Options are `human` and `mouse`.")}
        return(obj)
    })
    
    return(obj_list)
}


# Function for filtering list of Seurat objects by mitochondrial gene percentage
# Function first analyzes individual objects by cluster and removes clusters for which more than
    # a specified fraction of cells are at or above the acceptable threshold of mitochondrial gene percentage,
    # then removes all remaining cells that are at or above the acceptable threshold of mitochondrial gene percentage
    # (For example, if 90% of cells in a cluster have >5% mitochondrial gene content, the entire cluster is removed)
# Function prints four plots for each object:
    # UMAP, pre-filtering
    # Violin plot of mitochondrial gene percentage per cluster, pre-filtering
    # UMAP, post-filtering
    # Violin plot of mitochondrial gene percentage per cluster, post-filtering
# Input: a list of Seurat objects with a column in the metadata for the percent of 
    # genes in a cell that are mitochondrial genes (starting with "MT-" or "mt-")
# Output: a list of Seurat objects that have been filtered based on specified mitochondrial gene percentage 
    # threshold and cell majority threshold
# Arguments: 
    # obj_list = a list of Seurat objects 
    # mt.name = name of column in metadata that contains percent mitochondrial gene percent, default percent.mt
    # mt.pct.threshold = upper limit of mitochondrial gene percentage at or above which point cells are removed by filtering
            # default 5
    # majority.threshold = if a cluster contains this fraction of cells or more that are at or above the mt.pct.threshold value,
            # the entire cluster is removed by filtering
            # default 0.7, or 70%
filter_mt <- function(obj_list,
                      mt.name = "percent.mt",
                      mt.pct.threshold = 5,
                      majority.threshold = .7) {
     
    obj_list <- lapply(obj_list, function(obj){
        # For object name, use orig.ident
        obj_name <- obj$orig.ident[[1]]
        
        # Plot object pre-filtering
        p1 <- DimPlot(obj, reduction = "umap", label = TRUE) + NoLegend() + ggtitle(paste(obj_name, ": Pre-Filtering",sep=""))
        p2 <- VlnPlot(obj, features = c(mt.name),alpha = 0.3) + geom_hline(yintercept = 5) + ggtitle(paste(obj_name, ": Pre-Filtering",sep=""))
        print(plot_grid(p1, p2, ncol = 2))
        
        # Calculate how many cells are in each cluster
        cells_clusters <- as.data.frame(table(obj$seurat_clusters))
        colnames(cells_clusters) <- c("seurat_clusters","cells_per_cluster")

        # Count cells with percent.mt above threshold for each cluster, add as column
        mt_above_threshold <- as.data.frame(table(
        obj[, obj@meta.data[[mt.name]] > mt.pct.threshold]$seurat_clusters
        ))
        colnames(mt_above_threshold) <- c("seurat_clusters", "mt_above_threshold")
        cells_clusters <- merge(cells_clusters, mt_above_threshold, by = "seurat_clusters", all.x = TRUE)
        cells_clusters$mt_above_threshold[is.na(cells_clusters$mt_above_threshold)] <- 0
    
        # Calculate % of cells in each cluster that are above threshold
        cells_clusters$pct_above_threshold <- cells_clusters$mt_above_threshold/cells_clusters$cells_per_cluster

        # If % of cells is above the majority cutoff, remove the cluster from the object
        cells_clusters$above_majority_cutoff <- cells_clusters$pct_above_threshold > majority.threshold
        clusters_to_remove <- cells_clusters$seurat_clusters[cells_clusters$above_majority_cutoff == TRUE]
        obj_2 <- subset(x = obj, idents = clusters_to_remove, invert = TRUE)

        # Remove remaining cells over mt % threshold
        obj_2 <- obj_2[, obj_2@meta.data[[mt.name]] < mt.pct.threshold]

        # Plot object post-filtering
        p3 <- DimPlot(obj_2, reduction = "umap", label = TRUE) + NoLegend()+ ggtitle(paste(obj_name, ": Post-Filtering",sep=""))
        p4 <- VlnPlot(obj_2, features = c(mt.name), alpha = 0.3) + geom_hline(yintercept = 5)+ ggtitle(paste(obj_name, ": Post-Filtering",sep=""))
        print(plot_grid(p3, p4, ncol = 2))

        return(obj_2)
        })
    return(obj_list)
}



# Function for filtering a list of Seurat objects based on UMI count
# Function does the following:
    # Analyzes individual objects to calculate fraction of cells per cluster that are at or below the umi.low.threshold. 
    # Removes any cluster that contains a fraction of cells at or above the majority.threshold value that are at or below the umi.low.threshold value
        # (For example, if a cluster contains 90% of cells that are below 1,000 UMI count, the whole cluster is removed) 
    # Removes any remaining cells that are at or below the umi.low.threshold
    # Calculates top percentile of cells by UMI count
    # Removes any cells that are at or above the top percentile threshold (in order to remove possible doublets)
    # Prints four plots per object:
        # UMAP, pre-filtering
        # Violing plot of UMI count per cluster, pre-filtering
        # UMAP, post-filtering
        # Violing plot of UMI count per cluster, post-filtering
# Input: a list of Seurat objects that have a column in the metadata containing UMI counts
# Output: a list of Seurat objects that have been filtered to remove cells with low UMI counts or high UMI counts
# Arguments:
    # obj_list = a list of Seurat objects
    # umi.name = name of the column in the metadata that contains UMI count values, default "nCount_RNA"
    # umi.low.threshold = threshold of UMI count at or under which cells are removed, default 1,000
    # majority.threshold = if a cluster contains a fraction of cells below the umi.low.threshold that is at or above this number, the whole cluster is removed, default 0.7 or 70%
    # umi.percentile = the percentile of UMI count at or above which cells are removed, default 0.99

filter_umi <- function(obj_list,
                       umi.name = "nCount_RNA",
                       umi.low.threshold = 1000,
                       majority.threshold = .7,
                       umi.percentile = .99) {
    obj_list <- lapply(obj_list, function(obj){
       
        obj_name <- obj$orig.ident[[1]]
        
        # Plot object pre-filtering
        p1 <- DimPlot(obj, reduction = "umap", label = TRUE) + NoLegend() + ggtitle(paste(obj_name, ": Pre-Filtering",sep=""))
        p2 <- VlnPlot(obj, features = c(umi.name),alpha = 0.3) + geom_hline(yintercept = 5) + ggtitle(paste(obj_name, ": Pre-Filtering",sep=""))
        print(plot_grid(p1, p2, ncol = 2))

        # Calculate how many cells are in each cluster
        cells_clusters <- as.data.frame(table(obj$seurat_clusters))
        colnames(cells_clusters) <- c("seurat_clusters","cells_per_cluster")

        # Count cells with nCount_UMI below threshold for each cluster, add as column
        umi_below_threshold <- as.data.frame(table(
        obj[, obj@meta.data[[umi.name]] < umi.low.threshold]$seurat_clusters
        ))
        colnames(umi_below_threshold) <- c("seurat_clusters", "umi_below_threshold")
        cells_clusters <- merge(cells_clusters, umi_below_threshold, by = "seurat_clusters", all.x = TRUE)
        cells_clusters$umi_below_threshold[is.na(cells_clusters$umi_below_threshold)] <- 0

        # Calculate % of cells in each cluster that are above threshold
        cells_clusters$pct_below_threshold <- cells_clusters$umi_below_threshold/cells_clusters$cells_per_cluster

        # If % of cells is above the majority cutoff, remove the cluster from the object
        cells_clusters$above_majority_cutoff <- cells_clusters$pct_below_threshold > majority.threshold
        clusters_to_remove <- cells_clusters$seurat_clusters[cells_clusters$above_majority_cutoff == TRUE]
        obj_2 <- subset(x = obj, idents = clusters_to_remove, invert = TRUE)

        # Remove remaining cells under umi threshold 
        obj_2 <- obj_2[, obj_2@meta.data[[umi.name]] > umi.low.threshold]

        # Remove cells in top percentile of UMI count
        top_percentile <- quantile(obj_2$nCount_RNA, probs = umi.percentile)
        obj_2 <- obj_2[, obj_2@meta.data[[umi.name]] < top_percentile]

        # Plot object post-filtering
        p3 <- DimPlot(obj_2, reduction = "umap", label = TRUE) + NoLegend()+ ggtitle(paste(obj_name, ": Post-Filtering",sep=""))
        p4 <- VlnPlot(obj_2, features = c(umi.name), alpha = 0.3) + geom_hline(yintercept = 5)+ ggtitle(paste(obj_name, ": Post-Filtering",sep=""))
        print(plot_grid(p3, p4, ncol = 2))

        return(obj)
        })
        
    return(obj_list)
}

# Function apply_doubletfinder
# Function runs DoubletFinder on a list of Seurat objects
# It is recommended that, if using DoubletFinder, to not remove a top percentile of cells by UMI count,
    # Or at least assess results before removing a top percentile by UMI count
# Remember, DoubletFinder works best on datasets of non-homogenous cell types 
# Input: a list of Seurat objects
# Output: a list of Seurat objects which have been filtered such that cells labeled as doublets by DoubletFinder have been removed
# Arguments:
    # object_list : list of Seurat objects
    # nPCs : number of PCs to use, default 30
    # pN : pN value to use when running DoubletFinder, default 0.25
apply_doubletfinder <- function(object_list,
                               nPCs = 30,
                               pN = 0.25) {
    require(DoubletFinder)
    object_list <- lapply(object_list, function(object){
        object_name <- object$orig.ident[[1]]
        # DoubletFinder pipeline, no annotations
        sweep.res <- paramSweep(object, PCs = 1:nPCs, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        pK.set <- unique(sweep.stats$pK)[2]
        nExp_poi <- round(0.08*nrow(object@meta.data))
        object <- doubletFinder(object, 
                       PCs = 1:nPCs, 
                       pN = pN, 
                       pK = as.numeric(as.character(pK.set)),
                       nExp = nExp_poi, 
                       reuse.pANN = FALSE, 
                       sct = FALSE)
        # DoubletFinder gives a unqiue name to the classifications column in the metadata.
        # For consistency, copy classifications to a new column with the name "DF.classifications"
        object$DF.classifications <- object[[]][ncol(object[[]])] 
        object_filt <- subset(object, DF.classifications == "Singlet")

        # Print number of cells that were removed
        n_cells_pre <- ncol(object)
        n_cells_post <- ncol(object_filt)
        n_cells_removed <- n_cells_pre - n_cells_post
        print(object_name)
        print("Number of cells, pre-filtering:")
        print(n_cells_pre)
        print("Number of cells, post-filtering:")
        print(n_cells_post)
        print("Number of cells removed:")
        print(n_cells_removed)
        return(object_filt)
        })
    return(object_list)
}




# Function QC_objects
# Uses the above functions (add_pct_mt_rb, process_object, filter_mt, and filter_umi) 
    # to apply full QC pipeline to a list of seurat objects
# Input: A list of Seurat objects
# Output: A list of Seurat objects that has been filtered such that:
    # Cells with high mitochondrial gene content have been removed
    # Cells wtih low UMI count have been removed 
    # Cells with high UMI count (possible doublets) have been removed
# Arguments:
    # object_list : a list of Seurat objects
QC_objects <- function(object_list,
                      species = "human",
                      mt.name = "percent.mt",
                      rps.name = "percent.rps",
                      rpl.name = "percent.rpl",
                      rb.name = "percent.rb",
                      nfeatures = 5000, 
                      ndims = 10, 
                      resolution = 0.3, 
                      vars.to.regress = NULL,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000, 
                      selection.method = "vst",
                      mt.pct.threshold = 5,
                      mt.majority.threshold = .7,
                      umi.name = "nCount_RNA",
                      umi.low.threshold = 1000,
                      umi.majority.threshold = .7,
                      umi.percentile = .99,
		      doubletfinder = FALSE,
                      nPCs = 30,
                      pN = 0.25
                      )
    {
    if(doubletfinder == TRUE && umi.percentile < 1) {
        warning("You have chosen to run DoubletFinder to remove doublets AND to remove the top percentile of cells by UMI count. This may result in over-filtering. Double check to make sure this is what you want.")
        }
    message("Adding mitochondrial and ribosomal gene content to metadata")
    object_list <- add_pct_mt_rb(object_list,
                                 species = species,
                                 mt.name = mt.name,
                                 rps.name = rps.name,
                                 rpl.name = rpl.name,
                                 rb.name = rb.name)
    message("Successfully added mitochondrial and ribosomal gene content")
    message("Normalizing, scaling, finding dimensionality reduction, and clustering objects")
    object_list <- process_objects(object_list,
                                   nfeatures = nfeatures,
                                   ndims = ndims,
                                   resolution = resolution,
                                   vars.to.regress = vars.to.regress,
                                   normalization.method = normalization.method,
                                   scale.factor = scale.factor,
                                   selection.method = selection.method)
    message("Successfully processed objects")
    message("Filtering cells based on mitochondrial gene content")
    object_list <- filter_mt(object_list,
                            mt.name = mt.name,
                            mt.pct.threshold = mt.pct.threshold,
                            majority.threshold = mt.majority.threshold)
    message("Successfully filtered based on mitochondrial gene content")
    message("Re-processing objects")
    object_list <- process_objects(object_list,
                                   nfeatures = nfeatures,
                                   ndims = ndims,
                                   resolution = resolution,
                                   vars.to.regress = vars.to.regress,
                                   normalization.method = normalization.method,
                                   scale.factor = scale.factor,
                                   selection.method = selection.method)
    message("Filtering cells based on UMI count")
    object_list <- filter_umi(object_list,
                              umi.name = umi.name,
                              umi.low.threshold = umi.low.threshold,
                              majority.threshold = umi.majority.threshold,
                              umi.percentile = umi.percentile)
    message("Successfully filtered based on UMI count")
    if(doubletfinder == TRUE) {
        message("Applying DoubletFinder to remove doublets")
        object_list <- apply_doubletfinder(object_list,
                                           nPCs = nPCs,
                                           pN = pN)
        message("Successfully removed doublets")
     }
    message("Re-processing objects")
    object_list <- process_objects(object_list,
                                   nfeatures = nfeatures,
                                   ndims = ndims,
                                   resolution = resolution,
                                   vars.to.regress = vars.to.regress,
                                   normalization.method = normalization.method,
                                   scale.factor = scale.factor,
                                   selection.method = selection.method)
    message("Finished quality control of objects")
    return(object_list)
    }





