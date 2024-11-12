# Author: Ellen Bouchard
# Date Created: October 28 2024

# This file contains graphing functions for commonly used plots : 
# Horizontal Heatmap
# Volcano Plot

#Function for making volcano plots:

make_volcano = function( 
    de, 
    log_fc_cutoff = NULL, 
    percentile_fc_cutoff = NULL, 
    p_val_cutoff = NULL,  
    graph_title = 'add a title idiot', 
    overlap_metric = 17, 
    label_genes = TRUE, 
    label_specific_genes = FALSE,
    genes_to_label = c(),
    upcolor = 'Yellow', 
    downcolor = 'Blue', 
    midcolor = 'Gray', 
    labelcolor = 'Black', 
    remove_genes = c(),
    remove_labels = c(),
    overlap_metric_secondary = Inf, 
    visible_cutoffs = TRUE) {
    
    require(Seurat)
    require(ggplot2)
    require(ggrepel)
    if(!("gene" %in% colnames(de))) { print("Differential Expression Dataframe Must Have `gene` Column!")}
    
    de = subset(de,!(gene %in% remove_genes))
    
    up_de <- subset(de, avg_log2FC > 0)
    down_de <- subset(de, avg_log2FC < 0)

    if(is.null(percentile_fc_cutoff)){
        up_quant <- log_fc_cutoff
        down_quant <- -1 * log_fc_cutoff
        
    } else {
        up_quant <- quantile(up_de$avg_log2FC, c(percentile_fc_cutoff,1 - percentile_fc_cutoff))
        down_quant <- quantile(down_de$avg_log2FC, c(percentile_fc_cutoff,1 - percentilefc_cutoff)) 
        up_quant <- up_quant[[2]]
        down_quant <- down_quant[[1]]
    }

    de$reg = ""
    de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC > up_quant & de$avg_log2FC > 0] <- "UP"
    de$reg[de$p_val_adj < p_val_cutoff & de$avg_log2FC < down_quant & de$avg_log2FC < 0] <- "DOWN"
    de$name = de$gene
    de$name[de$reg == ""] <- ""
    de$name[de$name %in% remove_labels] <- ""
    if(label_specific_genes) {de <- de %>% mutate(name_specific = ifelse(gene  %in% genes_to_label & reg != "",gene, ""))}    

    min_p_val <- min(de$p_val_adj[de$p_val_adj > 0])
    de$p_val_adj[de$p_val_adj == 0] <- min_p_val*0.1
    
 
    plot = ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=reg, label=name)) + 
        geom_point(color = 'black', size = 2.5) + 
        theme_minimal() +
        scale_color_manual(breaks = c("DOWN", "", "UP"),values=c(downcolor, midcolor, upcolor)) + 
        geom_point(size = 1) + 
        ggtitle(graph_title) +
	theme(panel.grid = element_blank()) 
    
        if (label_genes) {plot = plot + ggrepel::geom_text_repel(aes(label = name), color = labelcolor, max.overlaps = overlap_metric, nudge_y = 1)}
        if (label_specific_genes) {plot = plot + geom_text_repel(aes(label = name_specific), color = labelcolor, max.overlaps = overlap_metric_secondary)}
	if (visible_cutoffs) {plot = plot + geom_vline(xintercept = up_quant, linetype=3) + geom_vline(xintercept = down_quant, linetype = 3) + geom_hline(yintercept = -log10(p_val_cutoff), linetype = 3)}
    
        return(plot)
    

}

#Function for making heatmap, from Dillon, altered to generate horizontal heatmaps:
horizontal_scheatmap = function(object, features = NULL, n_mostvar_feats = 50, average_within = FALSE, 
assay = 'RNA', slot = 'data', scaling_method = 'z_score', group.by = 'seurat_clusters', cellheight = NULL, cellwidth = NULL,
ident.order = unique(object[[group.by]][[group.by]]), show_row_names = TRUE, show_column_names = TRUE, 
gaps_after_which_samples = NULL, gaps_after_which_features = NULL, cluster_rows = TRUE, cluster_columns = FALSE, 
feature_annotations = NULL, sample_annotations = NULL, feature_annotation_colors = NULL, sample_annotation_colors = NULL, 
mode = 'box', colorvec = NULL, valuevec = c(-2, -1.5, 0, 1.5, 2), column_title = NULL, row_title = NULL, 
column_names_side = 'bottom', row_names_side = 'left', column_dend_side = 'top', row_dend_side = 'left', gridline_color = 'black', pdfsavepath = NULL, scale_bar_title = ' ',
column_colorvector = NULL, row_colorvector = NULL, rasterize = NULL, pctassay = 'RNA', pctslot = 'counts', dot_size_factor = 0.9, dotplot_outline_color = NA, adjust_for_viewport = TRUE) {
    
    # get viewport sizes
    
    ops = options()
    viewport_width = ops$repr.plot.width
    viewport_height = ops$repr.plot.height
    
    if (adjust_for_viewport) {
    viewport_width_adjust_factor = viewport_width/min(c(viewport_width, viewport_height))
    viewport_height_adjust_factor = viewport_height/min(c(viewport_width, viewport_height))
    } else {
        viewport_width_adjust_factor = 1
        viewport_height_adjust_factor = 1
    }
    
    print('----------')
    print(viewport_height_adjust_factor)
    print(viewport_width_adjust_factor)
    print('----------')
    # web sources: 
    # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
    
    # import packages
    suppressPackageStartupMessages({
    require(Seurat)
    require(ComplexHeatmap)
    require(circlize)
    require(viridis)
    })
    
    # perform initial checks
    
    if (!(class(object) == 'Seurat')) {
        stop('object is not of class "Seurat"')
    }
    if (!(assay %in% names(object@assays))) {
        stop(paste(assay, 'not found in assays of object'))
    }
    #if (!(slot %in% names(attributes(object@assays[[assay]])))) {
    #    stop(paste(slot, 'not found in slots of object'))
    #}
    
    print('available modes are dot and box')
    if (!(mode %in% c('dot', 'box'))) {
        stop('mode must be either dot or box')
    }
    print('available scaling methods are z_score, min_max, divbymax and none')
    if (!(scaling_method %in% c('z_score', 'min_max', 'none', 'divbymax'))) {
        stop('mode must be one of z_score, min_max, divbymax, or none')
    }
    
    
    if (mode == 'dot' & average_within == FALSE) {
        stop('dot mode is not usable on single cell-resolved data, please change teh average_within argument to TRUE')
    }
    
    if (length(colorvec) != length(valuevec) & !is.null(colorvec)) {
        stop('each element of the colorvec argument must have a corresponding element in the valuevec argument')
    }
    
    # generate necessary tables
    
    metadata = object@meta.data
    
    
    # if no features are provided, the top n most variable genes are assigned as features to be plotted
    
    if (is.null(features)) {
        print(paste('selecting the top', n_mostvar_feats, 'most variable genes to plot'))
        object = FindVariableFeatures(object = object, nfeatures = n_mostvar_feats, assay = assay)
        features = VariableFeatures(object, assay = assay)
    }
    
    # reorder idents of object prior to data extraction
   
    object[[group.by]][[1]] <- factor(x = object[[group.by]][[1]], levels = ident.order)
    
    # create expression matrix
    if (average_within) {
        
        data_mat = AverageExpression(object = object, assays = assay, layer = slot, group.by = group.by, features = features)[[assay]]
        data_mat = as.matrix(data_mat)
        data_mat = data_mat[features,,drop=FALSE]
        
    }else{
        data_mat = as.matrix(GetAssayData(object = object, assay = assay, slot = slot))
        data_mat = data_mat[features,]
    }  
    
    # create legend list
    
    lgd_list = list()
    
    
    if (mode == 'dot') {
        # this is all done assuming the data has been averaged within groups
        groups_to_pct = colnames(data_mat)
        features_to_pct = rownames(data_mat)
        pct_mat = as.matrix(GetAssayData(object = object, assay = pctassay, slot = pctslot))[features_to_pct,]
        
        pctlist = list()
        for (group_to_pct in groups_to_pct) {
            cells_of_group = rownames(metadata)[metadata[,group.by] == group_to_pct]
            pct_mat_group = pct_mat[,cells_of_group]
            pct_mat_group = pct_mat_group > 0
            pcts = rowSums(pct_mat_group)/ncol(pct_mat_group)
            pctlist[[group_to_pct]] = pcts
        }
        pctlist = as.data.frame(pctlist)
        rect_gp_fun = gpar(col = gridline_color, lwd = 0.5, type = 'none')
        
        # cell function for dotplot
        
        cell_function = function(j, i, x, y, width, height, fill) {
        
        
        
        
        #here
        grid.rect(x = x, y = y, width = width, height = height, 
            gp = gpar(col = gridline_color, fill = NA))
        if (is.null(cellheight)) {
            grid.circle(x = x, y = y, r = ((pctlist[i, j]) * 0.5 * min(unit.c(width * viewport_width_adjust_factor, height * viewport_height_adjust_factor))) * dot_size_factor, 
                    gp = gpar(fill = col_fun(data_mat[i, j]), col = dotplot_outline_color))
            
            }else{
            grid.circle(x = x, y = y, r = ((pctlist[i, j]) * 0.5 * min(unit.c(unit(cellwidth, 'mm'), unit(cellheight, 'mm')))) * dot_size_factor,
                    gp = gpar(fill = col_fun(data_mat[i, j]), col = dotplot_outline_color))
            }
            
            
            
        }
        
        if (is.null(cellheight)) {
            leg_unit_size = min(unit.c(unit(1/nrow(data_mat), "npc"), unit.c(unit(1/ncol(data_mat), "npc"))))
            print(leg_unit_size)
            leg_unit_size = unit(5, 'mm')
            lgd_list[['pct_legend']] = ComplexHeatmap::Legend(labels = c(0.25,0.75), title = "Percent Expressing",graphics = list(
                function(x, y, w, h) {
                    grid.circle(x = x, y = y, r = 0.25 * 0.5 * 0.9 * leg_unit_size,
                                                              gp = gpar(fill = "black"))},
                             function(x, y, w, h) {grid.circle(x = x, y = y, r = 0.75 * 0.5 * 0.9 * leg_unit_size,
                                                              gp = gpar(fill = "black"))})
    )
        }else{
        
        
        }
        
        
    }else{
      
        cell_function = NULL
        rect_gp_fun = gpar(col = gridline_color, lwd = 0.5, type = '1')
    }
    
    
    
    
    
    # check for presence of features
    if (!(all(features %in% rownames(data_mat)))) {
        stop('some features were not found in the appropriate assay/slot combination of the object')
    }
    
    
    # apply z score transformation
    
    
    
    if (scaling_method == 'z_score') {
        data_mat = t(scale(x = t(data_mat), center = TRUE, scale = TRUE))
    } else if (scaling_method == 'min_max') {
        row_mins = apply(X = data_mat, FUN = min, MARGIN = 1)
        row_maxs = apply(X = data_mat, FUN = max, MARGIN = 1)
        row_ranges = row_maxs - row_mins
        data_mat = (data_mat - row_mins)/row_ranges
        if (!identical(range(valuevec), c(0, 1))) {
             warning('It is reccomended that, when min-max scaling, one uses a color range going from 0 to 1, which is not the case')
        }
    } else if (scaling_method  == 'divbymax') {
            if (!identical(range(valuevec), c(0, 1))) {
             warning('It is reccomended that, when divbymax scaling, one uses a color range going from 0 to 1, which is not the case')
            }
        row_maxs = apply(X = data_mat, FUN = max, MARGIN = 1)
        data_mat = data_mat/row_maxs
    }else if (scaling_method == 'none') {
        print('no scaling method applied')
    }
    
    
    # if no color vec is set, create a default color vector that matches the length of the valuevec
    
    if (is.null(colorvec)) {
        colorvec = viridis::plasma(length(valuevec))
    }
    
    
    # create color scale (note, I an set -2 and 2 to the min and max of the matrix to change the behavior of color scaling)
    
    col_fun = colorRamp2(valuevec, colorvec, space = 'sRGB')
    
    
    # prepare metrics for image sizing (because complexheatmap is stupid)
    data_mat = t(data_mat)
    nr = nrow(data_mat)
    nc = ncol(data_mat)
    
    
    nc_adjust = nc/min(c(nc, nr))
    nr_adjust = nr/min(c(nc, nr))
    
    
    # assign row and column color vectors
    
    if (is.null(row_colorvector)) {row_colorvector = c(rep('black', nr))}
    if (is.null(column_colorvector)) {column_colorvector = c(rep('black', nc))}
    
    # some logic controlling cell height. Uses defaults if not supplied
    
    if (is.null(cellheight)) {height_measure = unit(nr, "null")} else {height_measure = nr*unit(cellheight, "mm")}
    if (is.null(cellwidth)) {width_measure = unit(nc, "null")} else {width_measure = nc*unit(cellwidth, "mm")}
    
            
    print('----------')
    print(height_measure)
    print(width_measure)
    print('----------')
    
    
    # remove gridline colors if not averaging (this is probably always )
    
    if (average_within == FALSE && !(is.na(gridline_color))) {
        warning('it is reccomended that, if not averaging within idents, to set gridline color to NA')
    }
    
    # rasterization logic
    
    if (is.null(rasterize)) {
        if (nr > 2000 | nc > 2000) {
            warning('rasterization is automatically set to true if matrix rows or columns exceed 2000, manually set if needed')
            rasterize = TRUE
        }else {
            rasterize = FALSE
        }
    }
    
    # create basic heatmap 
    heatmap_obj = ComplexHeatmap::Heatmap(matrix = data_mat, col = col_fun, na_col = 'black', border = TRUE, 
    column_title = column_title, row_title = row_title, cluster_rows = cluster_rows, cluster_columns = cluster_columns,column_names_side = column_names_side, row_names_side = row_names_side, 
    column_dend_side = column_dend_side, row_dend_side = row_dend_side, show_column_names = show_column_names, show_row_names = show_row_names, rect_gp = rect_gp_fun, 
    border_gp = gpar(col = "black", lty = 1, lwd = 3), name = scale_bar_title, row_names_gp = gpar(col = row_colorvector), 
    column_names_gp = gpar(col = column_colorvector), width = width_measure, height = height_measure, use_raster = rasterize, cell_fun = cell_function
)      
    
    
    # print(heatmap_obj)
    # print(lgd_list)
    return(draw(heatmap_obj, heatmap_legend_list = lgd_list))
#     if (!is.null(pdfsavepath)) {
#         pdf(pdfsavepath, width = , height = , bg = bg)
#         draw(heatmap_obj)
#         dev.off()
#     }   
#     return(heatmap_obj)
    
}


## make_DA_graph
## Simple bar  graph function for differential abundance of two variables in a Seurat object
# varX = the variable on the X axis
# varY = the variable used for the fill aesthetic of the bars
# Y axis is normalized proportion of varY per each varX group
# varX and varY must be names of columns in obj metadata
# Returns ggplot2 bar chart

make_DA_graph <- function(obj, varX, varY) {
     X_Y_counts <- as.data.frame(table(obj@meta.data[[varX]], obj@meta.data[[varY]]))
     colnames(X_Y_counts) <- c(varX, varY, "nCells")
     varYdf <- as.data.frame(table(obj@meta.data[[varY]]))
     colnames(varYdf) <- c(varY, "TotalCells")
     X_Y_counts <- left_join(X_Y_counts, varYdf, by = varY)
     X_Y_counts$Normalized <- X_Y_counts$nCells / X_Y_counts$TotalCells
    
    varX_sym <- sym(varX)
    varY_sym <- sym(varY)

    plot <- ggplot(X_Y_counts, aes(x = !!varX_sym, y = Normalized, fill = !!varY_sym)) +
            geom_bar(position = "fill", stat = "identity") +
            theme_minimal() +
            labs(x = varX, y = "Normalized Proportion", fill = varY) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    return(plot)
}

# Function to perform chi-square test, used for DA graphing function below
perform_chi_square_test <- function(a_count, b_count) {
  observed <- c(a_count, b_count)
  expected <- sum(observed) / 2
  chisq.test(c(a_count, b_count), p = c(0.5, 0.5))$p.value
}

## make_DA_graph_sig
## Bar graph function for differential abundance of two variables in a Seurat object,
# varX = the variable on the X axis
# varY = the variable used for the fill aesthetic of the bars; must have exactly two levels
# Y axis is normalized proportion of varY per each varX group
# varX and varY must be names of columns in obj metadata
# Returns ggplot2 bar chart

## Shows which categories are statistically significant by chi-square test (FDR p value)
## ONLY WORKS IF varY HAS EXACTLY TWO LEVELS, throws error if not

make_DA_graph_sig <- function(obj, varX, varY, pval_threshold = 0.01) {
    varY_unique <- unique(obj@meta.data[[varY]])
    if(length(varY_unique) != 2) {
        stop("Error: VarY must have exactly two levels!")
    }
    varY1 <- sym(varY_unique[1])
    varY2 <- sym(varY_unique[2])
    
     varX_sym <- sym(varX)
     varY_sym <- sym(varY)
     
     X_Y_counts <- as.data.frame(table(obj@meta.data[[varX]], obj@meta.data[[varY]]))
     colnames(X_Y_counts) <- c(varX, varY, "nCells")
     varYdf <- as.data.frame(table(obj@meta.data[[varY]]))
     colnames(varYdf) <- c(varY, "TotalCells")
     X_Y_counts <- left_join(X_Y_counts, varYdf, by = varY)
     X_Y_counts$Normalized <- X_Y_counts$nCells / X_Y_counts$TotalCells

    # To apply statistics, first make 'wide' dataframe,
    # Then use chi square test with FDR adjustment
    
    X_Y_counts_wide <- X_Y_counts %>%
        pivot_wider(id_cols = varX, names_from = varY, values_from = nCells) %>%
        rowwise() %>%
        mutate(p_value = perform_chi_square_test(!!varY1, !!varY2)) %>%
        mutate(fdr_p_value = p.adjust(p_value, method = "fdr"))
               
    sig_filtered <- X_Y_counts_wide %>%
                    filter(fdr_p_value < pval_threshold)

    X_Y_counts$Significant <- X_Y_counts[[varX]] %in% sig_filtered[[varX]]

    

    plot <- ggplot(X_Y_counts, aes(x = !!varX_sym, y = Normalized, fill = !!varY_sym)) +
            geom_bar(position = "fill", stat = "identity") +
            theme_minimal() +
            labs(x = varX, y = "Normalized Proportion", fill = varY) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
           geom_text(
                data = subset(X_Y_counts, Significant == TRUE),
                aes(label = "*", y = 1.01), # Adding stars above the bars
                size = 10,
                color = "black")
    
    return(plot)
} 

# GOterm_plot : Function for plotting GO terms from a dataframe of GO terms
# Returns dotplot with the following specifications:
# - X axis: Fold Enrichment
# - Y axis: Pathways, ordered by Fold Enrichment (highest to lowest)
# - Point color : scaled by FDR P Value; 
# - Point size : scaled by number of genes in the pathway 

# Arguments: 
# - GOterm_df : dataframe of GO terms; must have "Fold.Enrichment," "Enrichment.FDR," "nGenes," and "Pathway" columns
# - pval_threshold : function will only show pathways with an Enrichment.FDR value less than this threshold (default 0.05) 
# - fold_threshold : function will only show pathways with a Fold Enrichment value greater than this threshold (default 10) 
# - max_rows : default NULL; if set to a number, function will only show up to this number of pathways. If left NULL, function will show all significant pathways. 
# - color_low : color of dots representing MOST significant p vals (default red) 
# - color_high : color of dots representing LEAST significant p vals (default blue) 
# - size : size of dots with highest nGenes value (default 10) 
# - title : title of chart (default "GO Terms")
# - text_size : size of chart text (default 12) 

GOterm_plot <- function(GOterm_df, 
                        pval_threshold = 0.05, 
                        fold_threshold = 10, 
                        max_rows = NULL,
                        color_low = "red",
                        color_high = "blue",
                        size = 10,
                        title = "GO Terms",
                        text_size = 12) {

    df_filt <- GOterm_df %>% filter(Enrichment.FDR < pval_threshold, Fold.Enrichment > fold_threshold)
    df_filt$Pathway <- factor(df_filt$Pathway, levels = df_filt$Pathway[order(df_filt$Fold.Enrichment)])

    if(!is.null(max_rows)) {
        df_filt <- df_filt %>% slice(1:max_rows)
    }
    
    plot <- ggplot(df_filt, aes(x = Fold.Enrichment, y = Pathway)) +
          geom_point(aes(size = nGenes, color = Enrichment.FDR)) +
          scale_color_gradient(low = color_low, high = color_high) +
          scale_size_continuous(range = c(1, size)) +  # Adjust point size range as needed
          labs(x = "Fold Enrichment", y = NULL, color = "FDR P-value", size = "Number of Genes", title = title) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = text_size),
                panel.grid = element_blank(),
                plot.margin = margin(5, 5, 5, 100))
    return(plot)
    
} 
