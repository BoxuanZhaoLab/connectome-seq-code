# Functions for Purkinje Cell Analysis


# function 1: RunSeurat.Joe.v1, filter out unwanted data through mitochondria gene percentage, RNA features, RNA count and doublet percentage.
RunSeurat.Joe.v1 <- function(
    ObjList,
    MitoPercent,
    Filter_nFeature_RNA = c(100, 8000),
    Filter_nCount_RNA = c(500, 50000),
    nPCs = 15,
    DoubletPercent = 0.08
) {
  ObjList <- lapply(
    X = ObjList,
    FUN = function(Obj) {
      Obj[["percent_mt"]] <- PercentageFeatureSet(Obj, pattern = "^mt-")
      Obj <- subset(
        Obj,
        subset = 
          nFeature_RNA > Filter_nFeature_RNA[1] & 
          nFeature_RNA < Filter_nFeature_RNA[2] & 
          nCount_RNA > Filter_nCount_RNA[1] & 
          nCount_RNA < Filter_nCount_RNA[2] & 
          percent_mt < MitoPercent
      )
      return(Obj)
    }
  )
  return(ObjList)
}

# function 2: 

create_br_seurat <- function(data_path, # br for brain region
                             assay_type = "Antibody Capture",
                             synbar_features = c('SynBarpost', 'SynBarpre'),
                             default_assay = "RNA") {
  # Read in raw data  
  br.mrna.data <- Read10X(data_path)
  
  # Create the main Seurat object (Gene Expression)
  br_obj <- CreateSeuratObject(counts = br.mrna.data$`Gene Expression`)
  
  # Create a Synbar Assay
  synbar_data <- br.mrna.data[[assay_type]][synbar_features, ]
  synbar_data <- CreateAssayObject(counts = synbar_data)
  
  # Add Synbar
  br_obj[["Synbar"]] <- synbar_data
  
  # Set as default assay
  Seurat::DefaultAssay(br_obj) <- default_assay
  
  return(br_obj)
}


# function 3: 
create_harmony_integration <- function(br.list, 
                                       filter_param = 5, 
                                       cell_ids = NULL,
                                       nfeatures = 2000,
                                       npcs = 30,
                                       verbose = FALSE) {
  # Load required packages
  require(Seurat)
  require(harmony)
  
  # Step 1: Data filtering
  br.af.filtered <- RunSeurat.Joe.v1(br.list, filter_param)
  
  # Generate cell IDs if not provided
  if (is.null(cell_ids)) {
    cell_ids <- as.character(seq_along(br.af.filtered))
  } else if (length(cell_ids) != length(br.af.filtered)) {
    stop("cell_ids length must match filtered object count")
  }
  
  # Step 2: Object merging
  br.af <- merge(
    x = br.af.filtered[[1]],
    y = br.af.filtered[-1],
    add.cell.ids = cell_ids
  ) %>% JoinLayers()
  
  # Add batch information (extracted from merged cell names)
  br.af$batch <- sub("_.*", "", colnames(br.af))
  
  # Step 3: Normalization and feature selection
  br.af <- NormalizeData(br.af, verbose = verbose)
  
  # Batch-wise HVG selection
  batch_cells <- split(colnames(br.af), br.af$batch)
  var_features <- lapply(batch_cells, function(cells) {
    subset(br.af, cells = cells) %>% 
      FindVariableFeatures(
        selection.method = "vst",
        nfeatures = nfeatures,
        verbose = verbose
      ) %>% 
      VariableFeatures()
  }) %>% unlist() %>% unique()
  
  VariableFeatures(br.af) <- var_features
  
  # Step 4: Dimension reduction and batch correction
  br.af <- ScaleData(br.af, verbose = verbose) %>% 
    RunPCA(
      features = var_features,
      npcs = npcs,
      verbose = verbose
    ) %>% 
    RunHarmony(
      group.by.vars = "batch",
      verbose = verbose
    )
  
  return(br.af)
}



# Set consistent visual style
theme_paper <- theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90")
  )

# function4
analyze_celltype_pc <- function(seurat_obj, cn9_obj, celltype, con_level = list(positive = "Cseq-positive", negative = "Cseq-negative"),
                                output_dir = "pc_output_figures") {
  dir.create(output_dir, showWarnings = FALSE)
  celltype_dir <- output_dir
  
  # Subset and process cells of the specific type
  cells <- seurat_obj
  # Get DEGs
  de.markers <- FindMarkers(
    object         = cells,
    ident.1        = c("1","2","3","4","5","6"),       # can be numeric or character depending on your cluster IDs
    ident.2        = c("0"),              # (default) compare vs. all other cells
    only.pos       = F,             
    min.pct        = 0.25,             # expressed in at least 25% of cells
    logfc.threshold = 0.25             # minimum log2 fold-change threshold
  )
  criteria = list(log2FC = 1, pval = 20)
  top_genes <- de.markers %>%
    tibble::rownames_to_column("gene") %>%
    filter(avg_log2FC > criteria$log2FC & -log10(p_val_adj) > criteria$pval) %>%
    arrange(desc(-log10(p_val_adj)), desc(pct.1-pct.2), desc(avg_log2FC)) %>%
    head(10) %>%
    pull(gene)
  
  top_genes_30 <- de.markers %>%
    tibble::rownames_to_column("gene") %>%
    filter(avg_log2FC > criteria$log2FC & -log10(p_val_adj) > criteria$pval) %>%
    arrange(desc(-log10(p_val_adj)), desc(pct.1-pct.2), desc(avg_log2FC)) %>%
    head(30) %>%
    pull(gene)
  
  
  
  # Generate plots using your original plotting functions
  # 1. UMAP plots
  dimplot1 <- DimPlot(cells, 
                      reduction = "umap.scvi", 
                      label = FALSE,
                      pt.size = 0.5) +
    ggtitle(paste(celltype, "Subclusters")) +
    theme_paper
  
  size_vec <- ifelse(
    cells$connectivity == "Cseq-positive",
    2.5,          # bigger
    0.5           # smaller
  )
  
  dimplot2 <- DimPlot(cells, 
                      reduction = "umap.scvi", 
                      group.by = "connectivity",
                      pt.size = size_vec) +
    ggtitle("Connectivity Status") +
    theme_paper +
    scale_color_manual(values = c("Cseq-positive" = "#f26a11", 
                                  "Cseq-negative" = "#e6e6e6"))
  
  dimplot3 <- DimPlot(cells, 
                      reduction = "umap.scvi", 
                      group.by = "batch",
                      pt.size = 0.5) +
    ggtitle("Batch Integration") +
    theme_paper
  
  ggsave(file.path(celltype_dir, paste0(make.names(celltype), "_umap.pdf")), 
         dimplot1 + dimplot2 + dimplot3,
         width = 22, height = 7)
  
  # 2. Generate heatmap with cell number handling
  if(length(top_genes_30) > 0) {
    # Get expression matrix for top genes
    expr_mat <- GetAssayData(cells, slot = "data")[top_genes_30, ]
    
    # If there are too many cells, randomly sample 1000 cells per condition
    if(ncol(expr_mat) > 1000) {
      set.seed(42) # For reproducibility
      cells_per_condition <- 500
      
      # Sample cells from each condition
      connected_cells <- colnames(expr_mat)[cells$connectivity == "Cseq-positive"]
      nonconnected_cells <- colnames(expr_mat)[cells$connectivity == "Cseq-negative"]
      
      sampled_connected <- sample(connected_cells, 
                                  min(cells_per_condition, length(connected_cells)))
      sampled_nonconnected <- sample(nonconnected_cells, 
                                     min(cells_per_condition, length(nonconnected_cells)))
      
      # Combine sampled cells
      sampled_cells <- c(sampled_connected, sampled_nonconnected)
      
      # Subset expression matrix and scale
      expr_mat <- expr_mat[, sampled_cells]
      print(paste("Sampled", length(sampled_cells), "cells for heatmap visualization"))
    }
    
    # Scale the expression matrix
    expr_mat_scaled <- t(scale(t(expr_mat)))
    
    # Create annotation dataframe
    annotation_col <- data.frame(
      Connectivity = cells$connectivity[match(colnames(expr_mat), colnames(cells))],
      row.names = colnames(expr_mat)
    )
    
    annotation_col$Connectivity <- factor(
      annotation_col$Connectivity,
      levels = c("Cseq-positive","Cseq-negative"),      # original values
      labels = c("Positive-cluster","Negative-cluster")   # display labels
    )
    
    ann_colors <- list(
      Connectivity = setNames(
        c("#f26a11", "#e6e6e6"),                 # colours
        c("Positive-cluster", "Negative-cluster")# *new* labels
      )
    )
    
    pdf(file.path(celltype_dir, paste0(make.names(celltype), "_heatmap.pdf")), 
        width = 12, height = 6)
    
    
    pheatmap(expr_mat_scaled,
             annotation_col = annotation_col,
             show_colnames = FALSE,
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
             breaks = seq(-2, 2, length.out = 101),
             main = paste("Expression of Top Marker Genes in", celltype),
             clustering_method = "ward.D2",
             annotation_colors = ann_colors
    )
    dev.off()
  }
  
  # 3. Volcano plot
  deg_results_df <- de.markers %>%
    tibble::rownames_to_column("gene") %>%
    filter(-log10(p_val_adj) > 300) %>% 
    mutate(sig = abs(pct.1 - pct.2) > 0.5 & abs(avg_log2FC) > 2)
  
  volcano_plot <- ggplot(deg_results_df, 
                         aes(x = avg_log2FC, 
                             y = exp(10*abs(pct.1-pct.2))/10,
                             color = sig,
                             label = ifelse(gene %in% top_genes, gene, ""))) +
    geom_point(alpha = 0.6) +
    geom_text_repel(force = 10,
                    size = 6,
                    max.overlaps = Inf) +
    scale_color_manual(values = c("grey50", "#b4403e")) +
    theme_paper +
    labs(title = paste("Differential Expression in", celltype),
         x = "log2 Fold Change",
         y = "Scaled detection rate")
  
  ggsave(file.path(celltype_dir, paste0(make.names(celltype), "_volcano.pdf")), 
         volcano_plot, 
         width = 10, height = 8)
  
  # 4. Feature plots for top 6 genes
  feature_plots <- FeaturePlot(cells, 
                               features = top_genes, 
                               reduction = "umap.scvi",
                               ncol = 5,
                               order = TRUE) &
    theme_paper &
    scale_color_gradientn(colors = c("#e6e6e6", "#8176b2"))
  
  ggsave(file.path(celltype_dir, paste0(make.names(celltype), "_feature_plots.pdf")), 
         feature_plots, 
         width = 20, height = 8)
  
  # 5. Validation violin plots
  cells_cn9 <- cn9_obj
  plot_data <- data.frame()
  cells_connected <- subset(cells, connectivity == "Cseq-positive")
  cells_nonconnected <- subset(cells, connectivity == "Cseq-negative")
  
  for(gene in top_genes) {
    expr_data <- data.frame(
      gene = gene,
      condition = c(rep("Cseq-positive", ncol(cells_connected)),
                    rep("Cseq-negative", ncol(cells_nonconnected)),
                    rep("AAV1-labeled", ncol(cells_cn9))),
      
      expression = c(
        GetAssayData(cells_connected, slot = "data")[gene,],
        GetAssayData(cells_nonconnected, slot = "data")[gene,],
        GetAssayData(cells_cn9, slot = "data")[gene,]
      )
    )
    plot_data <- rbind(plot_data, expr_data)
  }
  
  plot_data$condition <- factor(plot_data$condition, levels = c("AAV1-labeled", "Cseq-positive", "Cseq-negative"))  # wanted order
  plot_data$gene <- factor(plot_data$gene, levels = top_genes)
  
  
  validation_plot <- ggplot(plot_data, 
                            aes(x = condition, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.9, fill = "white") +
    facet_wrap(~gene, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c( "#f26a11","#4292c9","#8176b2"),
                      breaks = c("AAV1-labeled", "Cseq-positive", "Cseq-negative")
    ) +
    theme_paper +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 24)
    ) +
    labs(title = paste("Expression Validation of Top Genes -", celltype),
         x = "Condition",
         y = "Expression Level")
  
  ggsave(file.path(celltype_dir, paste0(make.names(celltype), "_validation_violin.pdf")), 
         validation_plot, 
         width = 22, height = 10)
  
  # 6. Connectivity proportion plot
  connectivity_props <- table(
    Idents(cells),
    cells$connectivity) %>%
    as.data.frame() %>%
    tidyr::spread(Var2, Freq) %>%
    mutate(percent_connected = `Cseq-positive` / (`Cseq-positive` + `Cseq-negative`) * 100)
  
  prop_plot <- ggplot(connectivity_props, 
                      aes(x = reorder(Var1, -percent_connected), 
                          y = percent_connected)) +
    geom_bar(stat = "identity", fill = "#f26a11") +
    theme_paper +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      panel.grid.major.x = element_blank()
    ) +
    labs(x = "Subcluster", 
         y = "Percent Connected (%)",
         title = paste("Proportion of Connected Cells per Subcluster -", celltype)) +
    geom_text(aes(label = sprintf("%.1f%%", percent_connected)), 
              vjust = -0.5,
              size = 4)
  
  ggsave(file.path(celltype_dir, paste0(make.names(celltype), "_connectivity_proportion.pdf")), 
         prop_plot,
         width = 8, height = 6)
  
  # 7. Gene-wise violin plots
  vln_plot <- VlnPlot(cells, 
                      features = top_genes,
                      ncol = 5,
                      pt.size = 0) &
    theme_paper &
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  ggsave(file.path(celltype_dir, paste0(make.names(celltype), "_gene_violin.pdf")), 
         vln_plot,
         width = 25, height = 12)
  
  # 8. Two-condition violin plots
  plot_data_2 <- plot_data %>%                      # keep the original intact
    filter(condition != "AAV1-labeled") %>% 
    mutate(
      condition = factor(
        condition,
        levels = c("Cseq-positive", "Cseq-negative"),
        labels = c("Positive-cluster", "Negative-cluster")
      )
    )
  
  pal_two <- setNames(
    c("#f26a11", "#e6e6e6"),                       # colours
    c("Positive-cluster", "Negative-cluster")      # new names
  )
  
  
  
  two_condition_violin <- ggplot(plot_data_2 %>% 
                                   filter(condition != "AAV1-labeled"), 
                                 aes(x = condition, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.9, fill = "white") +
    facet_wrap(~gene, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = pal_two) +
    theme_paper +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 24)
    ) +
    labs(title = paste("Expression in Connected vs Non-connected", celltype, "Cells"),
         x = "Condition",
         y = "Expression Level")
  
  ggsave(file.path(celltype_dir, paste0(make.names(celltype), "_two_condition_violin.pdf")), 
         two_condition_violin,
         width = 22, height = 10)
  
  # Save results
  write.csv(de.markers, file.path(celltype_dir, paste0(make.names(celltype), "_DEGs.csv")))
  
  return(list(
    de.markers = de.markers,
    top_genes = top_genes,
    cells = cells
  ))
}
















