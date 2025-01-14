#' Generate a heatmap of scRNA-seq differential expression results
#' @author Ira Fleming, \email{ira.fleming@@cuanschutz.edu}
#' @param seurat_object Seurat object containing the single-cell data
#' @param markers_df Output dataframe from Seurat's FindAllMarkers() containing genes, logfold change, and p-values
#' @param n_markers_shown Number of positive markers per condition to show in heatmap
#' @param identity_comparison Column name in Seurat object metadata used to group cells e.g. "cell_type"
#' @param ranking_metric Column name to rank genes by (default: "avg_log2FC")
#' @param cluster Decision on whether genes should be clustered or just listed by gene rank
#' @param slot Which data slot to map
#' @return A ggplot object containing the heatmap
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @export

heatmarkers <- function(seurat_object,
                          markers_df,
                          n_markers_shown,
                          identity_comparison,
                          slot,
                          ranking_metric = "avg_log2FC",
                          cluster = F) {
  
  # Select top genes from differential expression results
  top_markers <- markers_df %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::top_n(n = n_markers_shown, 
          wt = .data[[ranking_metric]]) 
  
  top_markers <- top_markers[top_markers$gene == unique(top_markers$gene),]
  
  # Get expression data from Seurat object
  expr_data <- Seurat::FetchData(seurat_object, 
                        vars = top_markers$gene, 
                        layer = slot) %>%
    scale() %>%
    as.data.frame()
  
  # Add metadata to expression data
  expr_data$cell_id <- rownames(expr_data) 
  expr_data$group <- seurat_object[[identity_comparison]] %>% dplyr::pull()
  
  # Cluster rows (genes)
  gene_matrix <- t(expr_data[, top_markers$gene])
  gene_dist <- dist(gene_matrix)
  gene_clust <- hclust(gene_dist, method = "complete")
  gene_order <- gene_clust$labels[gene_clust$order]
  
  # Prepare plot data
  plot_data <- expr_data %>%
    tidyr::pivot_longer(-c(cell_id, group), 
                names_to = "gene", 
                values_to = "expression")
  
  # Order genes by clustering
  if(cluster == T){
    plot_data$gene <- factor(plot_data$gene, levels = gene_order)
  } else if(cluster == F){
    plot_data$gene <- factor(plot_data$gene, levels = rev(top_markers$gene))
  }
  colnames(plot_data) <- c("cell_id","group","gene","expression")

  max_abs <- max(abs(range(plot_data$expression)))
  value_limits <- c(-max_abs, max_abs)

  print(plot_data)
  # Create heatmap
  the_heatmap <- ggplot2::ggplot(plot_data, aes(x = cell_id, y = gene, fill = expression)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Scaled\nExpression", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 12),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      strip.text = element_text(size = 10, face = "bold")
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("")
  return(the_heatmap)
}

#heatmarkers <- create_heatmap(
#  seurat_object = sobj,
#  markers_df = sobj_de,
#  n_markers_shown = 15,
#  identity_comparison = "cell_type",
#  ranking_metric = "avg_log2FC",
#  cluster = F
#)

