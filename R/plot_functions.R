

# Function to create correlation heatmaps
plot_correlation_heatmap = function(data, title, cluster = FALSE, row_order = NULL,
                                    col_order = NULL, breaks) {
  pheatmap::pheatmap(mat = if (!is.null(row_order) && !is.null(col_order))
    data[row_order, col_order] else data,
    treeheight_row = 0, treeheight_col = 0, cluster_rows = cluster,
    cluster_cols = cluster, border_color = NA, legend = FALSE, show_colnames = FALSE,
    show_rownames = FALSE, breaks = breaks, color = cols, fontsize = 12, main = title)
}

# Function to create scatter plot with identity line
plot_correlation_scatter = function(data, title, lim_ax, point_color = "#1170aa") {
  ggplot(data, aes(x = Observed, y = Predicted)) +
    geom_point(color = point_color, shape = 1) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1) +
    labs(title = title, x = "Observed", y = "Predicted") +
    xlim(lim_ax) + ylim(lim_ax) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_text(size = 16))
}

plot_gene_network <- function(cor_mat, n = 50){
  mode(cor_mat) <- "numeric"
  C <- (cor_mat + t(cor_mat))/2
  diag(C) <- 1
  genes <- rownames(C)

  A <- C
  diag(A) <- 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_edges(g, E(g)[weight == 0])
  E(g)$w_abs <- abs(E(g)$weight)
  set.seed(42)
  wc <- cluster_walktrap(g, weights = E(g)$w_abs)
  membership <- membership(wc)
  V(g)$cluster <- factor(membership)
  V(g)$strength <- graph.strength(g, weights = E(g)$w_abs)
  V(g)$degree <- degree(g)


  dist_abs <- as.dist(1 - abs(C))
  hc <- hclust(dist_abs, method = "average")
  order_idx <- hc$order
  C_ord <- C[order_idx, order_idx]
  genes_ord <- rownames(C_ord)
  cluster_vec <- membership[genes_ord]

  E(g)$sign  <- ifelse(E(g)$weight > 0, "positive", "negative")
  E(g)$w_abs <- abs(E(g)$weight)

  g_pos <- g
  E(g_pos)$weight <- E(g)$w_abs
  E(g_pos)$sign   <- E(g)$sign
  E(g_pos)$w_abs  <- E(g)$w_abs
  V(g_pos)$cluster  <- V(g)$cluster
  V(g_pos)$strength <- V(g)$strength

  set.seed(123)
  lay <- ggraph::create_layout(g_pos, layout = "fr")

  n_labels <- min( n, igraph::gorder(g_pos))
  label_genes <- names(sort(igraph::graph.strength(g_pos, weights = E(g_pos)$w_abs),
                            decreasing = TRUE))[seq_len(n_labels)]

  lay$cluster  <- igraph::V(g_pos)$cluster[match(lay$name, igraph::V(g_pos)$name)]
  lay$strength <- igraph::V(g_pos)$strength[match(lay$name, igraph::V(g_pos)$name)]
  lay$label    <- ifelse(lay$name %in% label_genes, lay$name, NA_character_)

  clu_levels <- levels(igraph::V(g_pos)$cluster)
  clu_cols <- setNames(RColorBrewer::brewer.pal(max(3, length(clu_levels)) %/% 1 + 2, "Set3")[seq_along(clu_levels)], clu_levels)
  edge_cols <- c(positive = "#D73027", negative = "#4575B4")



  network_plot <- ggraph(lay) +
    geom_edge_link(aes(edge_alpha = w_abs, edge_width = w_abs, edge_color = sign)) +
    scale_edge_width(range = c(0.2, 0.9), guide = "none") +
    scale_edge_alpha(range = c(0.2, 0.9), guide = "none") +
    scale_edge_color_manual(values = edge_cols, name = "Edge sign") +
    geom_node_point(aes(color = cluster, size = strength)) +
    scale_color_manual(values = clu_cols, name = "Cluster") +
    scale_size_continuous(name = "Node strength", range = c(1, 3)) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3) +
    theme_void() +
    ggtitle("Gene co-expression network")  +
    theme(legend.position = "none")
  return(network_plot)
}


get_topPathways = function(gamma_mean, pathway_names, top_n, clean_names){
  k <- ncol(gamma_mean)
  plot_data <- data.frame()

  for (j in 1:k) {
    loadings <- gamma_mean[, j]
    final_idx <- order(abs(loadings), decreasing = TRUE)[1:top_n]
    pathways_display <- pathway_names[final_idx]
    # Clean names if requested
    if (clean_names) {
      pathways_display <- gsub("^REACTOME_|^GOMF_|^MIPS_|^IRIS_|^DMAP_", "",
                               pathways_display)
      pathways_display <- gsub("_", " ", pathways_display)
    }
    df_temp <- data.frame(Factor_num = j, Pathway = pathways_display,
                          Loading = gamma_mean[final_idx, j],
                          AbsLoading = abs(gamma_mean[final_idx, j]))
    plot_data <- rbind(plot_data, df_temp)
  }
  return(plot_data)
}


dotplot_top_pathways <- function(plot_data){
  ggplot(plot_data, aes(x = Loading, y = reorder(Pathway, AbsLoading),
                        color = Loading>0)) +
    geom_point(size=3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~ Factor_num, scales = "free_y", ncol = 2,
               labeller = labeller(Factor_num = function(x) paste0("Factor ", x))) +
    scale_color_manual(values = c("TRUE" = "brown1", "FALSE" = "royalblue3"),
                       labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                       name = "Loading") +
    #scale_size_continuous(range = c(2, 6)) +
    theme_minimal() +
    labs(x = "", y = "", size = "Absolute\nLoading", color = "Loading") +
    theme(axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 11, face = "bold"),
          legend.position = "none",
          plot.margin = margin(10, 10, 10, 10))  # Add left margin for labels
}








