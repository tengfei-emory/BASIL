

#' Plot correlation heatmap
#'
#' @param data Correlation matrix
#' @param title Plot title
#' @param cluster Whether to cluster rows/columns
#' @param row_order Optional row ordering
#' @param col_order Optional column ordering
#' @param breaks Color breaks
#'
#' @return A pheatmap object
#' @export
plot_correlation_heatmap = function(data, title, cluster = FALSE, row_order = NULL,
                                    col_order = NULL, breaks) {
  pheatmap::pheatmap(mat = if (!is.null(row_order) && !is.null(col_order))
    data[row_order, col_order] else data,
    treeheight_row = 0, treeheight_col = 0, cluster_rows = cluster,
    cluster_cols = cluster, border_color = NA, legend = FALSE, show_colnames = FALSE,
    show_rownames = FALSE, breaks = breaks, color = cols, fontsize = 12, main = title)
}

#' Plot gene co-expression network from correlation matrix
#'
#' @param cor_mat A correlation matrix (genes × genes) with row and
#'   column names indicating gene identifiers.
#' @param n Integer, the number ofgenes to label in the network. Genes are
#'   ranked by their total connection strength (sum of absolute correlations).
#'   Default is 50.
#'
#' @return A ggraph object showing the gene co-expression network with:
#'   \itemize{
#'     \item Node size proportional to connection strength
#'     \item Edge color indicating sign: red for positive, blue for negative correlations
#'     \item Edge width/transparency proportional to correlation magnitude
#'     \item Labels on the n most connected genes
#'   }
#'
#' @export
#'
plot_gene_network <- function(cor_mat, n = 50){
  mode(cor_mat) <- "numeric"
  C <- (cor_mat + t(cor_mat))/2
  diag(C) <- 1
  genes <- rownames(C)

  A <- C
  diag(A) <- 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE,
                                   diag = FALSE)
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
  clu_cols <- setNames(RColorBrewer::brewer.pal(max(3, length(clu_levels)) %/% 1 +
              2, "Set3")[seq_along(clu_levels)], clu_levels)
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

#' Extract top pathways for graphical matrix Gamma
#'
#' Identifies the pathways with the highest absolute loadings for each latent
#' factor from the Gamma matrix (pathway-level loadings). Optionally includes
#' credible intervals and cleans pathway names by removing common database prefixes.
#'
#' @param gamma_mean Posterior mean of graphical matriz Gamma (pathways × factors).
#'   Typically obtained from \code{compute_Gamma_samples}.
#' @param pathway_names Character vector of pathway names, corresponding to rows
#'   of \code{gamma_mean}. Should be the same length as \code{nrow(gamma_mean)}.
#' @param top_n Integer, the number of top pathways to extract per factor.
#' @param clean_names Logical, whether to clean pathway names by removing common
#'   database prefixes and replacing underscores with spaces.
#' @param ci_lower Optional numeric matrix (pathways × factors) of lower bounds
#'   of credible intervals for Gamma. If provided, will be included
#'   in the output. Default is NULL.
#' @param ci_upper Optional numeric matrix (pathways × factors) of upper bounds
#'   of credible intervals for Gamma. If provided, will be included
#'   in the output. Default is NULL.
#' @param max_length Optional integer specifying maximum length for pathway names.
#'   Names longer than this will be truncated. If NULL (default),
#'   no truncation is performed.
#'
#' @return A data frame with columns:
#'   \item{Factor_num}{Factor number (1 to k)}
#'   \item{Pathway}{Pathway name (cleaned if clean_names = TRUE)}
#'   \item{Loading}{Pathway value (can be positive or negative)}
#'   \item{AbsLoading}{Pathway absolute value  (used for ranking)}
#'   \item{CI_lower}{Lower bound of credible interval (if ci_lower provided)}
#'   \item{CI_upper}{Upper bound of credible interval (if ci_upper provided)}
#'
#' @export
#'
get_topPathways <- function(gamma_mean, pathway_names, top_n, clean_names = TRUE,
                            ci_lower = NULL, ci_upper = NULL, max_length = NULL){
  k <- ncol(gamma_mean)
  plot_data <- data.frame()
  for (j in 1:k) {
    loadings <- gamma_mean[, j]
    final_idx <- order(abs(loadings), decreasing = TRUE)[1:top_n]
    pathways_display <- pathway_names[final_idx]
    if (clean_names) {
      pathways_display <- gsub("^REACTOME_|^GOMF_|^MIPS_|^IRIS_|^DMAP_", "",
                               pathways_display)
      pathways_display <- gsub("_", " ", pathways_display)
    }
    # Truncate long names if max_length is specified
    if (!is.null(max_length)) {
      pathways_display <- ifelse(
        nchar(pathways_display) > max_length,
        paste0(substr(pathways_display, 1, max_length - 3), "..."),
        pathways_display
      )
    }

    df_temp <- data.frame(Factor_num = j, Pathway = pathways_display,
                          Loading = gamma_mean[final_idx, j],
                          AbsLoading = abs(gamma_mean[final_idx, j]))
    # Add credible intervals if provided
    if (!is.null(ci_lower) && !is.null(ci_upper)) {
      df_temp$CI_lower <- ci_lower[final_idx, j]
      df_temp$CI_upper <- ci_upper[final_idx, j]
    }
    plot_data <- rbind(plot_data, df_temp)
  }
  return(plot_data)
}


#' Create dot plot of top pathway by factor
#'
#' Visualizes the top pathways for each latent factor as a faceted dot plot,
#' with point color indicating the sign of the loading (positive or negative).
#' Optionally displays credible intervals as error bars if available in the data.
#'
#' @param plot_data A data frame with pathway loading information, typically
#'   produced by \code{get_topPathways}. Must contain columns:
#'   \describe{
#'     \item{Factor_num}{Factor identifier (for faceting)}
#'     \item{Pathway}{Pathway name (for y-axis)}
#'     \item{Loading}{Loading value (for x-axis and color)}
#'     \item{AbsLoading}{Absolute loading (for ordering pathways)}
#'     \item{CI_lower}{Optional: Lower bound of credible interval}
#'     \item{CI_upper}{Optional: Upper bound of credible interval}
#'   }
#'
#' @return A ggplot2 object showing a faceted dot plot with:
#'   \itemize{
#'     \item One panel per factor
#'     \item Pathways ordered by absolute loading (highest at top)
#'     \item Points colored by loading sign: red (positive), blue (negative)
#'     \item Horizontal lines showing credible intervals (if CI columns present)
#'     \item Vertical dashed line at zero
#'   }
#'
#' @export
#'
dotplot_top_pathways <- function(plot_data){
  p <- ggplot(plot_data, aes(x = Loading, y = reorder(Pathway, AbsLoading),
                             color = Loading>0))
  if ("CI_lower" %in% names(plot_data) && "CI_upper" %in% names(plot_data)) {
    p <- p + geom_linerange(aes(xmin = CI_lower, xmax = CI_upper),
                            size = 0.5, alpha = 0.6)
  }
  p <- p +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~ Factor_num, scales = "free_y", ncol = 2,
               labeller = labeller(Factor_num = function(x) paste0("Factor ", x))) +
    scale_color_manual(values = c("TRUE" = "brown1", "FALSE" = "royalblue3"),
                       labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                       name = "Loading") +
    theme_minimal() +
    labs(x = "", y = "", size = "Absolute\nLoading", color = "Loading") +
    theme(axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 11, face = "bold"),
          legend.position = "none",
          plot.margin = margin(10, 10, 10, 10))
  return(p)
}









