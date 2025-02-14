#' Calculate Community Centrality and Generate a GO Enrichment Network Plot
#'
#' @description
#' This function performs centrality analysis on a specified community within a graph,
#' filters GO terms from the resulting centrality measures, extracts the top-ranked GO
#' terms based on a selected metric, retrieves their ontology information from a specified
#' organism's annotation database, constructs a bipartite network, and visualizes the network
#' using ggraph.
#'
#' @param graph An igraph object representing the complete network.
#' @param community_obj A community detection object (e.g., returned by \code{cluster_louvain}) obtained from \code{graph}.
#' @param comm_index An integer specifying which community to analyze. Default is 1.
#' @param organism A character string specifying the organism. Allowed values are "Human" or "Mouse". Default is "Mouse".
#' @param TOP A numeric value indicating the number of top-ranked rows to extract for each centrality measure. Default is 10.
#'
#' @return A list with two components:
#' \describe{
#'   \item{plotData}{A ggplot object representing the network plot.}
#'   \item{tableData}{A data frame containing the filtered centrality measures for GO terms.}
#' }
#'
#' @details
#' The function first computes centrality metrics (degree, betweenness, closeness, and eigenvector centrality)
#' for the specified community. It then filters nodes corresponding to GO terms (using a pattern match on "GO")
#' and extracts the top \code{TOP} rows for each centrality metric. The function retrieves the ONTOLOGY information
#' for these GO terms from the specified organism's annotation database (using \code{AnnotationDbi::select}), constructs
#' a bipartite network, and visualizes the network using a manual layout generated from a bipartite layout. The edges
#' are displayed in black, and node aesthetics are mapped based on degree and class.
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggraph
#' @import tidygraph
#' @import igraph
#' @import wesanderson
#' @import gplots
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#'   # Assume 'my_graph' is an igraph object and 'comm' is obtained via cluster_louvain(my_graph)
#'   result <- calculate_community_centrality(graph = my_graph, community_obj = comm, comm_index = 1,
#'                                            organism = "Human", TOP = 10)
#'   print(result$plotData)
#'   head(result$tableData)
#' }
#'
#' @export
calculate_community_centrality <- function(graph, community_obj, comm_index = 1, organism = "Mouse", TOP = 10) {
  
  # Extract the list of communities from the community detection object
  comm_list <- communities(community_obj)
  
  # Load appropriate organism annotation package based on 'organism'
  OrgDb <- switch(organism,
                  "Human" = {
                    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
                      stop("Package 'org.Hs.eg.db' is required for Human data but not installed.")
                    }
                    org.Hs.eg.db
                  },
                  "Mouse" = {
                    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
                      stop("Package 'org.Mm.eg.db' is required for Mouse data but not installed.")
                    }
                    org.Mm.eg.db
                  },
                  stop("Unsupported organism. Please use 'Human' or 'Mouse'."))
  
  # Check if the specified community index is valid
  if (comm_index < 1 || comm_index > length(comm_list)) {
    stop("Invalid community index: ", comm_index)
  }
  
  # Create a subgraph from the nodes in the specified community
  subg <- induced_subgraph(graph, comm_list[[comm_index]])
  
  # Calculate centrality measures
  degree_cen <- degree(subg)
  betweenness_cen <- betweenness(subg)
  closeness_cen <- closeness(subg)
  eigenvector_cen <- eigen_centrality(subg)$vector
  
  # Compile the centrality measures into a data frame
  centrality_df <- data.frame(
    node = V(subg)$name,
    degree = degree_cen,
    betweenness = betweenness_cen,
    closeness = closeness_cen,
    eigenvector = eigenvector_cen,
    stringsAsFactors = FALSE
  )
  
  # Filter for GO terms from the centrality data
  genefunction <- centrality_df[grepl("GO", centrality_df$node), ]
  # Remove the node column so that only centrality measures remain for ranking
  genefunction$node <- NULL
  
  # Extract the top 'TOP' rows for each centrality metric
  top10_list <- apply(genefunction, 2, function(x) {
    idx <- order(x, decreasing = TRUE)
    top_idx <- idx[1:TOP]
    data.frame(row = rownames(genefunction)[top_idx],
               value = x[top_idx],
               stringsAsFactors = FALSE)
  })
  
  # For each centrality measure, extract the row names and add a column indicating the measure (sample)
  new_list <- lapply(names(top10_list), function(name) {
    df <- top10_list[[name]]
    data.frame(row = df$row, sample = name, stringsAsFactors = FALSE)
  })
  
  # Combine the results into one data frame
  combined_df <- do.call(rbind, new_list)
  
  # Retrieve ONTOLOGY information for the GO terms
  go_info <- AnnotationDbi::select(OrgDb,
                                   keys = combined_df$row,
                                   columns = c("ONTOLOGY"),
                                   keytype = "GO")
  mat1_common <- go_info[match(combined_df$row, go_info$GO), ]
  combined_df$ONTOLOGY <- mat1_common$ONTOLOGY
  
  # Prepare edge data
  edge_data <- combined_df[, c("row", "sample", "ONTOLOGY")] %>%
    dplyr::rename(from = sample, to = row) %>%
    dplyr::select(from, to)
  
  variable_info.sum <- combined_df[, c(1, 3)]
  
  # Build node data by pivoting edge data and merging ontology information
  node_data <- edge_data %>%
    dplyr::select(from, to) %>%
    tidyr::pivot_longer(cols = c(from, to),
                        names_to = "Class",
                        values_to = "node") %>%
    dplyr::mutate(
      class1 = variable_info.sum$ONTOLOGY[match(node, variable_info.sum$row)]
    ) %>%
    dplyr::select(node, class1) %>%
    dplyr::rename(Class = class1) %>%
    dplyr::distinct(node, .keep_all = TRUE) %>%
    dplyr::arrange(Class) %>%
    dplyr::mutate(true_name = node)
  
  node_data$Class[is.na(node_data$Class)] <- "Algorithm"
  
  code_vec <- setNames(as.numeric(factor(node_data$node, levels = node_data$node)),
                       node_data$node)
  
  edge_data2 <- data.frame(from = code_vec[edge_data$from],
                           to = code_vec[edge_data$to])
  
  # Create a tidygraph object
  total_graph <- tidygraph::tbl_graph(nodes = node_data,
                                      edges = edge_data2,
                                      directed = TRUE) %>%
    dplyr::mutate(Degree = centrality_degree(mode = 'all'))
  
  # Define color palettes
  pal <- wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")
  pal2 <- gplots::colorpanel(300, low = 'steelblue',
                             mid = "#F9F0D9",
                             high = 'coral4')
  
  g <- total_graph
  
  # Set bipartite mapping
  V(g)$type <- bipartite_mapping(g)$type
  
  # Generate a bipartite layout and adjust coordinates
  coords <- create_layout(g, layout = "bipartite") %>%
    dplyr::select(x, y)
  coords$y[coords$y == 0] <- 0.3
  coords <- coords %>%
    dplyr::mutate(theta = x / (max(x) + 1) * 2 * pi,
                  r = y + 1,
                  x = r * cos(theta),
                  y = r * sin(theta))
  
  my_graph <- create_layout(graph = g,
                            layout = "manual",
                            x = coords$x,
                            y = coords$y)
  
  # Set custom color values for node classes
  value.sum <- c("#436342", "#8f742f", "#d8995b",
                 "#ad5657", "#c6a58d", "#e4c97d", "#66a9b3",
                 "#62b2dc", "#765396", "#b5cf70")
  
  # Generate the network plot using ggraph
  plot_obj <- ggraph(my_graph, layout = "bipartite") +
    geom_edge_link(color = "black", edge_alpha = 1) +
    geom_node_point(
      aes(fill = Class,
          color = Class,
          size = Degree),
      shape = 21,
      show.legend = TRUE
    ) +
    scale_fill_manual(values = value.sum) +
    scale_color_manual(values = value.sum) +
    geom_node_text(
      aes(
        x = x * 1.03,
        y = y * 1.03,
        label = true_name,
        hjust = ifelse(Class %in% c("Algorithm"), 'inward', "outward"),
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        size = ifelse(Class %in% c("Algorithm"), 5, 5),
        colour = Class
      ),
      repel = FALSE,
      alpha = 1,
      show.legend = FALSE
    ) +
    guides(
      edge_width = guide_legend(title = "-log10(BH adjusted P value)",
                                override.aes = list(shape = NA)),
      edge_color = ggraph::guide_edge_colorbar(title = "Effect"),
      fill = guide_legend(
        title = "Class",
        override.aes = list(size = 7, linetype = "blank")
      ),
      size = guide_legend(title = "Degree", override.aes = list(linetype = 0))
    ) +
    ggraph::scale_edge_color_gradientn(colours = colorRampPalette(c("#5ca995", "#f8f1dd", "#ca5737"))(30)) +
    ggraph::scale_edge_width(range = c(0.1, 1)) +
    scale_size_continuous(range = c(1.5, 15)) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    ) +
    expand_limits(x = c(-3.2, 3.2), y = c(-3.2, 3.2))
  
  print(plot_obj)
  
  # Return the plot and the filtered gene function table
  result_list <- list(plotData = plot_obj, tableData = genefunction)
  return(result_list)
}
