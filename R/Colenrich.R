#' Colenrich: Integrative GO Enrichment and Community Detection
#'
#' @description
#' Performs GO enrichment analysis and community detection by integrating known gene functions from Gene Ontology (GO)
#' across multiple datasets. For each group (column) in the input data, the function retrieves GO annotations using
#' the appropriate organism annotation database, identifies the common GO terms across groups, and constructs a graph
#' based on these common terms. Louvain community detection is then applied to identify functionally connected modules.
#'
#' @param data A data frame in which each column represents a group containing gene identifiers.
#' @param ID A character string specifying the gene identifier type (e.g., "ENSEMBL"). Default is "ENSEMBL".
#' @param organism A character string specifying the organism. Allowed values are "Human" or "Mouse". Default is "Mouse".
#'
#' @return A list containing:
#' \describe{
#'   \item{merged_df}{A data frame of merged enrichment results containing common GO terms.}
#'   \item{community}{An igraph community object resulting from Louvain community detection.}
#'   \item{membership_df}{A data frame of graph nodes with their assigned community memberships.}
#'   \item{gene_function}{A data frame of nodes that correspond to GO terms.}
#' }
#'
#' @import AnnotationDbi
#' @importFrom igraph graph_from_data_frame cluster_louvain modularity
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' # Example input: a data frame with two groups of gene identifiers
#' gene_data <- data.frame(
#'   group1 = c("ENSG000001", "ENSG000002", "ENSG000003"),
#'   group2 = c("ENSG000002", "ENSG000003", "ENSG000004"),
#'   stringsAsFactors = FALSE
#' )
#' # For human data, ensure that the org.Hs.eg.db package is installed.
#' result <- Colenrich(gene_data, ID = "ENSEMBL", organism = "Human")
#' print(result$merged_df)
#' print(result$membership_df)
#' print(result$gene_function)
#' }
#'
#' @export
Colenrich <- function(data, ID = "ENSEMBL", organism = "Mouse", resolution = 1.1) {

  # Ensure necessary packages are available
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required but not installed.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required but not installed.")
  }

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

  enrichment_results1 <- list()
  enrichment_results2 <- list()

  # Process each column (group) in the input data frame
  for (col_name in colnames(data)) {
    message("Processing group: ", col_name)

    # Obtain gene list and remove NA values
    gene_list <- na.omit(data[[col_name]])

    if (length(gene_list) == 0) {
      message("Skipping ", col_name, " - No valid genes.")
      next
    }

    # Retrieve GO annotations using AnnotationDbi::select
    gene_annotations <- AnnotationDbi::select(OrgDb,
                                              keys = gene_list,
                                              columns = c("GO"),
                                              keytype = ID)

    if (!is.null(gene_annotations) && nrow(gene_annotations) > 0) {
      enrichment_results1[[col_name]] <- gene_annotations$GO
      enrichment_results2[[col_name]] <- gene_annotations[, 1:2]
    }
  }

  if (length(enrichment_results1) == 0) {
    stop("No enrichment results obtained. Please check your input data.")
  }

  # Identify common GO terms across groups
  common_elements <- Reduce(intersect, enrichment_results1)

  # Process enrichment results to extract rows corresponding to common GO terms
  sig_enrichment <- list()
  for (group in names(enrichment_results2)) {
    result_df <- as.data.frame(enrichment_results2[[group]])
    mat1_common <- result_df[result_df$GO %in% common_elements, ]
    sig_enrichment[[group]] <- mat1_common
  }

  # Merge the enrichment results from all groups
  merged_df <- do.call(rbind, sig_enrichment)

  # Construct a graph from the merged enrichment results and perform community detection
  g <- igraph::graph_from_data_frame(merged_df, directed = FALSE)
  community <- igraph::cluster_louvain(g, resolution = resolution)
  mod_value <- igraph::modularity(community)
  message("Modularity: ", mod_value)

  membership_result <- community$membership
  membership_df <- data.frame(node = community$names,
                              community = membership_result,
                              stringsAsFactors = FALSE)

  # Filter for nodes corresponding to GO terms
  gene_function <- membership_df[grepl("GO", membership_df$node), ]
  # 重新计算社区
  la <- igraph::cluster_louvain(g)
  # 绘图
  plot(la,
       g,
       vertex.size = 5,
       vertex.label.cex = 0.7,
       vertex.color = membership(la),
       layout = layout_with_fr(g),
       main = "GO Enrichment Network and Community Detection")
  # Return a list of results
  return(list(merged_df = merged_df,
              community = community,
              membership_df = membership_df,
              gene_function = gene_function,
              Object_of_drawing = g))
}
