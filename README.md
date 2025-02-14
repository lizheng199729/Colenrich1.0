# Colenrich

**Colenrich** is an R package that provides integrated Gene Ontology (GO) enrichment analysis and network visualization tools for high-throughput sequencing datasets. The package includes functions for combining datasets, performing community detection using the Louvain algorithm, calculating various centrality measures within communities, and generating customized network plots with ggraph. Colenrich leverages known gene functions to overcome the sparsity issues often encountered in overlapping gene signatures, facilitating drug mechanism studies and other bioinformatics applications.

## Features

- **GO Enrichment Analysis:** Integrates multiple datasets and retrieves common GO terms.
- **Community Detection:** Uses the Louvain algorithm to identify functionally connected modules.
- **Centrality Analysis:** Computes degree, betweenness, closeness, and eigenvector centralities.
- **Network Visualization:** Generates high-quality network plots using ggraph.
- **Built-in Example Data:** Includes sample data for demonstration purposes.

## Installation

You can install the development version of **Colenrich** directly from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install Colenrich from GitHub
devtools::install_github("YourGithubUsername/Colenrich")
```
## Usage
Basic GO Enrichment Analysis
Load the package and use the Colenrich function to perform GO enrichment analysis:
```r
#Import sample data
data("example_data", package = "Colenrich")
#Start calculating overlapping functions between different data sets and detect major communities using the Greedy community algorithm
#data: The input dataset, usually a data frame with each column representing a group of gene identifiers.
#ID = "ENSEMBL": Specifies that the gene identifiers are ENSEMBL IDs, guiding the function to use the appropriate annotation database.
#organism = "Human": Indicates that the data is from humans, ensuring the use of human-specific annotation databases (e.g., org.Hs.eg.db).
#resolution = 0.5: Adjusts the granularity of community detection; lower values yield fewer, larger communities, while higher values produce more, finer communities.
result <- Colenrich(data,ID = "ENSEMBL",organism = "Human",resolution = 0.5)
#Computationally interested communities with major features and visualizations
#result$Object_of_drawing: An igraph object used for constructing and visualizing the network.
#result$community: The community detection result (e.g., from cluster_louvain) derived from the network.
#comm_index = 1: Specifies that the analysis should focus on the first detected community.
#organism = "Human": Indicates that the human-specific annotation database (e.g., org.Hs.eg.db) should be used for retrieving GO annotations.
#TOP = 20: Instructs the function to extract the top 20 ranked rows for each centrality metric during gene filtering.
visualization <- calculate_community_centrality (result$Object_of_drawing,result$community,comm_index = 1, organism = "Human", TOP = 20)
