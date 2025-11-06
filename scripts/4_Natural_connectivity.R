# Load necessary packages
library(igraph)
library(ggplot2)
# === Step 1: Load edge list ===
edge_list <- read.csv("network_virome.csv")  # Replace with your file

# === Step 2: Convert to igraph object ===
g <- graph_from_data_frame(edge_list, directed = FALSE)

# === Step 3: Create binary adjacency matrix ===
adj_matrix <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
diag(adj_matrix) <- 0  # Remove self-loops

# Optional: Save adjacency matrix
write.csv(adj_matrix, "adjacency_matrix_as_micorbiome.csv", row.names = TRUE)

# === Step 4: Build igraph object from matrix for robustness analysis ===
g_adj <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# === Step 5: Define natural connectivity (normalized) ===
natural_connectivity <- function(graph) {
  if (vcount(graph) < 2) return(NA)
  A <- as.matrix(as_adjacency_matrix(graph, sparse = FALSE))
  eig_vals <- eigen(A, only.values = TRUE)$values
  nc_raw <- log(mean(exp(Re(eig_vals))))
  N <- vcount(graph)
  nc_norm <- nc_raw / (N - log(N))  # Normalized form
  return(nc_norm)
}

# === Step 6: Robustness test: hub-targeted removal ===
robustness_test <- function(graph, max_remove_prop = 0.8) {
  N <- vcount(graph)
  max_remove <- floor(N * max_remove_prop)
  
  ranked_nodes <- order(rank(betweenness(graph)) + rank(degree(graph)), decreasing = TRUE)
  proportions <- seq(1, max_remove) / N
  connectivity <- numeric(length(proportions))
  
  for (i in seq_along(proportions)) {
    remove_ids <- ranked_nodes[1:i]
    sub_g <- delete_vertices(graph, V(graph)[remove_ids])
    connectivity[i] <- natural_connectivity(sub_g)
  }
  
  return(data.frame(Proportion_Removed = proportions, Natural_Connectivity = connectivity))
}

# === Step 7: Run test and plot ===
set.seed(999)
as_robust_data <- robustness_test(g_adj)
