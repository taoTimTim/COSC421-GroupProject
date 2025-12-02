#This code shows clustrering and modularity without the coorinates of the eeg elecgtrodes.

library(igraph)
library(ggraph)
library(ggplot2)


# loading WPLI matrix


load_matrix <- function(path) {
  df <- read.csv(path, header = TRUE, check.names = FALSE)
  rownames(df) <- df[,1]
  df <- df[,-1]
  m <- as.matrix(df)
  m <- apply(m, 2, as.numeric)
  rownames(m) <- rownames(df)
  return(m)
}

# Keep top X% edges (density thresholding)

keep_top_density <- function(mat, dens) {
  diag(mat) <- 0
  mat[is.na(mat)] <- 0
  mat <- (mat + t(mat)) / 2
  
  n <- nrow(mat)
  max_edges <- n * (n - 1) / 2
  k <- round(dens * max_edges)
  
  ut <- mat[upper.tri(mat)]
  ut_sorted <- sort(ut, decreasing = TRUE)
  
  if (k < 1) k <- 1
  if (k > length(ut_sorted)) k <- length(ut_sorted)
  
  thr <- ut_sorted[k]
  
  out <- matrix(0, n, n)
  out[mat >= thr] <- mat[mat >= thr]
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  diag(out) <- 0
  
  return(out)
}

densities <- c(0.10, 0.15, 0.20, 0.25)


# path to files


paths <- list(
  alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
  alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
  alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
  beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
  beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
  beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)



#load matrixes and align the labels

mats <- lapply(paths, load_matrix)
labels <- rownames(mats[[1]])
mats <- lapply(mats, function(m) m[labels, labels])


#making the average of the graph

build_force_graph <- function(condition) {
  W <- mats[[condition]]
  
  W_list <- list()
  for (d in densities) {
    W_list[[as.character(d)]] <- keep_top_density(W, d)
  }
  
  # average adjacency matrix
  W_avg <- Reduce("+", W_list) / length(densities)
  
  # build graph
  g <- graph_from_adjacency_matrix(W_avg,
                                   mode = "undirected",
                                   weighted = TRUE,
                                   diag = FALSE)
  
  # clustering coefficient
  cl <- suppressWarnings(transitivity(g, type = "localundirected",
                                      weights = E(g)$weight))
  cl[is.na(cl)] <- 0
  V(g)$clust <- cl
  
  # community detection
  com <- cluster_louvain(g, weights = E(g)$weight)
  V(g)$community <- com$membership
  
  list(graph = g, modularity = modularity(com))
}

#plot graph

plot_force_graph <- function(condition) {
  res <- build_force_graph(condition)
  g <- res$graph
  mod <- res$modularity
  
  # rescale edge weights
  w <- E(g)$weight
  ws <- (w - min(w)) / (max(w) - min(w) + 1e-6)
  E(g)$ws <- 0.3 + 2.0 * ws
  
#   p <- ggraph(g, layout = "fr") +       # <-- this is with weighted edges
#     geom_edge_link(aes(width = ws),
#                    alpha = 0.25,
#                    colour = "grey40") +
#     geom_node_point(aes(size = clust,
#                         colour = factor(community))) +
#     geom_node_text(aes(label = name), 
#                    size = 3, 
#                    repel = TRUE) +
#     scale_size(range = c(2, 10)) +
#     guides(size = guide_legend("Clustering"),
#            colour = guide_legend("Module")) +
#     theme_void() +
#     ggtitle(paste0("Force-Directed Network: ", condition,
#                    "\nModularity = ", round(mod, 3)))
  
p <- ggraph(g, layout = "fr") +       # <-- this is without the weights of edges
  geom_edge_link(alpha = 0.15, colour = "grey50", width = 0.3) +   # << simplified edges
  geom_node_point(aes(size = clust, colour = factor(community))) +
  geom_node_text(aes(label = name), size = 3, repel = TRUE) +
  scale_size(range = c(2, 10)) +
  guides(size = guide_legend("Clustering"),
         colour = guide_legend("Module")) +
  theme_void() +
#   ggtitle(paste0("Force-Directed Network: ", condition,
#                  "\nModularity = ", round(mod, 3)))
ggtitle(paste0("Force-Directed Network: ", condition))


  print(p)
  return(p)
}

#graphs to run:

 #plot_force_graph("alpha_med1")
# plot_force_graph("alpha_med2")
# plot_force_graph("alpha_thinking")
# plot_force_graph("beta_med1")
# plot_force_graph("beta_med2")
# plot_force_graph("beta_thinking")

