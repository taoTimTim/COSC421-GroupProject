#This code shows clustrering and modularity with the coorinates of the eeg elecgtrodes.

library(igraph)
library(ggraph)
library(ggplot2)
library(gridExtra)


#loading the eeg cordinates 
coords <- read.csv("R/Q2/eeg_64coords.csv")



#load wpli matrix
load_matrix <- function(path) {
  df <- read.csv(path, header = TRUE, check.names = FALSE)
  rownames(df) <- df[,1]
  df <- df[,-1]
  m <- as.matrix(df)
  m <- apply(m, 2, as.numeric)
  rownames(m) <- rownames(df)
  return(m)
}


# threshold function
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


#load data
paths <- list(
  alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
  alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
  alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
  beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
  beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
  beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

mats <- lapply(paths, load_matrix)

labels <- rownames(mats[[1]])
mats <- lapply(mats, function(m) m[labels, labels])
coords <- coords[match(labels, coords$electrode), ]


#make the average graph
build_avg_graph <- function(condition) {
  W <- mats[[condition]]
  
  W_list <- list()
  for (d in densities) {
    W_list[[as.character(d)]] <- keep_top_density(W, d)
  }
  
  W_avg <- Reduce("+", W_list) / length(densities)
  
  g <- graph_from_adjacency_matrix(W_avg, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  V(g)$name <- labels
  V(g)$x <- coords$x
  V(g)$y <- coords$y
  
  # clustering coefficient
  cl <- suppressWarnings(transitivity(g, type = "localundirected", weights = E(g)$weight))
  cl[is.na(cl)] <- 0
  V(g)$clust <- cl
  
  # modularity community detection
  com <- cluster_louvain(g, weights = E(g)$weight)
  V(g)$community <- com$membership
  
  list(graph = g, modularity = modularity(com))
}


#plot graph
plot_condition <- function(condition, title_prefix = "") {
  res <- build_avg_graph(condition)
  g <- res$graph
  
  w <- E(g)$weight
  ws <- (w - min(w)) / (max(w) - min(w) + 1e-6)
  E(g)$ws <- 0.5 + 2.5 * ws
  
  p <- ggraph(g, layout = "manual", x = V(g)$x, y = V(g)$y) +
    geom_edge_link(alpha = 0.15, colour = "grey50", width = 0.3) +
    geom_node_point(aes(size = clust, colour = factor(community))) +
    geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
    scale_size(range = c(2, 7)) +
    guides(size = guide_legend("Clustering"),
           colour = guide_legend("Module")) +
    theme_void() +
    ggtitle(paste0(title_prefix, condition))
  
  return(p)
}


#save the graph 
save_plot_condition <- function(condition, prefix = "", width = 7, height = 6, dpi = 300) {
  p <- plot_condition(condition, prefix)
  
  safe_prefix <- gsub("[^A-Za-z0-9_-]", "", prefix)
  
  filename_png <- paste0(safe_prefix, condition, ".png") #getting the png file of the graph 
  filename_pdf <- paste0(safe_prefix, condition, ".pdf")
  
  ggsave(filename_png, p, width = width, height = height, dpi = dpi)
  ggsave(filename_pdf, p, width = width, height = height)
  
  message("Saved: ", filename_png, " and ", filename_pdf)
}


#save all plots

save_plot_condition("alpha_med1", "Alpha_")
#save_plot_condition("alpha_med2", "Alpha_")
#save_plot_condition("alpha_thinking", "Alpha_")

#save_plot_condition("beta_med1", "Beta_")
#save_plot_condition("beta_med2", "Beta_")
#save_plot_condition("beta_thinking", "Beta_")
