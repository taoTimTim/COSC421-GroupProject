library(igraph)

load_mat <- function(path) {
  df <- read.csv(path, header = TRUE, check.names = FALSE)
  rownames(df) <- df[,1]
  df <- df[,-1]
  m <- as.matrix(df)
  m <- apply(m, 2, as.numeric)
  rownames(m) <- rownames(df)
  m
}

keep_top <- function(mat, dens = 0.15) {
  diag(mat) <- 0
  mat[is.na(mat)] <- 0
  mat <- (mat + t(mat)) / 2
  n <- nrow(mat)
  max_edges <- n * (n - 1) / 2
  k <- max(1, round(dens * max_edges))
  ut <- mat[upper.tri(mat)]
  thr <- sort(ut, decreasing = TRUE)[k]
  out <- matrix(0, n, n)
  out[mat >= thr] <- mat[mat >= thr]
  diag(out) <- 0
  out
}

make_plot <- function(mat, name, out_path) {
  m <- keep_top(mat, 0.15)
  g <- graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE, diag = FALSE)
  s <- strength(g, weights = E(g)$weight)
  n_top <- max(1, ceiling(0.05 * length(s)))
  idx <- order(s, decreasing = TRUE)[seq_len(n_top)]
  cols <- rep("gray70", length(s))
  cols[idx] <- "red"
  png(out_path, width = 2000, height = 2000, res = 300)
  par(mar = c(0,0,0,0))
  plot(g,
       layout = layout_in_circle(g),
       vertex.label = V(g)$name,
       vertex.label.cex = 0.5,
       vertex.size = 6,
       vertex.color = cols,
       edge.width = sqrt(E(g)$weight) * 2,
       edge.color = "gray60",
       axes = FALSE,
       main = paste(name, "- 15%"))
  dev.off()
}

paths <- list(
  meditation_alpha_1 = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
  meditation_alpha_2 = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
  thinking_alpha     = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
  meditation_beta_1  = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
  meditation_beta_2  = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
  thinking_beta      = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

for (nm in names(paths)) {
  cat("plotting", nm, "\n")
  mat <- load_mat(paths[[nm]])
  out_file <- file.path("figures", paste0(gsub("_", "", nm), "_network_15.png"))
  make_plot(mat, nm, out_file)
}

cat("done\n")
