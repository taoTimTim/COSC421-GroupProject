library(igraph)

## load matrix from csv (first col = labels)
load_mat <- function(path) {
  df <- read.csv(path, header = TRUE, check.names = FALSE)

  # rownames come from first col
  rn <- df[,1]
  df <- df[,-1]

  m <- as.matrix(df)
  m <- apply(m, 2, as.numeric)   # make sure everything is numeric
  rownames(m) <- rn

  m
}

## keep top X% edges by weight (undirected)
keep_top <- function(mat, dens = 0.15) {

  diag(mat) <- 0
  mat[is.na(mat)] <- 0

  # symmetrize (wPLI should already be but just in case)
  mat <- (mat + t(mat)) * 0.5

  n <- nrow(mat)
  total_edges <- n * (n - 1) / 2
  k <- round(dens * total_edges)
  if (k < 1) k <- 1

  ut <- mat[upper.tri(mat)]
  thr <- sort(ut, decreasing = TRUE)[k]

  out <- matrix(0, n, n)
  keep <- mat >= thr
  out[keep] <- mat[keep]
  diag(out) <- 0

  out
}

## quick plotting wrapper
make_plot <- function(mat, name, out_path) {

  m <- keep_top(mat, 0.15)
  g <- graph_from_adjacency_matrix(m, mode="undirected",
                                   weighted=TRUE, diag=FALSE)

  s <- strength(g, weights = E(g)$weight)

  # pick roughly top 5%
  n_top <- ceiling(0.05 * length(s))
  top_idx <- order(s, decreasing=TRUE)[1:n_top]

  col_vec <- rep("gray70", length(s))
  col_vec[top_idx] <- "red"

  png(out_path, width=2000, height=2000, res=300)
  par(mar=c(0,0,2,0))

  plot(
    g,
    layout = layout_in_circle(g),
    vertex.label = V(g)$name,
    vertex.label.cex = 0.55,
    vertex.size = 6,
    vertex.color = col_vec,
    edge.width = sqrt(E(g)$weight) * 1.8,
    edge.color = "gray60",
    main = paste(name, "(15% threshold)")
  )

  dev.off()
}

## ------------------ main loop ------------------

paths <- list(
  meditation_alpha_1 = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
  meditation_alpha_2 = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
  thinking_alpha     = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
  meditation_beta_1  = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
  meditation_beta_2  = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
  thinking_beta      = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

for (nm in names(paths)) {
  cat("plotting:", nm, "\n")
  mat <- load_mat(paths[[nm]])


  out_file <- file.path("figures", paste0(gsub("_", "", nm), "_15pct.png"))

  make_plot(mat, nm, out_file)
}

cat("done.\n")
