#calculating the clustering and modularity for the bar plot

library(igraph)

# Load adjacency matrix from CSV -----------------------------
load_matrix <- function(path) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
    rownames(df) <- df[,1]
    df <- df[,-1]
    mat <- as.matrix(df)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rownames(df)
    return(mat)
}

# Network metrics (global) -----------------------------------
network_metrics <- function(g) {
    data.frame(
        mean_strength = mean(strength(g)),
        clustering = mean(transitivity(g, type="localundirected",
                                       weights = E(g)$weight), na.rm = TRUE),
        modularity = modularity(cluster_louvain(g, weights = E(g)$weight))
    )
}

# Keep top X% edges (density thresholding) -------------------
keep_top_density <- function(mat, dens) {
    diag(mat) <- 0
    mat[is.na(mat)] <- 0
    mat <- (mat + t(mat)) / 2

    n <- nrow(mat)
    max_edges <- n * (n - 1) / 2
    k <- round(dens * max_edges)

    ut <- mat[upper.tri(mat)]
    thr <- sort(ut, decreasing = TRUE)[k]

    mat_thr <- matrix(0, n, n)
    mat_thr[mat >= thr] <- mat[mat >= thr]
    diag(mat_thr) <- 0
    return(mat_thr)
}

# Densities to average across --------------------------------
densities <- c(0.10, 0.15, 0.20, 0.25)

# Paths to your six matrices ---------------------------------
paths <- list(
    alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
    alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
    alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
    beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
    beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
    beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

# Load matrices ----------------------------------------------
mats <- lapply(paths, load_matrix)

# Compute final metrics --------------------------------------
final_results <- data.frame(
    condition = character(),
    mean_strength = numeric(),
    clustering = numeric(),
    modularity = numeric(),
    stringsAsFactors = FALSE
)

for (name in names(mats)) {
    mat <- mats[[name]]
    metric_rows <- list()

    for (d in densities) {
        mat_d <- keep_top_density(mat, d)
        g_d <- graph_from_adjacency_matrix(mat_d, mode="undirected", weighted=TRUE)
        metric_rows[[as.character(d)]] <- network_metrics(g_d)
    }

    # Combine → average across densities
    M <- do.call(rbind, metric_rows)
    avg_row <- colMeans(M[, sapply(M, is.numeric)], na.rm = TRUE)


    final_results <- rbind(final_results,
                           data.frame(condition = name,
                                      t(avg_row),
                                      row.names = NULL))
}

# Print table 
cat("\n===== FINAL GLOBAL METRIC RESULTS (Average across densities) =====\n")

# Round only numeric columns
final_results_rounded <- final_results
numeric_cols <- sapply(final_results_rounded, is.numeric)
final_results_rounded[, numeric_cols] <- round(final_results_rounded[, numeric_cols], 4)

print(final_results_rounded)


# Base R barplot for modularity
barplot(
    height = final_results$modularity,
    names.arg = final_results$condition,
    las = 2,
    main = "Modularity Across Conditions",
    ylab = "Modularity",
    cex.names = 0.8
)


# ===============================
# BASE R BAR PLOTS + SAVING
# ===============================

# ----------- 1. Save Modularity Barplot -----------

png("modularity_barplot.png", width = 1400, height = 900, res = 150)
barplot(
    height = final_results$modularity,
    names.arg = final_results$condition,
    las = 2,
    main = "Modularity Across Conditions",
    ylab = "Modularity",
    col = "skyblue",
    cex.names = 0.8
)
dev.off()

pdf("modularity_barplot.pdf", width = 10, height = 6)
barplot(
    height = final_results$modularity,
    names.arg = final_results$condition,
    las = 2,
    main = "Modularity Across Conditions",
    ylab = "Modularity",
    col = "skyblue",
    cex.names = 0.8
)
dev.off()


# ----------- 2. Save Clustering Barplot -----------

png("clustering_barplot.png", width = 1400, height = 900, res = 150)
barplot(
    height = final_results$clustering,
    names.arg = final_results$condition,
    las = 2,
    main = "Clustering Across Conditions",
    ylab = "Clustering Coefficient",
    col = "palegreen3",
    cex.names = 0.8
)
dev.off()

pdf("clustering_barplot.pdf", width = 10, height = 6)
barplot(
    height = final_results$clustering,
    names.arg = final_results$condition,
    las = 2,
    main = "Clustering Across Conditions",
    ylab = "Clustering Coefficient",
    col = "palegreen3",
    cex.names = 0.8
)
dev.off()
