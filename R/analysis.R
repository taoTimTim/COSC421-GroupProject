# This script loads wpli matrices and turns them into brain networks.
# eeg gives fully connected matrices, so we need to remove weak edges or it's just noise
# Instead of picking one cutoff, we test multiple densities (10%, 15%, 20%, 25%).
# For each density we keep the strongest edges, build the graph, and get metrics,
# then we average the metrics so we're not relying on one arbitrary threshold (so the data is less likely to be skewed).
# this is standard in connectivity research and keeps things fair between edges
# Here are the metrics we care about (based on the project proposal): strength, clustering, modularity, betweenness

# 15% graph is used for the plots because it looks clean(ish) and still shows the structure

# Quick note on the "multiple densities" thing:
# Instead of picking one % of edges to keep (like only the top 5%), we try a few levels (10%, 15%, 20%, 25%)
# This lets us see if the results are consistent no matterwhat the cutoff is (cutoff for which edges we don't use/aren't good eough)
# If we only used one %, our data might be skewed and think there's an effect when it's just the threshold, so averaging across several thresholds
# makes the results way more reliable


library(igraph)

# load matrix from csv
load_matrix <- function(path) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
    
    # First column is channel labels, set row names
    rownames(df) <- df[,1]
    df <- df[,-1]

    mat <- as.matrix(df)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rownames(df)
    return(mat)
}

# network metrics
network_metrics <- function(g) {
    data.frame(
        nodes = gorder(g),
        edges = gsize(g),
        density = edge_density(g),
        mean_strength = mean(strength(g)),
        clustering = mean(transitivity(g, type = "localundirected", weights = E(g)$weight), na.rm = TRUE),
        modularity = modularity(cluster_louvain(g, weights = E(g)$weight)),
        avg_betw = mean(betweenness(g, normalized = TRUE))
    )
}

# keep top x% edges (density thresholding)
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

densities <- c(0.10, 0.15, 0.20, 0.25)

# file paths
paths <- list(
    alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
    alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
    alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
    beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
    beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
    beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

# load and build graphs
mats <- lapply(paths, load_matrix)

final_results <- list()

for (name in names(mats)) {
    mat <- mats[[name]]
    metric_rows <- list()

    for (d in densities) {
        mat_d <- keep_top_density(mat, d)
        g_d <- graph_from_adjacency_matrix(mat_d, mode="undirected", weighted=TRUE)
        metric_rows[[as.character(d)]] <- network_metrics(g_d)
    }

    M <- do.call(rbind, metric_rows)
    avg <- colMeans(M[, -1], na.rm = TRUE)
    final_results[[name]] <- avg

    cat("\n", name, "\n")
    print(round(avg, 4))
}

# plot all graphs in a simple layout
par(mfrow = c(2,3))

# plot network at 15% density
par(mfrow = c(2,3))
for (name in names(mats)) {
    mat_15 <- keep_top_density(mats[[name]], 0.15)
    g_15 <- graph_from_adjacency_matrix(mat_15, mode="undirected", weighted=TRUE)

    plot(
        g_15,
        main = name,
        layout = layout_in_circle,
        vertex.label = NA,
        vertex.size = 4,
        edge.width = sqrt(E(g_15)$weight) * 2,
        edge.color = "gray60"
    )
}
