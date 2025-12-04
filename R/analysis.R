



library(igraph)
load_matrix <- function(path) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
    rownames(df) <- df[,1]
    df <- df[,-1]
    mat <- as.matrix(df)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rownames(df)
    mat
}
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
    mat_thr
}
densities <- c(0.10, 0.15, 0.20, 0.25)
paths <- list(
    alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
    alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
    alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
    beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
    beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
    beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)
mins <- lapply(paths, load_matrix)
final_results <- list()
for (nm in names(mins)) {
    mat <- mins[[nm]]
    stuff <- list()
    for (d in densities) {
        mat_d <- keep_top_density(mat, d)
        g_d <- graph_from_adjacency_matrix(mat_d, mode="undirected", weighted=TRUE)
        stuff[[as.character(d)]] <- network_metrics(g_d)
    }
    M <- do.call(rbind, stuff)
    avg <- colMeans(M[, -1], na.rm = TRUE)
    final_results[[nm]] <- avg
    cat("\n", nm, "\n")
    print(round(avg, 4))
}
par(mfrow = c(2,3))
par(mfrow = c(2,3))
for (nm in names(mins)) {
    mat_15 <- keep_top_density(mins[[nm]], 0.15)
    g_15 <- graph_from_adjacency_matrix(mat_15, mode="undirected", weighted=TRUE)
    plot(g_15, main = nm, layout = layout_in_circle, vertex.label = NA, vertex.size = 4, edge.width = sqrt(E(g_15)$weight) * 2, edge.color = "gray60")
}
