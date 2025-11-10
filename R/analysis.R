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

# build the graph
matrix_to_graph <- function(mat) {
    diag(mat) <- 0
    mat[is.na(mat)] <- 0
    mat <- (mat + t(mat))/2 # force symmetry

    # NOTE: thresholding is commented out right now.
    # When thresholding was enabled, all six graphs ended up with identical 
    # node count, edge count, and mean degree. That makes interpretation unhelpful.
    # With no threshold, the graphs are extremely dense (almost fully connected),
    # so the edge structure is hard to visually interpret.
    # We need to decide an appropriate thresholding method for visualization and analysis
    # For now, I have kept all raw connectivity values, this is purely to get the data into networks
    
    # there are a TON of edges, so keep only stronger connections
    # thr <- quantile(mat[mat > 0], 0.95)  # top 5%
    # mat[mat < thr] <- 0

    graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)
}

# file pats
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
graphs <- lapply(mats, matrix_to_graph)

# quick metrics (raw)
metrics <- function(g) {
    data.frame(
        nodes = gorder(g),
        edges = gsize(g),
        mean_degree = mean(degree(g))
    )
}

results <- do.call(rbind, lapply(graphs, metrics))
rownames(results) <- names(graphs)
print(results)

# plot all graphs in a simple layout
par(mfrow = c(2,3))

for (name in names(graphs)) {
    plot(
        graphs[[name]],
        main = name,
        layout = layout_in_circle,
        vertex.label = NA,
        vertex.size = 4,
        edge.width = E(graphs[[name]])$weight * 2,
        edge.color = "gray"
    )
}
