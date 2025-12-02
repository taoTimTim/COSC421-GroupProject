library(igraph)

# load matrix from csv file
load_matrix <- function(path) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
    rownames(df) <- df[,1] # first column contains the names
    df <- df[,-1]

    mat <- as.matrix(df)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rownames(df)
    return(mat)
}

# define the frontal and posterior electrodes (basically figure out whether an electrode is part of the frontal or posterior part of the brain)
frontal <- c("Fp1","Fp2","AF3","AF4","F3","F4","F7","F8","Fz")
posterior <- c("P3","P4","P7","P8","Pz","PO3","PO4","O1","O2")

paths <- list(
    alpha_med1 = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
    alpha_med2 = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
    alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
    beta_med1 = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
    beta_med2 = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
    beta_thinking = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

# compute the frontal-posterior connectivity
fp_connectivity <- function(mat, frontal, posterior) {
    channels <- rownames(mat)

    frontal_ids <- which(channels %in% frontal)
    posterior_ids <- which(channels %in% posterior)

    # slice the matrix for these connections only, get rid of noise
    submat <- mat[frontal_ids, posterior_ids, drop = FALSE]

    return(mean(submat, na.rm = TRUE))
}


# run code for all groups
results <- data.frame(
    Condition = names(paths),
    FP_connectivity = NA
)

for (i in seq_along(paths)) {
    name <- names(paths)[i]
    mat <- load_matrix(paths[[i]])
    fp <- fp_connectivity(mat, frontal, posterior)
    results$FP_connectivity[i] <- fp
}

# print the results
cat("\n===== Frontal-Posterior Connectivity Results =====\n")
print(results)


# plotting cool graphs fr
plot_brain_regions <- function(mat, frontal, posterior, title = "Network") {
    g <- graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE)

    channels <- rownames(mat)

    # assigning colours
    col <- rep("gray70", length(channels))
    col[channels %in% frontal] <- "blue"
    col[channels %in% posterior] <- "red"

    # node sizes
    vs <- rep(5, length(channels))

    # frontal is blue, posterior is red, gray maps to other regions in the brain (central, temporal, etc.)
    plot(
        g,
        vertex.color = col,
        vertex.size = vs,
        vertex.label.cex = 0.6,
        vertex.label.color = "black",
        edge.width = sqrt(E(g)$weight) * 2,
        edge.color = "gray80",
        layout = layout_in_circle,
        main = title
    )
}

plot_fp_edges <- function(mat, frontal, posterior, title="") {

    channels <- rownames(mat)

    # get indices
    f_ids <- which(channels %in% frontal)
    p_ids <- which(channels %in% posterior)

    # FP-only matrix
    fp_mat <- matrix(0, nrow(mat), ncol(mat))
    rownames(fp_mat) <- channels
    colnames(fp_mat) <- channels

    fp_mat[f_ids, p_ids] <- mat[f_ids, p_ids]
    fp_mat[p_ids, f_ids] <- mat[p_ids, f_ids]

    g <- graph_from_adjacency_matrix(fp_mat, mode="undirected", weighted=TRUE)

    # Colors
    col <- rep("gray80", length(channels))
    col[f_ids] <- "blue" # frontal
    col[p_ids] <- "red" # posterior

    plot(
        g,
        vertex.color = col,
        vertex.size = 6,
        vertex.label = NA,
        edge.width = sqrt(E(g)$weight) * 4,
        edge.color = "purple",
        layout = layout_in_circle,
        main = ""
    )

    title (
        title,
        cex.main = 2.8,
        font.main = 2,
        line = -2
    )
}

# Layout: 2 rows, 3 columns
par(mfrow=c(2,3))

# Margins for each subplot
par(mar=c(1, 1, 4, 1)) # bottom, left, top, right

# Better spacing between plots
par(oma=c(1,1,1,1))

for (name in names(paths)) {
    mat <- load_matrix(paths[[name]])
    plot_fp_edges(mat, frontal, posterior, title=name)
}

