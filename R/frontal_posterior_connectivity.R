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

# compute the frontal-posterior connectivity
fp_connectivity <- function(mat, frontal, posterior, densities = c(0.1, 0.15, 0.2, 0.25)) {
    fp_vals <- c()

    for (d in densities) {
        mat_d <- keep_top_density(mat, d)

        channels <- rownames(mat_d)
        f_ids <- which(channels %in% frontal)
        p_ids <- which(channels %in% posterior)

        submat <- mat_d[f_ids, p_ids, drop = FALSE]

        val <- mean(submat[submat > 0], na.rm = TRUE)
        if (is.nan(val)) val <- 0
        fp_vals <- c(fp_vals, val)
    }

    return(mean(fp_vals, na.rm = TRUE))
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
# this function plots the entire network, all edges, not just FP
# Not needed for FP analysis, but useful for comparison
# To run it
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

# to see the brain regions plot
mat <- load_matrix(paths$beta_thinking) # switch the path if you want a different group plotted
plot_brain_regions(mat, frontal, posterior, title="")


# plot the edges raw (no thresholding)
# This shows the full connectivity pattern that naturally exists from the raw (averaged) data in the csv files
# For this particular question, it is better to use the raw data because the data is naturally being cut down by only considering the posteror and frontal electrodes
plot_fp_edges_raw <- function(mat, frontal, posterior, title="") {
    channels <- rownames(mat)
    f_ids <- which(channels %in% frontal)
    p_ids <- which(channels %in% posterior)

    fp_mat <- matrix(0, nrow(mat), ncol(mat))
    rownames(fp_mat) <- channels
    colnames(fp_mat) <- channels

    fp_mat[f_ids, p_ids] <- mat[f_ids, p_ids]
    fp_mat[p_ids, f_ids] <- t(mat[f_ids, p_ids])

    g <- graph_from_adjacency_matrix(fp_mat, mode="undirected", weighted = TRUE)

    # Node colors
    col <- rep("gray80", length(channels))
    col[f_ids] <- "dodgerblue3"
    col[p_ids] <- "firebrick2"

    plot(
        g,
        vertex.color = col,
        vertex.size = 6,
        vertex.label = NA,
        edge.width = sqrt(E(g)$weight) * 5,
        edge.color = "purple",
        layout = layout_in_circle,
        main = title
    )
}

# matrix thresholding helper
# Applies the same proportional density thresholding method used in analysis.R: sort all edges, keep only top % (10%, 15%, 20%, 25%)
# then extracts ONLY the frontal-posterior edges from the thresholded matrix

# this allows the FP plots to match the team's thresholding methodology used elsewhere in the project
threshold_fp_matrix <- function(mat, frontal, posterior, dens) {

    n <- nrow(mat)
    diag(mat) <- 0
    mat[is.na(mat)] <- 0
    mat <- (mat + t(mat)) / 2

    max_edges <- n * (n-1) / 2
    k <- round(dens * max_edges)

    ut <- mat[upper.tri(mat)]
    thr <- sort(ut, decreasing=TRUE)[k]

    thr_mat <- matrix(0, n, n)
    thr_mat[mat >= thr] <- mat[mat >= thr]
    diag(thr_mat) <- 0

    channels <- rownames(mat)
    f_ids <- which(channels %in% frontal)
    p_ids <- which(channels %in% posterior)

    fp_mat <- matrix(0, nrow(mat), ncol(mat))
    rownames(fp_mat) <- channels
    colnames(fp_mat) <- channels

    fp_mat[f_ids, p_ids] <- thr_mat[f_ids, p_ids]
    fp_mat[p_ids, f_ids] <- t(thr_mat[f_ids, p_ids])

    return(fp_mat)
}

# thresholded plot, same thresolding as done in analysis.R
# implements the multidensity threshold averaging
# Thresholds the matrix at 10%, 15%, 20%, 25%, thn extracts the FP edges from it, averages the FP adjacency matrices, then plots only the strongest edges

# this matches the methodology done in analysis.R
plot_fp_edges_multidensity <- function(mat, frontal, posterior, title="") {

    densities <- c(0.10, 0.15, 0.20, 0.25)

    fp_mats <- lapply(densities, function(d)
        threshold_fp_matrix(mat, frontal, posterior, d)
    )

    avg_fp_mat <- Reduce("+", fp_mats) / length(fp_mats)

    g <- graph_from_adjacency_matrix(avg_fp_mat, mode="undirected", weighted=TRUE)

    channels <- rownames(mat)
    f_ids <- which(channels %in% frontal)
    p_ids <- which(channels %in% posterior)

    col <- rep("gray80", length(channels))
    col[f_ids] <- "dodgerblue3"
    col[p_ids] <- "firebrick2"

    # highlight strongest FP edges
    weights <- E(g)$weight
    keep_thr <- sort(weights, decreasing=TRUE)[max(1, round(length(weights)*0.3))]

    E(g)$color <- ifelse(E(g)$weight >= keep_thr, "purple", NA)
    E(g)$width <- ifelse(E(g)$weight >= keep_thr, sqrt(E(g)$weight)*6, 0)

    plot(
        g,
        vertex.color = col,
        vertex.size = 6,
        vertex.label = NA,
        layout = layout_in_circle,
        main = title
    )
}


png("fp_raw_6panel.png", width=3000, height=2000, res=250)
par(mfrow=c(2,3))
par(mar=c(1,1,2,1))
par(oma=c(3,0,3,0)) 

for (name in names(paths)) {
    mat <- load_matrix(paths[[name]])
    plot_fp_edges_raw(mat, frontal, posterior, title=name)
}

mtext("Raw Frontal–Posterior Connectivity (All FP Edges, No Thresholding)",
      outer=TRUE, cex=2.2, font=2)

dev.off()

png("fp_thresholded_multidensity_6panel.png", width=3000, height=2000, res=250)
par(mfrow=c(2,3))
par(mar=c(1,1,2,1))
par(oma=c(3,0,3,0)) 

for (name in names(paths)) {
    mat <- load_matrix(paths[[name]])
    plot_fp_edges_multidensity(mat, frontal, posterior, title=name)
}

mtext("Frontal–Posterior Connectivity (Multi-Density Thresholding: 10–25%)",
      outer=TRUE, cex=2.2, font=2)

dev.off()