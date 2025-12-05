library(igraph)

source("R/analysis.R")
exists("keep_top_density")
exists("network_metrics")
exists("load_matrix")

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
fp_connectivity <- function(mat, frontal, posterior, densities = c(0.1, 0.15, 0.2, 0.25)) {
    fp_values <- c()

    for (d in densities) {
        threshold_matrix <- keep_top_density(mat, d)

        f_ids <- which(rownames(threshold_matrix) %in% frontal)
        p_ids <- which(rownames(threshold_matrix) %in% posterior)

        fp_sub <- threshold_matrix[f_ids, p_ids, drop = FALSE]
        values <- fp_sub[fp_sub > 0]

        if (length(values) == 0) {
            fp_values <- c(fp_values, 0)
        }
        else {
            fp_values <- c(fp_values, mean(values))
        }
    }

    return(mean(fp_values))
}


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
    colors <- rep("gray80", length(channels))
    colors[f_ids] <- "dodgerblue3"
    colors[p_ids] <- "firebrick2"

    plot(
        g,
        vertex.color = colors,
        vertex.size = 6,
        vertex.label = NA,
        edge.width = sqrt(E(g)$weight) * 4,
        edge.color = "purple",
        layout = layout_in_circle,
        main = title
    )
}

plot_fp_edges_thresholded <- function(mat, frontal, posterior, densities = c(0.10, 0.15, 0.20, 0.25), title = "") {
    fp_mats <- list()

    for (d in densities) {
        threshold <- keep_top_density(mat, d)

        f_ids <- which(rownames(threshold) %in% frontal)
        p_ids <- which(rownames(threshold) %in% posterior)

        fp <- matrix(0, nrow(threshold), ncol(threshold))
        rownames(fp) <- rownames(threshold)
        colnames(fp) <- colnames(threshold)

        fp[f_ids, p_ids] <- threshold[f_ids, p_ids]
        fp[p_ids, f_ids] <- t(threshold[f_ids, p_ids])

        fp_mats[[length(fp_mats) + 1]] <- fp
    }

    average_fp <- Reduce("+", fp_mats) / length(fp_mats)

    g <- graph_from_adjacency_matrix(average_fp, mode = "undirected", weighted = TRUE)

    channels <- rownames(average_fp)
    colors <- rep("gray80", length(channels))
    colors[channels %in% frontal] <- "dodgerblue3"
    colors[channels %in% posterior] <- "firebrick2"

    plot(
        g,
        vertex.color = colors,
        vertex.size = 6,
        vertex.label = NA,
        edge.width = sqrt(E(g)$weight) * 6,
        edge.color = "purple",
        layout = layout_in_circle,
        main = title
    )
}

results <- data.frame(Condition=names(paths), FP_connectivity=NA)

# Thresholded FP metric
for (i in seq_along(paths)) {
    mat <- load_matrix(paths[[i]])
    results$FP_connectivity[i] <- fp_connectivity(mat, frontal, posterior)
}

print(results)

png("fp_raw_6panel.png", width=3000, height=2000, res=250)
par(mfrow=c(2,3), mar=c(1,1,2,1))

for (name in names(paths)) {
    plot_fp_edges_raw(load_matrix(paths[[name]]), frontal, posterior, title=name)
}
dev.off()


png("fp_thresholded_6panel.png", width=3000, height=2000, res=250)
par(mfrow=c(2,3), mar=c(1,1,2,1))

for (name in names(paths)) {
    plot_fp_edges_thresholded(load_matrix(paths[[name]]), frontal, posterior, title=name)
}
dev.off()
