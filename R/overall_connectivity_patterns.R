library(igraph)
library(ggplot2)

#Load Matrix
load_matrix <- function(path) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
    rownames(df) <- df[,1]
    df <- df[,-1]
    mat <- as.matrix(df)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rownames(df)
    return(mat)
}

# density thresholding
keep_top_density <- function(mat, dens) {
    diag(mat) <- 0
    mat[is.na(mat)] <- 0
    mat <- (mat + t(mat)) / 2

    n <- nrow(mat)
    max_edges <- n * (n - 1) / 2
    k <- round(dens * max_edges)

    ut <- mat[upper.tri(mat)]
    thr <- sort(ut, decreasing = TRUE)[k]

    bin <- matrix(0, n, n)
    bin[mat >= thr] <- 1
    diag(bin) <- 0
    return(bin)
}

# Pearson similarity
network_vector <- function(bin_mat) {
    return(bin_mat[upper.tri(bin_mat)])
}

pearson_similarity <- function(mat1, mat2) {
    v1 <- network_vector(mat1)
    v2 <- network_vector(mat2)
    return(cor(v1, v2))
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

# Load all matrices
mats <- lapply(paths, load_matrix)


# store binary matrices for similarity 
binary_mats_by_thresh <- list()


# For each condition
for (name in names(mats)) {
  binary_mats_by_thresh[[name]] <- list() 
  for (d in densities) {
    bin_mat <- keep_top_density(mats[[name]], d)
    binary_mats_by_thresh[[name]][[as.character(d)]] <- bin_mat 
  }
}

# Similarity calculations
similarity_results <- data.frame(
    Pair = c("alpha_med1 vs alpha_thinking",
             "alpha_med2 vs alpha_thinking",
             "beta_med1 vs beta_thinking",
             "beta_med2 vs beta_thinking"),
    Avg_Similarity = NA
)

pairs <- list(
    c("alpha_med1", "alpha_thinking"),
    c("alpha_med2", "alpha_thinking"),
    c("beta_med1", "beta_thinking"),
    c("beta_med2", "beta_thinking")
)

for (i in seq_along(pairs)) {
    p <- pairs[[i]]
    sims <- c()

    for (d in densities) {
        bin1 <- binary_mats_by_thresh[[p[1]]][[as.character(d)]]
        bin2 <- binary_mats_by_thresh[[p[2]]][[as.character(d)]]
        sims <- c(sims, pearson_similarity(bin1, bin2))
    }

    similarity_results$Avg_Similarity[i] <- mean(sims)
}

print("average pearson similarities")
print(similarity_results)


# Similarity Bar Plot
png("Q3a_similarity_meditation_vs_thinking.png", width = 2000, height = 1500, res = 200)

ggplot(similarity_results, aes(x = Pair, y = Avg_Similarity)) +
    geom_bar(stat="identity", fill="darkorange") +
    theme_minimal() +
    labs(title="Meditation vs Thinking - Average Network Similarity",
         x="Pair",
         y="Pearson Similarity") +
    theme(axis.text.x = element_text(angle=30, hjust=1))

dev.off()