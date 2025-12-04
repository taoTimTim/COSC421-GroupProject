library(igraph)

load_matrix <- function(path) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
    rownames(df) <- df[,1] # first column contains the names
    df <- df[,-1]

    mat <- as.matrix(df)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rownames(df)
    return(mat)
}

frontal <- c("Fp1","Fp2","AF3","AF4","F3","F4","F7","F8","Fz")
posterior <- c("P3","P4","P7","P8","Pz","PO3","PO4","O1","O2")

region_connectivity <- function(mat, region) {
    ids <- which(rownames(mat) %in% region)
    submat <- mat[ids, ids, drop=FALSE]
    diag(submat) <- 0
    return(mean(submat, na.rm = TRUE))
}

paths <- list(
    alpha_med1 = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
    alpha_med2 = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
    alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
    beta_med1 = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
    beta_med2 = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
    beta_thinking = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)
results <- data.frame(
    Condition = names(paths),
    Frontal_Frontal = NA,
    Posterior_Posterior = NA
)

for (i in seq_along(paths)) {
    mat <- load_matrix(paths[[i]])
    results$Frontal_Frontal[i] <- region_connectivity(mat, frontal)
    results$Posterior_Posterior[i] <- region_connectivity(mat, posterior)
}

print(results)



pdf("Q4_region_connectivity_barlot.png", width = 10, height = 8)

mat <- t(as.matrix(results[, c("Frontal_Frontal", "Posterior_Posterior")]))
colnames(mat) <- results$Condition
rownames(mat) <- c("Frontal-Frontal", "Posterior-Posterior")

barplot(mat,
        beside = TRUE,           
        col = c("dodgerblue3","firebrick2"),
        ylim = c(0, max(mat) * 1.2),
        main = "Average Within-Region Connectivity (F-F vs P-P)",
        ylab = "Average Connectivity",
        xlab = "Condition",
        legend.text = TRUE,
        args.legend = list(x = "topright", bty="n"))


dev.off()
