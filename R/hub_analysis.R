library(igraph)
library(dplyr)

load_matrix <- function(path) {
    df <- read.csv(path, header = TRUE, check.names = FALSE)
    rownames(df) <- df[,1]
    df <- df[,-1]
    mat <- as.matrix(df)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rownames(df)
    mat
}

keep_top_density <- function(mat, dens) {
    diag(mat) <- 0
    mat[is.na(mat)] <- 0
    mat <- (mat + t(mat)) / 2
    n <- nrow(mat)
    max_edges <- n * (n - 1) / 2
    k <- max(1, round(dens * max_edges))
    ut <- mat[upper.tri(mat)]
    thr <- sort(ut, decreasing = TRUE)[k]
    out <- matrix(0, n, n)
    out[mat >= thr] <- mat[mat >= thr]
    diag(out) <- 0
    out
}

densities <- c(0.10, 0.15, 0.20, 0.25)

hub_analysis <- function(mat, name) {
    cat("\nprocessing", name, "\n")
    labs <- rownames(mat)
    cat("  electrodes:", length(labs), "\n")
    dm <- list()
    for (d in densities) {
        cat("  density:", d * 100, "%\n")
        m <- keep_top_density(mat, d)
        rownames(m) <- labs
        colnames(m) <- labs
        g <- graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE, diag = FALSE)
        dm[[as.character(d)]] <- list(
            degree = degree(g),
            strength = strength(g, weights = E(g)$weight),
            betweenness = betweenness(g, normalized = TRUE)
        )
    }
    if (length(dm) == 0) stop("no graphs")
    cat("  densities used:", paste(names(dm), collapse = ", "), "\n")
    grab <- function(k) {
        tmp <- do.call(cbind, lapply(dm, function(x) x[[k]]))
        rowMeans(tmp, na.rm = TRUE)
    }
    df <- data.frame(
        electrode = labs,
        degree = as.numeric(grab("degree")),
        strength = as.numeric(grab("strength")),
        betweenness = as.numeric(grab("betweenness")),
        stringsAsFactors = FALSE
    )
    df <- df[order(df$strength, decreasing = TRUE), ]
    rownames(df) <- NULL
    cat("  rows:", nrow(df), "\n")
    df
}

paths <- list(
    alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
    alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
    alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
    beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
    beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
    beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

cat("\nload mats\n")
mats <- lapply(paths, function(p) {
    cat(" ", p, "\n")
    load_matrix(p)
})
names(mats) <- names(paths)

cat("\nrun hubs\n")
hubs_by_condition <- list()
for (nm in names(mats)) {
    cat("\n--", nm, "--\n")
    hubs_by_condition[[nm]] <- hub_analysis(mats[[nm]], nm)
}

cat("\nalpha quick\n")
print(head(hubs_by_condition$alpha_thinking, 5))
print(head(hubs_by_condition$alpha_med1, 5))
print(head(hubs_by_condition$alpha_med2, 5))

cat("\nbeta quick\n")
print(head(hubs_by_condition$beta_thinking, 5))
print(head(hubs_by_condition$beta_med1, 5))
print(head(hubs_by_condition$beta_med2, 5))

cat("\nsaving all hubs\n")
all_hubs_list <- lapply(names(hubs_by_condition), function(n) {
    df <- hubs_by_condition[[n]]
    if (is.null(df)) return(NULL)
    df$condition <- n
    df
})
all_hubs_list <- all_hubs_list[!sapply(all_hubs_list, is.null)]
if (length(all_hubs_list) == 0) stop("no hub data")
all_hubs <- do.call("rbind", all_hubs_list)
rownames(all_hubs) <- NULL
all_hubs$band <- ifelse(grepl("alpha", all_hubs$condition), "alpha", "beta")
all_hubs$state <- ifelse(grepl("thinking", all_hubs$condition), "thinking", "meditation")
cat("  combined rows:", nrow(all_hubs), "\n")

region_map <- list(
    frontal = c("Fp1","Fp2","F3","F4","F7","F8","Fz","AF3","AF4","AF7","AF8"),
    central = c("C3","C4","Cz","CP1","CP2","CP3","CP4","CP5","CP6"),
    parietal = c("P3","P4","Pz","PO3","PO4","PO7","PO8"),
    occipital = c("O1","O2","Oz","Iz"),
    temporal = c("T3","T4","T5","T6","TP7","TP8","TP9","TP10")
)
get_region <- function(e) {
    for (r in names(region_map)) {
        if (e %in% region_map[[r]]) return(r)
    }
    "other"
}
all_hubs$region <- sapply(all_hubs$electrode, get_region)

dir.create("src/results/hubs", recursive = TRUE, showWarnings = FALSE)
write.csv(all_hubs, "src/results/hubs/all_hub_metrics_multi_density.csv", row.names = FALSE)
cat("  saved all_hub_metrics_multi_density.csv\n")

cat("\nexporting top5\n")
dir.create("src/results/hubs/top5", recursive = TRUE, showWarnings = FALSE)
top5_list <- list()
for (nm in names(hubs_by_condition)) {
    df <- hubs_by_condition[[nm]]
    if (is.null(df)) next
    t5 <- head(df, 5)
    t5$band <- ifelse(grepl("alpha", nm), "alpha", "beta")
    t5$state <- ifelse(grepl("thinking", nm), "thinking", "meditation")
    t5$region <- sapply(t5$electrode, get_region)
    top5_list[[nm]] <- t5
    write.csv(t5, paste0("src/results/hubs/top5/top5_", nm, ".csv"), row.names = FALSE)
}
cat("  saved top5 csvs\n")

cat("\nregion summary\n")
region_summary <- all_hubs %>%
    group_by(band, state, region) %>%
    summarise(
        avg_strength = mean(strength, na.rm = TRUE),
        avg_degree = mean(degree, na.rm = TRUE),
        avg_betweenness = mean(betweenness, na.rm = TRUE),
        hub_count = n(),
        .groups = "drop"
    ) %>%
    arrange(band, state, -avg_strength)
write.csv(region_summary, "src/results/hubs/region_summary.csv", row.names = FALSE)
print(region_summary)

cat("\nplots\n")
dir.create("src/results/hubs/plots", recursive = TRUE, showWarnings = FALSE)

pdf("src/results/hubs/plots/top5_hubs_per_condition.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))
for (nm in names(top5_list)) {
    df <- top5_list[[nm]]
    barplot(df$strength, names.arg = df$electrode, main = paste("Top 5 Hubs:", nm), ylab = "Strength", col = "steelblue", las = 2)
}
dev.off()
cat("  top5_hubs_per_condition.pdf\n")

pdf("src/results/hubs/plots/hub_region_distribution.pdf", width = 10, height = 6)
region_counts <- all_hubs %>% group_by(band, state, region) %>% summarise(count = n(), .groups = "drop")
regions_order <- c("frontal","central","parietal","occipital","temporal","other")
mk <- function(band_label) {
    dat <- region_counts %>% filter(band == band_label)
    regs <- regions_order[regions_order %in% dat$region]
    states <- c("meditation","thinking")
    m <- matrix(0, nrow = length(states), ncol = length(regs))
    rownames(m) <- states; colnames(m) <- regs
    for (i in seq_along(states)) {
        for (j in seq_along(regs)) {
            v <- dat$count[dat$state == states[i] & dat$region == regs[j]]
            if (length(v) == 0) v <- 0
            m[i, j] <- v
        }
    }
    m
}
alpha_mat <- mk("alpha"); beta_mat <- mk("beta")
par(mfrow = c(1, 2))
barplot(alpha_mat, beside = TRUE, col = c("red","blue"), main = "Alpha Band: Hub Distribution by Region", ylab = "Number of Hub Electrodes", las = 2)
legend("topright", legend = rownames(alpha_mat), fill = c("red","blue"))
barplot(beta_mat, beside = TRUE, col = c("red","blue"), main = "Beta Band: Hub Distribution by Region", ylab = "Number of Hub Electrodes", las = 2)
legend("topright", legend = rownames(beta_mat), fill = c("red","blue"))
dev.off()
cat("  hub_region_distribution.pdf\n")

pdf("src/results/hubs/plots/meditation_vs_thinking_strength.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
alpha_med <- all_hubs %>% filter(band == "alpha", state == "meditation") %>% arrange(-strength) %>% head(10)
alpha_think <- all_hubs %>% filter(band == "alpha", state == "thinking") %>% arrange(-strength) %>% head(10)
beta_med <- all_hubs %>% filter(band == "beta", state == "meditation") %>% arrange(-strength) %>% head(10)
beta_think <- all_hubs %>% filter(band == "beta", state == "thinking") %>% arrange(-strength) %>% head(10)
if (nrow(alpha_med) > 0 && nrow(alpha_think) > 0) {
    barplot(c(alpha_med$strength, alpha_think$strength), names.arg = c(alpha_med$electrode, alpha_think$electrode), col = c(rep("steelblue", nrow(alpha_med)), rep("coral", nrow(alpha_think))), main = "Alpha Band: Top 10 Hubs (Meditation vs Thinking)", ylab = "Strength", las = 2)
    legend("topright", c("Meditation", "Thinking"), fill = c("steelblue", "coral"))
}
if (nrow(beta_med) > 0 && nrow(beta_think) > 0) {
    barplot(c(beta_med$strength, beta_think$strength), names.arg = c(beta_med$electrode, beta_think$electrode), col = c(rep("steelblue", nrow(beta_med)), rep("coral", nrow(beta_think))), main = "Beta Band: Top 10 Hubs (Meditation vs Thinking)", ylab = "Strength", las = 2)
    legend("topright", c("Meditation", "Thinking"), fill = c("steelblue", "coral"))
}
dev.off()
cat("  meditation_vs_thinking_strength.pdf\n")

cat("\ndone\n")
