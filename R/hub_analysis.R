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

cond_band <- function(nm) ifelse(grepl("alpha", nm), "alpha", "beta")
cond_state <- function(nm) {
    if (grepl("med1", nm)) return("med1")
    if (grepl("med2", nm)) return("med2")
    if (grepl("thinking", nm)) return("thinking")
    "unknown"
}
cond_meta <- data.frame(
    condition = names(paths),
    band = sapply(names(paths), cond_band),
    state = sapply(names(paths), cond_state),
    stringsAsFactors = FALSE
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

cat("\nTop 5 hubs by condition (band + state)\n")
for (nm in names(mats)) {
    meta <- cond_meta[cond_meta$condition == nm, ]
    cat("\n", nm, " (band: ", meta$band, ", state: ", meta$state, ")\n", sep = "")
    t5 <- head(hubs_by_condition[[nm]], 5)
    t5$condition <- nm
    t5$band <- meta$band
    t5$state <- meta$state
    t5 <- t5[, c("condition", "band", "state", "electrode", "degree", "strength", "betweenness")]
    print(t5)
}

cat("\nsaving all hubs\n")
all_hubs_list <- lapply(names(hubs_by_condition), function(n) {
    df <- hubs_by_condition[[n]]
    if (is.null(df)) return(NULL)
    meta <- cond_meta[cond_meta$condition == n, ]
    df$condition <- n
    df$band <- meta$band
    df$state <- meta$state
    df
})
all_hubs_list <- all_hubs_list[!sapply(all_hubs_list, is.null)]
if (length(all_hubs_list) == 0) stop("no hub data")
all_hubs <- do.call("rbind", all_hubs_list)
rownames(all_hubs) <- NULL
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
    meta <- cond_meta[cond_meta$condition == nm, ]
    t5 <- head(df, 5)
    t5$band <- meta$band
    t5$state <- meta$state
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

# precompute band max for consistent y axis within each band
band_max <- list(alpha = 0, beta = 0)
for (nm in names(top5_list)) {
    df <- top5_list[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    # same band same y so folks can eyeball easy
    b <- cond_band(nm)
    band_max[[b]] <- max(band_max[[b]], max(df$strength, na.rm = TRUE))
}
if (!is.finite(band_max$alpha) || band_max$alpha <= 0) band_max$alpha <- 1
if (!is.finite(band_max$beta) || band_max$beta <= 0) band_max$beta <- 1

pdf("src/results/hubs/plots/top5_hubs_per_condition.pdf", width = 14, height = 10)
par(mfrow = c(2, 3))
for (nm in names(paths)) {
    df <- top5_list[[nm]]
    meta <- cond_meta[cond_meta$condition == nm, ]
    if (is.null(df) || nrow(df) == 0) {
        plot.new(); title(main = paste("No data:", nm))
    } else {
        ylim_max <- if (meta$band == "alpha") band_max$alpha else band_max$beta
        barplot(df$strength,
                names.arg = df$electrode,
                main = paste("Top 5 Hubs:", nm, "(", meta$band, meta$state, ")"),
                ylab = "Strength",
                col = "steelblue",
                las = 2,
                ylim = c(0, ylim_max * 1.05))
    }
}
dev.off()
cat("  top5_hubs_per_condition.pdf\n")

pdf("src/results/hubs/plots/hub_region_distribution.pdf", width = 10, height = 6)
region_counts <- all_hubs %>% group_by(band, state, region) %>% summarise(count = n(), .groups = "drop")
regions_order <- c("frontal","central","parietal","occipital","temporal","other")
mk <- function(band_label) {
    dat <- region_counts %>% filter(band == band_label)
    regs <- regions_order[regions_order %in% dat$region]
    states <- c("med1","med2","thinking")
    states <- states[states %in% dat$state]
    if (length(states) == 0 || length(regs) == 0) return(matrix(0, nrow = 0, ncol = 0))
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
if (length(alpha_mat)) {
    barplot(alpha_mat, beside = TRUE, col = c("steelblue","skyblue3","coral")[seq_len(nrow(alpha_mat))], main = "Alpha Band: Hub Distribution by Region", ylab = "Number of Hub Electrodes", las = 2)
    legend("topright", legend = rownames(alpha_mat), fill = c("steelblue","skyblue3","coral")[seq_len(nrow(alpha_mat))])
} else {
    plot.new(); title(main = "Alpha: no data")
}
if (length(beta_mat)) {
    barplot(beta_mat, beside = TRUE, col = c("steelblue","skyblue3","coral")[seq_len(nrow(beta_mat))], main = "Beta Band: Hub Distribution by Region", ylab = "Number of Hub Electrodes", las = 2)
    legend("topright", legend = rownames(beta_mat), fill = c("steelblue","skyblue3","coral")[seq_len(nrow(beta_mat))])
} else {
    plot.new(); title(main = "Beta: no data")
}
dev.off()
cat("  hub_region_distribution.pdf\n")

pdf("src/results/hubs/plots/meditation_vs_thinking_strength.pdf", width = 12, height = 6)
par(mfrow = c(2, 2), mar = c(8, 4, 4, 1))

alpha_med1 <- all_hubs %>% filter(band == "alpha", state == "med1") %>% arrange(-strength) %>% head(10)
alpha_med2 <- all_hubs %>% filter(band == "alpha", state == "med2") %>% arrange(-strength) %>% head(10)
alpha_think <- all_hubs %>% filter(band == "alpha", state == "thinking") %>% arrange(-strength) %>% head(10)
beta_med1 <- all_hubs %>% filter(band == "beta", state == "med1") %>% arrange(-strength) %>% head(10)
beta_med2 <- all_hubs %>% filter(band == "beta", state == "med2") %>% arrange(-strength) %>% head(10)
beta_think <- all_hubs %>% filter(band == "beta", state == "thinking") %>% arrange(-strength) %>% head(10)

# use a consistent y-axis per band for easier comparisons
band_compare_max <- list(
    alpha = max(c(alpha_med1$strength, alpha_med2$strength, alpha_think$strength), na.rm = TRUE),
    beta  = max(c(beta_med1$strength, beta_med2$strength, beta_think$strength), na.rm = TRUE)
)
# dont mix med1/med2 just line em up fair
if (!is.finite(band_compare_max$alpha) || band_compare_max$alpha <= 0) band_compare_max$alpha <- 1
if (!is.finite(band_compare_max$beta)  || band_compare_max$beta <= 0)  band_compare_max$beta  <- 1

plot_state_vs_thinking <- function(band_label, state_label, med_df, think_df, ylim_max) {
    vals <- c(med_df$strength, think_df$strength)
    labs <- c(paste0(med_df$electrode, " (", state_label, ")"),
              paste0(think_df$electrode, " (thinking)"))
    cols <- c(rep("steelblue", nrow(med_df)), rep("coral", nrow(think_df)))
    if (length(vals) == 0 || all(!is.finite(vals))) {
        plot.new(); title(main = paste(band_label, state_label, "vs thinking: no data"))
        return()
    }
    barplot(vals,
            names.arg = labs,
            col = cols,
            main = paste0(band_label, " band: ", state_label, " vs thinking"),
            ylab = "Strength",
            las = 2,
            ylim = c(0, ylim_max * 1.05))
    legend("topright", c(state_label, "thinking"), fill = c("steelblue", "coral"))
}

plot_state_vs_thinking("Alpha", "med1", alpha_med1, alpha_think, band_compare_max$alpha)
plot_state_vs_thinking("Alpha", "med2", alpha_med2, alpha_think, band_compare_max$alpha)
plot_state_vs_thinking("Beta", "med1", beta_med1, beta_think, band_compare_max$beta)
plot_state_vs_thinking("Beta", "med2", beta_med2, beta_think, band_compare_max$beta)

dev.off()
cat("  meditation_vs_thinking_strength.pdf\n")

cat("\ndone\n")
