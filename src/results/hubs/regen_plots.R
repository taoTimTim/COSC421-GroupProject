# Regenerate PNG plots from CSV outputs
# Produces: top5_hubs_per_condition.png, hub_region_distribution.png, meditation_vs_thinking_strength.png

# Read top5 CSVs
top5_dir <- "src/results/hubs/top5"
top5_files <- list.files(top5_dir, pattern = "top5_.*\\.csv$", full.names = TRUE)

top5_list <- list()
for (f in top5_files) {
  name <- gsub("^.*/top5_|\\.csv$", "", f)
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  top5_list[[name]] <- df
}

# 1) Top 5 hubs per condition (2x3)
png("src/results/hubs/plots/top5_hubs_per_condition.png", width = 1400, height = 1000, res = 150)
par(mfrow = c(2,3), mar = c(8,4,4,2))
for (n in names(top5_list)) {
  df <- top5_list[[n]]
  if (!is.null(df) && nrow(df) > 0) {
    barplot(df$strength,
            names.arg = df$electrode,
            main = paste("Top 5 Hubs:", n),
            ylab = "Strength",
            col = "steelblue",
            las = 2)
  } else {
    plot.new(); title(main = paste("No data:", n))
  }
}
dev.off()

# 2) Hub region distribution (alpha and beta side-by-side)
region_summary <- read.csv("src/results/hubs/region_summary.csv", stringsAsFactors = FALSE)
regions <- unique(region_summary$region)
# Ensure consistent ordering
regions_order <- c("frontal","central","parietal","occipital","temporal","other")
regions <- regions_order[regions_order %in% regions]

# Build matrices for alpha and beta with rows = states (meditation, thinking), cols = regions
make_mat <- function(band) {
  dat <- region_summary[region_summary$band == band, ]
  states <- c("meditation", "thinking")
  m <- matrix(0, nrow = length(states), ncol = length(regions))
  rownames(m) <- states
  colnames(m) <- regions
  for (i in seq_along(states)) {
    s <- states[i]
    for (j in seq_along(regions)) {
      r <- regions[j]
      v <- dat$count[dat$state == s & dat$region == r]
      if (length(v) == 0) v <- 0
      m[i,j] <- v
    }
  }
  return(m)
}

# But region_summary may not have 'count' column (it had hub_count). Try both
if ("count" %in% names(region_summary)) {
  # nothing
} else if ("hub_count" %in% names(region_summary)) {
  names(region_summary)[names(region_summary) == "hub_count"] <- "count"
}

alpha_mat <- make_mat("alpha")
beta_mat <- make_mat("beta")

png("src/results/hubs/plots/hub_region_distribution.png", width = 1200, height = 600, res = 150)
par(mfrow = c(1,2), mar = c(6,4,4,2))
barplot(alpha_mat, beside = TRUE, col = c("red","blue"), main = "Alpha Band: Hub Distribution by Region", ylab = "Number of Hub Electrodes", las = 2)
legend("topright", legend = rownames(alpha_mat), fill = c("red","blue"))
barplot(beta_mat, beside = TRUE, col = c("red","blue"), main = "Beta Band: Hub Distribution by Region", ylab = "Number of Hub Electrodes", las = 2)
legend("topright", legend = rownames(beta_mat), fill = c("red","blue"))
dev.off()

# 3) Meditation vs Thinking strength comparison (top10 each)
all_hubs <- read.csv("src/results/hubs/all_hub_metrics_15pct.csv", stringsAsFactors = FALSE)

png("src/results/hubs/plots/meditation_vs_thinking_strength.png", width = 1200, height = 600, res = 150)
par(mfrow = c(1,2), mar = c(8,4,4,2))
# Alpha
alpha_med <- subset(all_hubs, band == "alpha" & state == "meditation")[order(-all_hubs$strength[all_hubs$band=="alpha" & all_hubs$state=="meditation"]), ]
alpha_med <- head(alpha_med, 10)
alpha_think <- subset(all_hubs, band == "alpha" & state == "thinking")[order(-all_hubs$strength[all_hubs$band=="alpha" & all_hubs$state=="thinking"]), ]
alpha_think <- head(alpha_think, 10)
if (nrow(alpha_med) > 0 && nrow(alpha_think) > 0) {
  barplot(c(alpha_med$strength, alpha_think$strength), names.arg = c(alpha_med$electrode, alpha_think$electrode), col = c(rep("steelblue", nrow(alpha_med)), rep("coral", nrow(alpha_think))), main = "Alpha Band: Top 10 Hubs (Meditation vs Thinking)", ylab = "Strength", las = 2)
  legend("topright", c("Meditation", "Thinking"), fill = c("steelblue", "coral"))
} else {
  plot.new(); title(main = "Alpha: insufficient data")
}

# Beta
beta_med <- subset(all_hubs, band == "beta" & state == "meditation")[order(-all_hubs$strength[all_hubs$band=="beta" & all_hubs$state=="meditation"]), ]
beta_med <- head(beta_med, 10)
beta_think <- subset(all_hubs, band == "beta" & state == "thinking")[order(-all_hubs$strength[all_hubs$band=="beta" & all_hubs$state=="thinking"]), ]
beta_think <- head(beta_think, 10)
if (nrow(beta_med) > 0 && nrow(beta_think) > 0) {
  barplot(c(beta_med$strength, beta_think$strength), names.arg = c(beta_med$electrode, beta_think$electrode), col = c(rep("steelblue", nrow(beta_med)), rep("coral", nrow(beta_think))), main = "Beta Band: Top 10 Hubs (Meditation vs Thinking)", ylab = "Strength", las = 2)
  legend("topright", c("Meditation", "Thinking"), fill = c("steelblue", "coral"))
} else {
  plot.new(); title(main = "Beta: insufficient data")
}
dev.off()

cat("PNG regeneration completed\n")
