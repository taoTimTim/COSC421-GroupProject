# Hub Electrode Analysis
# Research Question 1: Which electrodes act as hubs during meditation compared to thinking?
# This script identifies hub electrodes using degree, strength, and betweenness centrality

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

# keep top x% edges (density thresholding)
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

# Hub analysis function - compute centrality measures at 15% density
hub_analysis <- function(mat, name) {
    cat("\nProcessing:", name, "\n")
    
    # Apply 15% density thresholding
    mat_15 <- keep_top_density(mat, 0.15)
    
    # Store electrode names BEFORE any operations
    electrode_names <- rownames(mat_15)
    cat("  Electrodes found:", length(electrode_names), "\n")
    
    # Ensure matrix has proper dimnames
    colnames(mat_15) <- electrode_names
    
    # Convert to graph - handle errors better
    g_15 <- tryCatch({
        graph_from_adjacency_matrix(mat_15, mode="undirected", weighted=TRUE, diag=FALSE)
    }, error = function(e) {
        cat("  ERROR creating graph:", e$message, "\n")
        return(NULL)
    })
    
    # Check if graph creation failed
    if (is.null(g_15)) {
        cat("  Failed to create graph!\n")
        return(NULL)
    }
    
    # Get vertex count
    cat("  Graph vertices:", length(V(g_15)), "\n")
    
    # Calculate centrality measures
    degree_cent <- degree(g_15)
    strength_cent <- strength(g_15, weights = E(g_15)$weight)
    betweenness_cent <- betweenness(g_15, normalized = TRUE)
    
    cat("  Degree values:", length(degree_cent), "\n")
    cat("  Strength values:", length(strength_cent), "\n")
    
    # Create hub dataframe using the stored electrode names (not vertex names)
    hub_df <- data.frame(
        electrode = electrode_names,
        degree = as.numeric(degree_cent),
        strength = as.numeric(strength_cent),
        betweenness = as.numeric(betweenness_cent),
        stringsAsFactors = FALSE
    )
    
    # Sort by strength
    hub_df <- hub_df[order(hub_df$strength, decreasing = TRUE), ]
    rownames(hub_df) <- NULL
    
    cat("  Created hub dataframe with", nrow(hub_df), "rows\n")
    
    return(hub_df)
}

# File paths
paths <- list(
    alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
    alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
    alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
    beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
    beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
    beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

# Load matrices
cat("\n===== LOADING MATRICES =====\n")
mats <- list()
for (path_name in names(paths)) {
    tryCatch({
        cat("Loading:", path_name, "from", paths[[path_name]], "\n")
        mats[[path_name]] <- load_matrix(paths[[path_name]])
        cat("  ✓ Success - dimensions:", dim(mats[[path_name]]), "\n")
    }, error = function(e) {
        cat("  ✗ ERROR loading", path_name, ":", e$message, "\n")
    })
}

# ===== COMPUTE HUBS FOR ALL CONDITIONS =====
cat("\n===== HUB ELECTRODES ANALYSIS (15% DENSITY) =====\n")

hubs_by_condition <- list()

for (name in names(mats)) {
    cat("\n---\n")
    result <- tryCatch({
        hub_analysis(mats[[name]], name)
    }, error = function(e) {
        cat("  ✗ CRITICAL ERROR in hub_analysis():", e$message, "\n")
        cat("     Stack trace:", paste(e$call), "\n")
        return(NULL)
    })
    
    hubs_by_condition[[name]] <- result
    
    if (!is.null(result)) {
        cat("✓ Successfully created hub dataframe for", name, "\n")
        cat("\n", name, " - Top 10 Hub Electrodes:\n")
        print(head(result, 10))
    } else {
        cat("✗ FAILED to create hub dataframe for", name, "\n")
    }
}

# ===== ALPHA BAND COMPARISON =====
cat("\n\n===== ALPHA BAND COMPARISON =====\n")
cat("\nAlpha Thinking - Top 5 hubs:\n")
print(head(hubs_by_condition$alpha_thinking, 5))
cat("\nAlpha Meditation 1 - Top 5 hubs:\n")
print(head(hubs_by_condition$alpha_med1, 5))
cat("\nAlpha Meditation 2 - Top 5 hubs:\n")
print(head(hubs_by_condition$alpha_med2, 5))

# ===== BETA BAND COMPARISON =====
cat("\n\n===== BETA BAND COMPARISON =====\n")
cat("\nBeta Thinking - Top 5 hubs:\n")
print(head(hubs_by_condition$beta_thinking, 5))
cat("\nBeta Meditation 1 - Top 5 hubs:\n")
print(head(hubs_by_condition$beta_med1, 5))
cat("\nBeta Meditation 2 - Top 5 hubs:\n")
print(head(hubs_by_condition$beta_med2, 5))

# ===== 1. SAVE ALL HUB RESULTS INTO A SINGLE CSV =====
cat("\n\n===== SAVING HUB RESULTS =====\n")

# Check if we have any valid hub data
valid_hubs <- sum(!sapply(hubs_by_condition, is.null))
cat("Valid hub dataframes:", valid_hubs, "/", length(hubs_by_condition), "\n")

if (valid_hubs == 0) {
    cat("✗ ERROR: No valid hub dataframes to combine!\n")
    stop("Cannot proceed without hub data")
}

# Create results directory if it doesn't exist
tryCatch({
    dir.create("src/results/hubs", recursive = TRUE, showWarnings = FALSE)
    cat("✓ Created directory: src/results/hubs\n")
}, error = function(e) {
    cat("✗ ERROR creating directory:", e$message, "\n")
})

# Combine all hub data into one table - using do.call with rbind properly
tryCatch({
    all_hubs_list <- lapply(names(hubs_by_condition), function(n) {
        df <- hubs_by_condition[[n]]
        if (is.null(df)) {
            cat("  Skipping NULL dataframe:", n, "\n")
            return(NULL)
        }
        df$condition <- n
        return(df)
    })
    
    # Remove NULL entries
    all_hubs_list <- all_hubs_list[!sapply(all_hubs_list, is.null)]
    
    if (length(all_hubs_list) == 0) {
        stop("No valid dataframes to combine after filtering NULLs")
    }
    
    all_hubs <- do.call("rbind", all_hubs_list)
    rownames(all_hubs) <- NULL
    
    cat("✓ Combined hub data:", nrow(all_hubs), "total rows\n")
    cat("  Conditions included:", paste(unique(all_hubs$condition), collapse=", "), "\n")
}, error = function(e) {
    cat("✗ ERROR combining hub data:", e$message, "\n")
    stop("Cannot continue without combined hub data")
})

# ===== 2. ADD BAND AND STATE LABELS =====
tryCatch({
    all_hubs$band <- ifelse(grepl("alpha", all_hubs$condition), "alpha", "beta")
    all_hubs$state <- ifelse(grepl("thinking", all_hubs$condition), "thinking", "meditation")
    cat("✓ Added band and state labels\n")
}, error = function(e) {
    cat("✗ ERROR adding labels:", e$message, "\n")
})

# ===== 3. ADD ELECTRODE REGION LABELS =====
region_map <- list(
    frontal = c("Fp1","Fp2","F3","F4","F7","F8","Fz","AF3","AF4","AF7","AF8"),
    central = c("C3","C4","Cz","CP1","CP2","CP3","CP4","CP5","CP6"),
    parietal = c("P3","P4","Pz","PO3","PO4","PO7","PO8"),
    occipital = c("O1","O2","Oz","Iz"),
    temporal = c("T3","T4","T5","T6","TP7","TP8","TP9","TP10")
)

get_region <- function(elec) {
    for (region in names(region_map)) {
        if (elec %in% region_map[[region]]) return(region)
    }
    return("other")
}

tryCatch({
    all_hubs$region <- sapply(all_hubs$electrode, get_region)
    cat("✓ Added electrode region labels\n")
    cat("  Regions found:", paste(unique(all_hubs$region), collapse=", "), "\n")
}, error = function(e) {
    cat("✗ ERROR adding regions:", e$message, "\n")
})

# Save combined hub table
tryCatch({
    write.csv(all_hubs,
              "src/results/hubs/all_hub_metrics_15pct.csv",
              row.names = FALSE)
    cat("✓ Saved: all_hub_metrics_15pct.csv\n")
}, error = function(e) {
    cat("✗ ERROR saving all_hub_metrics_15pct.csv:", e$message, "\n")
})

# ===== 4. EXPORT TOP 5 HUBS PER CONDITION =====
cat("\n===== EXPORTING TOP 5 HUBS =====\n")

top5_list <- list()
for (name in names(hubs_by_condition)) {
    df <- hubs_by_condition[[name]]
    if (!is.null(df)) {
        top5_df <- head(df, 5)
        top5_df$band <- ifelse(grepl("alpha", name), "alpha", "beta")
        top5_df$state <- ifelse(grepl("thinking", name), "thinking", "meditation")
        top5_df$region <- sapply(top5_df$electrode, get_region)
        top5_list[[name]] <- top5_df
    } else {
        cat("  Skipping NULL dataframe:", name, "\n")
        top5_list[[name]] <- NULL
    }
}

dir.create("src/results/hubs/top5", recursive = TRUE, showWarnings = FALSE)

tryCatch({
    for (name in names(top5_list)) {
        if (!is.null(top5_list[[name]])) {
            write.csv(top5_list[[name]],
                      paste0("src/results/hubs/top5/top5_", name, ".csv"),
                      row.names = FALSE)
        }
    }
    cat("✓ Saved: top5_*.csv for each condition\n")
}, error = function(e) {
    cat("✗ ERROR saving top5 CSVs:", e$message, "\n")
})

# ===== 5. REGION SUMMARY ANALYSIS =====
cat("\n===== REGION SUMMARY ANALYSIS =====\n")

tryCatch({
    library(dplyr)
    cat("✓ dplyr library loaded\n")
}, error = function(e) {
    cat("✗ ERROR loading dplyr:", e$message, "\n")
    stop("Cannot proceed without dplyr")
})

tryCatch({
    region_summary <- all_hubs %>%
        group_by(band, state, region) %>%
        summarise(
            avg_strength = mean(strength, na.rm = TRUE),
            avg_degree = mean(degree, na.rm = TRUE),
            avg_betweenness = mean(betweenness, na.rm = TRUE),
            hub_count = n(),
            .groups = 'drop'
        ) %>%
        arrange(band, state, -avg_strength)
    
    cat("✓ Created region summary\n")
    cat("  Summary rows:", nrow(region_summary), "\n")
    
    write.csv(region_summary,
              "src/results/hubs/region_summary.csv",
              row.names = FALSE)
    
    cat("✓ Saved: region_summary.csv\n")
    cat("\nRegion Summary:\n")
    print(region_summary)
}, error = function(e) {
    cat("✗ ERROR creating region summary:", e$message, "\n")
})

# ===== 6. VISUALIZATIONS =====
cat("\n===== GENERATING VISUALIZATIONS =====\n")

# Create output directory for plots
dir.create("src/results/hubs/plots", recursive = TRUE, showWarnings = FALSE)

# Plot 1: Top 5 hubs per condition (strength barplot)
tryCatch({
    pdf("src/results/hubs/plots/top5_hubs_per_condition.pdf", width = 14, height = 10)
    par(mfrow = c(2, 3))
    
    for (name in names(top5_list)) {
        df <- top5_list[[name]]
        if (!is.null(df)) {
            barplot(df$strength,
                    names.arg = df$electrode,
                    main = paste("Top 5 Hubs:", name),
                    ylab = "Strength",
                    col = "steelblue",
                    las = 2)
        }
    }
    
    dev.off()
    cat("✓ Saved: top5_hubs_per_condition.pdf\n")
}, error = function(e) {
    cat("✗ ERROR creating top5_hubs_per_condition.pdf:", e$message, "\n")
})

# Plot 2: Region distribution of hubs
tryCatch({
    pdf("src/results/hubs/plots/hub_region_distribution.pdf", width = 10, height = 6)
    
    region_counts <- all_hubs %>%
        group_by(band, state, region) %>%
        summarise(count = n(), .groups = 'drop')
    
    par(mfrow = c(1, 2))
    
    # Alpha band
    alpha_data <- region_counts %>% filter(band == "alpha")
    if (nrow(alpha_data) > 0) {
        plot(alpha_data$region, alpha_data$count,
             main = "Alpha Band: Hub Distribution by Region",
             xlab = "Region", ylab = "Number of Hub Electrodes",
             type = "h", lwd = 3, col = c("red", "blue")[as.factor(alpha_data$state)])
    }
    
    # Beta band
    beta_data <- region_counts %>% filter(band == "beta")
    if (nrow(beta_data) > 0) {
        plot(beta_data$region, beta_data$count,
             main = "Beta Band: Hub Distribution by Region",
             xlab = "Region", ylab = "Number of Hub Electrodes",
             type = "h", lwd = 3, col = c("red", "blue")[as.factor(beta_data$state)])
    }
    
    dev.off()
    cat("✓ Saved: hub_region_distribution.pdf\n")
}, error = function(e) {
    cat("✗ ERROR creating hub_region_distribution.pdf:", e$message, "\n")
})

# Plot 3: Strength comparison - Meditation vs Thinking
tryCatch({
    pdf("src/results/hubs/plots/meditation_vs_thinking_strength.pdf", width = 12, height = 6)
    par(mfrow = c(1, 2))
    
    # Alpha
    alpha_med <- all_hubs %>% filter(band == "alpha", state == "meditation") %>% arrange(-strength) %>% head(10)
    alpha_think <- all_hubs %>% filter(band == "alpha", state == "thinking") %>% arrange(-strength) %>% head(10)
    
    if (nrow(alpha_med) > 0 && nrow(alpha_think) > 0) {
        barplot(c(alpha_med$strength, alpha_think$strength),
                names.arg = c(alpha_med$electrode, alpha_think$electrode),
                col = c(rep("steelblue", nrow(alpha_med)), rep("coral", nrow(alpha_think))),
                main = "Alpha Band: Top 10 Hubs (Meditation vs Thinking)",
                ylab = "Strength",
                las = 2)
        legend("topright", c("Meditation", "Thinking"), fill = c("steelblue", "coral"))
    }
    
    # Beta
    beta_med <- all_hubs %>% filter(band == "beta", state == "meditation") %>% arrange(-strength) %>% head(10)
    beta_think <- all_hubs %>% filter(band == "beta", state == "thinking") %>% arrange(-strength) %>% head(10)
    
    if (nrow(beta_med) > 0 && nrow(beta_think) > 0) {
        barplot(c(beta_med$strength, beta_think$strength),
                names.arg = c(beta_med$electrode, beta_think$electrode),
                col = c(rep("steelblue", nrow(beta_med)), rep("coral", nrow(beta_think))),
                main = "Beta Band: Top 10 Hubs (Meditation vs Thinking)",
                ylab = "Strength",
                las = 2)
        legend("topright", c("Meditation", "Thinking"), fill = c("steelblue", "coral"))
    }
    
    dev.off()
    cat("✓ Saved: meditation_vs_thinking_strength.pdf\n")
}, error = function(e) {
    cat("✗ ERROR creating meditation_vs_thinking_strength.pdf:", e$message, "\n")
})

cat("\n===== ALL RESULTS SAVED =====\n")
cat("Location: src/results/hubs/\n")
cat("Files created:\n")
cat("  - all_hub_metrics_15pct.csv (complete table)\n")
cat("  - top5/*.csv (top 5 hubs per condition)\n")
cat("  - region_summary.csv (region statistics)\n")
cat("  - plots/ (PDF visualizations)\n")
