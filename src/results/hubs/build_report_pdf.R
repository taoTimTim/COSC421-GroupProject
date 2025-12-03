# Build a combined PDF report embedding text content and PNG figures
out_pdf <- "src/results/hubs/hub_report_combined.pdf"
png_dir <- "src/results/hubs/plots"

# helper to wrap text into multiple lines
wrap_lines <- function(text, width=90) {
  paste(strwrap(text, width=width), collapse="\n")
}

# Check for png package
if (!requireNamespace("png", quietly = TRUE)) {
  cat("Package 'png' is required but not installed. Please install it: install.packages('png')\n")
  quit(status=1)
}

library(png)
library(grid)

pdf(out_pdf, width=11, height=8.5)
on.exit(dev.off(), add = TRUE)

# Title page
plot.new()
title(main = "Hub Electrode Analysis - Research Question 1", cex.main = 1.6)
text(x=0.5, y=0.6, labels = paste("Generated:", Sys.Date()), cex=1.0)
text(x=0.5, y=0.5, labels = "Source: src/results/hubs/", cex=0.9)

# Short summary page
plot.new()
summary_text <- c(
  "Summary:",
  "This report summarizes the hub electrode analysis comparing meditation and thinking states for Alpha and Beta bands.",
  "Centrality measures used: Degree, Strength (weighted), Betweenness (normalized).",
  "Thresholding for visualization: top 15% strongest edges (proportional density)."
)
text(x=0, y=1, labels = wrap_lines(paste(summary_text, collapse="\n\n"), width=80), adj = c(0,1), cex=0.9)

# Region summary table page
region_file <- "src/results/hubs/region_summary.csv"
if (file.exists(region_file)) {
  reg <- read.csv(region_file, stringsAsFactors = FALSE)
  # try to render as table using gridExtra if available
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    grid.newpage()
    grid.text("Region Summary", gp=gpar(fontsize=14))
    grid.table(reg)
  } else {
    plot.new()
    txt <- capture.output(print(reg))
    text(x=0, y=1, labels = paste(txt, collapse="\n"), adj=c(0,1), cex=0.7)
  }
}

# Insert each PNG figure on its own page
png_files <- c(
  file.path(png_dir, "top5_hubs_per_condition.png"),
  file.path(png_dir, "hub_region_distribution.png"),
  file.path(png_dir, "meditation_vs_thinking_strength.png")
)
for (p in png_files) {
  if (file.exists(p)) {
    img <- readPNG(p)
    grid.newpage()
    grid.raster(img)
  } else {
    plot.new(); text(0.5,0.5, paste("Missing image:", p))
  }
}

# Top5 tables pages
top5_dir <- "src/results/hubs/top5"
if (dir.exists(top5_dir)) {
  files <- list.files(top5_dir, pattern = "top5_.*\\.csv$", full.names = TRUE)
  for (f in files) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    title_h <- basename(f)
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      grid.newpage()
      grid.text(title_h, gp=gpar(fontsize=12))
      grid.table(df)
    } else {
      plot.new()
      txt <- capture.output(print(df))
      text(x=0, y=1, labels = paste(txt, collapse="\n"), adj=c(0,1), cex=0.8)
    }
  }
}

dev.off()
cat("Created", out_pdf, "\n")
