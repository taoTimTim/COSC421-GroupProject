library(ggplot2)
library(reshape2)

# function to load a matrix from CSV
load_matrix <- function(path) {
  df <- read.csv(path, header = TRUE, check.names = FALSE)
  rownames(df) <- df[,1]  
  df <- df[,-1]
  mat <- as.matrix(df)
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- rownames(df)
  return(mat)
}

# File paths for averaged matrices
paths <- list(
  alpha_med1     = "src/results/averages/alpha_med1breath/alpha_med1breath_wpli_average_sub1-2.csv",
  alpha_med2     = "src/results/averages/alpha_med2/alpha_med2_wpli_average_sub1-2.csv",
  alpha_thinking = "src/results/averages/alpha_thinking/alpha_thinking_wpli_average_sub1-2.csv",
  beta_med1      = "src/results/averages/beta_med1breath/beta_med1breath_wpli_average_sub1-2.csv",
  beta_med2      = "src/results/averages/beta_med2/beta_med2_wpli_average_sub1-2.csv",
  beta_thinking  = "src/results/averages/beta_thinking/beta_thinking_wpli_average_sub1-2.csv"
)

# create a heatmap for a matrix and save it as PNG
plot_heatmap <- function(mat, title, filename) {

  # edge list , pair of electrodes and their connectivity value
  edge_list_df <- melt(mat)
  colnames(edge_list_df) <- c("Electrode1", "Electrode2", "Connectivity")
  
  # create the heatmap
  p <- ggplot(edge_list_df, aes(x = Electrode1, y = Electrode2, fill = Connectivity)) +
    geom_tile() + scale_fill_gradient(low = "white", high = "red") + theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 5),
      axis.title = element_blank()
    ) + labs(title = title)
  
  # save as PNG
  ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
}


# loop over all conditions
for (name in names(paths)) {
  mat <- load_matrix(paths[[name]])
  plot_heatmap(mat, title = name, filename = paste0(name, "_heatmap.png"))
}
