# Install and load necessary packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)
library(RColorBrewer)

# Get ColorBrewer Set1 colors
set1_colors <- brewer.pal(9, "Set1")

# Sorted row means plot
sorted_row_means <- sort(rowMeans(full_matrix), decreasing = TRUE)
row_means_df <- data.frame(Mean = sorted_row_means)
row_means_df$Row <- seq_len(nrow(row_means_df))

ggplot(row_means_df, aes(x = Row, y = Mean)) +
  geom_point(color = set1_colors[1], size = 1) +
  labs(title = "Sorted Row Densities", x = "Row", y = "Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Sorted column densities plot with a different color
sorted_col_densities <- sort(colMeans(full_matrix), decreasing = TRUE)
col_densities_df <- data.frame(Density = sorted_col_densities)

ggplot(col_densities_df, aes(x = seq_along(Density), y = Density)) +
  geom_point(color = set1_colors[2], size = 1) +
  labs(title = "Sorted Column Densities", x = "Column", y = "Density") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
