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


#######################
# Calculate densities for the matrices separately
original_densities <- rowMeans(full_matrix)
pruned_densities <- rowMeans(pruned_matrix)

# Merge the matrices based on the common languages
merged_densities <- merge(original_densities, pruned_densities, by = "row.names", all = TRUE)

# Rename the columns
colnames(merged_densities) <- c("Language", "Original_Density", "Pruned_Density")

# Create the densities data frame
densities_df <- data.frame(merged_densities)

# Sort the data frame by original density
densities_df <- densities_df[order(densities_df$Original_Density, decreasing = TRUE), ]

# Plot densities before and after pruning
ggplot(densities_df, aes(x = reorder(Language, -Original_Density), y = Original_Density)) +
  geom_point(aes(color = "Original Density"), size = 1) +
  geom_point(data = densities_df[!is.na(densities_df$Pruned_Density), ],
             aes(x = reorder(Language, -Original_Density), y = Pruned_Density, color = "Pruned Density"), size = 1) +
  labs(title = "Density Before and After Pruning", x = "Language", y = "Density") +
  scale_color_manual(values = c("Original Density" = set1_colors[1], "Pruned Density" = set1_colors[2])) +
  theme_minimal() +
  theme(axis.text.x = element_blank())


