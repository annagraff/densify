install.packages("ggplot2")  # Install the package if you haven't
library(ggplot2)  # Load the package

# Sorted row means plot
sorted_row_means <- sort(rowMeans(full_matrix), decreasing = TRUE)
row_means_df <- data.frame(Mean = sorted_row_means)
ggplot(row_means_df, aes(x = reorder(Mean, Mean), y = Mean)) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  labs(title = "Sorted Row Means", x = "Row", y = "Mean") +
  theme_minimal()

# Sorted column means plot
sorted_col_means <- sort(colMeans(full_matrix), decreasing = TRUE)
col_means_df <- data.frame(Mean = sorted_col_means)
ggplot(col_means_df, aes(x = reorder(Mean, Mean), y = Mean)) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  labs(title = "Sorted Column Means", x = "Column", y = "Mean") +
  theme_minimal()
