# Step 1: Set the dimensions of the matrix (rows and columns)
num_rows <- 1084
num_cols <- 284

# Step 2: Set the desired proportions of 0s and 1s (e.g., 70% 0s and 30% 1s)
proportion_1s <- 0.7747
proportion_0s <- 1-proportion_1s

# Step 3: Calculate the number of 0s and 1s based on the desired proportions
num_zeros <- round(num_rows * num_cols * proportion_0s)
num_ones <- num_rows * num_cols - num_zeros

# Step 4: Generate random binary values (0 or 1) with the specified proportions
random_binary_values <- sample(c(rep(0, num_zeros), rep(1, num_ones)))

# Step 5: Reshape the vector of binary values into a matrix
full_matrix <- matrix(random_binary_values, nrow = num_rows, ncol = num_cols)

# Display the standard deviation
sqrt(var(rowSums(full_matrix)/ncol(full_matrix))+var(colSums(full_matrix)/nrow(full_matrix)))
