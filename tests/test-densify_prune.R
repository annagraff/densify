# Disclaimer: This is a first version of a test case generated with ChatGPT to help Marc learn the syntax. ChatGPT has limitations in accurately understanding" complex formats and dimensions, especially in real-world scenarios.


context("Test densify_prune function")

test_that("It correctly prunes and densifies the matrix", {
  # Create sample data for testing
  original_data <- data.frame(
    Glottocode = c("A", "B", "C", "D"),
    Var1 = c(1, 2, 3, 4),
    Var2 = c(5, 6, 7, 8)
  )

  documentation <- data.frame(
    iteration = 1:3,
    removed_tax = c("NA", "B", "NA"),
    removed_var = c("NA", "Var1", "NA")
  )

  # Test densify_prune function with an optimum of 2
  pruned_densified <- densify_prune(original_data, documentation, optimum = 2)

  # Check if the result is a data frame
  expect_is(pruned_densified, "data.frame")

  # Check if the resulting data frame has the correct dimensions
  expect_equal(dim(pruned_densified), c(2, 1))

  # Check if the row names are as expected
  expect_equal(rownames(pruned_densified), c("A", "C"))

  # Check if the column names are as expected
  expect_equal(colnames(pruned_densified), "Var2")
})

