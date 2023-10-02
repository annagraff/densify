# Disclaimer: This is a first version of a test case generated with ChatGPT to help Marc learn the syntax. ChatGPT has limitations in accurately understanding" complex formats and dimensions, especially in real-world scenarios.


context("Test build_flat_taxonomy_matrix function")

test_that("It correctly flattens a taxonomy", {
  # Create sample data for testing
  id <- c("A", "B", "C", "D", "E", "F")
  parent_id <- c(NA, "A", "B", "B", "A", "F")

  # Test the build_flat_taxonomy_matrix function
  flat_taxonomy <- build_flat_taxonomy_matrix(id, parent_id)

  # Check if the result is a data frame
  expect_is(flat_taxonomy, "data.frame")

  # Check if the number of columns is correct (number of levels in the taxonomy)
  expect_equal(ncol(flat_taxonomy), 3)  # Adjust the expected value as needed

  # Check if the levels are correctly filled
  expect_equal(flat_taxonomy$level1, c("A", "B", "B", "A", "F"))
  expect_equal(flat_taxonomy$level2, c("B", "A", "F", NA, NA))
  expect_equal(flat_taxonomy$level3, c("C", "E", NA, NA, NA))

  # Add more checks as needed based on your specific expectations

})

