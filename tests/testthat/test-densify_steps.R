# Not sure we need the sample data in the end, but maybe it makes things easier.

context("Test densify_steps function with incorrect input values")

test_that("It handles incorrect input values", {
  # Create sample data for testing
  original_data <- data.frame(
    Taxon1 = c(1, 2, 3, NA),
    Taxon2 = c(0, 1, 1, 0),
    Taxon3 = c(1, 0, 1, 1),
    Taxon4 = c(0, 0, 0, 0)
  )

  taxonomy <- data.frame(
    id = as.character(1:4),
    parent_id = c("FamilyA", "FamilyB", "FamilyA", "FamilyC")
  )

  # Test the densify_steps function with negative max_steps
  expect_error(densify_steps(original_data, max_steps = -1))

  # Test the densify_steps function with zero max_steps
  expect_error(densify_steps(original_data, max_steps = 0))

  # Test the densify_steps function with non-integer max_steps
  expect_error(densify_steps(original_data, max_steps = 2.5))

  # Test the densify_steps function with non-integer variability_threshold
  expect_error(densify_steps(original_data, variability_threshold = 2.5))

  # Test the densify_steps function with mean_type not in the list
  expect_error(densify_steps(original_data, mean_type = "median"))

  # Test the densify_steps function with taxonomy_weight out of range
  expect_error(densify_steps(original_data, use_taxonomy = TRUE, taxonomy_matrix = taxonomy, taxonomy_weight = 1.5))

  # Test the densify_steps function with coding_weight out of range
  expect_error(densify_steps(original_data, coding_weight = -0.5))

  # Add more test cases as needed for other incorrect input values

})

