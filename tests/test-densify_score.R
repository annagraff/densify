# Not sure we need the sample data in the end, but maybe it makes things easier.

context("Test densify_score function with unexpected numerical inputs")

test_that("It handles unexpected numerical inputs", {
  # Create sample data for testing
  iteration_log <- data.frame(
    prop_coded_data = c(0.5, 0.6, 0.7, 0.8),
    available_data_points = c(400, 300, 200, 100),  # Decreasing values
    worst_tax_abs_coding_density = c(0.1, 0.2, 0.3, 0.4),
    worst_var_abs_coding_density = c(0.2, 0.3, 0.4, 0.5),
    taxonomic_index = c(40, 30, 20, 10)  # Decreasing values
  )

  # Test the densify_score function with non-numeric exponents
  expect_error(densify_score(iteration_log, exponent_prop_coded_data = "invalid"))

  # Test the densify_score function with non-numeric exponents_lowest_taxon_coding_score
  expect_error(densify_score(iteration_log, exponent_lowest_taxon_coding_score = "invalid"))

  # Test the densify_score function with non-numeric exponents_lowest_variable_coding_score
  expect_error(densify_score(iteration_log, exponent_lowest_variable_coding_score = "invalid"))

  # Test the densify_score function with non-numeric exponents_taxonomic_diversity
  expect_error(densify_score(iteration_log, exponent_taxonomic_diversity = "invalid"))

  # Add more test cases as needed for other unexpected numerical inputs

})

