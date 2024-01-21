
test_that("It correctly flattens a taxonomy and controls arguments", {
  # Test case with valid input data
  id        <- c("A", "B", "C", "D", "E", "F")
  parent_id <- c( NA, "A", "B", "B", "A", "B")
  newick <- ape::read.tree(text = "((C, D, F)B, E)A;")

  expected_output <- structure(list(
    id = c("A", "B", "C", "D", "E", "F"), 
    level1 = c("A", "A", "A", "A", "A", "A"), 
    level2 = c(NA, "B", "B", "B", "E", "B"), 
    level3 = c(NA, NA, "C", "D", NA, "F")), 
    class = "data.frame", 
    row.names = c(NA, -6L)
  )

  # Test the build_flat_taxonomy_matrix function
  expect_identical(dplyr::arrange(as_flat_taxonomy_matrix(id = id, parent_id = parent_id), id), expected_output)
  expect_identical(dplyr::arrange(as_flat_taxonomy_matrix(newick), id), expected_output)

  # Test cases with empty id and parent_id vectors
  empty_taxonomy <- as_flat_taxonomy_matrix(id = character(), parent_id = character())

  # Check if the result is a data frame
  expect_s3_class(empty_taxonomy, "data.frame")

  # Check if the data frame is empty
  expect_true(nrow(empty_taxonomy) == 0L)

  # Test cases with id and parent_id of different sizes
  id_mismatched <- c("A", "B", "C", "D")
  parent_id_mismatched <- c(NA, "A", "B")

  expect_error(as_flat_taxonomy_matrix(id = id_mismatched, parent_id = parent_id_mismatched))

  # Add more test cases as needed to cover out-of-bounds numbers and other scenarios

})
