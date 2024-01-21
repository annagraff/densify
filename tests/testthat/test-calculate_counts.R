test_that("calculate_counts returns correct counts", {
  # Test case 1: Factor vector
  x1 <- factor(c("A", "B", "B", "C", "C", "C", "D", "D", "D", "D"))
  expect_equal(calculate_counts(x1), c(1L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L, 4L))

  # Test case 2: Character vector
  x2 <- c("apple", "orange", "apple", "banana", "banana", "banana")
  expect_equal(calculate_counts(x2), c(2L, 1L, 2L, 3L, 3L, 3L))


  # Test case 3: Character vector with NAs
  x3 <- c("apple", "orange", NA, "apple", "banana", NA, "banana", "banana")
  expect_equal(calculate_counts(x3), c(2L, 1L, 1L, 2L, 3L, 1L, 3L, 3L))


  # Test case 4: Empty vector
  x4 <- numeric(0)
  expect_equal(calculate_counts(x4), integer(0))
})