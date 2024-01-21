test_that("calculate_variability returns correct frequencies", {
  # Test case 1: Numeric vector with repeated values
  var1 <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4)
  expect_equal(calculate_variability(list(var1), 2), 3)


  # Test case 2: Factor vector
  var2 <- factor(c("A", "B", "B", "C", "C", "C", "D", "D", "D", "D"))
  expect_equal(calculate_variability(list(var2), 2), 3)

  # Test case 3: Integer vector with missing values
  var3 <- c(1, 2, NA, 3, 3, 3, 4, 4, 4, 4)
  expect_equal(calculate_variability(list(var3), 2), 3)

  # Test case 4: Empty vector
  var4 <- numeric(0)
  expect_equal( calculate_variability(list(var4), 2), 0)

  # Test case 5: List of multiple vectors
  var5_1 <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4)
  var5_2 <- factor(c("A", "B", "B", "C", "C", "C", "D", "D", "D", "D"))
  expect_equal(calculate_variability(list(var5_1, var5_2), 2), c(3, 3))
})