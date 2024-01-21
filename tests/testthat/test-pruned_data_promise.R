test_that("make_pruned_data_factory returns a function", {
  data <- data.frame(a = 1:10, b = 11:20, c = 21:30)
  factory <- make_pruned_data_factory(data)

  expect_type(factory, "closure")
})

test_that("print.pruned_data_promise prints correct information", {
  data <- data.frame(a = 1:10, b = 11:20, c = 21:30)
  factory <- make_pruned_data_factory(data)
  promise <- factory(c(1, 3, 5), c("a", "c"))

  expect_snapshot(print(promise))
})

test_that("as.data.frame.pruned_data_promise returns correct data frame", {
  data <- data.frame(a = 1:10, b = 11:20, c = 21:30)
  factory <- make_pruned_data_factory(data)
  promise <- factory(c(1, 3, 5), c("a", "c"))

  df <- as.data.frame(promise)
  expect_equal(df, vctrs::vec_slice(data[c("a", "c")], c(1, 3, 5)))
})
