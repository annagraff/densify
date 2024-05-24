test_that("init_pruning_state initializes correctly", {
  data <- data.frame(
    id = c("taxon1", "taxon2", "taxon3", "taxon4"),
    var1 = c("a", "a", "a", "a"),
    var2 = c(NA, "b", "a", "b"),
    var3 = c("b", NA, "b", "a")
  )

  taxonomy <- cbind(c(1, 1, 2, 2))
  rownames(taxonomy) <- c("taxon1", "taxon2", "taxon3", "taxon4")

  vars <- c(2L, 3L, 4L)
  names(vars) <- c("var1", "var2", "var3")

  ids <- c("taxon1", "taxon2", "taxon3", "taxon4")

  state <- init_pruning_state(
    data,
    vars,
    ids,
    taxonomy,
    scoring_fn = rowSums,
    density_mean_weights = list(coding = 1, taxonomy = 0)
  )


  expect_type(state, "list")

  # basic data shape
  expect_true(is.matrix(state$matrix))
  expect_true(is.matrix(state$taxonomy))
  expect_identical(nrow(state$matrix), nrow(data))
  expect_identical(nrow(state$matrix), nrow(state$taxonomy))
  expect_s3_class(state$data_levels, "data.frame")
  expect_identical(nrow(state$matrix), nrow(state$data_levels))
  expect_identical(ncol(state$matrix), ncol(state$data_levels))

  # indicator matrix
  expect_equal(state$matrix, matrix(c(1L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L), ncol = 3, byrow = TRUE, dimnames = list(ids, names(vars))))

  # indices
  expect_type(state$indices, "list")
  expect_type(state$indices$rows, "integer")
  expect_identical(length(state$indices$rows), nrow(state$matrix))
  expect_type(state$indices$cols, "integer")
  expect_identical(length(state$indices$cols), ncol(data))
  expect_type(state$indices$vars, "integer")
  expect_identical(length(state$indices$vars), ncol(state$matrix))

  expect_type(state$params, "list")
  expect_type(state$weights, "list")

  # row weights and scores
  expect_type(state$weights$rows, "list")
  expect_true(is.matrix(state$weights$rows$weights))
  expect_identical(nrow(state$weights$rows$weights), nrow(state$matrix))
  expect_true(is.numeric(state$weights$rows$scores))
  expect_identical(length(state$weights$rows$scores), nrow(state$matrix))

  # column weights and scores
  expect_type(state$weights$cols, "list")
  expect_true(is.matrix(state$weights$cols$weights))
  expect_identical(nrow(state$weights$cols$weights), ncol(state$matrix))
  expect_true(is.numeric(state$weights$cols$scores))
  expect_identical(length(state$weights$cols$scores), ncol(state$matrix))
})


test_that("prune_indices prunes rows and columns correctly", {
 data <- data.frame(
    id = c("taxon1", "taxon2", "taxon3", "taxon4"),
    var1 = c("a", "a", "a", "a"),
    var2 = c(NA, "b", "a", "b"),
    var3 = c("b", NA, "b", "a")
  )

  taxonomy <- cbind(c(1, 1, 2, 2))
  rownames(taxonomy) <- c("taxon1", "taxon2", "taxon3", "taxon4")

  vars <- c(2L, 3L, 4L)

  ids <- c("taxon1", "taxon2", "taxon3", "taxon4")

  state <- init_pruning_state(
    data,
    vars,
    ids,
    taxonomy,
    scoring_fn = rowSums,
    density_mean_weights = list(coding = 1, taxonomy = 0)
  )

  # Prune rows and columns
  indices <- data.frame(axis = c(1, 2, 1), index = c(2, 1, 3))
  state <- prune_indices(state, indices)

  expect_equal(state$matrix, matrix(c(0L, 1L,  1L, 1L), ncol = 2, byrow = TRUE, dimnames = list(c("taxon1", "taxon4"), c("var2", "var3"))))
  expect_equal(state$taxonomy, taxonomy[c(1L, 4L), , drop = FALSE])
  expect_equal(state$data_levels, `rownames<-`(encode_unique_levels(data[vars])[c(1, 4), c(2, 3)], NULL))
  expect_equal(state$indices$rows, c(1, 4))
  expect_equal(state$indices$cols, c(1L, 3L, 4L))
  expect_equal(state$indices$vars, c(2L, 3L))
})



test_that("prune_non_informative_data prunes correctly", {
  data <- data.frame(
    id = c("taxon1", "taxon2", "taxon3", "taxon4"),
    var1 = c(NA, "a", "a", "a"),
    var2 = c(NA, NA, NA, NA),
    var3 = c(NA, NA, "b", "a")
  )

  taxonomy <- data.frame(
    taxon = c("taxon1", "taxon2", "taxon3", "taxon4"),
    level = c(1, 1, 2, 2)
  )

  vars <- c(2L, 3L, 4L)
  ids <- c("taxon1", "taxon2", "taxon3", "taxon4")

  state <- init_pruning_state(
    data,
    vars,
    ids,
    taxonomy,
    scoring_fn = rowSums,
    density_mean_weights = list(coding = 1, taxonomy = 0)
  )

  # Prune non-informative data
  changes <- NULL
  new_state <- prune_non_informative_data(state, check_taxa = TRUE, check_vars = TRUE, min_variability = NA, .changes = setter(changes))

  expect_equal(rownames(new_state$matrix), c("taxon2", "taxon3", "taxon4"))
  expect_equal(colnames(new_state$matrix), c("var1", "var3"))
  expect_equal(new_state$indices$rows, c(2L, 3L, 4L))
  expect_equal(new_state$indices$cols, c(1L, 2L,  4L))
  expect_equal(new_state$indices$vars, c(2L, 3L))
})



