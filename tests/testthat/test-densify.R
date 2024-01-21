# Test case 1: Basic densify with minimal arguments
test_that("densify with minimal arguments", {
  result <- densify(WALS[1:10, 1:10], cols = tidyselect::everything())
  expect_s3_class(result, "densify_result")
})

# Test case 2: Densify with taxonomy
test_that("densify with taxonomy", {
  result <- densify(WALS[1:10, 1:10], cols = -Glottocode, taxonomy = glottolog_languoids, taxon_id = "Glottocode")
  expect_s3_class(result, "densify_result")
})

# Test case 3: Densify with custom scoring method
test_that("densify with custom scoring method", {
  result <- densify(WALS[1:10, 1:10], cols = tidyselect::everything(), scoring = "arithmetic")
  expect_s3_class(result, "densify_result")
})

# Test case 4: Densify with minimum variability threshold
test_that("densify with minimum variability threshold", {
  result <- densify(WALS[1:10, 1:10], cols = tidyselect::everything(), min_variability = 2)
  expect_s3_class(result, "densify_result")
})

# Test case 5: Densify without taxonomic diversity consideration
test_that("densify without taxonomic diversity consideration", {
  result <- densify(WALS[1:10, 1:10], cols = -Glottocode, taxonomy = glottolog_languoids, taxon_id = "Glottocode", consider_taxonomic_diversity = FALSE)
  expect_s3_class(result, "densify_result")
})

# Test case 6: Densify without cols arguments 
test_that("densify without cols specified (should show a warning)", {
  expect_snapshot(result <- densify(WALS[1:10, 1:10]))
  expect_s3_class(result, "densify_result")
})

# Test case 7: Densify without taxon_id 
test_that("densify without cols specified (should show a warning and guess Glottocode)", {
  expect_snapshot(result <- densify(WALS[1:10, 1:10], taxonomy = glottolog_languoids))
  expect_s3_class(result, "densify_result")
})


# Test case 8: Densify with invalid data (should result in an error)
test_that("densify with invalid data (should result in an error)", {
  expect_error(densify("invalid_data"))
})

# Test case 9: Densify with invalid taxon_id (should result in an error)
test_that("densify with invalid taxon_id (should result in an error)", {
  expect_error(densify(WALS, cols = -Glottocode, taxon_id = "InvalidColumn"))
})

# Test case 10: Densify with invalid taxonomy (should result in an error)
test_that("densify with invalid taxonomy (should result in an error)", {
  expect_error(densify(WALS, cols = -Glottocode, taxonomy = list()))
})


# Test case 11: Densify with mismatched taxonomies (should result in an error)
test_that("densify with invalid taxon_id (should result in an error)", {
  WALS$Glottocode[c(1, 4, 45)] <- c("fake1234556", "fake3414", "fake34521")

  expect_error(densify(WALS, cols = -Glottocode, taxonomy = glottolog_languoids, taxon_id = "Glottocode"))
})