# Define a test case
test_that("densify() generates correct statistics", {
  result <- densify(WALS[1:10, 1:100], -Glottocode, taxonomy = glottolog_languoids, taxon_id = "Glottocode")
  taxonomy <- as_flat_taxonomy_matrix(glottolog_languoids)
  data <- lapply(result$data, function(x) {
    x <- as.data.frame(x)
    g <- x$Glottocode
    x <- as.matrix(!is.na(x[names(x) != "Glottocode"]))
    rownames(x) <- g

    x
  })  
    
  # helper function that convers non-finite values to NA
  nonfinite_to_na <- function(x) ifelse(is.finite(x), x, NA_real_)
  min0 <- function(x) {
    if (length(x) == 0) return(NA_real_)
    nonfinite_to_na(min(x))
  }
  max0 <- function(x) {
    if (length(x) == 0) return(NA_real_)
    nonfinite_to_na(max(x))
  }

  # number of data points
  expect_equal(result$n_rows, sapply(data, nrow))
  expect_equal(result$n_cols, sapply(data, ncol))
  expect_equal(result$n_data_points, sapply(data, sum))
  expect_equal(result$coding_density, nonfinite_to_na(result$n_data_points/(result$n_rows*result$n_cols)))

  # row stats
  expect_equal(result$row_coding_density_min, sapply(data, function(x) min0(rowSums(x)/ncol(x))))
  expect_equal(result$row_coding_density_median, sapply(data, function(x) median(rowSums(x)/ncol(x))))
  expect_equal(result$row_coding_density_max, sapply(data, function(x) max0(rowSums(x)/ncol(x))))

  # col stats
  expect_equal(result$col_coding_density_min, sapply(data, function(x) min0(colSums(x)/nrow(x))))
  expect_equal(result$col_coding_density_median, sapply(data, function(x) median(colSums(x)/nrow(x))))
  expect_equal(result$col_coding_density_max, sapply(data, function(x) max0(colSums(x)/nrow(x))))

  # taxonomic index
  taxonomic_index <- sapply(data, function(x) {
    ii <- match(rownames(x), taxonomy$id)
    counts <- table(taxonomy$level1[ii])

    props <- counts/sum(counts)
    -sum(props * log(props))
  })
  expect_equal(result$taxonomic_index, taxonomic_index)
})
