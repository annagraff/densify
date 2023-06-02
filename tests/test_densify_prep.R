# Create a test data frame with factor columns
testData <- data.frame(
  A = as.factor(c(NA, "value1", "value2")),
  B = as.factor(c("value3", "value4", NA))
)

# Call the function
result <- densify_prep(testData)

# Validate the output
expectedOutput <- matrix(c(0, 1, 1, 1, 1, 0), nrow = 3)
identical(result, expectedOutput)
