context("solve_slope")

X <- as.matrix(
  read.table(system.file("testdata", "testData_X.txt", package = "glmSlope"), 
  header = FALSE,
  colClasses = "numeric")
)
y <- as.matrix(
  read.table(system.file("testdata", "testData_y.txt", package = "glmSlope"), 
  header = FALSE,
  colClasses = "numeric")
)
lambda <- as.matrix(
  read.table(system.file("testdata", "testData_lambda.txt", package = "glmSlope"), 
  header = FALSE, 
  colClasses = "numeric")
)

info_linear <- solve_slope(X, y, lambda, model = 'linear')

test_that("slope linear solver - results consistency", {
  expect_equal(
    as.vector(info$w), 
    c(-0.78334637, 0.02170683, -1.13560148, 0, 0.02170683, 0, 0, 0, 0, 0)
  )
})


