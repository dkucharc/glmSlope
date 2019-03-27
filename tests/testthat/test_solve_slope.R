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
    as.vector(info_linear$w), 
    c(-0.78334637, 0.02170683, -1.13560148, 0, 0.02170683, 0, 0, 0, 0, 0)
  )
  expect_equal(info_linear$status, 1)
})


y_bin <- y
y_med <- median(y_bin)
ind <- as.vector(y_bin > y_med)
y_bin[ind] <- 1
y_bin[!ind] <- -1

info_logistic <- solve_slope(X, y_bin, lambda, model = 'logistic')
test_that("slope logistc solver - results consistency", {
  expect_equal(
    as.vector(info_logistic$w), 
    c(-0.5246236, -0.5246236, -0.2851559, 0, 0, 0, 0, 0, 0, 0),
    tolerance = 1e-6
  )
  expect_equal(info_logistic$status, 1)
})



