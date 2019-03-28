context("solve_slope")

X <- as.matrix(
  read.table(system.file("testdata", "testData_X.txt", package = "glmSLOPE"), 
  header = FALSE,
  colClasses = "numeric")
)
y <- as.matrix(
  read.table(system.file("testdata", "testData_y.txt", package = "glmSLOPE"), 
  header = FALSE,
  colClasses = "numeric")
)
lambda <- as.matrix(
  read.table(system.file("testdata", "testData_lambda.txt", package = "glmSLOPE"), 
  header = FALSE, 
  colClasses = "numeric")
)

X <- c(0.53766714, 1.833885, -2.2588469, 0.86217332, 0.31876524, -1.3076883, -0.43359202, 0.34262447, 3.5783969, 2.769437, -1.3498869, 3.0349235, 0.72540422, -0.063054873, 0.7147429, -0.20496606, -0.12414435, 1.4896976, 1.4090345, 1.4171924, 0.67149713, -1.2074869, 0.71723865, 1.6302353, 0.48889377, 1.034693, 0.72688513, -0.30344092, 0.29387147, -0.7872828, 0.88839563, -1.1470701, -1.0688705, -0.80949869, -2.9442842, 1.4383803, 0.32519054, -0.75492832, 1.3702985, -1.7115164, -0.10224245, -0.24144704, 0.31920674, 0.3128586, -0.86487992, -0.030051296, -0.16487902, 0.62770729, 1.0932657, 1.109273)
dim(X) <- c(5, 10)

y <- c(1.0734014, -5.3021346, 1.096639, -0.39124089, -0.92884291)
dim(y) <- c(5,1)

lambda <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)

info_linear <- solve_slope(X, y, lambda, model = 'linear')

test_that("slope linear solver - results consistency", {
  expect_equal(
    as.vector(info_linear$w), 
    c(-0.78334637, 0.02170683, -1.13560148, 0, 0.02170683, 0, 0, 0, 0, 0)
  )
  expect_equal(info_linear$status, 1)
})

X
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



