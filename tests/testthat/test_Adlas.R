context("Adlas")

X <- as.matrix(read.table(system.file("testdata", "testData_X.txt",
                                package = "glmSlope"), header=FALSE,colClasses="numeric"))
y <- as.matrix(read.table(system.file("testdata", "testData_y.txt",
                                      package = "glmSlope"), header=FALSE,colClasses="numeric"))
lambda <- as.matrix(read.table(system.file("testdata", "testData_lambda.txt",
                                      package = "glmSlope"), header=FALSE,colClasses="numeric"))

info <- Adlas(X,y,lambda)

test_that("adlas results consistency", {
  expect_equal(info$L, 42.68131)
  expect_equal(info$objPrimal, 1.962045, tolerance=1e-7)
})