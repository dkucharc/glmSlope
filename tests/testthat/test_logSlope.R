context("LogisticSlope")

X <- as.matrix(read.table(system.file("testdata", "testData_X.txt",
                                package = "glmSLOPE"), header=FALSE,colClasses="numeric"))
y <- as.matrix(read.table(system.file("testdata", "testData_y.txt",
                                      package = "glmSLOPE"), header=FALSE,colClasses="numeric"))
lambda <- as.matrix(read.table(system.file("testdata", "testData_lambda.txt",
                                      package = "glmSLOPE"), header=FALSE,colClasses="numeric"))

y_med <- median(y)
ind <- as.vector(y>y_med)
y[ind] <- 1
y[!ind] <- -1

info <- LogisticSlope(X,y,lambda)

test_that("logistic slope results consistency", {
  expect_equal(as.vector(info$w), c(-0.5246236, -0.5246236, -0.2851559, 0, 0, 0, 0, 0, 0, 0), tolerance=1e-7)
})