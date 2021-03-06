context("Adlas")

X <- as.matrix(read.table(system.file("testdata", "testData_X.txt",
                                package = "glmSLOPE"), header=FALSE,colClasses="numeric"))
y <- as.matrix(read.table(system.file("testdata", "testData_y.txt",
                                      package = "glmSLOPE"), header=FALSE,colClasses="numeric"))
lambda <- as.matrix(read.table(system.file("testdata", "testData_lambda.txt",
                                      package = "glmSLOPE"), header=FALSE,colClasses="numeric"))

info <- Adlas(X,y,lambda)

test_that("adlas results consistency", {
  expect_equal(as.vector(info$w), c(-0.78334637,
                                    0.02170683,
                                    -1.13560148,
                                    0.00000000,
                                    0.02170683,
                                    0.00000000,
                                    0.00000000,
                                    0.00000000,
                                    0.00000000,
                                    0.00000000))
})


