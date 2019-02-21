X      <- as.matrix(read.table('data/testData_A.txt'     ,header=FALSE,colClasses="numeric"))
y      <- as.matrix(read.table('data/testData_b.txt'     ,header=FALSE,colClasses="numeric"))
lambda <- as.matrix(read.table('data/testData_lambda.txt',header=FALSE,colClasses="numeric"))

colnames(X) <- NULL
colnames(y) <- NULL
colnames(lambda) <- NULL

info <- Adlas(X,y,lambda)
