#' Sorted L1 parameters estimation solver for the generalized linear models
#'
#' Solves the sorted L1 penalized probelns for the following regression models:
#' \itemize{
#'   \item{linear}{
#'     \deqn{ 
#'       \arg\!\min_{w} \frac{1}{2}\Vert Xw - y \Vert_2^2 + \sum_{i=1}^p \lambda_i |w|_{(i)},
#'     }
#'   }
#'   \item{logistic}{
#'     \deqn{
#'       \arg\!\min_{w}\sum_{i=1}^{n}\left(\log{\left(1+\exp{\left(-y_ix_i^Tw\right)}\right)}\right) +
#'       \sum_{i=1}^p \lambda_i |w|_{(i)},
#'     }
#'   }
#'}
#' where \eqn{X} is an \eqn{n\times p} matrix, \eqn{y\in R^n} (linear) or \eqn{y\in \{0,1\}^n} (logistic) depending on the model selection,
#' and \eqn{|w|_{(i)}} denotes the \eqn{i}-th largest entry in \eqn{|w|}.
#' @param X an \eqn{n}-by-\eqn{p} matrix
#' @param y a vector of length \eqn{n}
#' @param lambda vector of length \eqn{p}, sorted in decreasing order
#' @param model a description of the regression model. Supported models: \code{"linear"} (default) and \code{"logistic"}
#' @param initial initial guess for \eqn{w}
#' @param max_iter maximum number of iterations in the gradient descent
#' @param grad_iter number of iterations between gradient updates
#' @param opt_iter number of iterations between checks for optimality
#' @param tol_infeas tolerance for infeasibility
#' @param tol_rel_gap tolerance for relative gap between primal and dual
#'  problems
#'
#' @return The solution vector \eqn{w}
#'
#' @details This optimization problem is convex and is solved using an
#' accelerated proximal gradient descent method.
#' @examples 
#' # Linear model
#' # Test data
#' X <- c(0.53766714, 1.833885, -2.2588469, 0.86217332, 0.31876524, -1.3076883, -0.43359202, 0.34262447, 3.5783969, 2.769437, -1.3498869, 3.0349235, 0.72540422, -0.063054873, 0.7147429, -0.20496606, -0.12414435, 1.4896976, 1.4090345, 1.4171924, 0.67149713, -1.2074869, 0.71723865, 1.6302353, 0.48889377, 1.034693, 0.72688513, -0.30344092, 0.29387147, -0.7872828, 0.88839563, -1.1470701, -1.0688705, -0.80949869, -2.9442842, 1.4383803, 0.32519054, -0.75492832, 1.3702985, -1.7115164, -0.10224245, -0.24144704, 0.31920674, 0.3128586, -0.86487992, -0.030051296, -0.16487902, 0.62770729, 1.0932657, 1.109273)
#' dim(X) <- c(5, 10)
#' y <- c(1.0734014, -5.3021346, 1.096639, -0.39124089, -0.92884291)
#' dim(y) <- c(5,1)
#' y.bin <- c(1, -1, 1, -1, -1)
#' dim(y.bin) <- c(5,1)
#' lambda <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
#' # Estimate parameters for linear model
#' solve_slope(X, y, lambda, model = 'linear')
#' # Estimate parameters for logistic model
#' solve_slope(X, y.bin, lambda, model = 'logistic')
#' 
#' @references M. Bogdan et al. (2015) \emph{SLOPE--Adaptive variable selection via convex optimization}, \url{http://dx.doi.org/10.1214/15-AOAS842} 
#' 
#' @export
# Adapted from SLOPE_solver.R from the SLOPE: Sorted L1 Penalized Estimation (SLOPE) package
solve_slope <- function(X, y, lambda, model = c("linear", "logistic"), initial = NULL, max_iter = 10000, grad_iter = 20, opt_iter = 1,
                        tol_infeas = 1e-6, tol_rel_gap = 1e-6) {
  model <- match.arg(model)
  # -------------------------------------------------------------
  # Start times
  # -------------------------------------------------------------
  t0 <- proc.time()[3]

  # -------------------------------------------------------------
  # Ensure that lambda is non-increasing
  # -------------------------------------------------------------
  n <- length(lambda)
  if ((n > 1) && any(lambda[2:n] > lambda[1:n - 1])) {
    stop("Lambda must be non-increasing.")
  }

  # -------------------------------------------------------------
  # Initialize
  # -------------------------------------------------------------
  # Get problem dimension
  n <- ncol(X)

  # Get initial lower bound on the Lipschitz constant
  rndKind <- RNGkind()
  rndState <- tryCatch({
    .Random.seed
  }, error = function(m) {})
  set.seed(0, kind = "Mersenne-Twister", normal.kind = "Inversion")
  w <- matrix(rnorm(n), c(n, 1))
  w <- w / sqrt(sum(w^2))
  w <- t(X) %*% (X %*% w)
  L <- sqrt(sum(w^2))
  if (!is.null(rndState)) {
    set.seed(rndState[-1], kind = rndKind[1], normal.kind = rndKind[2])
  }

  # Set constants
  STATUS_RUNNING <- 0
  STATUS_OPTIMAL <- 1
  STATUS_ITERATIONS <- 2

  # Initialize parameters and iterates
  w.init <- if (is.null(initial)) rep(0, n) else initial
  t <- 1
  eta <- 2
  lambda <- matrix(lambda, nrow = length(lambda))
  Y <- diag(as.vector(y), length(y))
  YX <- Y %*% X
  y <- matrix(y, nrow = length(y))
  w <- w.init
  v <- w
  Xw <- X %*% w
  f.prev <- Inf
  iter <- 0
  status <- STATUS_RUNNING
  # Logistic


  # Deal with Lasso case
  lasso_mode <- (length(lambda) == 1)
  if (lasso_mode) {
    prox_func <- function(x, lambda) {
      return(sign(x) * pmax(abs(x) - lambda, 0))
    }
  } else {
    prox_func <- function(v1, v2) {
      return(prox_sorted_L1(v1, v2))
    }
  }

  # -------------------------------------------------------------
  # Main loop
  # -------------------------------------------------------------
  while (TRUE) {
    # Compute the gradient at f(v)
    if (model == "linear") {
      if ((iter %% grad_iter) == 0) {
        r <- (X %*% v) - y
      } else {
        r <- (Xw + ((t.prev - 1) / t) * (Xw - Xw.prev)) - y
      }
      g <- as.double(crossprod(X, r))
      f <- as.double(crossprod(r)) / 2
    } else if (model == "logistic") {
      r <- 1 / (1 + exp(YX %*% v))
      g <- -as.double(crossprod(YX, r))
      log1mr <- log(1 - r)
      f <- -sum(log1mr)
    } else {
      stop("Not supported model")
    }

    # Increment iteration count
    iter <- iter + 1

    # Check optimality conditions
    if ((iter %% opt_iter) == 0) {
      # Compute 'dual', check infeasibility and gap
      if (lasso_mode) {
        infeas <- max(norm(g, "I") - lambda, 0)
        objPrimal <- f + lambda * norm(v, "1")
      } else {
        gs <- sort(abs(g), decreasing = TRUE)
        vs <- sort(abs(v), decreasing = TRUE)
        infeas <- max(max(cumsum(gs - lambda)), 0)
        # Compute primal and dual objective
        objPrimal <- f + as.double(crossprod(lambda, vs))
      }
      if (model == "linear") {
        objDual <- -f - as.double(crossprod(r, y))
      } else if (model == "logistic") {
        objDual <- as.double(crossprod(r - 1, log1mr)) - as.double(crossprod(r, log(r)))
      } else {
        (
          stop("Not supported model")
        )
      }

      # Check primal-dual gap
      if ((abs(objPrimal - objDual) / max(1, objPrimal) < tol_rel_gap) &&
        (infeas < tol_infeas * lambda[[1]])) {
        status <- STATUS_OPTIMAL
      }
    }

    # Stopping criteria
    if ((status == 0) && (iter >= max_iter)) {
      status <- STATUS_ITERATIONS
    }

    if (status != 0) {
      break
    }

    # Keep copies of previous values
    Xw.prev <- Xw
    w.prev <- w
    f.prev <- f
    t.prev <- t

    # Lipschitz search
    while (TRUE) { # Compute prox mapping
      w <- prox_func(v - (1 / L) * g, lambda / L)
      d <- w - v

      Xw <- X %*% w
      if (model == "linear") {
        r <- Xw - y
        f <- as.double(crossprod(r)) / 2
      } else if (model == "logistic") {
        r <- 1 / (1 + exp(YX %*% w))
        f <- -sum(log(1 - r))
      } else {
        stop("Not supported model")
      }
      q <- f.prev + as.double(crossprod(d, g)) + (L / 2) * as.double(crossprod(d))

      if (q < f * (1 - 1e-12)) {
        L <- L * eta
      } else {
        break
      }
    } # Lipschitz search

    # Update
    t <- (1 + sqrt(1 + 4 * t^2)) / 2
    v <- w + ((t.prev - 1) / t) * (w - w.prev)
  } # While (TRUE)


  # Information structure
  info <- c()
  info$runtime <- proc.time()[3] - t0
  info$objPrimal <- objPrimal
  info$objDual <- objDual
  info$infeas <- infeas
  info$status <- status

  info$w <- v
  info$L <- L



  return(info)
} # Function Adlas
