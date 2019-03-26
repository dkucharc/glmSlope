#' Sorted L1 parameters estimation solver for the generalized linear models
#'
#' Solves the sorted L1 penalized probelns for the following regression models:
#' \itemize{
#'   \item{linear}{
#'   \deqn{ \arg\!\min_{w} \frac{1}{2}\Vert Xw - y \Vert_2^2 +
#'       \sum_{i=1}^p \lambda_i |w|_{(i)},}
#'   }
#'   \item{logistic}{
#'   
#'  \deqn{
#'  \arg\!\min_{w}\sum_{i=1}^{n}\left(\log{\left(1+\exp{\left(-y_ix_i^Tw\right)}\right)}\right) +
#'       \sum_{i=1}^p \lambda_i |w|_{(i)},}
#'       }
#' }
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
  STATUS_MSG <- c("Optimal", "Iteration limit reached")

  # Initialize parameters and iterates
  wInit = if (is.null(initial)) rep(0, n) else initial
  t <- 1
  eta <- 2
  lambda <- matrix(lambda, nrow = length(lambda))
  y <- matrix(y, nrow = length(y))
  w <- wInit
  v <- w
  Xw <- X %*% w
  fPrev <- Inf
  iter <- 0
  status <- STATUS_RUNNING

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
    if ((iter %% grad_iter) == 0) # Includes first iterations
    {
      r <- (X %*% v) - y
      g <- t(X) %*% r
      f <- as.double(crossprod(r)) / 2
    }
    else {
      r <- (Xw + ((tPrev - 1) / t) * (Xw - XwPrev)) - y
      g <- t(X) %*% r
      f <- as.double(crossprod(r)) / 2
    }

    # Increment iteration count
    iter <- iter + 1

    # Check optimality conditions
    if ((iter %% opt_iter) == 0) { # Compute 'dual', check infeasibility and gap
      if (lasso_mode) {
        infeas <- max(norm(g, "I") - lambda, 0)
        objPrimal <- f + lambda * norm(v, "1")
        objDual <- -f - as.double(crossprod(r, y))
      }
      else {
        gs <- sort(abs(g), decreasing = TRUE)
        ys <- sort(abs(v), decreasing = TRUE)
        infeas <- max(max(cumsum(gs - lambda)), 0)

        # Compute primal and dual objective
        objPrimal <- f + as.double(crossprod(lambda, ys))
        objDual <- -f - as.double(crossprod(r, y))
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
    XwPrev <- Xw
    wPrev <- w
    fPrev <- f
    tPrev <- t

    # Lipschitz search
    while (TRUE) { # Compute prox mapping
      w <- prox_func(v - (1 / L) * g, lambda / L)
      d <- w - v

      Xw <- X %*% w
      r <- Xw - y
      f <- as.double(crossprod(r)) / 2
      q <- fPrev + as.double(crossprod(d, g)) + (L / 2) * as.double(crossprod(d))

      if (q >= f * (1 - 1e-12)) {
        break
      } else {
        L <- L * eta
      }
    } # Lipschitz search

    # Update
    t <- (1 + sqrt(1 + 4 * t^2)) / 2
    v <- w + ((tPrev - 1) / t) * (w - wPrev)
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
