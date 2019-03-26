#' Sorted L1 solver
#'
#' Solves the sorted L1 penalized regression problem: given a matrix \eqn{X},
#' a vector \eqn{y}, and a decreasing vector \eqn{\lambda}, find the vector
#' \eqn{w} minimizing
#' \deqn{\frac{1}{2}\Vert Xw - y \Vert_2^2 +
#'       \sum_{i=1}^p \lambda_i |w|_{(i)}.}
#'
#' @param X an \eqn{n}-by-\eqn{p} matrix
#' @param y a vector of length \eqn{n}
#' @param lambda vector of length \eqn{p}, sorted in decreasing order
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

solve_slope <- function(X, y, lambda, initial = NULL, max_iter = 10000, grad_iter = 20, opt_iter = 1,
                        tol_infeas = 1e-6, tol_rel_gap = 1e-6) {
  # Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes

  # This file is part of SLOPE Toolbox version 1.0.
  #
  #    The SLOPE Toolbox is free software: you can redistribute it
  #    and/or  modify it under the terms of the GNU General Public License
  #    as published by the Free Software Foundation, either version 3 of
  #    the License, or (at your option) any later version.
  #
  #    The SLOPE Toolbox is distributed in the hope that it will
  #    be useful, but WITHOUT ANY WARRANTY; without even the implied
  #    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  #    See the GNU General Public License for more details.
  #
  #    You should have received a copy of the GNU General Public License
  #    along with the SLOPE Toolbox. If not, see
  #    <http://www.gnu.org/licenses/>.

  # -------------------------------------------------------------
  # Start times
  # -------------------------------------------------------------
  t0 <- proc.time()[3]

  # -------------------------------------------------------------
  # Define function for retrieving option fields with defaults
  # -------------------------------------------------------------
  getDefaultField <- function(options, name, default) {
    if (!is.null(options[[name]])) {
      return(options[[name]])
    }
    else {
      return(default)
    }
  }

  # -------------------------------------------------------------
  # Parse parameters
  # -------------------------------------------------------------
  options <- list()
  iterations <- getDefaultField(options, "iterations", 10000)
  verbosity <- getDefaultField(options, "verbosity", 1)
  optimIter <- getDefaultField(options, "optimIter", 1)
  gradIter <- getDefaultField(options, "gradIter", 20)
  tolInfeas <- getDefaultField(options, "tolInfeas", 1e-6)
  tolRelGap <- getDefaultField(options, "tolRelGap", 1e-6)
  wInit <- getDefaultField(options, "wInit", vector())

  # -------------------------------------------------------------
  # Ensure that lambda is non-increasing
  # -------------------------------------------------------------
  n <- length(lambda)
  matrix(lambda, c(1, n))
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
  if (length(wInit) == 0) wInit <- matrix(0, n, 1)
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
  Aprods <- 2
  ATprods <- 1

  # Deal with Lasso case
  modeLasso <- (length(lambda) == 1)
  if (modeLasso) {
    proxFunction <- function(x, lambda) {
      return(sign(x) * pmax(abs(x) - lambda, 0))
    }
  } else {
    proxFunction <- function(v1, v2) {
      return(prox_sorted_L1(v1, v2))
    }
  }

  if (verbosity > 0) {
    printf <- function(...) invisible(cat(sprintf(...)))
    printf("%5s  %9s   %9s  %9s  %9s\n", "Iter", "||r||_2", "Gap", "Infeas.", "Rel. gap")
  }

  # -------------------------------------------------------------
  # Main loop
  # -------------------------------------------------------------
  while (TRUE) {
    # Compute the gradient at f(v)
    if ((iter %% gradIter) == 0) # Includes first iterations
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
    if ((iter %% optimIter) == 0) { # Compute 'dual', check infeasibility and gap
      if (modeLasso) {
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

      # Format string
      if (verbosity > 0) {
        str <- sprintf("   %9.2e  %9.2e  %9.2e", objPrimal - objDual, infeas / lambda[[1]], abs(objPrimal - objDual) / max(1, objPrimal))
      }

      # Check primal-dual gap
      if ((abs(objPrimal - objDual) / max(1, objPrimal) < tolRelGap) &&
        (infeas < tolInfeas * lambda[[1]])) {
        status <- STATUS_OPTIMAL
      }
    }
    else {
      str <- ""
    }

    if (verbosity > 0) {
      if ((verbosity == 2) ||
        ((verbosity == 1) && ((iter %% optimIter) == 0))) {
        printf("%5d  %9.2e%s\n", iter, f, str)
      }
    }

    # Stopping criteria
    if ((status == 0) && (iter >= iterations)) {
      status <- STATUS_ITERATIONS
    }

    if (status != 0) {
      if (verbosity > 0) {
        printf("Exiting with status %d -- %s\n", status, STATUS_MSG[[status]])
      }
      break
    }

    # Keep copies of previous values
    XwPrev <- Xw
    wPrev <- w
    fPrev <- f
    tPrev <- t

    # Lipschitz search
    while (TRUE) { # Compute prox mapping
      w <- proxFunction(v - (1 / L) * g, lambda / L)
      d <- w - v

      Xw <- X %*% w
      r <- Xw - y
      f <- as.double(crossprod(r)) / 2
      q <- fPrev + as.double(crossprod(d, g)) + (L / 2) * as.double(crossprod(d))

      Aprods <- Aprods + 1

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
  info$Aprods <- Aprods + ceiling(iter / gradIter)
  info$ATprods <- ATprods + iter
  info$objPrimal <- objPrimal
  info$objDual <- objDual
  info$infeas <- infeas
  info$status <- status

  info$w <- v
  info$L <- L



  return(info)
} # Function Adlas
