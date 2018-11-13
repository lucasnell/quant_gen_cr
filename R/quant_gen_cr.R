#' Multiple trials for a given set of parameter values.
#'
#' @param n_trials Number of trials to perform.
#' @param N0 Vector containing the starting abundance(s) for each resource.
#' @param P0 Vector containing the starting abundance(s) for each consumer.
#' @param V0 Matrix containing the starting trait value(s) for each resource.
#'     Each row refers to a particular resource.
#' @param U0 Matrix containing the starting trait value(s) for each consumer.
#'     Each row refers to a particular consumer.
#' @param r Base growth rate for resources.
#' @param a Baseline competition coefficient for resources.
#' @param f Cost of defensive traits for resources.
#' @param b Coefficient for attack rate of consumers on resources.
#' @param c Base growth rate for consumers.
#' @param m Mortality rate for consumers.
#' @param g Cost of offensive traits for consumers.
#' @param etaN Non-additive effect of increasing 2+ traits for resources.
#' @param etaP Non-additive effect of increasing 2+ traits for consumers.
#' @param sig2N Additive genetic variance for resources.
#' @param sig2P Additive genetic variance for consumers.
#' @param delta SD for log-normal distribution generating perturbation after the
#'     starting time steps.
#' @param delta2 SD for log-normal distribution generating perturbation at the
#'     halfway point in the primary time steps.
#' @param start_t Number of starting time steps.
#' @param max_t Number of primary time steps.
#' @param density_threshold Threshold for extinction. Defaults to `1e-5`.
#' @param precision Threshold for considering rows of traits to be the same.
#'     Defaults to `1e-2`.
#' @param n_cores Number of cores to use for processing.
#'     Ignored if not compiled using OpenMP. Defaults to `1`.
#' @param show_progress Boolean for whether to show a progress bar. Defaults to `TRUE`.
#' @param return_every Integer specifing the number of iterations (after the starting
#'     iterations) at which to output values of abundances and trait values.
#'     Values <= 0 mean that only final values are output.
#'     NOTE: A high value of `max_t` and a low value of `return_every` will result in
#'     a massive object being output.
#'     Defaults to `0`.
#'
#' @return A list containing the following components:
#' * `resource`: Info for resource(s).
#' * `consumer`: Info for consumer(s).
#'
#' Both `resource` and `consumer` have the same structure:
#'
#' * If `return_every <= 0`, a list of matrices (one for each trial),
#'   each matrix with the following columns (in order):
#'     * The final fitness
#'     * The final selection pressure
#'     * The final trait values (1 column for each trait)
#' * If `return_every > 0`, a list with the following items:
#'     * `final`: Matrix of final values for each trial (just like if
#'       `return_every <= 0`).
#'     * `past_abund`: A list of matrices (one for each trial).
#'       Each matrix contains species abundances through time,
#'       where rows are species and columns are time.
#'     * `past_traits`: A list of 3D arrays (one for each rial).
#'       Each array contains species traits through time,
#'       where rows are species, columns are traits, "slices" are time.
#'
#'
#' @examples
#' n_trials <- 50
#' n <- 20  # number of resources
#' p <- 20  # number of consumers
#' q <- 3   # number of traits
#' precision <- 0.1
#'
#' max_t <- 2e+06
#' start_t <- 1000
#'
#' # turns evolution "on" with TRUE\n
#' evoN <- TRUE
#' evoP <- TRUE
#'
#' # additive genetic variances:
#' sig2N <- 1e-10
#' sig2P <- 1e-10
#'
#' if (evoN) sig2N <- 0.05
#' if (evoP) sig2P <- 0.05
#'
#' # For resource...
#' r <- 0.05  # growth rate
#' a <- 2     # baseline competition coefficient
#' f <- 0.1   # cost of defensive traits
#' b <- 0.05  # coefficient for attack rate
#' # For consumer...
#' cc <- 1    # growth rate
#' m <- 0.01  # mortality rate
#' g <- 0.01  # cost of offense
#'
#' # Non-additive tradeoffs for increasing >1 traits
#' etaN <- 0.3
#' etaP <- 0.2
#'
#' initN <- 0.4
#' initP <- 0.1
#' initV <- 0.1
#' initU <- 1
#'
#' delta_vec <- c(0.001, 0.01, 0.1, 0.3)
#' delta <- delta_vec[1]
#' delta2 <- 0.3
#'
#' N0 <- matrix(initN, nrow = n, ncol = 1)
#' P0 <- matrix(initP, nrow = p, ncol = 1)
#' V0 <- matrix(initV, nrow = n, ncol = q)
#' U0 <- matrix(initU, nrow = p, ncol = q)
#'
#' # Added to make testing/building take less time:
#' max_t <- max_t / 100
#' n_trials <- n_trials / 10
#'
#' # Max # cores available
#' n_cores <- 1
#' # You could also do this:
#' # n_cores <- parallel::detectCores()
#'
#' quantgen <- quantgen_trials(n_trials, N0, P0, V0, U0,
#'                             r, a, f, b, cc, m, g, etaN, etaP,
#'                             sig2N, sig2P, delta, delta2, start_t, max_t,
#'                             precision = precision,
#'                             n_cores = n_cores,
#'                             show_progress = FALSE)
#'
#' quantgen
#'
#' @export
#'
quantgen_trials <- function(n_trials, N0, P0, V0, U0, r, a, f, b, c, m, g, etaN,
                            etaP, sig2N, sig2P, delta, delta2, start_t, max_t,
                            density_threshold = 1e-5,
                            precision = 1e-2,
                            n_cores = 1,
                            show_progress = TRUE,
                            return_every = 0) {

    if (return_every < 0) return_every <- 0

    call_ <- match.call()
    # So it doesn't show the whole function if using do.call:
    if (call_[1] != as.call(quote(quantgen_trials()))) {
        call_[1] <- as.call(quote(quantgen_trials()))
    }

    qg_obj <- quantgen_trials_(n_trials, N0, P0, V0, U0, r, a, f, b, c, m, g, etaN,
                               etaP, sig2N, sig2P, delta, delta2, start_t, max_t,
                               density_threshold, precision,
                               n_cores, show_progress, return_every)

    class(qg_obj) <- "quantgen"

    qg_obj$call <- call_

    return(qg_obj)

}



#' Print a `quantgen` object.
#'
#' @param x an object of class \code{quantgen}.
#' @param digits the number of digits to be printed.
#' @param ... arguments passed to and from other methods.
#'
#' @export
#' @noRd
#'
print.quantgen <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    cat("< quantgen object >\n")
    if (is.null(x$resource$past_traits)) {
        cat("  - past values not stored\n")
        cat(sprintf("  - # trials: %i\n", length(x$resource)))
    } else {
        cat(sprintf("  - past values stored every %i timesteps\n", x$call$return_every))
        cat(sprintf("  - # trials = %i\n", length(x$resource$final)))
    }

    invisible(NULL)

}


