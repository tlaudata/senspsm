#' Function to perform sensitivity analysis of calibrated confounders
#' This function calculate ATT, SE , and lambda and gamma values for a set of conditional probabilities for U
#'
#' @param data a dataframe
#' @param treatment (`character`) Name of treatment variable
#' @param indvar (`character`) Vector of names of covariate variables
#' @param outvar (`character`) Name of the outcome variable
#' @param categorical (`boolean`) TRUE if categorical, FALSE if not
#' @param p00 (`numeric`) Value of the probability U=1 when T=0 and Y=0
#' @param p01 (`numeric`) Value of the probability U=1 when for T=0 and Y=1
#' @param p10 (`numeric`) Value of the probability U=1 when for T=1 and Y=0
#' @param p11 (`numeric`) Value of the probability U=1 when for T=1 and Y=1
#' @param R (`numeric`) Number of simulations
#' @param workers (`numeric`) Number of workers for parallel processing
#' @param digits (`numeric`) Number of digits
#' @param distance (`character`) Type of the distance; "logit" by default
#' @param caliperval (`numeric`) Value of the caliper
#' @importFrom magrittr "%>%"
#'
#' @return a data frame with ATT, SE, and lambda and gamma values
#' @export
#' @examples
#' # Example 1
#' 
#' # Import data
#' data(lalonde, package = "Matching")
#' library(magrittr)
#' 
#' data <- dplyr::mutate(lalonde,
#'   black = factor(black),
#'   hisp = factor(hisp),
#'   married = factor(married)
#' )
#' 
#' xvars <- c("age", "educ", "black", "hisp", "married", "re75", "re74")
#' outcome <- "re78"
#' treatment <- "treat"
#' 
#' # Example 1 (Sensitivity using the effect of 'calibrated' confounders)
#' sensconf(data = data, treatment = treatment, indvar = xvars, outvar = outcome, categorical = FALSE, R = 5)
#' 
#' # Example 2 (Sensitivity analysis to characterize the killer confounders)
#' killconf(data = data, treatment = treatment, indvar = xvars, outvar = outcome, categorical = FALSE, d = 0.3, s = 0.3, R = 5)
sensconf <- function(data, treatment, indvar, outvar, categorical, p00 = 0.5, p01 = 0.5, p10 = 0.5, p11 = 0.5, R = 1000, workers = 6, digits = 3,
                     distance = "logit", caliperval = 0.1) {
  if ((!identical(sort(unique(data[[treatment]])), as.integer(c(0, 1)))) && (!identical(sort(unique(data[[treatment]])), as.numeric(c(0, 1))))) stop("Treatment variable should be binary")
  if (categorical == TRUE && ((!identical(sort(unique(data[[outvar]])), as.integer(c(0, 1)))) && (!identical(sort(unique(data[[outvar]])), as.numeric(c(0, 1)))))) stop("Outcome variable should be binary")
  if ((p00 <= 0) | (p01 <= 0) | (p10 <= 0) | (p11 <= 0)) stop("Probabilities cannot be set to zero or negative")

  if (categorical == TRUE) {
    data <- data %>%
      dplyr::mutate(ycatnum := !!rlang::sym(outvar)) %>%
      dplyr::mutate(ycat = factor(ycatnum))
  } else {
    thresval <- mean(data[[outvar]])
    data <- data %>%
      dplyr::mutate(ycatnum := ifelse(!!rlang::sym(outvar) < thresval, 0, 1)) %>%
      dplyr::mutate(ycat = factor(ycatnum))
  }

  future::plan(future::multicore, workers = workers)

  temp <- rsample::bootstraps(data, times = R)
  temp$results <- furrr::future_map(temp$splits, ~ sensbas(
    trial = .x, indvar = indvar, outvar = outvar, treatment = treatment,
    p00 = p00, p01 = p01, p10 = p10, p11 = p11, distance = distance, caliperval = caliperval
  ),
  .options = furrr::furrr_options(seed = TRUE)
  )

  ressim <- unlist(lapply(temp$results, unlist))
  ressim <- matrix(ressim, ncol = 4, byrow = TRUE)

  ATT <- round(mean(ressim[, 1]), digits = digits)

  se <- round(sqrt(sum(ressim[, 2]^2) / R + (R + 1) / (R * (R - 1)) * sum((ressim[, 1] - ATT)^2)), digits = digits)

  gamma <- round(mean(ressim[, 3]), digits = digits)

  lambda <- round(mean(ressim[, 4]), digits = digits)

  output_res <- data.frame(ATT = ATT, se = se, gamma = gamma, lambda = lambda)

  return(output_res)
}

