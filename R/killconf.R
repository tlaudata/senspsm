#' Function to perform sensitivity analysis of 'killer' confounders
#' This function takes a data frame, and d and s values as arguments, as defined by Ichino to calculate ATT, SE and ranges for lambda and gamma
#'
#' @param data a dataframe
#' @param treatment (`character`) Name of treatment variable
#' @param indvar (`character`) Vector of names of covariate variables
#' @param outvar (`character`) Name of the outcome variable
#' @param categorical (`boolean`) TRUE if categorical, FALSE if not
#' @param d (`numeric`) Value of the outcome effect of U in the absence of treatment
#' @param s (`numeric`) Value of the effect of U on the treatment selection
#' @param pu (`numeric`) Probability of U=1
#' @param R (`numeric`) Number of simulations
#' @param workers (`numeric`) Number of workers for parallel processing
#' @param digits (`numeric`) Number of digits
#' @param distance (`character`) Type of the distance; "logit" by default
#' @param caliperval (`numeric`) Value of the caliper
#' @importFrom magrittr "%>%"
#'
#' @return a data frame with ATT, SE and ranges for lambda and gamma
#' @export
killconf <- function(data, treatment, indvar, outvar, categorical = TRUE, d, s, pu = 0.4, R = 1000, workers = 6, digits = 3, distance = "logit", caliperval = 0.1) {
  if ((!identical(sort(unique(data[[treatment]])), as.integer(c(0, 1)))) && (!identical(sort(unique(data[[treatment]])), as.numeric(c(0, 1))))) stop("Treatment variable should be binary")
  if (categorical == TRUE && ((!identical(sort(unique(data[[outvar]])), as.integer(c(0, 1)))) && (!identical(sort(unique(data[[outvar]])), as.numeric(c(0, 1)))))) stop("Outcome variable should be binary")

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

  evalue <- data %>%
    dplyr::group_by(!!rlang::sym(treatment)) %>%
    dplyr::count(ycat) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::select(prop) %>%
    dplyr::pull()

  fvalue <- data %>%
    dplyr::count(!!rlang::sym(treatment), ycat) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::select(prop) %>%
    dplyr::pull()

  p00 <- (pu - d * fvalue[2] - (fvalue[4] + fvalue[3]) * (s + d * evalue[2]) / (evalue[3] + evalue[4])) / (fvalue[1] + fvalue[2] + (evalue[1] + evalue[2]) / (evalue[3] + evalue[4]) * (fvalue[4] + fvalue[3]))
  p01 <- d + p00
  p11 <- (s + p00 * (evalue[1]) + (d + p00) * evalue[2]) / (evalue[3] + evalue[4])
  p10 <- p11
  dres <- p01 - p00
  p1dot <- p10 * evalue[3] + p11 * evalue[4]
  p0dot <- p00 * evalue[1] + p01 * evalue[2]
  sres <- p1dot - p0dot

  future::plan(future::multisession, workers = workers)

  tempdat <- rsample::bootstraps(data, times = R)

  tempdat$results <- furrr::future_map(tempdat$splits, ~ sensbas(
    trial = .x, treatment = treatment, indvar = indvar,
    outvar = outvar, p00 = p00, p01 = p01, p10 = p10, p11 = p11, distance = distance,
    caliperval = caliperval
  ),
  .options = furrr::furrr_options(seed = TRUE)
  )

  ressim <- unlist(lapply(tempdat$results, unlist))
  ressim <- matrix(ressim, ncol = 4, byrow = TRUE)

  ATT <- mean(ressim[, 1])

  se <- sqrt(sum(ressim[, 2]^2) / R + (R + 1) / (R * (R - 1)) * sum((ressim[, 1] - ATT)^2))

  min_gamma <- round(min(ressim[, 3]), digits = digits)
  max_gamma <- round(max(ressim[, 3]), digits = digits)
  gamma_res <- paste0("[", min_gamma, " ; ", max_gamma, "]")

  min_lambda <- round(min(ressim[, 4]), digits = digits)
  max_lambda <- round(max(ressim[, 4]), digits = digits)
  lambda_res <- paste0("[", min_lambda, " ; ", max_lambda, "]")

  return(data.frame(d = d, s = s, att = paste0(round(ATT, digits = digits), " (", round(se, digits = digits), ")"), gamma = gamma_res, lambda = lambda_res))
}
