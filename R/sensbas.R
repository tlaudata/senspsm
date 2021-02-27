#' Basic function to perform sensitivity analysis in parallel
#' This function takes bootstrapped data frame and fraction of U=1 by treatment/outcome to calculate ATT and values of lambda and gamma
#'
#' @param trial a rsample bootstrapped dataframe
#' @param treatment (`character`) Name of treatment variable
#' @param indvar (`character`) Vector of names of covariate variables
#' @param outvar (`character`) Name of the outcome variable
#' @param p00 (`numeric`) Value of the probability U=1 when T=0 and Y=0
#' @param p01 (`numeric`) Value of the probability U=1 when for T=0 and Y=1
#' @param p10 (`numeric`) Value of the probability U=1 when for T=1 and Y=0
#' @param p11 (`numeric`) Value of the probability U=1 when for T=1 and Y=1
#' @param caliperval (`numeric`) Value of the caliper
#' @param distance (`character`) Type of the distance; "logit" by default
#' @importFrom magrittr "%>%""
#'
#' @return a list containing ATT value, SE, lambda and gamma values
#' @export
sensbas <- function(trial, treatment, indvar, outvar, p00 = 0.5, p01 = 0.5, p10 = 0.5, p11 = 0.5, caliperval = 0.1, distance = "logit") {
  formraw <- paste(treatment, "~", paste(indvar, collapse = "+"))
  formmod <- as.formula(paste(formraw, "+u", collapse = "+"))
  formout <- as.formula(paste("ycat ~", paste(c(indvar, "u"), collapse = "+")))
  out <- try({
    df <- rsample::analysis(trial) %>%
      dplyr::mutate(pval = dplyr::case_when(
        !!rlang::sym(treatment) == "0" & ycat == "0" ~ p00,
        !!rlang::sym(treatment) == "0" & ycat == "1" ~ p01,
        !!rlang::sym(treatment) == "1" & ycat == "0" ~ p10,
        !!rlang::sym(treatment) == "1" & ycat == "1" ~ p11
      ))


    b <- df$pval %>%
      furrr::future_map(~ rbinom(1, 1, .), .options = furrr::furrr_options(seed = TRUE)) %>%
      unlist()

    df$u <- b

    df <- df %>%
      dplyr::mutate(u = factor(u))

    # PS matching
    mod_match <- MatchIt::matchit(formmod,
      method = "nearest", data = df, distance = distance, caliper = caliperval
    )

    df <- df %>%
      dplyr::mutate(id = as.integer(dplyr::row_number()))

    temp_match <- mod_match$match.matrix %>%
      data.frame(stringsAsFactors = FALSE)
    temp_match <- temp_match %>%
      tibble::rownames_to_column() %>%
      na.omit()
    temp <- temp_match %>%
      setNames(c("1", "0")) %>%
      dplyr::mutate(id_pairs = dplyr::row_number()) %>%
      tidyr::gather(key = treatment, value = id, -id_pairs) %>%
      dplyr::mutate(across(everything(), ~ as.integer(.)))

    df2 <- df %>%
      dplyr::inner_join(temp, by = c("id")) %>%
      dplyr::mutate(id = factor(id))

    mod_dat <- df2 %>%
      dplyr::select(!!rlang::sym("ycatnum"), !!rlang::sym(treatment), !!rlang::sym("id_pairs")) %>%
      na.omit() %>%
      dplyr::mutate(id_pairs = as.numeric(id_pairs)) %>%
      dplyr::mutate(treatment := factor(!!rlang::sym(treatment))) %>%
      dplyr::mutate(id_pairs = factor(id_pairs))

    att <- mean(mod_dat$ycatnum[mod_dat[[treatment]] == 1]) - mean(mod_dat$ycatnum[mod_dat[[treatment]] == 0])

    nt <- length(mod_dat$ycatnum[mod_dat[[treatment]] == 1])
    seval <- sqrt(var(mod_dat$ycatnum[mod_dat[[treatment]] == 1]) / nt +
      var(mod_dat$ycatnum[mod_dat[[treatment]] == 0]) / nt)

    modgamma <- stats::glm(formout,
      data = df %>% dplyr::filter(!!rlang::sym(treatment) == 0),
      family = binomial(link = distance)
    )

    gamma <- as.numeric(exp(tail(coef(modgamma), n = 1)))

    modlambda <- stats::glm(formmod,
      data = df,
      family = binomial(link = distance)
    )

    lambda <- as.numeric(exp(tail(coef(modlambda), n = 1)))

    list(att, seval, gamma, lambda)
  })

  return(out)
}

