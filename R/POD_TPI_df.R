#' My Exported Function
#'
#' This function does something useful.
#'
#' @param x A parameter description.
#' @return The result of the function.
#' @export
#'
POD_TPI_df <- function(df, W, MaxDose,
                       p_T, epsilon = c(0.05, 0.05),
                       niter = 1000,
                       type_dose_decision = "mTPI",
                       type_p_prior = "uniform",
                       type_t_model = "pwuniform",
                       type_suspension = 2,
                       type_safety = 2,
                       q1 = 1, q2 = 0.15){

  source("helper_functions.R")
  source("TITE-setup.R")

  dose_vec <- df$dose
  event_time_vec <- df$event_time
  event_vec <- df$event
  n_d <- sum(event_vec == 1, na.rm = TRUE)
  m_d <- sum(event_vec == 0, na.rm = TRUE)
  r_d <- sum(is.na(event_vec))

  ## For the nth element in the event_vec, if it's equal to one or NA, the nth element in the event_time_vec should be less than W and if it's equal to 0, the nth element in the event_time_vec should equal to 28.
  validate_lists <- function(list1, list2, W) {
    if (length(list1) != length(list2)) {
      stop("event_vec and event_time_vec must have equal length!")
    }

    is_valid <- TRUE
    for (i in seq_along(list1)) {
      if (is.na(list1[[i]]) || list1[[i]] == 1) {
        if (list2[[i]] >= W) {
          is_valid <- FALSE
          break
        }
      } else if (list1[[i]] == 0) {
        if (list2[[i]] != W) {
          is_valid <- FALSE
          break
        }
      }
    }

    return(is_valid)
  }

  stopifnot(validate_lists(event_vec, event_time_vec, W))

  if(length(unique(dose_vec[is.na(event_vec)])) > 1) stop("Dose level must be identical within pending patients!")

  return(POD_TPI_decision(W, MaxDose,
                          dose_vec, event_time_vec, event_vec,
                          n_d, m_d, r_d, p_T, epsilon = c(0.05, 0.05),
                          niter = 1000,
                          type_dose_decision = "mTPI",
                          type_p_prior = "uniform",
                          type_t_model = "pwuniform",
                          type_suspension = 2,
                          type_safety = 2,
                          q1 = 1, q2 = 0.15))
}



