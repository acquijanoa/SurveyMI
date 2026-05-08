#' Generate multinomial bootstrap replicate weights by stratum
#'
#' Builds bootstrap replicate weights for a complex survey design using
#' stratum-specific multinomial PSU resampling with size `n_h - 1` and equal
#' selection probabilities, then rescales by `n_h / (n_h - 1)`.
#'
#' @param design A survey design object containing strata, PSU identifiers, and
#' base sampling weights.
#' @param B Integer number of bootstrap replicates to generate.
#' @param normalize If `TRUE` (default), each replicate column is divided by its
#'   unweighted sample mean so replicate weights have mean one, following the
#'   Kim-style bootstrap (divide each replicate by its unweighted mean).
#'
#' @return A data frame with one row per original observation, including stratum
#'   and PSU identifiers plus replicate-weight columns `Wt_1` to `Wt_B`.
#'
#' @importFrom dplyr `%>%` all_of arrange bind_rows filter group_by left_join n_distinct
#'   select starts_with summarise
#' @export
#' @examples
#' design <- svydesign(
#'   id = ~cluster_id,
#'   strata = ~strata_id,
#'   data = mydata,
#'   weights = ~weight_var,
#'   nest = TRUE
#' )
#'
#' w_boot <- multinom_weights(design, B = 200)
#'

multinom_weights <- function(design, B = 50, normalize = TRUE) {
   ## create a dataframe that contains strata, psu and sampling weights
   df <- data.frame(
     `strata` = unname(design$strata),
     `psu` = unname(design$cluster),
     `weights` = unname(design$allprob),
     `row_id` = seq_len(nrow(design$variables))
   )

   # Count the number of PSUs within strata
   psu_count <- df %>%
     group_by(.data$strata) %>%
     summarise(n_psus = n_distinct(.data$psu), .groups = "drop")

   # Accumulate weighted rows across all strata
   strata_out <- vector("list", length = nrow(psu_count))

   for (s in seq_len(nrow(psu_count))) {
     i <- psu_count$strata[s]
     nh <- psu_count$n_psus[s]

     if (nh <= 1) {
       stop(paste0("Stratum ", i, " has <= 1 PSU; multinomial bootstrap requires nh > 1."))
     }

     df_i <- df %>% filter(.data$strata == i)
     psu_ids <- unique(df_i[["psu"]])

     # Draw B random vectors from a multinomial distribution and rescale
     kh <- nh / (nh - 1)
     nh_star <- as.data.frame(
       kh * rmultinom(B, size = nh - 1, prob = rep(1 / nh, nh))
     )
     colnames(nh_star) <- paste0("r_star", seq_len(B))
     nh_star$psu <- psu_ids

     df_i <- left_join(df_i, nh_star, by = "psu")

     for (k in seq_len(B)) {
       wt_name <- paste0("Wt_", k)
       df_i[[wt_name]] <- df_i[["weights"]] * df_i[[paste0("r_star", k)]]
     }

     strata_out[[s]] <- df_i %>% select(-starts_with("r_star"))
   }

   out <- bind_rows(strata_out) %>%
     arrange(.data$row_id) %>%
     select(-all_of(c("weights", "row_id")))

   if (normalize) {
     for (k in seq_len(B)) {
       wn <- paste0("Wt_", k)
       mw <- mean(out[[wn]], na.rm = TRUE)
       if (!is.finite(mw) || mw <= 0) {
         stop("Cannot normalize `", wn, "`: non-positive or non-finite mean.")
       }
       out[[wn]] <- out[[wn]] / mw
     }
   }

   out
}

