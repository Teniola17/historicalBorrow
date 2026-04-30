#' Robustify a Time-Dependent Prior by mixing with a vague component
#'
#' Similar to [robustify_map()], creates a robust TDP (RTDP) by mixing the
#' time-dependent prior's theta_star with a vague prior. Protects against
#' prior-data conflict when the current data are inconsistent with the
#' temporal trend.
#'
#' @param tdp_obj A `tdp_prior` object from [fit_tdp()].
#' @param weight Scalar in `(0, 1)`. Weight assigned to the TDP component.
#'   Default `0.8`.
#' @param vague_prior Named list specifying the vague component:
#'   \describe{
#'     \item{`dist`}{`"normal"` (only option currently).}
#'     \item{`mean`}{Mean of the vague normal. Default `0`.}
#'     \item{`sd`}{SD of the vague normal. Default `10` (binary/count) or
#'       `1000` (continuous).}
#'   }
#'
#' @return An object of classes `c("rtdp_prior", "tdp_prior")`. All slots from
#'   `tdp_prior` are preserved, plus:
#'   \describe{
#'     \item{`mixture_rtdp`}{Extended mixture with vague component appended.}
#'     \item{`rtdp_weight`}{The TDP weight supplied.}
#'     \item{`vague_prior`}{The vague prior specification used.}
#'     \item{`ess_prior_rtdp`}{ESS of the RTDP (`weight * ess_prior`).}
#'   }
#'
#' @seealso [fit_tdp()], [fit_borrowing_model()]
#'
#' @export
robustify_tdp <- function(tdp_obj,
                          weight      = 0.8,
                          vague_prior = list()) {

  if (!inherits(tdp_obj, "tdp_prior"))
    stop("`tdp_obj` must be a `tdp_prior` object from fit_tdp().", call. = FALSE)

  checkmate::assert_number(weight, lower = 1e-6, upper = 1 - 1e-6,
    .var.name = "weight")

  vague_prior <- .fill_vague_defaults(vague_prior, tdp_obj$outcome_type)

  mix_tdp <- tdp_obj$mixture

  # Scale TDP component weights by the TDP weight
  scaled_weights <- mix_tdp$weights * weight

  # Append vague component
  new_weights <- c(scaled_weights, 1 - weight)
  new_means   <- c(mix_tdp$means,  vague_prior$mean)
  new_sds     <- c(mix_tdp$sds,    vague_prior$sd)
  new_K       <- mix_tdp$K + 1L

  mixture_rtdp <- list(
    weights = new_weights / sum(new_weights),
    means   = new_means,
    sds     = new_sds,
    K       = new_K
  )

  ess_rtdp <- weight * tdp_obj$ess_prior

  obj <- tdp_obj
  obj$mixture_rtdp   <- mixture_rtdp
  obj$rtdp_weight    <- weight
  obj$vague_prior    <- vague_prior
  obj$ess_prior_rtdp <- round(ess_rtdp, 2)

  class(obj) <- c("rtdp_prior", "tdp_prior")
  obj
}
