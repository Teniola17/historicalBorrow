#' Robustify a MAP prior by mixing with a vague component
#'
#' Constructs a robust MAP (RMAP) prior by forming a weighted mixture of the
#' MAP prior and a user-specified vague (non-informative) prior. The vague
#' component protects against prior-data conflict: if the current trial data
#' are inconsistent with the historical data, the posterior will naturally
#' down-weight the MAP component.
#'
#' @param map_obj A `map_prior` object from [fit_map()].
#' @param weight Scalar in `(0, 1)`. Weight assigned to the MAP component.
#'   The vague component receives weight `1 - weight`. Default `0.8`.
#' @param vague_prior Named list specifying the vague component:
#'   \describe{
#'     \item{`dist`}{`"normal"` (only option currently).}
#'     \item{`mean`}{Mean of the vague normal. Default `0`.}
#'     \item{`sd`}{SD of the vague normal. Default `10` (binary/count) or
#'       `1000` (continuous).}
#'   }
#'
#' @return An object of classes `c("rmap_prior", "map_prior")`. All slots of
#'   the original `map_prior` are preserved, and the following are added:
#'   \describe{
#'     \item{`mixture_rmap`}{Extended `list(weights, means, sds, K)` with the
#'       vague component appended as the final element.}
#'     \item{`rmap_weight`}{The MAP weight supplied.}
#'     \item{`vague_prior`}{The vague prior specification used.}
#'     \item{`ess_prior_rmap`}{Approximate ESS of the RMAP prior
#'       (`weight * ess_prior`).}
#'   }
#'
#' @seealso [fit_map()], [fit_borrowing_model()]
#'
#' @examples
#' \dontrun{
#' pd  <- prepare_data(hist_data, curr_data, "binary")
#' mp  <- fit_map(pd)
#' rmp <- robustify_map(mp, weight = 0.8)
#' print(rmp)
#' plot(rmp)
#' }
#'
#' @export
robustify_map <- function(map_obj,
                          weight      = 0.8,
                          vague_prior = list()) {

  if (!inherits(map_obj, "map_prior"))
    stop("`map_obj` must be a `map_prior` object from fit_map().", call. = FALSE)

  checkmate::assert_number(weight, lower = 1e-6, upper = 1 - 1e-6,
    .var.name = "weight")

  vague_prior <- .fill_vague_defaults(vague_prior, map_obj$outcome_type)

  mix_map <- map_obj$mixture

  # Scale MAP component weights by the MAP weight
  scaled_weights <- mix_map$weights * weight

  # Append the single vague component
  new_weights <- c(scaled_weights, 1 - weight)
  new_means   <- c(mix_map$means,  vague_prior$mean)
  new_sds     <- c(mix_map$sds,    vague_prior$sd)
  new_K       <- mix_map$K + 1L

  mixture_rmap <- list(
    weights = new_weights / sum(new_weights),  # re-normalise for numerical safety
    means   = new_means,
    sds     = new_sds,
    K       = new_K
  )

  ess_rmap <- weight * map_obj$ess_prior

  obj <- map_obj
  obj$mixture_rmap  <- mixture_rmap
  obj$rmap_weight   <- weight
  obj$vague_prior   <- vague_prior
  obj$ess_prior_rmap <- round(ess_rmap, 2)

  class(obj) <- c("rmap_prior", "map_prior")
  obj
}

# ── Vague prior defaults ───────────────────────────────────────────────────────

.fill_vague_defaults <- function(vp, outcome_type) {
  default_sd <- switch(outcome_type,
    binary     = 10,
    continuous = 1000,
    count      = 10
  )
  list(
    dist = vp$dist %||% "normal",
    mean = vp$mean %||% 0,
    sd   = vp$sd   %||% default_sd
  )
}
