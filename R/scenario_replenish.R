#' Scenario-driven replenishment under joint demand uncertainty
#'
#' @description
#' Demand for several products often rises and falls together in the same
#' period, so resampling each product separately can misstate joint
#' stockout risk. `scenario_replenish()` builds *joint* demand scenarios
#' for several stock-keeping units (SKUs) and chooses a replenishment
#' plan that globally minimizes the empirical conditional value-at-risk
#' (CVaR) of stockout plus holding loss over integer plans.
#'
#' The procedure:
#'
#' 1. **Align histories by period.** Every SKU's history is restricted to
#'    the periods shared by all SKUs, so that each scenario resamples the
#'    same periods for every SKU.
#' 2. **Draw joint bootstrap scenarios.** For each scenario, one vector
#'    of period indices is drawn with replacement and applied to every
#'    SKU, preserving cross-SKU co-movement.
#' 3. **Fit and predict natively.** In each scenario, each SKU's
#'    resampled history is fitted with the same Gaussian engine used by
#'    [fit.infer()] ([stats::glm()] with `family = "gaussian"`), and
#'    demand is predicted from that SKU's own future features. Negative
#'    predictions are truncated to zero. No precomputed demand is ever
#'    substituted for these fits and predictions.
#' 4. **Optimize exactly.** The empirical CVaR of total stockout plus
#'    holding loss is minimized globally over all integer replenishment
#'    plans respecting each SKU's `capacity` and the total `budget`. The
#'    CVaR retains the fractional weight of scenarios sitting exactly at
#'    the tail boundary rather than averaging a top-`ceil()` count of
#'    scenario losses.
#'
#' If any step fails---invalid input, a failed or rank-deficient fit, a
#' non-finite prediction, or a candidate space larger than
#' `max_combinations`---the function errors and returns nothing rather
#' than delivering a partial result.
#'
#' @param data A data frame with one row per SKU and period, containing
#'   the `sku` and `period` columns, the response column named in
#'   `formula`, and every predictor column used in `formula`.
#' @param formula A two-sided formula whose left-hand side is a bare
#'   numeric column of `data` (e.g. `demand ~ price + promo`). The same
#'   formula is fitted per SKU with `family = "gaussian"`, as in
#'   [fit.infer()].
#' @param future A data frame with exactly one row per SKU, containing
#'   the `sku` column and the future value of every predictor in
#'   `formula`.
#' @param sku,period Unquoted names of the SKU and period columns in
#'   `data`; `future` must use the same `sku` column name.
#' @param reps Number of joint bootstrap scenarios.
#' @param unit_cost Budget cost per unit ordered.
#' @param holding_cost Loss per unit left over after demand is met.
#' @param stockout_cost Loss per unit of unmet demand.
#' @param capacity Maximum units to order per SKU; must be
#'   whole-numbered.
#' @param budget Total budget; the sum of `unit_cost * order` across SKUs
#'   must not exceed it.
#' @param alpha CVaR level in `(0, 1)`; the mean loss over the worst
#'   `1 - alpha` tail mass of scenarios is minimized.
#' @param max_combinations Safety limit on the number of candidate
#'   integer plans; the function errors if enumeration would exceed it.
#' @param x An object of class `scenario_replenishment`, i.e. output from
#'   `scenario_replenish()`.
#' @param ... Arguments passed to methods (currently ignored).
#'
#' @details
#' `unit_cost`, `holding_cost`, `stockout_cost`, and `capacity` may each
#' be a single number recycled to every SKU or a named numeric vector
#' with one element per SKU.
#'
#' @return An object of class `scenario_replenishment`, a list with
#'   components:
#'
#'   * `plan`: tibble with the optimal integer `order` per SKU, its
#'     `unit_cost`, `spend`, `capacity`, and `capacity_effective`
#'     (the capacity tightened by the largest scenario demand).
#'   * `budget`: tibble with `budget`, `spent`, and `remaining`.
#'   * `objective`: tibble with `alpha`, the number of `scenarios`, the
#'     empirical `cvar` and `var` (value-at-risk), the `mean_loss`, the
#'     `tail_mass` in scenario-equivalents, the counts of scenarios above
#'     and at the VaR, and the fractional `boundary_weight` entering the
#'     CVaR.
#'   * `tail`: tibble of tail scenarios with their CVaR `weight` (1 for
#'     scenarios above the VaR, fractional at the boundary).
#'   * `losses`: tibble of per-scenario `stockout_loss`, `holding_loss`,
#'     and total `loss` under the chosen plan.
#'   * `demand`: tibble of joint demand scenarios with the raw Gaussian
#'     `prediction` and the `demand` truncated at zero per scenario and
#'     SKU.
#'   * `fits`: tibble of the native per-scenario fits (`term`,
#'     `estimate`) per SKU, mirroring [fit.infer()] output.
#'   * `population`: tibble of the aligned sampling population that
#'     scenarios were resampled from.
#'   * `scenario_index`: tibble of the resampled period indices per
#'     scenario, shared across SKUs.
#'   * `meta`: a list of inputs needed to interpret the plan.
#'
#' @section Scope:
#' All reported risk quantities are *empirical* quantities of the
#' regression-bootstrap scenarios resampled from the aligned history.
#' They describe resampling risk only; they are not a guarantee of future
#' demand coverage or of population-level risk.
#'
#' @includeRmd man-roxygen/seeds.Rmd
#'
#' @examples
#' set.seed(1)
#'
#' # two SKUs whose demands move together period by period
#' history <- expand.grid(period = 1:24, sku = c("a", "b"))
#' history$demand <- ifelse(history$sku == "a", 20, 35) +
#'   2 * history$period +
#'   rep(rnorm(24, sd = 3), times = 2) +
#'   rnorm(48, sd = 1)
#'
#' future <- data.frame(sku = c("a", "b"), period = 25)
#'
#' scenario_replenish(
#'   history,
#'   demand ~ period,
#'   future,
#'   sku = sku,
#'   period = period,
#'   reps = 100,
#'   unit_cost = 2,
#'   holding_cost = 1,
#'   stockout_cost = 8,
#'   capacity = 200,
#'   budget = 800,
#'   alpha = 0.9
#' )
#'
#' @export
scenario_replenish <- function(
  data,
  formula,
  future,
  sku,
  period,
  reps = 1000,
  unit_cost,
  holding_cost,
  stockout_cost,
  capacity,
  budget,
  alpha = 0.95,
  max_combinations = 1e6
) {
  sku_col <- rlang::as_name(rlang::enquo(sku))
  period_col <- rlang::as_name(rlang::enquo(period))

  check_scenario_args(
    data = data,
    formula = formula,
    future = future,
    sku_col = sku_col,
    period_col = period_col,
    reps = reps,
    alpha = alpha,
    budget = budget,
    max_combinations = max_combinations
  )

  response <- as.character(formula[[2L]])
  rhs_vars <- all.vars(formula[[3L]])

  # align histories to the periods shared by every SKU --------------------
  aligned <- align_panel(data, sku_col, period_col, response, rhs_vars)
  sku_ids <- aligned$sku_ids
  sku_keys <- as.character(sku_ids)
  n_sku <- length(sku_ids)
  n_periods <- length(aligned$periods)

  factored <- factorize_histories(aligned$histories, aligned$cat_vars)
  future_std <- standardize_future(
    future,
    sku_col,
    sku_ids,
    rhs_vars,
    aligned$num_vars,
    aligned$cat_vars,
    factored$levels
  )

  # per-SKU cost and capacity parameters -----------------------------------
  unit_cost <- recycle_sku_param(unit_cost, n_sku, sku_keys, "unit_cost")
  holding_cost <- recycle_sku_param(
    holding_cost, n_sku, sku_keys, "holding_cost"
  )
  stockout_cost <- recycle_sku_param(
    stockout_cost, n_sku, sku_keys, "stockout_cost"
  )
  capacity <- recycle_sku_param(
    capacity, n_sku, sku_keys, "capacity",
    whole = TRUE
  )

  # joint bootstrap scenarios: one shared index vector per scenario -------
  reps <- as.integer(reps)
  if (n_periods * reps >= .Machine$integer.max) {
    cli_abort("`reps` times the number of aligned periods is too large.")
  }
  scenario_idx <- matrix(
    sample.int(n_periods, n_periods * reps, replace = TRUE),
    nrow = n_periods,
    ncol = reps
  )
  scenarios <- fit_predict_scenarios(
    factored$histories,
    future_std,
    formula,
    scenario_idx,
    sku_ids
  )

  # exact global optimization over integer plans ---------------------------
  capacity_eff <- pmin(capacity, ceiling(apply(scenarios$demand, 2L, max)))
  if (sum(capacity_eff + 1) * reps > 1e8) {
    cli_abort(
      c(
        "The loss matrices required for optimization would be too large.",
        i = "Reduce `capacity` or `reps`.",
        i = "No partial result is returned."
      )
    )
  }
  loss_mats <- build_loss_mats(
    scenarios$demand, capacity_eff, holding_cost, stockout_cost
  )
  plans <- enumerate_plans(capacity_eff, unit_cost, budget, max_combinations)
  optimum <- minimize_empirical_cvar(plans, loss_mats, alpha, reps)
  best <- plans[optimum$index, ]

  # per-scenario losses and tail risk under the optimal plan ---------------
  losses <- scenario_losses(
    scenarios$demand, best, holding_cost, stockout_cost
  )
  tail_info <- cvar_empirical(losses$loss, alpha)

  # assemble the complete result; any earlier failure returns nothing -----
  spent <- sum(best * unit_cost)
  at_or_above <- losses$loss >= tail_info$var
  tail <- losses[at_or_above, c("scenario", "loss")]
  tail$weight <- ifelse(
    losses$loss[at_or_above] > tail_info$var,
    1,
    tail_info$boundary_weight / tail_info$n_at
  )
  tail <- tail[order(-tail$loss, tail$scenario), ]

  out <- list(
    plan = tibble::tibble(
      sku = sku_ids,
      order = as.integer(best),
      unit_cost = unname(unit_cost),
      spend = unname(best * unit_cost),
      capacity = unname(capacity),
      capacity_effective = unname(capacity_eff)
    ),
    budget = tibble::tibble(
      budget = budget,
      spent = spent,
      remaining = budget - spent
    ),
    objective = tibble::tibble(
      alpha = alpha,
      scenarios = reps,
      cvar = tail_info$cvar,
      var = tail_info$var,
      mean_loss = mean(losses$loss),
      tail_mass = tail_info$tail_mass,
      n_tail_above = tail_info$n_above,
      n_tail_at = tail_info$n_at,
      boundary_weight = tail_info$boundary_weight
    ),
    tail = tibble::as_tibble(tail),
    losses = losses,
    demand = tibble::tibble(
      scenario = rep(seq_len(reps), times = n_sku),
      sku = rep(sku_ids, each = reps),
      prediction = as.vector(scenarios$prediction),
      demand = as.vector(scenarios$demand)
    ),
    fits = scenarios$fits,
    population = dplyr::bind_rows(aligned$histories),
    scenario_index = tibble::tibble(
      scenario = rep(seq_len(reps), each = n_periods),
      draw = rep(seq_len(n_periods), times = reps),
      period_index = as.vector(scenario_idx),
      period = aligned$periods[as.vector(scenario_idx)]
    ),
    meta = list(
      formula = formula,
      response = response,
      sku = sku_col,
      period = period_col,
      reps = reps,
      alpha = alpha,
      n_periods = n_periods,
      holding_cost = unname(holding_cost),
      stockout_cost = unname(stockout_cost)
    )
  )

  class(out) <- "scenario_replenishment"
  out
}

#' @rdname scenario_replenish
#' @export
print.scenario_replenishment <- function(x, ...) {
  obj <- x$objective
  cat("# Scenario-driven replenishment\n")
  cat(sprintf(
    "SKUs: %d | aligned periods: %d | scenarios: %d | alpha: %g\n",
    nrow(x$plan), x$meta$n_periods, x$meta$reps, obj$alpha
  ))
  cat(sprintf(
    "Budget: %s | spent: %s | remaining: %s\n",
    fmt_num(x$budget$budget),
    fmt_num(x$budget$spent),
    fmt_num(x$budget$remaining)
  ))
  cat(sprintf(
    "Empirical CVaR of stockout + holding loss: %s (VaR: %s, mean: %s)\n",
    fmt_num(obj$cvar), fmt_num(obj$var), fmt_num(obj$mean_loss)
  ))
  cat("\nReplenishment plan (integer units):\n")
  print(x$plan)
  cat(sprintf(
    paste0(
      "\nTail risk: %s scenario-equivalents = %d scenario(s) above the ",
      "VaR at full weight + %d at the VaR sharing fractional weight %s.\n"
    ),
    fmt_num(obj$tail_mass),
    obj$n_tail_above,
    obj$n_tail_at,
    fmt_num(obj$boundary_weight)
  ))
  cat("CVaR keeps the fractional tail-boundary weight (no top-ceil ",
      "averaging).\n", sep = "")
  cat(sprintf(
    paste0(
      "Note: empirical risk over the %d regression-bootstrap scenarios ",
      "only; not a\n  guarantee of future demand coverage or ",
      "population-level risk.\n"
    ),
    x$meta$reps
  ))
  invisible(x)
}

# Input validation ---------------------------------------------------------

check_scenario_args <- function(
  data,
  formula,
  future,
  sku_col,
  period_col,
  reps,
  alpha,
  budget,
  max_combinations,
  call = caller_env()
) {
  check_type(data, is.data.frame, call = call)
  check_type(future, is.data.frame, call = call)

  if (!inherits(formula, "formula") || length(formula) != 3L) {
    cli_abort(
      "`formula` must be a two-sided formula, e.g. `demand ~ price`.",
      call = call
    )
  }
  if (!rlang::is_symbol(formula[[2L]])) {
    cli_abort(
      "The response of `formula` must be a bare column name in `data`, \\
       not a transformation.",
      call = call
    )
  }

  response <- as.character(formula[[2L]])
  rhs_vars <- all.vars(formula[[3L]])

  if (identical(sku_col, period_col)) {
    cli_abort("`sku` and `period` must be different columns.", call = call)
  }
  if (response %in% c(sku_col, period_col)) {
    cli_abort(
      "The response of `formula` must be distinct from the `sku` and \\
       `period` columns.",
      call = call
    )
  }
  if (sku_col %in% rhs_vars) {
    cli_abort(
      "The `sku` column cannot be a predictor: each SKU is fitted \\
       separately.",
      call = call
    )
  }
  if (response %in% rhs_vars) {
    cli_abort(
      "The response column cannot also be a predictor in `formula`.",
      call = call
    )
  }

  missing_data <- setdiff(
    c(sku_col, period_col, response, rhs_vars),
    names(data)
  )
  if (length(missing_data) > 0L) {
    cli_abort(
      "{qty(missing_data)}Column{?s} {.field {missing_data}} not found \\
       in `data`.",
      call = call
    )
  }
  missing_future <- setdiff(c(sku_col, rhs_vars), names(future))
  if (length(missing_future) > 0L) {
    cli_abort(
      "{qty(missing_future)}Column{?s} {.field {missing_future}} not \\
       found in `future`.",
      call = call
    )
  }

  if (!is_single_number(reps, min_val = 1) || reps != floor(reps)) {
    cli_abort(
      "`reps` must be a single whole number of at least 1.",
      call = call
    )
  }
  if (
    !is_single_number(
      alpha,
      min_val = 0,
      max_val = 1,
      include_min_val = FALSE,
      include_max_val = FALSE
    )
  ) {
    cli_abort(
      "`alpha` must be a single number strictly between 0 and 1.",
      call = call
    )
  }
  if (!is_single_number(budget, min_val = 0)) {
    cli_abort("`budget` must be a single non-negative number.", call = call)
  }
  if (!is_single_number(max_combinations, min_val = 1)) {
    cli_abort(
      "`max_combinations` must be a single number of at least 1.",
      call = call
    )
  }

  invisible(TRUE)
}

check_complete_columns <- function(df, cols, arg, call = caller_env()) {
  for (cl in cols) {
    x <- df[[cl]]
    if (anyNA(x)) {
      cli_abort(
        "The {.field {cl}} column of `{arg}` contains missing values; \\
         scenario replenishment requires complete data.",
        call = call
      )
    }
    if (is.numeric(x) && any(!is.finite(x))) {
      cli_abort(
        "The {.field {cl}} column of `{arg}` contains non-finite values.",
        call = call
      )
    }
  }

  invisible(TRUE)
}

# Panel alignment ----------------------------------------------------------

align_panel <- function(
  data,
  sku_col,
  period_col,
  response,
  rhs_vars,
  call = caller_env()
) {
  used_cols <- c(sku_col, period_col, response, rhs_vars)

  if (anyDuplicated(data[c(sku_col, period_col)])) {
    cli_abort(
      "`data` must have at most one row per `sku`-`period` combination.",
      call = call
    )
  }

  check_complete_columns(data, used_cols, "data")

  if (!is.numeric(data[[response]])) {
    cli_abort(
      "The response column {.field {response}} must be numeric for the \\
       Gaussian fit.",
      call = call
    )
  }
  data[[response]] <- as.numeric(data[[response]])

  num_vars <- character(0)
  cat_vars <- character(0)
  for (v in rhs_vars) {
    x <- data[[v]]
    if (is.numeric(x) || inherits(x, c("Date", "POSIXt"))) {
      data[[v]] <- as.numeric(x)
      num_vars <- c(num_vars, v)
    } else if (is.character(x) || is.logical(x) || is.factor(x)) {
      data[[v]] <- as.character(x)
      cat_vars <- c(cat_vars, v)
    } else {
      cli_abort(
        "The {.field {v}} column of `data` has an unsupported type for \\
         the Gaussian fit.",
        call = call
      )
    }
  }

  sku_ids <- unique(data[[sku_col]])
  histories <- split(data[used_cols], data[[sku_col]])
  histories <- histories[as.character(sku_ids)]

  # intersect periods across SKUs, preserving the period column's type
  common <- NULL
  for (df in histories) {
    ps <- df[[period_col]]
    common <- if (is.null(common)) ps else common[common %in% ps]
  }
  common <- sort(unique(common))

  if (length(common) < 2L) {
    cli_abort(
      "The SKUs share fewer than 2 common periods; cannot align \\
       histories for the joint bootstrap.",
      call = call
    )
  }

  n_all <- length(Reduce(union, lapply(histories, `[[`, period_col)))
  n_dropped <- n_all - length(common)
  if (n_dropped > 0L) {
    cli_inform(
      "{n_dropped} period{?s} dropped to align all SKUs to shared \\
       periods."
    )
  }

  histories <- lapply(histories, function(df) {
    df <- df[df[[period_col]] %in% common, , drop = FALSE]
    df[order(df[[period_col]]), , drop = FALSE]
  })

  list(
    histories = histories,
    periods = common,
    sku_ids = sku_ids,
    num_vars = num_vars,
    cat_vars = cat_vars
  )
}

# Give every categorical predictor fixed levels from that SKU's own
# history, so resamples keep a stable design matrix and genuinely missing
# levels surface as aliased (NA) coefficients rather than silent shifts.
factorize_histories <- function(histories, cat_vars) {
  levels_map <- vector("list", length(histories))
  names(levels_map) <- names(histories)
  out <- vector("list", length(histories))
  names(out) <- names(histories)

  for (s in names(histories)) {
    df <- histories[[s]]
    levels_s <- list()
    for (v in cat_vars) {
      lev <- sort(unique(df[[v]]))
      df[[v]] <- factor(df[[v]], levels = lev)
      levels_s[[v]] <- lev
    }
    levels_map[[s]] <- levels_s
    out[[s]] <- df
  }

  list(histories = out, levels = levels_map)
}

standardize_future <- function(
  future,
  sku_col,
  sku_ids,
  rhs_vars,
  num_vars,
  cat_vars,
  levels_map,
  call = caller_env()
) {
  keys <- as.character(sku_ids)
  future_skus <- as.character(future[[sku_col]])

  if (anyDuplicated(future_skus)) {
    cli_abort("`future` must have exactly one row per SKU.", call = call)
  }
  if (!setequal(future_skus, keys) || length(future_skus) != length(keys)) {
    cli_abort(
      "`future` must have exactly one row per SKU in `data`.",
      call = call
    )
  }
  future <- future[match(keys, future_skus), , drop = FALSE]

  check_complete_columns(future, rhs_vars, "future")

  out <- vector("list", length(rhs_vars))
  names(out) <- rhs_vars
  for (v in rhs_vars) {
    x <- future[[v]]
    if (v %in% num_vars) {
      if (!(is.numeric(x) || inherits(x, c("Date", "POSIXt")))) {
        cli_abort(
          "The {.field {v}} column of `future` must be numeric, \\
           matching `data`.",
          call = call
        )
      }
      out[[v]] <- as.numeric(x)
    } else {
      if (!(is.character(x) || is.logical(x) || is.factor(x))) {
        cli_abort(
          "The {.field {v}} column of `future` must be categorical, \\
           matching `data`.",
          call = call
        )
      }
      out[[v]] <- as.character(x)
    }
  }

  # validate each future categorical value against that SKU's own levels
  for (v in cat_vars) {
    vals <- out[[v]]
    for (i in seq_along(vals)) {
      if (!vals[i] %in% levels_map[[keys[i]]][[v]]) {
        cli_abort(
          "The future value {.val {vals[i]}} of {.field {v}} for SKU \\
           {.val {keys[i]}} was not observed in that SKU's history.",
          call = call
        )
      }
    }
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

# Per-SKU scalar-or-named-vector parameters --------------------------------

recycle_sku_param <- function(
  x,
  n,
  keys,
  arg,
  whole = FALSE,
  call = caller_env()
) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    cli_abort("`{arg}` must be a finite numeric vector.", call = call)
  }
  if (length(x) == 1L) {
    out <- rep(as.numeric(x), n)
  } else if (length(x) == n) {
    if (is.null(names(x)) || !setequal(names(x), keys)) {
      cli_abort(
        "When `{arg}` has more than one element it must be a named \\
         vector whose names are the SKUs.",
        call = call
      )
    }
    out <- as.numeric(x[keys])
  } else {
    cli_abort(
      "`{arg}` must be a single number or a named vector with one \\
       element per SKU.",
      call = call
    )
  }
  if (any(out < 0)) {
    cli_abort("`{arg}` must be non-negative.", call = call)
  }
  if (whole && any(out != floor(out))) {
    cli_abort(
      "`{arg}` must be whole-numbered (integer units).",
      call = call
    )
  }

  out
}

# Scenario generation --------------------------------------------------------

# Fit each SKU's resampled history with the same Gaussian engine as
# `fit.infer()` and predict from that SKU's own future features.
fit_predict_scenarios <- function(
  histories,
  future,
  formula,
  scenario_idx,
  sku_ids,
  call = caller_env()
) {
  reps <- ncol(scenario_idx)
  keys <- names(histories)
  n_sku <- length(histories)
  prediction <- matrix(NA_real_, nrow = reps, ncol = n_sku)
  fits <- vector("list", reps * n_sku)
  i <- 0L

  for (r in seq_len(reps)) {
    idx <- scenario_idx[, r]
    for (s in seq_len(n_sku)) {
      i <- i + 1L
      resampled <- histories[[s]][idx, , drop = FALSE]
      model <- tryCatch(
        stats::glm(formula, data = resampled, family = stats::gaussian),
        error = function(e) {
          cli_abort(
            "The Gaussian fit failed in scenario {r} for SKU \\
             {.val {keys[s]}}: {conditionMessage(e)}",
            call = call
          )
        }
      )
      tidied <- tidy_glm_coefs(model)
      if (any(is.na(tidied$estimate))) {
        cli_abort(
          c(
            "Scenario {r} for SKU {.val {keys[s]}} produced a \\
             rank-deficient Gaussian fit (aliased coefficients).",
            i = "This usually means a resampled history lost a \\
             categorical level or has collinear columns.",
            i = "No partial result is returned."
          ),
          call = call
        )
      }
      pred <- tryCatch(
        stats::predict(
          model,
          newdata = future[s, , drop = FALSE],
          type = "response"
        ),
        error = function(e) {
          cli_abort(
            "Prediction failed in scenario {r} for SKU \\
             {.val {keys[s]}}: {conditionMessage(e)}",
            call = call
          )
        }
      )
      pred <- unname(pred)[1L]
      if (!is.finite(pred)) {
        cli_abort(
          "Scenario {r} for SKU {.val {keys[s]}} produced a non-finite \\
           prediction; no partial result is returned.",
          call = call
        )
      }
      prediction[r, s] <- pred
      fits[[i]] <- dplyr::mutate(
        tidied,
        scenario = r,
        sku = sku_ids[s],
        .before = 1L
      )
    }
  }

  list(
    fits = dplyr::bind_rows(fits),
    prediction = prediction,
    demand = pmax(prediction, 0)
  )
}

# Optimization ---------------------------------------------------------------

# Loss of ordering q units when demand is d: holding cost on leftovers
# plus stockout cost on unmet demand. One matrix per SKU with one row per
# candidate order quantity and one column per scenario.
build_loss_mats <- function(demand, max_units, holding, stockout) {
  mats <- vector("list", ncol(demand))
  for (s in seq_len(ncol(demand))) {
    qs <- seq_len(max_units[s] + 1) - 1
    diff <- outer(qs, demand[, s], `-`)
    mats[[s]] <- holding[s] * pmax(diff, 0) + stockout[s] * pmax(-diff, 0)
  }

  mats
}

# Enumerate all integer plans within capacity and budget, aborting before
# allocation if the candidate space exceeds `limit`.
enumerate_plans <- function(
  max_units,
  costs,
  budget,
  limit,
  call = caller_env()
) {
  tol <- 1e-09 * max(1, abs(budget))
  plans <- matrix(integer(0), nrow = 1L, ncol = 0L)
  spent <- 0

  for (s in seq_along(max_units)) {
    q_max <- max_units[s]
    if (costs[s] > 0) {
      q_max <- min(q_max, floor((budget + tol) / costs[s]))
    }
    n_q <- q_max + 1
    if (nrow(plans) * n_q > limit) {
      cli_abort(
        c(
          "Enumerating feasible plans requires more than \\
           {format(limit, scientific = FALSE)} candidate combinations.",
          i = "Reduce `capacity` or `budget`, or raise \\
           `max_combinations`.",
          i = "No partial result is returned."
        ),
        call = call
      )
    }
    idx <- rep(seq_len(nrow(plans)), each = n_q)
    q_new <- rep(seq_len(n_q) - 1L, times = nrow(plans))
    spent <- spent[idx] + q_new * costs[s]
    keep <- spent <= budget + tol
    plans <- cbind(plans[idx[keep], , drop = FALSE], q_new[keep])
    spent <- spent[keep]
  }

  plans
}

# Position of the empirical VaR in the sorted losses, robust to floating
# point error in `alpha * n`.
tail_k <- function(alpha, n) {
  k <- ceiling(alpha * n - 1e-09)
  max(1L, min(as.integer(n), as.integer(k)))
}

# Empirical CVaR with fractional weight on the boundary scenario at the
# VaR, not a top-ceil average.
cvar_value <- function(x, k, alpha_n, tail_mass) {
  v <- sort.int(x, partial = k)[k]
  above <- x > v
  boundary <- max(0, sum(!above) - alpha_n)
  (sum(x[above]) + boundary * v) / tail_mass
}

cvar_empirical <- function(loss, alpha) {
  n <- length(loss)
  k <- tail_k(alpha, n)
  alpha_n <- alpha * n
  tail_mass <- n - alpha_n
  v <- sort.int(loss, partial = k)[k]
  above <- loss > v
  n_le <- n - sum(above)
  boundary <- max(0, n_le - alpha_n)

  list(
    var = v,
    cvar = (sum(loss[above]) + boundary * v) / tail_mass,
    tail_mass = tail_mass,
    n_above = sum(above),
    n_at = sum(loss == v),
    boundary_weight = boundary
  )
}

# Exact global minimization of the empirical CVaR over all enumerated
# plans, evaluating candidates in chunks to bound memory.
minimize_empirical_cvar <- function(plans, loss_mats, alpha, reps) {
  n <- nrow(plans)
  k <- tail_k(alpha, reps)
  alpha_n <- alpha * reps
  tail_mass <- reps - alpha_n
  chunk_size <- max(1L, min(4096L, floor(2e7 / reps)))
  best_value <- Inf
  best_index <- NA_integer_

  for (start in seq(1L, n, by = chunk_size)) {
    rows <- start:min(n, start + chunk_size - 1L)
    total <- matrix(0, nrow = length(rows), ncol = reps)
    for (s in seq_along(loss_mats)) {
      total <- total + loss_mats[[s]][plans[rows, s] + 1L, , drop = FALSE]
    }
    values <- apply(total, 1L, cvar_value, k, alpha_n, tail_mass)
    chunk_min <- min(values)
    if (chunk_min < best_value) {
      best_value <- chunk_min
      best_index <- rows[which.min(values)]
    }
  }

  list(index = best_index, value = best_value)
}

scenario_losses <- function(demand, plan, holding, stockout) {
  stockout_loss <- 0
  holding_loss <- 0
  for (s in seq_len(ncol(demand))) {
    d <- demand[, s]
    q <- plan[s]
    stockout_loss <- stockout_loss + stockout[s] * pmax(d - q, 0)
    holding_loss <- holding_loss + holding[s] * pmax(q - d, 0)
  }

  tibble::tibble(
    scenario = seq_len(nrow(demand)),
    stockout_loss = stockout_loss,
    holding_loss = holding_loss,
    loss = stockout_loss + holding_loss
  )
}

fmt_num <- function(x) {
  format(signif(x, 6), trim = TRUE, scientific = FALSE)
}
