# Helper data -------------------------------------------------------------

sr_panel <- function(n = 10, seed = 11) {
  set.seed(seed)
  data.frame(
    sku = rep(c("a", "b"), each = n),
    period = rep(seq_len(n), 2),
    demand = c(
      10 + seq_len(n) + rnorm(n, sd = 1),
      20 + 2 * seq_len(n) + rnorm(n, sd = 2)
    ),
    x = rep(seq_len(n), 2) + rnorm(2 * n, sd = 0.2)
  )
}

sr_future <- data.frame(sku = c("a", "b"), x = c(11, 11))

# Independent CVaR via the Rockafellar-Uryasev minimization, for
# cross-checking the fractional-weight implementation.
ru_cvar <- function(loss, alpha) {
  n <- length(loss)
  f <- function(eta) eta + sum(pmax(loss - eta, 0)) / (n * (1 - alpha))
  stats::optimize(f, interval = range(loss))$objective
}

# Structure and reproducibility --------------------------------------------

test_that("scenario_replenish returns a complete, well-formed object", {
  set.seed(1)
  res <- scenario_replenish(
    sr_panel(), demand ~ x, sr_future,
    sku = sku, period = period,
    reps = 50, unit_cost = 2, holding_cost = 1, stockout_cost = 5,
    capacity = 30, budget = 100, alpha = 0.9
  )

  expect_s3_class(res, "scenario_replenishment")
  expect_named(
    res,
    c(
      "plan", "budget", "objective", "tail", "losses", "demand", "fits",
      "population", "scenario_index", "meta"
    )
  )

  # plan: integer orders within capacity and budget
  expect_equal(nrow(res$plan), 2)
  expect_true(all(res$plan$order == floor(res$plan$order)))
  expect_true(all(res$plan$order >= 0))
  expect_true(all(res$plan$order <= res$plan$capacity))
  expect_lte(res$budget$spent, res$budget$budget + 1e-9)
  expect_equal(res$budget$spent, sum(res$plan$spend))
  expect_equal(res$budget$remaining, res$budget$budget - res$budget$spent)

  # capacity_effective is tightened by the largest scenario demand
  max_d <- tapply(res$demand$demand, res$demand$sku, max)
  expect_equal(
    res$plan$capacity_effective,
    pmin(res$plan$capacity, ceiling(max_d[as.character(res$plan$sku)])),
    ignore_attr = TRUE
  )

  # one row per scenario x sku in demand and fits
  expect_equal(nrow(res$demand), 50 * 2)
  expect_true(all(res$demand$demand >= 0))
  expect_true(all(res$demand$demand == pmax(res$demand$prediction, 0)))
  expect_equal(unique(res$fits$term), c("intercept", "x"))
  expect_equal(nrow(res$fits), 50 * 2 * 2)

  # population is aligned: same periods for every sku
  periods <- split(res$population$period, res$population$sku)
  expect_equal(periods[[1]], periods[[2]])

  # scenario_index: one shared draw vector per scenario
  expect_equal(nrow(res$scenario_index), 50 * 10)
  expect_true(all(res$scenario_index$period_index %in% seq_len(10)))

  # losses are consistent with the objective
  expect_equal(
    res$losses$loss,
    res$losses$stockout_loss + res$losses$holding_loss
  )
  expect_equal(res$objective$cvar, cvar_empirical(res$losses$loss, 0.9)$cvar)
  expect_equal(nrow(res$tail), res$objective$n_tail_above +
    res$objective$n_tail_at)
  expect_equal(sum(res$tail$weight), res$objective$tail_mass)
})

test_that("scenario_replenish is reproducible given a seed", {
  args <- list(
    sr_panel(), demand ~ x, sr_future,
    sku = "sku", period = "period",
    reps = 30, unit_cost = 2, holding_cost = 1, stockout_cost = 5,
    capacity = 30, budget = 100
  )
  set.seed(42)
  res1 <- do.call(scenario_replenish, args)
  set.seed(42)
  res2 <- do.call(scenario_replenish, args)

  expect_identical(res1$demand, res2$demand)
  expect_identical(res1$plan, res2$plan)
  expect_identical(res1$objective, res2$objective)
  expect_identical(res1$scenario_index, res2$scenario_index)
})

# Joint bootstrap and native fits -------------------------------------------

test_that("scenarios share period indices and use native glm fits", {
  set.seed(99)
  reps <- 5
  res <- scenario_replenish(
    sr_panel(n = 8), demand ~ x, sr_future,
    sku = sku, period = period,
    reps = reps, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
    capacity = 30, budget = 100
  )

  idx <- matrix(res$scenario_index$period_index, nrow = 8, ncol = reps)

  for (r in seq_len(reps)) {
    for (s in c("a", "b")) {
      hist_s <- as.data.frame(res$population[res$population$sku == s, ])
      # the same index vector resamples both SKUs in this scenario
      m <- stats::glm(
        demand ~ x,
        data = hist_s[idx[, r], ],
        family = stats::gaussian
      )
      p <- stats::predict(
        m,
        newdata = sr_future[sr_future$sku == s, ],
        type = "response"
      )

      got <- res$demand[res$demand$scenario == r & res$demand$sku == s, ]
      expect_equal(got$prediction, unname(p))
      expect_equal(got$demand, max(0, unname(p)))

      cf <- res$fits[res$fits$scenario == r & res$fits$sku == s, ]
      expect_equal(
        cf$estimate[match(c("intercept", "x"), cf$term)],
        unname(stats::coef(m))
      )
    }
  }
})

test_that("negative predictions are truncated to zero demand", {
  set.seed(1)
  n <- 6
  declining <- data.frame(
    sku = "a",
    period = seq_len(n),
    demand = 100 - 20 * seq_len(n),
    x = seq_len(n)
  )
  future <- data.frame(sku = "a", x = 10)

  res <- scenario_replenish(
    declining, demand ~ x, future,
    sku = sku, period = period,
    reps = 10, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
    capacity = 50, budget = 100
  )

  expect_true(all(res$demand$prediction < 0))
  expect_true(all(res$demand$demand == 0))
  # with zero demand in every scenario, ordering nothing is optimal
  expect_equal(res$plan$order, 0L)
  expect_equal(res$objective$cvar, 0)
})

# Empirical CVaR with fractional tail weight ---------------------------------

test_that("cvar_empirical keeps the fractional tail-boundary weight", {
  loss <- 1:10

  # alpha * n is not an integer: the boundary scenario gets weight 0.5
  res <- cvar_empirical(loss, 0.75)
  expect_equal(res$var, 8)
  expect_equal(res$cvar, (9 + 10 + 0.5 * 8) / 2.5)
  expect_equal(res$boundary_weight, 0.5)
  # a top-ceil average would give mean(8:10) = 9 instead
  expect_true(res$cvar != mean(8:10))

  # alpha * n is an integer: boundary weight is zero
  res_int <- cvar_empirical(loss, 0.8)
  expect_equal(res_int$cvar, mean(9:10))
  expect_equal(res_int$boundary_weight, 0)

  # ties at the VaR share the boundary weight
  res_ties <- cvar_empirical(c(5, 5, 5, 1, 2), 0.6)
  expect_equal(res_ties$cvar, 5)
  expect_equal(res_ties$n_at, 3)

  # constant losses
  expect_equal(cvar_empirical(rep(7, 20), 0.9)$cvar, 7)

  # cross-check against the Rockafellar-Uryasev minimization
  set.seed(3)
  l <- rnorm(200, mean = 10, sd = 5)
  expect_equal(cvar_empirical(l, 0.9)$cvar, ru_cvar(l, 0.9), tolerance = 1e-4)
  expect_equal(cvar_empirical(l, 0.77)$cvar, ru_cvar(l, 0.77), tolerance = 1e-4)
})

# Optimization ---------------------------------------------------------------

test_that("the optimizer globally minimizes empirical CVaR", {
  set.seed(7)
  alpha <- 0.8
  res <- scenario_replenish(
    sr_panel(), demand ~ x, sr_future,
    sku = sku, period = period,
    reps = 80, unit_cost = c(a = 1, b = 1), holding_cost = c(a = 1, b = 2),
    stockout_cost = c(a = 5, b = 3), capacity = c(a = 4, b = 4),
    budget = 5, alpha = alpha
  )

  d <- matrix(res$demand$demand, nrow = 80, ncol = 2)
  h <- c(1, 2)
  p <- c(5, 3)

  # independent brute force over every feasible integer plan
  candidates <- expand.grid(q1 = 0:4, q2 = 0:4)
  candidates <- candidates[candidates$q1 + candidates$q2 <= 5 + 1e-9, ]
  cand_cvar <- apply(candidates, 1, function(q) {
    loss <- h[1] * pmax(q[1] - d[, 1], 0) + p[1] * pmax(d[, 1] - q[1], 0) +
      h[2] * pmax(q[2] - d[, 2], 0) + p[2] * pmax(d[, 2] - q[2], 0)
    ru_cvar(loss, alpha)
  })

  # the reported objective matches the best feasible candidate
  expect_equal(
    res$objective$cvar, min(cand_cvar),
    tolerance = 1e-4
  )

  # the reported plan is feasible and achieves the reported objective
  q <- res$plan$order
  expect_lte(sum(q), 5 + 1e-9)
  loss_q <- h[1] * pmax(q[1] - d[, 1], 0) + p[1] * pmax(d[, 1] - q[1], 0) +
    h[2] * pmax(q[2] - d[, 2], 0) + p[2] * pmax(d[, 2] - q[2], 0)
  expect_equal(ru_cvar(loss_q, alpha), res$objective$cvar, tolerance = 1e-4)
})

test_that("budget and capacity constraints are respected", {
  set.seed(5)
  res <- scenario_replenish(
    sr_panel(), demand ~ x, sr_future,
    sku = sku, period = period,
    reps = 20, unit_cost = c(a = 2, b = 3), holding_cost = 1,
    stockout_cost = 5, capacity = c(a = 10, b = 10), budget = 0
  )
  expect_equal(res$plan$order, c(0L, 0L))
  expect_equal(res$budget$spent, 0)

  set.seed(5)
  res_cap <- scenario_replenish(
    sr_panel(), demand ~ x, sr_future,
    sku = sku, period = period,
    reps = 20, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
    capacity = c(a = 0, b = 10), budget = 100
  )
  expect_equal(res_cap$plan$order[res_cap$plan$sku == "a"], 0L)
})

test_that("named per-SKU parameter vectors are matched by name", {
  set.seed(8)
  res <- scenario_replenish(
    sr_panel(), demand ~ x, sr_future,
    sku = sku, period = period,
    reps = 20, unit_cost = c(b = 3, a = 2), holding_cost = 1,
    stockout_cost = 5, capacity = 10, budget = 100
  )
  expect_equal(res$plan$unit_cost[res$plan$sku == "a"], 2)
  expect_equal(res$plan$unit_cost[res$plan$sku == "b"], 3)
})

test_that("enumeration beyond max_combinations fails without a result", {
  out <- tryCatch(
    scenario_replenish(
      sr_panel(), demand ~ x, sr_future,
      sku = sku, period = period,
      reps = 10, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
      capacity = 10, budget = 100, max_combinations = 5
    ),
    error = function(e) e
  )
  expect_s3_class(out, "error")
  expect_match(conditionMessage(out), "candidate combinations")
})

# Alignment ------------------------------------------------------------------

test_that("histories are aligned to shared periods with a message", {
  panel <- sr_panel()
  panel <- panel[!(panel$sku == "b" & panel$period == 5), ]

  expect_message(
    res <- scenario_replenish(
      panel, demand ~ x, sr_future,
      sku = sku, period = period,
      reps = 10, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
      capacity = 10, budget = 100
    ),
    "1 period dropped"
  )
  expect_false(5 %in% res$population$period)
  expect_equal(length(unique(res$population$period)), 9)
  expect_true(all(res$scenario_index$period != 5))
})

test_that("fewer than two shared periods is an error", {
  panel <- sr_panel()
  panel <- panel[!(panel$sku == "b" & panel$period %in% 2:10), ]

  expect_error(
    scenario_replenish(
      panel, demand ~ x, sr_future,
      sku = sku, period = period,
      reps = 10, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
      capacity = 10, budget = 100
    ),
    "fewer than 2 common periods"
  )
})

test_that("duplicated sku-period rows are an error", {
  panel <- rbind(sr_panel(), sr_panel()[1, ])
  expect_error(
    scenario_replenish(
      panel, demand ~ x, sr_future,
      sku = sku, period = period,
      reps = 10, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
      capacity = 10, budget = 100
    ),
    "at most one row per"
  )
})

# Input validation -----------------------------------------------------------

test_that("invalid inputs error and return no partial result", {
  base <- list(
    data = sr_panel(), formula = demand ~ x, future = sr_future,
    sku = "sku", period = "period", reps = 10, unit_cost = 1,
    holding_cost = 1, stockout_cost = 5, capacity = 10, budget = 100
  )
  run <- function(...) {
    args <- base
    dots <- list(...)
    args[names(dots)] <- dots
    tryCatch(
      do.call(scenario_replenish, args),
      error = function(e) conditionMessage(e)
    )
  }

  # missing values in data or future
  bad <- sr_panel()
  bad$demand[3] <- NA
  expect_match(run(data = bad), "missing values")
  bad_f <- sr_future
  bad_f$x[1] <- NA
  expect_match(run(future = bad_f), "missing values")

  # non-numeric response
  bad2 <- sr_panel()
  bad2$demand <- as.character(bad2$demand)
  expect_match(run(data = bad2), "must be numeric")

  # transformed response
  expect_match(run(formula = log(demand) ~ x), "bare column name")

  # sku as predictor
  expect_match(run(formula = demand ~ sku), "cannot be a predictor")

  # bad scalars
  expect_match(run(alpha = 1), "between 0 and 1")
  expect_match(run(alpha = 0), "between 0 and 1")
  expect_match(run(reps = 0), "at least 1")
  expect_match(run(reps = 2.5), "whole number")
  expect_match(run(budget = -1), "non-negative")

  # parameter vectors
  expect_match(run(unit_cost = c(1, 2)), "named vector")
  expect_match(run(unit_cost = c(a = 1, c = 2)), "named vector")
  expect_match(run(capacity = 2.5), "whole-numbered")
  expect_match(run(holding_cost = -1), "non-negative")

  # future with missing or extra SKU
  expect_match(run(future = sr_future[1, ]), "exactly one row per SKU")
  extra <- rbind(sr_future, data.frame(sku = "c", x = 1))
  expect_match(run(future = extra), "exactly one row per SKU")

  # future categorical level not observed in that SKU's history
  panel_cat <- sr_panel()
  panel_cat$promo <- rep(c("x", "y"), length.out = nrow(panel_cat))
  future_bad <- data.frame(sku = c("a", "b"), x = 1, promo = "z")
  expect_match(
    run(
      data = panel_cat, formula = demand ~ x + promo,
      future = future_bad
    ),
    "not observed"
  )

  # every failure above produced an error message, not a result
  expect_type(run(alpha = 1), "character")
})

test_that("rank-deficient scenario fits error instead of returning NA", {
  set.seed(1)
  panel <- sr_panel()
  panel$x2 <- 2 * panel$x
  future <- data.frame(sku = c("a", "b"), x = 11, x2 = 22)

  expect_error(
    scenario_replenish(
      panel, demand ~ x + x2, future,
      sku = sku, period = period,
      reps = 10, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
      capacity = 10, budget = 100
    ),
    "rank-deficient"
  )
})

test_that("a resample losing a categorical level errors cleanly", {
  panel <- data.frame(
    sku = rep("a", 3),
    period = 1:3,
    demand = c(10, 12, 11),
    promo = c("x", "y", "y")
  )
  future <- data.frame(sku = "a", promo = "x")

  # resamples of size 3 frequently miss one of the two promo levels
  set.seed(2)
  expect_error(
    scenario_replenish(
      panel, demand ~ promo, future,
      sku = sku, period = period,
      reps = 30, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
      capacity = 5, budget = 10
    ),
    "Gaussian fit failed|rank-deficient"
  )
})

test_that("categorical predictors work when levels are observed", {
  set.seed(13)
  n <- 12
  panel <- data.frame(
    sku = rep(c("a", "b"), each = n),
    period = rep(seq_len(n), 2),
    demand = c(
      10 + seq_len(n) + rnorm(n, sd = 0.5),
      20 + seq_len(n) + rnorm(n, sd = 0.5)
    ),
    promo = rep(c("x", "y"), length.out = 2 * n)
  )
  future <- data.frame(sku = c("a", "b"), promo = c("x", "y"))

  res <- scenario_replenish(
    panel, demand ~ promo, future,
    sku = sku, period = period,
    reps = 20, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
    capacity = 40, budget = 100
  )
  expect_s3_class(res, "scenario_replenishment")
  expect_true(all(c("intercept", "promoy") %in% res$fits$term))
})

test_that("intercept-only formulas and numeric SKU ids work", {
  set.seed(21)
  panel <- data.frame(
    sku = rep(1:2, each = 8),
    period = rep(1:8, 2),
    demand = c(10 + rnorm(8), 20 + rnorm(8))
  )
  future <- data.frame(sku = 1:2)

  res <- scenario_replenish(
    panel, demand ~ 1, future,
    sku = sku, period = period,
    reps = 10, unit_cost = 1, holding_cost = 1, stockout_cost = 5,
    capacity = 30, budget = 100
  )
  expect_s3_class(res, "scenario_replenishment")
  expect_equal(unique(res$fits$term), "intercept")
  expect_equal(res$plan$sku, 1:2)
})

# Print method -----------------------------------------------------------------

test_that("print method summarizes the plan and tail risk honestly", {
  set.seed(1)
  res <- scenario_replenish(
    sr_panel(), demand ~ x, sr_future,
    sku = sku, period = period,
    reps = 40, unit_cost = 2, holding_cost = 1, stockout_cost = 5,
    capacity = 30, budget = 100, alpha = 0.9
  )

  expect_output(print(res), "Scenario-driven replenishment")
  expect_output(print(res), "Empirical CVaR")
  expect_output(print(res), "Replenishment plan")
  expect_output(print(res), "fractional weight")
  expect_output(print(res), "no top-ceil")
  expect_output(print(res), "empirical risk over the 40")
  expect_output(print(res), "not a")
  expect_output(print(res), "population-level risk")
})
