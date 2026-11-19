# Scenario-driven replenishment under joint demand uncertainty

Demand for several products often rises and falls together in the same
period, so resampling each product separately can misstate joint
stockout risk. `scenario_replenish()` builds *joint* demand scenarios
for several stock-keeping units (SKUs) and chooses a replenishment plan
that globally minimizes the empirical conditional value-at-risk (CVaR)
of stockout plus holding loss over integer plans.

The procedure:

1.  **Align histories by period.** Every SKU's history is restricted to
    the periods shared by all SKUs, so that each scenario resamples the
    same periods for every SKU.

2.  **Draw joint bootstrap scenarios.** For each scenario, one vector of
    period indices is drawn with replacement and applied to every SKU,
    preserving cross-SKU co-movement.

3.  **Fit and predict natively.** In each scenario, each SKU's resampled
    history is fitted with the same Gaussian engine used by
    [`fit.infer()`](https://infer.tidymodels.org/dev/reference/fit.infer.md)
    ([`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
    `family = "gaussian"`), and demand is predicted from that SKU's own
    future features. Negative predictions are truncated to zero. No
    precomputed demand is ever substituted for these fits and
    predictions.

4.  **Optimize exactly.** The empirical CVaR of total stockout plus
    holding loss is minimized globally over all integer replenishment
    plans respecting each SKU's `capacity` and the total `budget`. The
    CVaR retains the fractional weight of scenarios sitting exactly at
    the tail boundary rather than averaging a top-`ceil()` count of
    scenario losses.

If any step fails—invalid input, a failed or rank-deficient fit, a
non-finite prediction, or a candidate space larger than
`max_combinations`—the function errors and returns nothing rather than
delivering a partial result.

## Usage

``` r
scenario_replenish(
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
  max_combinations = 1e+06
)

# S3 method for class 'scenario_replenishment'
print(x, ...)
```

## Arguments

- data:

  A data frame with one row per SKU and period, containing the `sku` and
  `period` columns, the response column named in `formula`, and every
  predictor column used in `formula`.

- formula:

  A two-sided formula whose left-hand side is a bare numeric column of
  `data` (e.g. `demand ~ price + promo`). The same formula is fitted per
  SKU with `family = "gaussian"`, as in
  [`fit.infer()`](https://infer.tidymodels.org/dev/reference/fit.infer.md).

- future:

  A data frame with exactly one row per SKU, containing the `sku` column
  and the future value of every predictor in `formula`.

- sku, period:

  Unquoted names of the SKU and period columns in `data`; `future` must
  use the same `sku` column name.

- reps:

  Number of joint bootstrap scenarios.

- unit_cost:

  Budget cost per unit ordered.

- holding_cost:

  Loss per unit left over after demand is met.

- stockout_cost:

  Loss per unit of unmet demand.

- capacity:

  Maximum units to order per SKU; must be whole-numbered.

- budget:

  Total budget; the sum of `unit_cost * order` across SKUs must not
  exceed it.

- alpha:

  CVaR level in `(0, 1)`; the mean loss over the worst `1 - alpha` tail
  mass of scenarios is minimized.

- max_combinations:

  Safety limit on the number of candidate integer plans; the function
  errors if enumeration would exceed it.

- x:

  An object of class `scenario_replenishment`, i.e. output from
  `scenario_replenish()`.

- ...:

  Arguments passed to methods (currently ignored).

## Value

An object of class `scenario_replenishment`, a list with components:

- `plan`: tibble with the optimal integer `order` per SKU, its
  `unit_cost`, `spend`, `capacity`, and `capacity_effective` (the
  capacity tightened by the largest scenario demand).

- `budget`: tibble with `budget`, `spent`, and `remaining`.

- `objective`: tibble with `alpha`, the number of `scenarios`, the
  empirical `cvar` and `var` (value-at-risk), the `mean_loss`, the
  `tail_mass` in scenario-equivalents, the counts of scenarios above and
  at the VaR, and the fractional `boundary_weight` entering the CVaR.

- `tail`: tibble of tail scenarios with their CVaR `weight` (1 for
  scenarios above the VaR, fractional at the boundary).

- `losses`: tibble of per-scenario `stockout_loss`, `holding_loss`, and
  total `loss` under the chosen plan.

- `demand`: tibble of joint demand scenarios with the raw Gaussian
  `prediction` and the `demand` truncated at zero per scenario and SKU.

- `fits`: tibble of the native per-scenario fits (`term`, `estimate`)
  per SKU, mirroring
  [`fit.infer()`](https://infer.tidymodels.org/dev/reference/fit.infer.md)
  output.

- `population`: tibble of the aligned sampling population that scenarios
  were resampled from.

- `scenario_index`: tibble of the resampled period indices per scenario,
  shared across SKUs.

- `meta`: a list of inputs needed to interpret the plan.

## Details

`unit_cost`, `holding_cost`, `stockout_cost`, and `capacity` may each be
a single number recycled to every SKU or a named numeric vector with one
element per SKU.

## Scope

All reported risk quantities are *empirical* quantities of the
regression-bootstrap scenarios resampled from the aligned history. They
describe resampling risk only; they are not a guarantee of future demand
coverage or of population-level risk.

## Reproducibility

When using the infer package for research, or in other cases when exact
reproducibility is a priority, be sure the set the seed for R’s random
number generator. infer will respect the random seed specified in the
[`set.seed()`](https://rdrr.io/r/base/Random.html) function, returning
the same result when
[`generate()`](https://infer.tidymodels.org/dev/reference/generate.md)ing
data given an identical seed. For instance, we can calculate the
difference in mean `age` by `college` degree status using the `gss`
dataset from 10 versions of the `gss` resampled with permutation using
the following code.

    set.seed(1)

    gss |>
      specify(age ~ college) |>
      hypothesize(null = "independence") |>
      generate(reps = 5, type = "permute") |>
      calculate("diff in means", order = c("degree", "no degree"))

    ## Response: age (numeric)
    ## Explanatory: college (factor)
    ## Null Hypothesis: indepe...
    ## # A tibble: 5 x 2
    ##   replicate   stat
    ##       <int>  <dbl>
    ## 1         1 -0.531
    ## 2         2 -2.35
    ## 3         3  0.764
    ## 4         4  0.280
    ## 5         5  0.350

Setting the seed to the same value again and rerunning the same code
will produce the same result.

    # set the seed
    set.seed(1)

    gss |>
      specify(age ~ college) |>
      hypothesize(null = "independence") |>
      generate(reps = 5, type = "permute") |>
      calculate("diff in means", order = c("degree", "no degree"))

    ## Response: age (numeric)
    ## Explanatory: college (factor)
    ## Null Hypothesis: indepe...
    ## # A tibble: 5 x 2
    ##   replicate   stat
    ##       <int>  <dbl>
    ## 1         1 -0.531
    ## 2         2 -2.35
    ## 3         3  0.764
    ## 4         4  0.280
    ## 5         5  0.350

Please keep this in mind when writing infer code that utilizes
resampling with
[`generate()`](https://infer.tidymodels.org/dev/reference/generate.md).

## Examples

``` r
set.seed(1)

# two SKUs whose demands move together period by period
history <- expand.grid(period = 1:24, sku = c("a", "b"))
history$demand <- ifelse(history$sku == "a", 20, 35) +
  2 * history$period +
  rep(rnorm(24, sd = 3), times = 2) +
  rnorm(48, sd = 1)

future <- data.frame(sku = c("a", "b"), period = 25)

scenario_replenish(
  history,
  demand ~ period,
  future,
  sku = sku,
  period = period,
  reps = 100,
  unit_cost = 2,
  holding_cost = 1,
  stockout_cost = 8,
  capacity = 200,
  budget = 800,
  alpha = 0.9
)
#> # Scenario-driven replenishment
#> SKUs: 2 | aligned periods: 24 | scenarios: 100 | alpha: 0.9
#> Budget: 800 | spent: 322 | remaining: 478
#> Empirical CVaR of stockout + holding loss: 9.67372 (VaR: 7.35715, mean: 4.87807)
#> 
#> Replenishment plan (integer units):
#> # A tibble: 2 × 6
#>   sku   order unit_cost spend capacity capacity_effective
#>   <fct> <int>     <dbl> <dbl>    <dbl>              <dbl>
#> 1 a        73         2   146      200                 74
#> 2 b        88         2   176      200                 89
#> 
#> Tail risk: 10 scenario-equivalents = 10 scenario(s) above the VaR at full weight + 1 at the VaR sharing fractional weight 0.
#> CVaR keeps the fractional tail-boundary weight (no top-ceil averaging).
#> Note: empirical risk over the 100 regression-bootstrap scenarios only; not a
#>   guarantee of future demand coverage or population-level risk.
```
