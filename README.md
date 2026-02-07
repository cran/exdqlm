exdqlm — Extended Dynamic Quantile Linear Models
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/AntonioAPDL/exdqlm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AntonioAPDL/exdqlm/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

`exdqlm` fits **dynamic quantile** state-space models using the
**extended Asymmetric Laplace** (exAL) error family. It targets
time-series applications where the goal is to model **conditional
quantiles** (not only means) while keeping a familiar state-space
specification (design/evolution matrices and a state vector). In
**v0.3.0**, we add **optional C++ bridges** for speed (Kalman filter and
samplers) and **ELBO monitoring** for the variational routine. Defaults
remain **pure-R** for easy installation.

> **Terminology.** We say **exAL** for the extended Asymmetric Laplace
> family used here. It extends the standard AL by adding a **skewness**
> parameter; AL is the special case with zero skewness.

## Installation

CRAN (when available):

``` r
install.packages("exdqlm")
```

Development (GitHub):

``` r
# install.packages("pak")
pak::pak("AntonioAPDL/exdqlm")
```

## Quick start (≤ 10 lines)

Local-level model at a **single quantile** (the median). We fix
**scale** and **skewness** to keep it fast and stable for CRAN; we keep
the **pure-R** path.

``` r
set.seed(1)
library(exdqlm)

T      <- 120
state  <- cumsum(rnorm(T, sd = 0.2))
y      <- state + rnorm(T, sd = 1.0)

model  <- list(FF = matrix(1), GG = matrix(1), m0 = 0, C0 = 100)
options(exdqlm.use_cpp_kf = FALSE, exdqlm.use_cpp_samplers = FALSE)

fit <- exdqlmISVB(
  y = y, p0 = 0.5, model = model, df = 0.98, dim.df = 1,
  fix.sigma = TRUE, sig.init = 1.0,
  fix.gamma = TRUE, gam.init = 0.0
)
#> ISVB converged: 2 iterations, 0.514 seconds

tail(fit$diagnostics$elbo, 3)
#> [1] -113.62048  -67.45699
```

## Core concepts (at a glance)

- **State-space skeleton**: *design* (`FF`) and *evolution* (`GG`)
  matrices with a prior for the **state vector** (`m0`, `C0`).
- **Quantile of interest**: `p0` (e.g., `0.1`, `0.5`, `0.9`).
- **exAL errors**: controlled by **scale** and **skewness**; fixing them
  often stabilizes small examples.
- **Discount factors**: `df` and `dim.df` control evolution per block
  (e.g., trend vs seasonality).
- **ELBO**: recorded at `fit$diagnostics$elbo` (weakly monotone up to
  importance-sampling noise).

## What’s new in v0.3.0

- **C++ Kalman filter bridge** (drop-in): set
  `options(exdqlm.use_cpp_kf = TRUE)` (default **FALSE**; examples keep
  it OFF).
- **ELBO monitoring** for IS-VB: per-iteration values in
  `fit$diagnostics$elbo`.
- **Optional C++ samplers** for posterior draws: set
  `options(exdqlm.use_cpp_samplers = TRUE)` (default **FALSE**).
- **Portability**: OpenMP is **optional and gated**; package runs
  serially when unavailable.

> Keep both options **FALSE** in examples/CI. Enable them locally if
> your toolchain supports compiled code.

### Runtime options (summary)

| Option                    | Default | Effect                           | Use when…                                |
|---------------------------|:-------:|----------------------------------|------------------------------------------|
| `exdqlm.use_cpp_kf`       |  FALSE  | C++ Kalman filter bridge         | you have compilers/OpenMP and want speed |
| `exdqlm.use_cpp_samplers` |  FALSE  | C++ samplers for posterior draws | same as above; keep OFF on CRAN/examples |

Set with:

``` r
options(exdqlm.use_cpp_kf = TRUE)
options(exdqlm.use_cpp_samplers = TRUE)
```

## Minimal examples (CRAN-safe)

### 1) Single-quantile fit on built-in data (tiny slice)

Trend + seasonality + one regressor (`nino34`). **Note**: `FF` for the
regressor is `1 × T`. Combine components **pairwise**.

``` r
set.seed(2)
T <- 150
y <- log(BTflow[seq_len(T)])
x <- nino34[seq_len(T)]

trend.comp <- polytrendMod(order = 1, m0 = 0, C0 = 1)
seas.comp  <- seasMod(p = 12, h = 1, C0 = diag(1, 2))

# 1-d regressor block (explicit 1 x T design)
reg.comp <- list(m0 = 0, C0 = 1, FF = matrix(x, nrow = 1), GG = matrix(1))

# combine pairwise
base.mod <- combineMods(trend.comp, seas.comp)
model    <- combineMods(base.mod, reg.comp)

# one discount per block: (trend, seasonal[2-d], reg)
df     <- c(1.00, 0.98, 1.00)
dim.df <- c(1,       2,   1)

options(exdqlm.use_cpp_kf = FALSE, exdqlm.use_cpp_samplers = FALSE)

fit <- exdqlmISVB(
  y = y, p0 = 0.5, model = model,
  df = df, dim.df = dim.df,
  fix.sigma = TRUE, sig.init = 0.2,
  fix.gamma = TRUE, gam.init = 0.0
)
#> ISVB converged: 2 iterations, 0.547 seconds

# quick checks
tail(fit$diagnostics$elbo, 2)
#> [1] -1078.8072  -934.9032
dim(fit$theta.out$sm)  # state-dimension x time
#> [1]   4 150
```

### 2) exAL helper sanity check (CDF ↔ quantile)

``` r
set.seed(3)
x      <- seq(-2, 2, length.out = 5)
p0     <- 0.25
mu     <- 0
sigma  <- 1
gamma  <- 0.0

# CDF then invert with QF — should approximately return x
cdf_vals <- pexal(x,  p0 = p0, mu = mu, sigma = sigma, gamma = gamma)
x_back   <- qexal(cdf_vals, p0 = p0, mu = mu, sigma = sigma, gamma = gamma)

round(cbind(x, x_back), 4)
#>       x x_back
#> [1,] -2     -2
#> [2,] -1     -1
#> [3,]  0      0
#> [4,]  1      1
#> [5,]  2      2

# A few random draws
rexal(5, p0 = p0, mu = mu, sigma = sigma, gamma = gamma)
#> [1] -0.5296664  5.4402490  0.7934288  0.4376783  2.5354967
```

> **CRAN-safety.** All examples set a seed, use tiny data, finish in a
> few seconds, and keep the pure-R path by default.

## FAQ / Troubleshooting

- **It runs slowly.** Use short series (≤ 200), fix
  **scale**/**skewness**, and keep discount factors near but below one
  (≈ 0.96–0.99). Enable C++ bridges only if your toolchain supports
  them.

- **ELBO dips slightly—bug?** Small downward blips are expected from
  importance-sampling noise. Look for an overall upward trend; if not,
  simplify the model or adjust variance/discounts.

- **OpenMP not available.** That’s fine. It is optional. Everything runs
  serially; examples here use the pure-R path.

- **Numerical stability tips.** Avoid extremely tight `C0`; start with
  moderate priors (e.g., `C0` around 1–100 for simple models), and fix
  **scale**/**skewness** for initial runs.

## How to cite

Barata, R., Prado, R., & Sansó, B. (2022). *Fast inference for
time-varying quantiles via flexible dynamic models with application to
the characterization of atmospheric rivers*. **Annals of Applied
Statistics**, 16(1), 247–271. <https://doi.org/10.1214/21-AOAS1497>

## License

MIT © The authors. See `LICENSE`.

## Getting help

Open an issue: <https://github.com/AntonioAPDL/exdqlm/issues>
