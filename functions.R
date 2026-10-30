## Functions shared by the two applications of the paper
## "On bootstrap in generalized estimating equations".

library(glmtoolbox)
library(aplore3)      # polypharm data
library(nlme)         # Milk data
library(gamlss)       # worm plots
library(pbivnorm)     # bivariate normal CDF (NORTA calibration)
library(Matrix)       # nearPD
library(ggplot2)
library(parallel)
library(pbapply)
pbapply::pboptions(type = "txt")

## ---- Model fitting ----
## The data always contain the columns `cid` (cluster) and `wave` (integer
## time); `waves` defines the lags used by the correlation estimator. Under
## independence `waves` is irrelevant, and glmgee() fails with `waves` in
## unbalanced data, so the model is then fitted without it.
fit_gee <- function(spec, data) {
  fit <- try(suppressWarnings(
    if (identical(spec$corstr, "Independence"))
      glmtoolbox::glmgee(spec$formula, id = cid, family = spec$family, corstr = spec$corstr, data = data)
    else glmtoolbox::glmgee(spec$formula, id = cid, waves = wave, family = spec$family,
                            corstr = spec$corstr, data = data)), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$phi)) return(NULL)
  fit
}

coef_names <- function(fit) {
  cf <- coef(fit)
  nm <- rownames(as.matrix(cf))
  if (is.null(nm)) nm <- names(cf)
  if (is.null(nm)) nm <- paste0("beta", seq_along(cf) - 1)
  nm
}

## beta (p), phi, rho1, rho2, converged (0/1)
summarize_fit <- function(fit, p) {
  if (is.null(fit)) return(c(rep(NA_real_, p + 3), 0))
  c(as.numeric(coef(fit)), fit$phi, fit$corr[1, 2], fit$corr[1, 3],
    as.numeric(isTRUE(fit$converged)))
}

## Fitted means in the row order of `data`
mean_from <- function(spec, data, beta) {
  X <- model.matrix(delete.response(terms(spec$formula)), data = data)
  as.vector(spec$linkinv(X %*% beta))
}

var_fun <- function(dist, mu) if (dist == "bernoulli") mu * (1 - mu) else mu^2

## ---- Data-generating process ----
## Target correlation matrix of each cluster: rho_lags[h] at lag h = |t_j - t_l|.
## rho_lags = fit$corr[1, -1] already contains the AR(2) Yule-Walker extension
## computed by glmtoolbox.
target_matrices <- function(idx, wave, rho_lags) {
  lapply(idx, function(ii) {
    w <- wave[ii]; L <- abs(outer(w, w, "-"))
    R <- matrix(0, length(ii), length(ii))
    pos <- L > 0 & L <= length(rho_lags)
    R[pos] <- rho_lags[L[pos]]
    diag(R) <- 1
    R
  })
}

## NORTA calibration (Emrich and Piedmonte, 1991; Cario and Nelson, 1997).
## For each pair (j, l), find the latent correlation r such that
##   Phi2(z_j, z_l; r) = p_j p_l + rho * sqrt(p_j q_j p_l q_l),
## with z = qnorm(p) and y = 1{Z > qnorm(1 - p)}. Phi2 is increasing in r, so a
## vectorized bisection solves all pairs at once. Targets outside the Frechet
## bounds are set at the bound and counted.
calibrate_norta <- function(mu, idx, target, n_iter = 50) {
  pairs <- do.call(rbind, lapply(seq_along(idx), function(i) {
    R <- target[[i]]
    if (nrow(R) < 2) return(NULL)
    jl <- which(upper.tri(R) & R != 0, arr.ind = TRUE)
    if (!nrow(jl)) return(NULL)
    cbind(i, jl[, 1], jl[, 2], mu[idx[[i]][jl[, 1]]], mu[idx[[i]][jl[, 2]]], R[jl])
  }))
  Rlat <- lapply(target, function(R) diag(nrow(R)))
  if (is.null(pairs)) return(list(R = Rlat, diag = c(n_pairs = 0, n_bound = 0, max_error = 0)))

  p1 <- pmin(pmax(pairs[, 4], 1e-10), 1 - 1e-10)
  p2 <- pmin(pmax(pairs[, 5], 1e-10), 1 - 1e-10)
  rho <- pairs[, 6]
  z1 <- qnorm(p1); z2 <- qnorm(p2)
  sdp <- sqrt(p1 * (1 - p1) * p2 * (1 - p2))
  target_p11 <- p1 * p2 + rho * sdp

  lim <- 1 - 1e-7
  lo <- rep(-lim, length(rho)); hi <- rep(lim, length(rho))
  at_bound <- target_p11 > pbivnorm(z1, z2, hi) | target_p11 < pbivnorm(z1, z2, lo)
  for (it in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    above <- pbivnorm(z1, z2, mid) > target_p11
    hi[above] <- mid[above]; lo[!above] <- mid[!above]
  }
  r <- (lo + hi) / 2
  rho_obtained <- (pbivnorm(z1, z2, r) - p1 * p2) / sdp

  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1]; j <- pairs[k, 2]; l <- pairs[k, 3]
    Rlat[[i]][j, l] <- Rlat[[i]][l, j] <- r[k]
  }
  list(R = Rlat,
       diag = c(n_pairs = length(r), n_bound = sum(at_bound),
                max_error = if (all(at_bound)) NA else max(abs(rho_obtained - rho)[!at_bound])))
}

## Latent matrices and their Cholesky factors. In the gamma model with small
## phi (large shape 1/phi) the normal-to-gamma transformation is nearly linear
## and no calibration is needed (checked with the oracle estimators).
prepare_dgp <- function(dist, mu, phi, rho_lags, idx, wave) {
  target <- target_matrices(idx, wave, rho_lags)
  if (dist == "bernoulli") {
    cal <- calibrate_norta(mu, idx, target)
    lat <- cal$R; dg <- cal$diag
  } else {
    lat <- target; dg <- c(n_pairs = NA, n_bound = NA, max_error = NA)
  }
  n_npd <- 0
  U <- lapply(lat, function(R) {
    U <- tryCatch(chol(R), error = function(e) NULL)
    if (is.null(U)) {
      n_npd <<- n_npd + 1
      U <- chol(as.matrix(Matrix::nearPD(R, corr = TRUE)$mat))
    }
    U
  })
  list(dist = dist, mu = mu, phi = phi, idx = idx, U = U,
       threshold = if (dist == "bernoulli") qnorm(1 - mu) else NULL,
       diag = c(dg, n_nearpd = n_npd, n_clusters = length(idx)))
}

generate_y <- function(dgp) {
  z <- numeric(length(dgp$mu))
  for (i in seq_along(dgp$idx)) {
    ii <- dgp$idx[[i]]
    z[ii] <- crossprod(dgp$U[[i]], rnorm(length(ii)))
  }
  if (dgp$dist == "bernoulli") {
    as.numeric(z > dgp$threshold)
  } else {
    qgamma(pnorm(z), shape = 1 / dgp$phi, scale = dgp$mu * dgp$phi)
  }
}

## Oracle estimators (true mu and phi): they measure what the generation
## delivers, without any estimation. Lags are defined by `wave`.
rho_oracle <- function(dist, y, mu, phi, idx, wave, lags = 1:2) {
  r <- (y - mu) / sqrt(phi * var_fun(dist, mu))
  sapply(lags, function(h) {
    mean(unlist(lapply(idx, function(ii) {
      w <- wave[ii]; m <- match(w + h, w); ok <- !is.na(m)
      r[ii[ok]] * r[ii[m[ok]]]
    })))
  })
}
phi_oracle <- function(dist, y, mu) mean((y - mu)^2 / var_fun(dist, mu))

## ---- Residuals, bootstrap and envelope ----
## Quantile residual (randomized for binary responses)
quantile_residual <- function(dist, y, mu, phi) {
  if (dist == "bernoulli") {
    a <- ifelse(y == 1, 1 - mu, 0); b <- ifelse(y == 1, 1, 1 - mu)
    u <- runif(length(y), a, b)
  } else {
    u <- pgamma(y, shape = 1 / phi, scale = mu * phi)
  }
  qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
}

## Cluster bootstrap sample
resample_clusters <- function(data, idx) {
  s <- sample.int(length(idx), length(idx), replace = TRUE)
  d <- data[unlist(idx[s], use.names = FALSE), , drop = FALSE]
  d$cid <- rep(seq_along(s), lengths(idx[s]))
  d
}

## Grouped jackknife with random assignment of clusters to J groups
## (J = NA or J >= number of clusters: complete jackknife)
jackknife_groups <- function(n_cl, J) {
  if (is.na(J) || J >= n_cl) return(seq_len(n_cl))
  sample(rep_len(seq_len(J), n_cl))
}

## BCa interval (Efron and Tibshirani, 1993, ch. 14), with ties counted as 1/2
## in z0, the proportion restricted to [1/(2B), 1 - 1/(2B)], and quantiles by
## the type 6 rule. It reduces to the percentile interval when z0 = a = 0.
BCa_interval <- function(boot, theta, jack, conf = 0.95) {
  boot <- boot[is.finite(boot)]; jack <- jack[is.finite(jack)]
  B <- length(boot)
  if (B < 10 || length(jack) < 3 || !is.finite(theta)) return(c(NA_real_, NA_real_))
  prop <- (sum(boot < theta) + 0.5 * sum(boot == theta)) / B
  z0 <- qnorm(min(max(prop, 1 / (2 * B)), 1 - 1 / (2 * B)))
  d <- mean(jack) - jack
  den <- 6 * sum(d^2)^1.5
  a <- if (den > 0) sum(d^3) / den else 0
  zc <- qnorm(1 - (1 - conf) / 2)
  zz <- z0 + c(-zc, zc)
  if (any(1 - a * zz <= 0)) return(c(NA_real_, NA_real_))
  p <- pnorm(z0 + zz / (1 - a * zz))
  unname(quantile(boot, p, type = 6))
}

## Pointwise 95% band, proportion of points outside it, and global Monte Carlo
## rank test: T = number of positions at which an ordered residual vector is
## among the k most extreme values of {observed, B simulated}; the p-value is
## the proportion of the B + 1 vectors with T at least as large as the observed.
envelope_test <- function(r_obs, sim_ord, level = 0.95) {
  sim_ord <- sim_ord[stats::complete.cases(sim_ord), , drop = FALSE]
  B <- nrow(sim_ord); al <- (1 - level) / 2
  r_ord <- sort(r_obs)
  lower <- apply(sim_ord, 2, quantile, probs = al)
  upper <- apply(sim_ord, 2, quantile, probs = 1 - al)
  med   <- apply(sim_ord, 2, quantile, probs = 0.5)
  pct_outside <- 100 * mean(r_ord < lower | r_ord > upper)

  M <- rbind(r_ord, sim_ord)
  k <- max(1, floor(al * (B + 1)))
  rk <- apply(M, 2, rank, ties.method = "random")
  extreme <- rk <= k | rk >= (B + 2 - k)
  T_stat <- rowSums(extreme)
  p_value <- mean(T_stat >= T_stat[1])
  list(lower = lower, upper = upper, med = med, r_ord = r_ord,
       pct_outside = pct_outside, p_value = p_value, B = B)
}

## ---- Parallel execution and cache ----
## Each task k has its own random number stream (L'Ecuyer-CMRG with
## nextRNGStream), so the results do not depend on the number of cores. With
## `ckpt`, each task saves its result as soon as it finishes, so an interrupted
## run resumes without loss.
EXPORT <- c("fit_gee", "coef_names", "summarize_fit", "mean_from", "var_fun",
            "target_matrices", "calibrate_norta", "prepare_dgp", "generate_y",
            "rho_oracle", "phi_oracle", "quantile_residual",
            "resample_clusters", "jackknife_groups", "BCa_interval",
            "envelope_test", "COVERAGE_METRICS")

run_parallel <- function(n_tasks, FUN, ctx, seed, n_cores = N_CORES, ckpt = NULL) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  s <- .Random.seed; seeds <- vector("list", n_tasks)
  for (k in seq_len(n_tasks)) { seeds[[k]] <- s; s <- parallel::nextRNGStream(s) }
  if (!is.null(ckpt)) dir.create(ckpt, recursive = TRUE, showWarnings = FALSE)

  task <- function(k) {
    f <- if (is.null(.ckpt)) NULL else file.path(.ckpt, sprintf("task_%05d.rds", k))
    if (!is.null(f) && file.exists(f)) return(readRDS(f))
    assign(".Random.seed", .seeds[[k]], envir = globalenv())
    r <- .FUN(k, .ctx)
    if (!is.null(f)) saveRDS(r, f)
    r
  }
  environment(task) <- globalenv()

  if (n_cores > 1 && n_tasks > 1) {
    cl <- parallel::makeCluster(min(n_cores, n_tasks))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressMessages({ library(glmtoolbox); library(pbivnorm); library(Matrix) })
      RNGkind("L'Ecuyer-CMRG"); NULL
    })
    parallel::clusterExport(cl, EXPORT, envir = globalenv())
    .FUN <- FUN; .ctx <- ctx; .seeds <- seeds; .ckpt <- ckpt
    parallel::clusterExport(cl, c(".FUN", ".ctx", ".seeds", ".ckpt"), envir = environment())
    res <- pbapply::pblapply(seq_len(n_tasks), task, cl = cl)
  } else {
    assign(".FUN", FUN, envir = globalenv()); assign(".ctx", ctx, envir = globalenv())
    assign(".seeds", seeds, envir = globalenv()); assign(".ckpt", ckpt, envir = globalenv())
    on.exit(rm(list = c(".FUN", ".ctx", ".seeds", ".ckpt"), envir = globalenv()), add = TRUE)
    res <- pbapply::pblapply(seq_len(n_tasks), task)
  }
  res
}

## Cache: each step is stored under a key built from CACHE_TAG, the Monte Carlo
## sizes and the fitted model (formula, correlation structure and estimates).
hash_key <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  saveRDS(x, f, version = 3); substr(unname(tools::md5sum(f)), 1, 8)
}
model_id <- function(obj) list(formula = paste(deparse(obj$spec$formula), collapse = ""),
                               corstr = obj$spec$corstr, n = nrow(obj$data),
                               beta = signif(obj$beta, 8), phi = signif(obj$phi, 8),
                               rho = signif(obj$rho_lags, 8))
step_model <- function(name) {
  pref <- sub("_.*", "", name)
  if (pref %in% c("ber", "gam") && exists(pref, envir = globalenv())) model_id(get(pref, envir = globalenv())) else NULL
}
step_key <- function(name, params, model = step_model(name)) hash_key(list(CACHE_TAG, params, model))
checkpoint_dir <- function(name, params, model = step_model(name))
  file.path(DIR_CACHE, sprintf("ckpt_%s_%s", name, step_key(name, params, model)))

with_cache <- function(name, params, expr, model = step_model(name)) {
  f <- file.path(DIR_CACHE, sprintf("%s_%s.rds", name, step_key(name, params, model)))
  if (USE_CACHE && file.exists(f)) {
    cat(sprintf("[cache] %s\n", basename(f)))
    return(readRDS(f))
  }
  t0 <- Sys.time()
  val <- force(expr)
  attr(val, "minutes") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(val, f)
  cat(sprintf("[%s] %.1f min\n", name, attr(val, "minutes")))
  val
}

## ---- Simulation tasks ----
task_boot_app <- function(k, ctx) {
  d <- resample_clusters(ctx$data, ctx$idx)
  summarize_fit(fit_gee(ctx$spec, d), ctx$p)
}

task_jack_app <- function(k, ctx) {
  out <- unlist(ctx$idx[ctx$groups == k], use.names = FALSE)
  summarize_fit(fit_gee(ctx$spec, ctx$data[-out, , drop = FALSE]), ctx$p)
}

## One envelope replicate (parametric bootstrap with refit)
envelope_replicate <- function(spec, data, dgp) {
  y <- generate_y(dgp)
  d <- data; d[[spec$resp]] <- y
  fit <- fit_gee(spec, d)
  if (is.null(fit)) return(NULL)
  mu <- mean_from(spec, d, as.numeric(coef(fit)))
  sort(quantile_residual(spec$dist, y, mu, fit$phi))
}
EXPORT <- c(EXPORT, "envelope_replicate")

task_envelope <- function(k, ctx) envelope_replicate(ctx$spec, ctx$data, ctx$dgp)

## Residuals (in the original order) of a sample simulated from the fitted model and refitted
residual_replicate <- function(spec, data, dgp) {
  y <- generate_y(dgp)
  d <- data; d[[spec$resp]] <- y
  fit <- fit_gee(spec, d)
  if (is.null(fit)) return(NULL)
  quantile_residual(spec$dist, y, mean_from(spec, d, as.numeric(coef(fit))), fit$phi)
}
EXPORT <- c(EXPORT, "residual_replicate")
task_residuals <- function(k, ctx) residual_replicate(ctx$spec, ctx$data, ctx$dgp)

## Type I error / power: generate from the true model, fit spec_fit, build a
## complete envelope with refitting and return the proportion outside the band
## and the p-value of the global test.
task_type1 <- function(k, ctx) {
  spec <- ctx$spec_fit
  y <- generate_y(ctx$dgp_true)
  d <- ctx$data; d[[spec$resp]] <- y
  fit <- fit_gee(spec, d)
  if (is.null(fit)) return(c(pct_outside = NA, p_value = NA, failed_inner = NA))
  mu <- mean_from(spec, d, as.numeric(coef(fit)))
  r_obs <- quantile_residual(spec$dist, y, mu, fit$phi)
  phi_dgp <- if (spec$dist == "bernoulli") 1 else fit$phi
  dgp_fit <- prepare_dgp(spec$dist, mu, phi_dgp, as.numeric(fit$corr[1, -1]), ctx$idx, ctx$wave)
  sims <- lapply(seq_len(ctx$B_inner), function(b) envelope_replicate(spec, d, dgp_fit))
  failed <- sum(vapply(sims, is.null, logical(1)))
  sims <- do.call(rbind, sims[!vapply(sims, is.null, logical(1))])
  if (is.null(sims) || nrow(sims) < 5) return(c(pct_outside = NA, p_value = NA, failed_inner = failed))
  et <- envelope_test(r_obs, sims)
  c(pct_outside = et$pct_outside, p_value = et$p_value, failed_inner = failed)
}

## Simulation study: estimates, oracle estimators and Wald coverage
task_sim <- function(k, ctx) {
  y <- generate_y(ctx$dgp)
  d <- ctx$data; d[[ctx$spec$resp]] <- y
  orc <- c(rho_oracle(ctx$spec$dist, y, ctx$dgp$mu, ctx$truth$phi, ctx$idx, ctx$wave),
           phi_oracle(ctx$spec$dist, y, ctx$dgp$mu))
  fit <- fit_gee(ctx$spec, d)
  est <- summarize_fit(fit, ctx$p)
  cov <- rep(NA_real_, 3 * ctx$p)
  if (!is.null(fit)) {
    b <- as.numeric(coef(fit)); z <- qnorm(0.975)
    types <- c("robust", "model", "bias-corrected")
    for (t in seq_along(types)) {
      se <- tryCatch(sqrt(diag(vcov(fit, type = types[t]))), error = function(e) rep(NA_real_, ctx$p))
      cov[(t - 1) * ctx$p + seq_len(ctx$p)] <- as.numeric(abs(b - ctx$truth$beta) <= z * se)
    }
  }
  c(est, orc, cov)
}

## Nested bootstrap for the coverage of the BCa, percentile and basic intervals
## of all parameters, with the robust Wald interval for beta in the same
## samples. Only converged fits are used.
COVERAGE_METRICS <- c("hat", "bca", "pct", "basic", "wald", "width_bca", "width_pct", "width_wald",
                      "boot_bias", "boot_se", "n_boot", "jack_ok")
task_coverage <- function(k, ctx) {
  names_par <- c(make.names(ctx$names), "phi", "rho1", "rho2"); p <- length(ctx$names); q <- p + 3
  out_na <- setNames(rep(NA_real_, q * length(COVERAGE_METRICS) + 1),
                     c(paste(rep(names_par, each = length(COVERAGE_METRICS)), COVERAGE_METRICS, sep = "_"), "converged"))
  out_na["converged"] <- 0
  y <- generate_y(ctx$dgp)
  d <- ctx$data; d[[ctx$spec$resp]] <- y
  fit <- fit_gee(ctx$spec, d)
  if (is.null(fit) || !isTRUE(fit$converged)) return(out_na)
  pars <- function(f) c(as.numeric(coef(f)), f$phi, f$corr[1, 2], f$corr[1, 3])
  ok <- function(f) !is.null(f) && isTRUE(f$converged)
  th <- pars(fit)
  se_rob <- tryCatch(sqrt(diag(vcov(fit, type = "robust"))), error = function(e) rep(NA_real_, p))

  boot <- matrix(NA_real_, ctx$R_boot, q)
  for (r in seq_len(ctx$R_boot)) {
    f <- fit_gee(ctx$spec, resample_clusters(d, ctx$idx))
    if (ok(f)) boot[r, ] <- pars(f)
  }
  g <- jackknife_groups(length(ctx$idx), ctx$J_cov)
  jack <- matrix(NA_real_, max(g), q)
  for (j in seq_len(max(g))) {
    f <- fit_gee(ctx$spec, d[-unlist(ctx$idx[g == j], use.names = FALSE), , drop = FALSE])
    if (ok(f)) jack[j, ] <- pars(f)
  }

  truth <- c(ctx$truth$beta, ctx$truth$phi, ctx$truth$rho1, ctx$truth$rho2)
  z <- qnorm(0.975)
  out <- out_na; out["converged"] <- 1
  for (j in seq_len(q)) {
    bj <- boot[is.finite(boot[, j]), j]
    ci_bca <- BCa_interval(bj, th[j], jack[, j])
    ci_pct <- if (length(bj) >= 10) unname(quantile(bj, c(0.025, 0.975), type = 6)) else c(NA, NA)
    ci_bas <- 2 * th[j] - rev(ci_pct)
    ci_wal <- if (j <= p) th[j] + c(-1, 1) * z * se_rob[j] else c(NA, NA)
    covers <- function(ci) if (any(is.na(ci))) NA_real_ else as.numeric(truth[j] >= ci[1] & truth[j] <= ci[2])
    out[paste(names_par[j], COVERAGE_METRICS, sep = "_")] <-
      c(th[j], covers(ci_bca), covers(ci_pct), covers(ci_bas), covers(ci_wal),
        diff(ci_bca), diff(ci_pct), diff(ci_wal),
        mean(bj) - th[j], sd(bj), length(bj), sum(is.finite(jack[, j])))
  }
  out
}

## ---- Tables and figures ----
## Absolute bias with MCSE and 95% interval; relative bias when |theta| > 0.05
## and always for phi and rho
bias_table <- function(est, truth, names_par) {
  est <- as.matrix(est)
  n  <- colSums(is.finite(est))
  m  <- colMeans(est, na.rm = TRUE)
  sdv <- apply(est, 2, sd, na.rm = TRUE)
  mcse <- sdv / sqrt(n); bias <- m - truth
  data.frame(Parameter = names_par, True = round(truth, 5), Mean = round(m, 5),
             Bias = round(bias, 5), MCSE = round(mcse, 5),
             CI95_bias = sprintf("[%.5f; %.5f]", bias - 1.96 * mcse, bias + 1.96 * mcse),
             Significant = ifelse(abs(bias) > 1.96 * mcse, "yes", "no"),
             Rel_bias_pct = ifelse(abs(truth) > 0.05 | names_par %in% c("phi", "rho1", "rho2"), round(100 * bias / abs(truth), 2), NA),
             Empirical_SD = round(sdv, 5), n = n, row.names = NULL)
}

coverage_table <- function(M, row_names) {
  n <- colSums(is.finite(M)); p <- colMeans(M, na.rm = TRUE)
  data.frame(Item = row_names, Coverage_pct = round(100 * p, 1),
             MCSE_pp = round(100 * sqrt(p * (1 - p) / n), 1), n = n, row.names = NULL)
}

plot_bootstrap <- function(boot, jack, obs, labels, nrow, ncol) {
  k <- ncol(boot)
  layout(matrix(c(seq_len(k), rep(k + 1, nrow * ncol - k)), nrow, ncol, byrow = TRUE)); on.exit(layout(1), add = TRUE)
  op <- par(mar = c(4, 4, 2, 1) + 0.1); on.exit(par(op))
  for (i in seq_len(k)) {
    x <- boot[is.finite(boot[, i]), i]
    h <- hist(x, plot = FALSE); dk <- density(x)
    ymax <- max(h$density, dk$y) * 1.4
    hist(x, freq = FALSE, col = "gray", border = "black", main = "",
         xlab = labels[[i]], ylab = "Density", ylim = c(0, ymax))
    lines(dk, col = "blue", lwd = 2)
    m_x <- mean(x); s_x <- sd(x)
    curve(dnorm(x, m_x, s_x), add = TRUE, col = "magenta", lwd = 2, lty = 2)
    abline(v = obs[i], col = "black", lwd = 2, lty = 2)
    abline(v = median(x), col = "green", lwd = 2)
    ci <- BCa_interval(x, obs[i], jack[, i])
    if (!anyNA(ci)) segments(ci[1], ymax * 0.05, ci[2], ymax * 0.05, lwd = 3)
    box()
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", c("Kernel Density", "Adjusted Normal", "Observed Value", "Bootstrap Median", "CI 95% BCa"),
         col = c("blue", "magenta", "black", "green", "black"),
         lty = c(1, 2, 2, 1, 1), lwd = c(2, 2, 2, 2, 4), bty = "n")
}

plot_envelope <- function(et, ylab) {
  n <- length(et$r_ord)
  df <- data.frame(theo = qnorm(ppoints(n)), r = et$r_ord, lower = et$lower, upper = et$upper, med = et$med)
  df$col <- ifelse(df$r >= df$lower & df$r <= df$upper, "Inside", "Outside")
  ggplot(df, aes(theo, r)) +
    geom_point(aes(color = col), size = 0.8) +
    geom_line(aes(y = med), linetype = "dashed") +
    geom_line(aes(y = lower)) + geom_line(aes(y = upper)) +
    scale_color_manual(values = c(Inside = "black", Outside = "red")) +
    guides(color = "none") +
    labs(x = "Theoretical Quantile N(0,1)", y = ylab) +
    theme_test()
}

## Worm plots in groups of equal size defined by the ordered values of x
wp_groups <- function(r, x, groups = 8, ylim_all = 0.8, label = quote(hat(mu))) {
  g <- ceiling(groups * rank(x, ties.method = "first") / length(x))
  op <- par(mfrow = c(ceiling(groups / 2), 2), mar = c(4, 4, 2, 1)); on.exit(par(op))
  for (k in seq_len(groups)) {
    sel <- g == k
    ttl <- bquote(.(label) %in% group("[", list(.(formatC(min(x[sel]), format = "f", digits = 3)), .(formatC(max(x[sel]), format = "f", digits = 3))), "]"))
    if (sum(sel) >= 5) {
      wp(resid = r[sel], ylim.all = ylim_all)
      title(main = ttl)
    } else {
      plot.new(); text(0.5, 0.5, paste(ttl, "with few observations"))
    }
  }
}

## Worm plots by groups of values of a variable (e.g. year); `groups` is a
## named list of values; if NULL, each distinct value is a group.
wp_by_value <- function(r, v, ylim_all = 0.8, label = "Group", groups = NULL) {
  if (is.null(groups)) { u <- sort(unique(v)); groups <- setNames(as.list(u), u) }
  op <- par(mfrow = c(ceiling(length(groups) / 2), 2), mar = c(4, 4, 2, 1)); on.exit(par(op))
  for (nm in names(groups)) {
    sel <- v %in% groups[[nm]]
    wp(resid = r[sel], ylim.all = ylim_all)
    title(main = paste(label, nm))
  }
  if (length(groups) %% 2 == 1) plot.new()
}

plot_simulation <- function(est, truth, labels, nrow, ncol) {
  k <- ncol(est)
  op <- par(mfrow = c(nrow, ncol), mar = c(4.6, 4.4, 1.5, 1), cex.lab = 1.5, cex.axis = 1.25); on.exit(par(op))
  for (i in seq_len(k)) {
    x <- est[is.finite(est[, i]), i]
    h <- hist(x, plot = FALSE); dk <- density(x)
    ymax <- max(h$density, dk$y) * 1.3
    hist(x, freq = FALSE, col = "#E07A5F", border = "white", main = "",
         xlab = labels[[i]], ylab = "Density", ylim = c(0, ymax))
    lines(dk, col = "#3D405B", lwd = 2.5)
    abline(v = truth[i], lwd = 2.5, lty = 2)
  }
  plot.new()
  legend("center", c("Kernel Density", "True Value"), col = c("#3D405B", "black"),
         lty = c(1, 2), lwd = 2.5, bty = "n", cex = 1.3)
}

save_pdf <- function(file, drawing, width = 8, height = 8) {
  grDevices::cairo_pdf(file.path(DIR_FIG, file), width = width, height = height)
  on.exit(dev.off())
  drawing()
}

## Heat plots of the quantile residuals:
## (a) quantile residual against the fitted value and (b) worm plot, as density plots of the
## residuals of R randomizations (R = 1 for continuous responses), with smooth trend curves.
## The worm plot bands assume independent residuals with known parameters.
heat_plot <- function(dist, y, mu, phi, R = 1, xlab = "Fitted value", ylim_wp = NULL) {
  n <- length(y)
  rq_matrix <- matrix(replicate(R, quantile_residual(dist, y, mu, phi)), nrow = n)
  rq <- c(rq_matrix)
  rq_theor <- qnorm(ppoints(n))
  rq_cente <- apply(rq_matrix, 2, sort) - rq_theor
  se <- sqrt(ppoints(n) * (1 - ppoints(n)) / n) / dnorm(rq_theor)
  conf_upp <- qnorm(0.975) * se
  conf_low <- qnorm(0.025) * se
  if (is.null(ylim_wp)) {
    m <- max(abs(range(rq_cente[6:(n - 5), ])), conf_upp[round(0.01 * n)])
    ylim_wp <- c(-1, 1) * ceiling(10 * 1.15 * m) / 10
  }
  ylim_r <- c(-1, 1) * max(3.5, ceiling(max(abs(rq))))
  ramp <- colorRampPalette(c("white", viridisLite::viridis(11, option = "inferno", direction = -1)))
  trend_line <- function(x, y) {
    lines(x, y, col = "white", lwd = 4)
    lines(x, y, col = "black", lwd = 1.8)
  }
  op <- par(mfrow = c(1, 2), mar = c(3.9, 4.0, 1.6, 0.6), mgp = c(2.6, 0.6, 0), las = 1,
            cex.lab = 1, cex.axis = 0.85, tcl = -0.3)
  on.exit(par(op))

  ## (a) quantile residual against the fitted value, with the smooth mean residual
  smoothScatter(rq ~ rep(mu, R), xlab = xlab, ylab = "Quantile residual", colramp = ramp,
                nbin = 200, nrpoints = 50, col = "grey30", pch = 19, cex = 0.3, ylim = ylim_r,
                useRaster = TRUE)
  abline(h = 0, lwd = 1.4, lty = 2)
  abline(h = c(-2, 2), lwd = 1, lty = 3, col = "grey40")
  mu_rep <- rep(mu, R)
  grid <- seq(min(mu), max(mu), length.out = 200)
  trend_line(grid, predict(mgcv::gam(rq ~ s(mu_rep, k = 10)), data.frame(mu_rep = grid)))
  u <- par("usr"); w <- diff(u[1:2]); h <- diff(u[3:4])
  x0 <- u[2] - 0.30 * w; x1 <- u[2] - 0.04 * w; y0 <- u[4] - 0.09 * h; y1 <- u[4] - 0.05 * h
  rasterImage(as.raster(matrix(ramp(100), nrow = 1)), x0, y0, x1, y1, interpolate = TRUE)
  rect(x0, y0, x1, y1, border = "grey40", lwd = 0.5)
  text(c(x0, x1), y0, c("low", "high"), pos = 1, cex = 0.6, offset = 0.25)
  text((x0 + x1) / 2, y1, "density", pos = 3, cex = 0.6, offset = 0.2)
  mtext("(a)", side = 3, adj = 0, line = 0.3, font = 2)

  ## (b) worm plot, with a smooth curve of the (median, if R > 1) worm
  smoothScatter(rq_cente ~ rep(rq_theor, R), xlab = "Theoretical quantile",
                ylab = "Centered quantile residual", colramp = ramp, nbin = 200, nrpoints = 0,
                ylim = ylim_wp, useRaster = TRUE)
  worm <- if (R > 1) apply(rq_cente, 1, median) else rq_cente[, 1]
  trend_line(rq_theor, predict(mgcv::gam(worm ~ s(rq_theor, k = 10))))
  abline(h = 0, lwd = 1.4, lty = 2)
  lines(rq_theor, conf_upp, lwd = 1.4, lty = 2)
  lines(rq_theor, conf_low, lwd = 1.4, lty = 2)
  mtext("(b)", side = 3, adj = 0, line = 0.3, font = 2)
  invisible(c(pct_outside_wp = mean(rq_cente > conf_upp | rq_cente < conf_low)))
}

## Moment bucket plot (De Bastiani et al., 2022) adapted to GEE: transformed moment skewness and
## kurtosis of the quantile residuals (continuous responses), with
##  - a cluster bootstrap cloud (clusters resampled) for the uncertainty of the observed point;
##  - a reference cloud from samples simulated from the fitted model (Gaussian copula) and
##    refitted, replacing the Jarque-Bera region, which assumes independent residuals;
##  - a Monte Carlo p-value from the Mahalanobis distances of the B + 1 points.
## r_obs: observed residuals; R_sim: N x B matrix of simulated residuals (original order);
## idx: list of the rows of each cluster; rows: subset of observations.
bucket_stats <- function(r_obs, R_sim, idx, rows = seq_along(r_obs), n_boot = 499, seed = 2028) {
  tsk <- function(r) { m <- gamlss.dist::momentSK(r); c(x = m$trans.mom.skew, y = m$trans.mom.kurt) }
  obs <- tsk(r_obs[rows])
  sim <- t(apply(R_sim[rows, , drop = FALSE], 2, tsk))
  cl <- lapply(idx, intersect, rows); cl <- cl[lengths(cl) > 0]
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  boot <- t(replicate(n_boot, tsk(r_obs[unlist(cl[sample.int(length(cl), replace = TRUE)])])))
  ## Monte Carlo test: Mahalanobis distances of the observed and simulated points
  all <- rbind(obs, sim)
  d <- mahalanobis(all, colMeans(all), cov(all))
  p <- mean(d >= d[1])
  ## 95% region of the simulated points: ellipse with radius given by their 95% distance quantile
  m <- colMeans(sim); S <- cov(sim)
  q95 <- quantile(mahalanobis(sim, m, S), 0.95, names = FALSE)
  th <- seq(0, 2 * pi, length.out = 200)
  ell <- t(m + t(chol(S)) %*% rbind(cos(th), sin(th)) * sqrt(q95))
  ## Jarque-Bera 95% region for independent residuals, for comparison
  n <- length(rows)
  u <- seq(-0.99, 0.99, length.out = 400); g1 <- u / (1 - abs(u))
  g2 <- (24 / n) * (qchisq(0.95, 2) - (n / 6) * g1^2)
  ok <- g2 >= 0; g1 <- g1[ok]; g2 <- sqrt(g2[ok])
  jb <- data.frame(x = c(g1, rev(g1)), y = c(g2, -rev(g2)))
  jb <- data.frame(x = jb$x / (1 + abs(jb$x)), y = jb$y / (1 + abs(jb$y)))
  list(obs = obs, sim = as.data.frame(sim), boot = as.data.frame(boot), p = p,
       ell = data.frame(x = ell[, 1], y = ell[, 2]), jb = jb, n = n,
       moments = gamlss.dist::momentSK(r_obs[rows])[c("mom.skew", "excess.mom.kurt")])
}

## ggplot of a bucket_stats object on the gray background of gamlss.ggplots
bucket_plot <- function(b, title = NULL, lim = NULL, base_size = 8) {
  g <- gamlss.ggplots::moment_gray_both() +
    ggplot2::geom_polygon(data = b$jb, ggplot2::aes(x, y), fill = "#56B4E9", alpha = 0.35, inherit.aes = FALSE) +
    ggplot2::geom_point(data = b$sim, ggplot2::aes(x, y), colour = "grey40", size = 0.35, alpha = 0.6, inherit.aes = FALSE) +
    ggplot2::geom_path(data = b$ell, ggplot2::aes(x, y), colour = "black", linetype = 2, linewidth = 0.4, inherit.aes = FALSE) +
    ggplot2::geom_point(data = b$boot, ggplot2::aes(x, y), colour = "#D55E00", size = 0.35, alpha = 0.45, inherit.aes = FALSE) +
    ggplot2::annotate("point", x = b$obs["x"], y = b$obs["y"], shape = 21, size = 2.2, fill = "black", colour = "white") +
    ggplot2::labs(x = "transformed moment skewness", y = "transformed moment excess kurtosis", title = title) +
    ggplot2::theme(legend.position = "none", text = ggplot2::element_text(size = base_size),
                   axis.text = ggplot2::element_text(size = base_size - 1),
                   plot.title = ggplot2::element_text(face = "bold", size = base_size + 1))
  if (!is.null(lim)) g <- g + ggplot2::coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim))
  g
}


## ---- Model selection ----
structures <- c(Independence = "Independence", Exchangeable = "Exchangeable",
                "AR(1)" = "AR-M-dependent(1)", "AR(2)" = "AR-M-dependent(2)", "AR(3)" = "AR-M-dependent(3)")
criteria <- c("CIC", "QIC", "GHYC", "RJC", "AGPC", "SGPC")
criteria_table <- function(form, family, data) {
  t(sapply(structures, function(cs) {
    m <- tryCatch(if (cs == "Independence") glmtoolbox::glmgee(form, id = cid, family = family, corstr = cs, data = data)
                  else glmtoolbox::glmgee(form, id = cid, waves = wave, family = family, corstr = cs, data = data),
                  error = function(e) NULL)
    if (is.null(m)) return(setNames(rep(NA_real_, length(criteria)), criteria))
    sapply(criteria, function(g) tryCatch(as.numeric(do.call(g, list(m, verbose = FALSE)))[1], error = function(e) NA_real_))
  }))
}

nested_tests <- function(formulas, pairs, family, data) {
  fits <- lapply(formulas, function(f) glmtoolbox::glmgee(f, id = cid, waves = wave, family = family,
                                                           corstr = "AR-M-dependent(2)", data = data))
  do.call(rbind, lapply(pairs, function(pr) {
    out <- c(Comparison = paste(pr, collapse = " vs "))
    for (tp in c("wald", "score")) {
      r <- as.data.frame(anova(fits[[pr[1]]], fits[[pr[2]]], test = tp, verbose = FALSE))
      out[paste0(tp, c("_chi2", "_df", "_p"))] <- c(round(r[nrow(r), 1], 2), r[nrow(r), 2], signif(r[nrow(r), ncol(r)], 3))
    }
    out
  }))
}

## ---- Summaries of the simulation studies ----
type1_summary <- function(lst) do.call(rbind, lapply(lst, function(d) {
  rej <- d$p_value < 0.05; n <- sum(is.finite(rej))
  data.frame(Scenario = d$Scenario[1], n = n,
             Rejection_pct = round(100 * mean(rej, na.rm = TRUE), 1),
             MCSE_pp = round(100 * sqrt(mean(rej, na.rm = TRUE) * (1 - mean(rej, na.rm = TRUE)) / n), 1),
             Mean_pct_outside = round(mean(d$pct_outside, na.rm = TRUE), 2),
             Failed_inner = sum(d$failed_inner, na.rm = TRUE))
}))

size_omitted <- function(term, X, scale) sd(resid(lm(term ~ X))) / scale

## Coverage of all parameters and comparison of the bootstrap bias and
## standard error with the true bias and standard deviation (simulation study)
coverage_summary <- function(C, truth_vec, est_sim, names_beta) {
  C <- C[C[, "converged"] == 1, , drop = FALSE]
  nm <- c(make.names(names_beta), "phi", "rho1", "rho2"); lab <- c(names_beta, "phi", "rho1", "rho2")
  p <- length(names_beta)
  cov <- do.call(rbind, lapply(seq_along(nm), function(j) {
    v <- nm[j]; met <- c("bca", "pct", "basic", if (j <= p) "wald")
    tc <- coverage_table(C[, paste0(v, "_", met), drop = FALSE], c("BCa", "Percentile", "Basic", if (j <= p) "Robust Wald"))
    width <- c(mean(C[, paste0(v, "_width_bca")], na.rm = TRUE), rep(mean(C[, paste0(v, "_width_pct")], na.rm = TRUE), 2),
               if (j <= p) mean(C[, paste0(v, "_width_wald")], na.rm = TRUE))
    data.frame(Parameter = lab[j], Method = tc$Item, Coverage_pct = tc$Coverage_pct, MCSE_pp = tc$MCSE_pp,
               Mean_width = signif(width, 4), n = tc$n)
  }))
  bias <- do.call(rbind, lapply(seq_along(nm), function(j) {
    v <- nm[j]; true_bias <- est_sim[, j] - truth_vec[j]
    data.frame(Parameter = lab[j],
               Boot_bias = signif(mean(C[, paste0(v, "_boot_bias")], na.rm = TRUE), 4),
               True_bias = signif(mean(true_bias, na.rm = TRUE), 4),
               Boot_SE = signif(mean(C[, paste0(v, "_boot_se")], na.rm = TRUE), 4),
               True_SD = signif(sd(est_sim[, j], na.rm = TRUE), 4),
               Valid_resamples = round(mean(C[, paste0(v, "_n_boot")], na.rm = TRUE), 1))
  }))
  list(coverage = cov, bias = bias, n_valid = nrow(C))
}
## Paired comparison of bootstrap and Wald intervals in the same samples
paired_comparison <- function(C, names_beta) {
  C <- C[C[, "converged"] == 1, , drop = FALSE]
  do.call(rbind, lapply(seq_along(names_beta), function(i) {
    v <- make.names(names_beta[i]); w <- C[, paste0(v, "_wald")]
    do.call(rbind, lapply(c("bca", "pct", "basic"), function(m) {
      x <- C[, paste0(v, "_", m)]; ok <- is.finite(x) & is.finite(w); d <- x[ok] - w[ok]
      b <- sum(x[ok] == 1 & w[ok] == 0); cc <- sum(x[ok] == 0 & w[ok] == 1)
      width <- C[ok, paste0(v, if (m == "bca") "_width_bca" else "_width_pct")]
      data.frame(Parameter = names_beta[i], Method = m, Diff_coverage_pp = round(100 * mean(d), 2),
                 Paired_MCSE_pp = round(100 * sd(d) / sqrt(length(d)), 2), Only_bootstrap = b, Only_Wald = cc,
                 p_McNemar = if (b + cc > 0) signif(binom.test(b, b + cc)$p.value, 3) else NA,
                 Width_ratio = round(mean(width) / mean(C[ok, paste0(v, "_width_wald")]), 3))
    }))
  }))
}
