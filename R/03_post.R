# =============================================================================
# 03_post.R
# 後処理（段階間の相関、得点価値の補正、中立環境、得点価値、リーグ平均が 100 の尺度）
# =============================================================================

# ---- 段階間の相関 --------------------------------------------------------------
# 各段階の打者効果の事後分布から、「データだけが示す値 yy」と「その精度 tau」を逆算する。
# 三振の段階は t 分布なので、標準偏差を正規分布の近似に直す（sqrt(12 / 10) 倍）。
corr_inputs <- function(season, ver) {
  specs <- prod_spec(season, ver == "pac")
  df_s <- make_season_prod(season)
  n_bat <- max(df_s$batter_idx)
  b <- map(set_names(corr_stages), \(s) {
    f <- readRDS(prod_path(season, ver, s))
    d <- unclass(f$draws(c("sigma[1]", paste0("z[", seq_len(n_bat), "]")),
                         format = "draws_matrix"))
    dev <- d[, 1] * d[, -1]
    list(m = colMeans(dev), v = apply(dev, 2, var),
         sig = mean(d[, 1]) * if (isTRUE(specs[[s]]$use_t == 2L)) sqrt(12 / 10) else 1)
  })
  m <- sapply(b, \(x) x$m)
  v <- sapply(b, \(x) x$v)
  sig <- map_dbl(b, \(x) x$sig)
  tau <- pmax(1 / v - matrix(1 / sig^2, n_bat, length(corr_stages), byrow = TRUE), 0)
  list(season = season, sig = sig, tau = tau, yy = ifelse(tau > 0, (m / v) / tau, 0),
       use = which(tabulate(df_s$batter_idx, n_bat) >= 200))
}

to_omega4 <- function(par) {
  l <- diag(4)
  l[lower.tri(l)] <- par
  l <- l / sqrt(rowSums(l^2))
  om <- l %*% t(l)
  dimnames(om) <- list(corr_stages, corr_stages)
  om
}

# 相関行列の周辺尤度（符号を反転）
nll_one <- function(par, inp) {
  s_mat <- diag(inp$sig) %*% to_omega4(par) %*% diag(inp$sig)
  s <- 0
  for (i in inp$use) {
    k <- inp$tau[i, ] > 0
    if (!any(k)) next
    vv <- s_mat[k, k, drop = FALSE] + diag(1 / inp$tau[i, k], sum(k))
    s <- s + 0.5 * (as.numeric(determinant(vv)$modulus) +
                      sum(inp$yy[i, k] * solve(vv, inp$yy[i, k])))
  }
  s
}

# シーズンごとの相関行列と、全シーズン共通の相関行列（指標に使うのは共通の方）
estimate_omega <- function(ver, years = seasons) {
  inputs <- map(years, \(s) { message(ver, " ", s); corr_inputs(s, ver) })
  each <- map(inputs, \(inp) to_omega4(optim(rep(0, 6), nll_one, inp = inp, method = "BFGS")$par))
  pooled <- to_omega4(optim(rep(0, 6), \(par) sum(map_dbl(inputs, \(inp) nll_one(par, inp))),
                            method = "BFGS")$par)
  list(each = set_names(each, years), pooled = pooled)
}

# ---- 得点価値の補正 --------------------------------------------------------------
# 目的変数 resid は、打席の得点価値から「同じイベント・同じ塁アウト状況の平均」を引いた値
make_value_data <- function(df_s, rv_s) {
  df_s |>
    left_join(select(rv_s, game_id, at_bat_number, rv_pa), by = c("game_id", "at_bat_number")) |>
    mutate(
      state = paste(bases, outs),
      on1 = substr(bases, 1, 1) != "0",
      on3 = substr(bases, 3, 3) != "0",
      sg = factor(case_when(
        bases == "000" ~ "empty",
        outs == 2 ~ "on_2out",
        on1 ~ "dp_opp",
        on3 ~ "r3_lt2",
        .default = "other_lt2"
      ), levels = c("empty", "dp_opp", "r3_lt2", "other_lt2", "on_2out")),
      spray_field = atan2(hit_coord_x - 125.42, 198.27 - hit_coord_y) * 180 / pi,
      hp1b_c = avg_hp_to_1b - mean(avg_hp_to_1b)
    ) |>
    mutate(rv_state = mean(rv_pa), n_cell = n(), .by = c(event, state)) |>
    mutate(rv_event = mean(rv_pa), .by = event) |>
    mutate(
      resid = rv_pa - if_else(n_cell >= 30, rv_state, rv_event),
      use = !is.na(launch_speed) & !is.na(launch_angle) & !is.na(spray_field) &
        abs(spray_field) <= 55 & !(event == "out" & events == "fielders_choice")
    )
}

fit_value <- function(d) bam(val_formula, data = d, discrete = TRUE, nthreads = 4)

# 各打球を、5つの塁アウト状況カテゴリに置き、リーグ全体の割合で重み付けして平均する
value_star <- function(g, d, pi_g) {
  lv <- levels(d$sg)
  m <- sapply(lv, \(l) as.numeric(predict(g, newdata = mutate(d, sg = factor(l, levels = lv)))))
  as.numeric(m %*% pi_g[lv])
}

# 打者 × イベント（アウト、単打、二塁打）ごとの値。件数に応じて 0 に向けて縮小する。
# tau は打者間の差の標準偏差（シーズンごとにデータから決まる）。
value_season <- function(season) {
  df_s <- make_season_prod(season)
  vd <- make_value_data(df_s, filter(rv_all, year == season))
  pi_g <- prop.table(table(vd$sg))
  map(c(out = "out", x1b = "x1b", x2b = "x2b"), \(ev) {
    d <- vd |> filter(event == ev, use)
    d$vs <- value_star(fit_value(d), d, pi_g)
    by <- d |> summarise(n = n(), delta = mean(vs), s2 = var(vs), .by = batter_idx)
    big <- by$n >= 30
    tau2 <- max(var(by$delta[big]) - mean(by$s2[big] / by$n[big]), 0)
    by |>
      mutate(s2 = coalesce(s2, mean(s2, na.rm = TRUE)),
             delta_hat = delta * tau2 / (tau2 + s2 / n),
             event = ev, tau = sqrt(tau2))
  }) |>
    list_rbind() |>
    mutate(season = season)
}

# ---- 1シーズン・1つの版の指標 ----------------------------------------------------
# 1. 段階ごとの事後サンプルを取り出す
# 2. 四球・三振・本塁打・長打を、6シーズン共通の相関 omega で縮小し直す（縮小し直した後の分布から引き直す）
# 3. Statcast データを使う段階を、実際の件数に合わせて補正する（係数は推定誤差つき）
# 4. 中立環境（リーグの実際の打席 n_ctx 個）で9カテゴリの確率を出す
# 5. 得点価値を掛け、得点価値の補正を加え、リーグ平均が 100 の尺度に直す
# rv_all、val_all はグローバルのものを使う。
post_season <- function(season, ver, omega, n_ctx = 3000, n_draw = 400) {
  out_path <- file.path(Sys.getenv("DRC_DATA_DIR"),
                        sprintf("result_%d_%s_v3.rds", season, ver))
  if (file.exists(out_path)) return(readRDS(out_path))
  set.seed(20261009)
  df_s <- make_season_prod(season)
  rv_s <- rv_all |> filter(year == season)
  r_pa <- rv_s$r_pa[1]
  specs <- prod_spec(season, ver == "pac")
  n_bat <- max(df_s$batter_idx)
  attr_cols <- c("x_logpa", "x_pitcher_bat")
  ball_stages <- c("hr", "hit", "xbh")
  k4 <- length(corr_stages)

  ctx <- df_s |>
    slice_sample(n = n_ctx) |>
    select(pitcher_idx, pitch_hand, venue_name, park_idx, fld_team_idx, x_temp, x_pos_pitcher)
  ps_levels <- levels(factor(paste(df_s$venue_name, df_s$bat_side)))
  side_prob <- df_s |>
    summarise(q = mean(bat_side == pitch_hand), .by = c(batter_idx, pitch_hand)) |>
    pivot_wider(names_from = pitch_hand, values_from = q, names_prefix = "q_same_vs_")
  batter_attr <- df_s |>
    summarise(batter_id = first(batter_id), batter_name = first(batter_name), n_pa = n(),
              share_l = mean(bat_side == "L"), x_logpa = first(x_logpa),
              x_pitcher_bat = first(x_pitcher_bat), .by = batter_idx) |>
    arrange(batter_idx) |>
    left_join(side_prob, by = "batter_idx") |>
    mutate(q_same_vs_L = coalesce(q_same_vs_L, as.numeric(share_l >= 0.5)),
           q_same_vs_R = coalesce(q_same_vs_R, as.numeric(share_l < 0.5)))

  get_stage <- function(st) {
    sp <- specs[[st]]
    f <- readRDS(prod_path(season, ver, st))
    n_level <- map_int(sp$g, \(g) max(df_s[[g]]))
    d <- unclass(f$draws(c("alpha", "beta", "sigma", "z"), format = "draws_matrix"))
    d <- d[sample(nrow(d), n_draw), ]
    beta <- d[, paste0("beta[", seq_along(sp$x), "]"), drop = FALSE]
    colnames(beta) <- sp$x
    sigma <- d[, paste0("sigma[", seq_along(sp$g), "]"), drop = FALSE]
    z <- d[, grep("^z\\[", colnames(d)), drop = FALSE]
    off <- c(0, cumsum(n_level))
    ctx_x <- setdiff(sp$x, attr_cols)
    eta_side <- map(c(L = "L", R = "R"), \(s) {
      xc <- cbind(
        x_bat_left = rep(as.numeric(s == "L"), n_ctx),
        x_same_hand = as.numeric(ctx$pitch_hand == s),
        x_temp = ctx$x_temp,
        x_pos_pitcher = ctx$x_pos_pitcher
      )[, ctx_x, drop = FALSE]
      eta <- beta[, ctx_x, drop = FALSE] %*% t(xc)
      for (g in seq_along(sp$g)[-1]) {
        idx <- switch(sp$g[g],
                      pitcher_idx = ctx$pitcher_idx,
                      park_idx = ctx$park_idx,
                      fld_team_idx = ctx$fld_team_idx,
                      park_side_idx = match(paste(ctx$venue_name, s), ps_levels)
        )
        eff <- sigma[, g] * z[, off[g] + replace(idx, is.na(idx), 1L), drop = FALSE]
        eff[, is.na(idx)] <- 0
        eta <- eta + eff
      }
      eta
    })
    attr_x <- intersect(sp$x, attr_cols)
    attr_term <- if (length(attr_x) > 0) {
      beta[, attr_x, drop = FALSE] %*% t(as.matrix(batter_attr[attr_x]))
    } else {
      matrix(0, n_draw, n_bat)
    }
    b_dev <- sigma[, 1] * z[, seq_len(n_level[1]), drop = FALSE]
    list(alpha = d[, "alpha"], b_dev = b_dev, attr = attr_term, b_batter = b_dev + attr_term,
         sig1 = mean(sigma[, 1]), eta = eta_side)
  }
  st_list <- map(set_names(names(stage_def)), get_stage)

  # (1) 相関による縮小し直し：縮小し直した後の分布からサンプルを取り直す
  m <- sapply(corr_stages, \(s) colMeans(st_list[[s]]$b_dev))
  v <- sapply(corr_stages, \(s) apply(st_list[[s]]$b_dev, 2, var))
  sig <- map_dbl(corr_stages, \(s) {
    st_list[[s]]$sig1 * if (isTRUE(specs[[s]]$use_t == 2L)) sqrt(12 / 10) else 1
  })
  tau <- pmax(1 / v - matrix(1 / sig^2, n_bat, k4, byrow = TRUE), 0)
  sinv <- solve(diag(sig) %*% omega %*% diag(sig))
  rhs <- (m / v) * (tau > 0)
  new_b <- array(0, c(n_draw, n_bat, k4))
  for (i in seq_len(n_bat)) {
    v_i <- solve(diag(tau[i, ], k4) + sinv)
    mu <- as.numeric(v_i %*% rhs[i, ])
    new_b[, i, ] <- sweep(matrix(rnorm(n_draw * k4), n_draw) %*% chol(v_i), 2, mu, "+")
  }
  for (j in seq_len(k4)) {
    s <- corr_stages[j]
    st_list[[s]]$b_batter <- new_b[, , j] + st_list[[s]]$attr
  }

  # (2) Statcast データを使う段階の補正係数：推定誤差つきで引き直す
  cal <- map(set_names(ball_stages), \(st) {
    d <- readRDS(file.path(fit_dir, sprintf("prod_%d_ystar_%s.rds", season, st)))
    b <- d |>
      summarise(n = n(), k = sum(y), x_raw = mean(ystar), .by = batter_idx) |>
      filter(n >= 30, x_raw > 0, x_raw < 1)
    fit <- glm(cbind(k, n - k) ~ qlogis(x_raw), family = quasibinomial(), data = b)
    list(coef = coef(fit), draws = MASS::mvrnorm(n_draw, coef(fit), vcov(fit)))
  })

  stage_p <- function(st, i, side) {
    x <- st_list[[st]]
    lin <- x$alpha + x$b_batter[, i] + x$eta[[side]]
    if (st %in% ball_stages) lin <- cal[[st]]$draws[, 1] + cal[[st]]$draws[, 2] * lin
    plogis(lin)
  }
  cat_probs_side <- function(i, side) {
    p <- map(set_names(names(stage_def)), \(st) stage_p(st, i, side))
    r4 <- (1 - p$hbp) * (1 - p$bb) * (1 - p$so)
    r6 <- r4 * (1 - p$hr) * (1 - p$roe)
    hit <- r6 * p$hit
    list(
      hbp = p$hbp, bb = (1 - p$hbp) * p$bb,
      so = (1 - p$hbp) * (1 - p$bb) * p$so,
      hr = r4 * p$hr, roe = r4 * (1 - p$hr) * p$roe,
      x1b = hit * (1 - p$xbh), x2b = hit * p$xbh * (1 - p$x3b),
      x3b = hit * p$xbh * p$x3b, out = r6 * (1 - p$hit)
    )
  }
  p_arr <- array(NA_real_, c(n_draw, n_bat, length(cats)), dimnames = list(NULL, NULL, cats))
  for (i in seq_len(n_bat)) {
    w_l <- if_else(ctx$pitch_hand == "L",
                   batter_attr$q_same_vs_L[i], 1 - batter_attr$q_same_vs_R[i])
    acc <- matrix(0, n_draw, length(cats))
    if (any(w_l > 0)) {
      acc <- acc + vapply(cat_probs_side(i, "L"), \(x) as.numeric(x %*% w_l), numeric(n_draw))
    }
    if (any(w_l < 1)) {
      acc <- acc + vapply(cat_probs_side(i, "R"), \(x) as.numeric(x %*% (1 - w_l)), numeric(n_draw))
    }
    p_arr[, i, ] <- acc / n_ctx
  }

  w <- df_s |>
    left_join(select(rv_s, game_id, at_bat_number, rv_pa), by = c("game_id", "at_bat_number")) |>
    summarise(w = mean(rv_pa), .by = event) |>
    mutate(event = factor(event, cats)) |>
    arrange(event) |>
    pull(w)
  theta_nv <- matrix(matrix(p_arr, ncol = length(cats)) %*% w, n_draw, n_bat)

  # (3) 得点価値の補正：打者ごとの値の不確実性つきで加える
  theta <- theta_nv
  v_s <- val_all |> filter(season == .env$season)
  for (ev in c("out", "x1b", "x2b")) {
    ve <- filter(v_s, event == ev)
    dh <- numeric(n_bat)
    ds <- rep(ve$tau[1], n_bat)
    dh[ve$batter_idx] <- ve$delta_hat
    ds[ve$batter_idx] <- sqrt(ve$tau^2 * (ve$s2 / ve$n) / (ve$tau^2 + ve$s2 / ve$n))
    dk <- sweep(sweep(matrix(rnorm(n_draw * n_bat), n_draw), 2, ds, "*"), 2, dh, "+")
    theta <- theta + p_arr[, , ev] * dk
  }

  wt <- batter_attr$n_pa * (1 - batter_attr$x_pitcher_bat)
  to_plus <- function(th) 100 * (1 + (th - as.numeric(th %*% wt) / sum(wt)) / r_pa)
  plus_mat <- to_plus(theta)
  plus_nv <- to_plus(theta_nv)
  freq <- unclass(table(factor(df_s$batter_idx, seq_len(n_bat)), factor(df_s$event, cats)))
  theta_raw <- as.numeric((freq / rowSums(freq)) %*% w)
  theta_raw_lg <- sum(theta_raw * wt) / sum(wt)

  res <- batter_attr |>
    transmute(
      season = season, batter_idx, batter_id, batter_name, n_pa,
      is_pitcher = x_pitcher_bat == 1,
      plus_raw = 100 * (1 + (theta_raw - theta_raw_lg) / r_pa),
      plus = colMeans(plus_mat),
      plus_sd = apply(plus_mat, 2, sd),
      plus_q05 = apply(plus_mat, 2, quantile, 0.05),
      plus_q95 = apply(plus_mat, 2, quantile, 0.95),
      plus_no_value = colMeans(plus_nv),
      raa = colMeans(theta - as.numeric(theta %*% wt) / sum(wt)) * n_pa
    )
  out <- list(res = res, p_mean = apply(p_arr, c(2, 3), mean), freq = freq / rowSums(freq),
              w = set_names(w, cats), r_pa = r_pa, cal = map(cal, \(x) x$coef))
  saveRDS(out, out_path)
  out
}
