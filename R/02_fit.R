# =============================================================================
# 02_fit.R
# 段階ごとのフィット（結果を使う段階と、Statcast データを使う段階）
# =============================================================================

# mod_bin（stage_binary_v3.stan）と mod_frac（stage_frac_v1.stan）は、
# 実行スクリプトの側でコンパイルしてグローバルに置く。

# ---- Stan に渡すデータ ----------------------------------------------------------
# 同じランダム効果・固定効果の組み合わせの打席をまとめて、行数を減らす
build_stan_data <- function(df_stage, x_cols, g_cols, n_level) {
  df_agg <- df_stage |>
    summarise(n = n(), y = sum(y), .by = all_of(c(g_cols, x_cols)))
  list(
    df_agg = df_agg,
    data = list(
      N = nrow(df_agg), K = length(x_cols), X = as.matrix(df_agg[x_cols]),
      G = length(g_cols), n_level = n_level, idx = as.matrix(df_agg[g_cols]),
      n = df_agg$n, y = df_agg$y
    )
  )
}

# ウォームアップの短縮：パラメータごとの事後分散を Fisher 情報量から概算し、
# 質量行列の初期値として渡す。並びは alpha, beta, log(sigma), z。
init_inv_metric <- function(df_agg, x_cols, g_cols, n_level, sigma0) {
  p0 <- sum(df_agg$y) / sum(df_agg$n)
  w <- p0 * (1 - p0)
  n_tot <- sum(df_agg$n)
  v_alpha <- sum(sigma0^2 / n_level) + 1 / (w * n_tot)
  v_beta <- map_dbl(x_cols, \(x) {
    xv <- df_agg[[x]]
    1 / (w * sum(df_agg$n * (xv - weighted.mean(xv, df_agg$n))^2))
  })
  v_log_sigma <- 1 / n_level
  v_z <- map(seq_along(g_cols), \(g) {
    n_lv <- numeric(n_level[g])
    s <- rowsum(df_agg$n, df_agg[[g_cols[g]]])
    n_lv[as.integer(rownames(s))] <- s[, 1]
    1 / (1 + sigma0[g]^2 * w * n_lv)
  })
  c(v_alpha, v_beta, v_log_sigma, unlist(v_z))
}

# ---- 段階の仕様 ----------------------------------------------------------------
# type  : "bin" は結果を使う段階、"ball" は Statcast データを使う段階
# x     : 固定効果、g : ランダム効果（1つ目は必ず打者）、sigma0 : 質量行列の概算に使う標準偏差
# use_t : 打者の事前分布。0 = 正規、2 = t（自由度12で固定）
prod_spec <- function(season, use_pa) {
  common <- c(if (season == 2021) "x_pitcher_bat", "x_bat_left", "x_same_hand",
              "x_pos_pitcher", if (use_pa) "x_logpa")
  bp <- c("batter_idx", "pitcher_idx")
  list(
    hbp = list(type = "bin", x = common, g = bp, sigma0 = c(0.4, 0.3), use_t = 0L),
    bb  = list(type = "bin", x = c(common, "x_temp"), g = c(bp, "park_idx"),
               sigma0 = c(0.4, 0.3, 0.1), use_t = 0L),
    so  = list(type = "bin", x = c(common, "x_temp"), g = c(bp, "park_idx"),
               sigma0 = c(0.4, 0.3, 0.1), use_t = 2L),
    hr  = list(type = "ball", x = common, g = bp, sigma0 = c(0.5, 0.2),
               formula = gam_hr_formula),
    roe = list(type = "bin", x = common, g = c(bp, "fld_team_idx", "park_idx"),
               sigma0 = c(0.2, 0.1, 0.1, 0.1), use_t = 0L),
    hit = list(type = "ball", x = common, g = bp, sigma0 = c(0.15, 0.1),
               formula = gam_hit_formula),
    xbh = list(type = "ball", x = common, g = bp, sigma0 = c(0.3, 0.15),
               formula = gam_hit_formula),
    x3b = list(type = "bin", x = common, g = c(bp, "fld_team_idx", "park_side_idx"),
               sigma0 = c(0.3, 0.1, 0.1, 0.2), use_t = 0L)
  )
}

prod_path <- function(season, ver, stage) {
  file.path(fit_dir, sprintf("prod_%d_%s_%s.rds", season, ver, stage))
}

# ---- 結果を使う段階 -----------------------------------------------------------
fit_bin_prod <- function(df_s, stage, spec, path) {
  if (file.exists(path)) return(invisible(path))
  n_level <- map_int(spec$g, \(g) max(df_s[[g]]))
  sd <- build_stan_data(make_stage(df_s, stage), spec$x, spec$g, n_level)
  im <- init_inv_metric(sd$df_agg, spec$x, spec$g, n_level, spec$sigma0)
  fit <- mod_bin$sample(
    data = c(sd$data, list(use_t = spec$use_t, nu_fixed = 12)),
    seed = 20261008, chains = 4, parallel_chains = 4,
    iter_warmup = 500, iter_sampling = 500, init = 0.1, refresh = 250,
    metric = "diag_e", inv_metric = im, output_dir = fit_dir
  )
  fit$save_object(path)
  invisible(path)
}

# ---- Statcast データを使う段階 -----------------------------------------------------------
# 学習済みの打球単位のモデルで、新しい打球の「中立環境での確率」と「実際の環境での確率」を出す。
# env は中立環境（球場、守備側チーム、気温の組）。
neutral_env <- function(gam_h, new, env) {
  lab <- sapply(gam_h$smooth, \(s) s$label)
  lab_park <- lab[grepl("park_f", lab)]
  lab_fld <- lab[grepl("fld_f", lab)]
  vars <- setdiff(all.vars(formula(gam_h)), "y")
  new_min <- as.data.frame(new[, vars])
  tm <- predict(gam_h, newdata = new_min, type = "terms")
  eta_base <- coef(gam_h)[["(Intercept)"]] +
    rowSums(tm[, setdiff(lab, c(lab_park, lab_fld)), drop = FALSE])

  parks <- levels(new$park_f)
  grid <- new_min[rep(seq_len(nrow(new_min)), times = length(parks)), ]
  grid$park_f <- factor(rep(parks, each = nrow(new_min)), levels = parks)
  g_mat <- matrix(predict(gam_h, newdata = grid, type = "terms", terms = lab_park)[, 1],
                  nrow = nrow(new_min), dimnames = list(NULL, parks))
  rm(grid)

  ctx_term <- coef(gam_h)[["x_temp"]] * env$x_temp
  if (length(lab_fld) > 0) {
    teams <- levels(new$fld_f)
    one <- new_min[rep(1, length(teams)), ]
    one$fld_f <- factor(teams, levels = teams)
    u_fld <- set_names(predict(gam_h, newdata = one, type = "terms", terms = lab_fld)[, 1],
                       teams)
    ctx_term <- ctx_term + u_fld[env$fld]
  }
  eta <- eta_base + g_mat[, env$park] +
    matrix(ctx_term, nrow(new_min), nrow(env), byrow = TRUE)
  list(p_neutral = rowMeans(plogis(eta)),
       p_env = as.numeric(predict(gam_h, newdata = new_min, type = "response")))
}

# 打球ごとの中立化した確率。試合単位の交差検証で、その打球を学習に使っていないモデルから出す。
# 対象外の打球（計測なしなど）は実際の結果（0/1）をそのまま使う。
ball_ystar <- function(df_s, stage, gam_formula, path, n_fold = 5) {
  if (file.exists(path)) return(readRDS(path))
  set.seed(20261008)
  pool <- make_stage(df_s, stage) |>
    mutate(
      spray_field = atan2(hit_coord_x - 125.42, 198.27 - hit_coord_y) * 180 / pi,
      spray_pull = if_else(bat_side == "L", spray_field, -spray_field),
      park_f = factor(venue_name),
      fld_f = factor(fld_team),
      ok = !is.na(launch_speed) & !is.na(launch_angle) & !is.na(spray_field)
    )
  pool$in_domain <- if (stage == "hr") {
    pool$ok & pool$launch_speed >= 90 & between(pool$launch_angle, 10, 50) &
      abs(pool$spray_field) <= 50
  } else {
    pool$ok & abs(pool$spray_field) <= 55
  }
  games <- unique(pool$game_id)
  fold_of <- set_names(sample(rep(seq_len(n_fold), length.out = length(games))), games)
  pool$cf <- fold_of[as.character(pool$game_id)]
  env <- df_s |> slice_sample(n = 300) |>
    transmute(park = as.character(venue_name), fld = as.character(fld_team), x_temp)
  pool$ystar <- as.numeric(pool$y)
  pool$p_env <- as.numeric(pool$y)
  for (k in seq_len(n_fold)) {
    gam_k <- bam(gam_formula, data = filter(pool, in_domain, cf != k),
                 family = binomial(), method = "fREML", discrete = TRUE, nthreads = 4)
    ii <- which(pool$in_domain & pool$cf == k)
    res <- neutral_env(gam_k, pool[ii, ], env)
    pool$ystar[ii] <- res$p_neutral
    pool$p_env[ii] <- res$p_env
  }
  out <- pool |>
    select(game_id, at_bat_number, batter_idx, pitcher_idx, starts_with("x_"),
           y, in_domain, ystar, p_env)
  saveRDS(out, path)
  out
}

# 打者の段階モデル。0/1 の代わりに確率を目的変数にし、情報量の比 phi で尤度を重み付けする。
# 固定効果の値は、保存済みのファイルではなく、その時点のシーズンデータから取る。
fit_ball_prod <- function(d, df_s, spec, path) {
  if (file.exists(path)) return(invisible(path))
  d <- d |>
    select(game_id, at_bat_number, batter_idx, pitcher_idx, y, in_domain, ystar, p_env) |>
    left_join(select(df_s, game_id, at_bat_number, all_of(spec$x)),
              by = c("game_id", "at_bat_number"))
  stopifnot(!anyNA(d[spec$x]))
  n_level <- map_int(spec$g, \(g) max(df_s[[g]]))
  df_agg <- d |>
    summarise(n = n(), y = sum(ystar), .by = all_of(c(spec$g, spec$x)))
  im <- init_inv_metric(df_agg, spec$x, spec$g, n_level, spec$sigma0)
  fit <- mod_frac$sample(
    data = list(N = nrow(df_agg), K = length(spec$x), X = as.matrix(df_agg[spec$x]),
                G = length(spec$g), n_level = n_level, idx = as.matrix(df_agg[spec$g]),
                n = df_agg$n, ystar = df_agg$y, phi = var(d$y) / var(d$ystar)),
    seed = 20261008, chains = 4, parallel_chains = 4,
    iter_warmup = 500, iter_sampling = 500, init = 0.1, refresh = 250,
    metric = "diag_e", inv_metric = im, output_dir = fit_dir
  )
  fit$save_object(path)
  invisible(path)
}

# ---- 1シーズン分（2つの版 × 8段階） -----------------------------------------------
# "pac" = 打席数（上限つき）を入れた版、"nopa" = 入れない版。
# 段階ごとに保存し、保存済みの段階は飛ばす。
run_season <- function(season) {
  df_s <- make_season_prod(season)
  for (ver in c("pac", "nopa")) {
    specs <- prod_spec(season, ver == "pac")
    for (st in names(specs)) {
      message(sprintf("=== %d / %s / %s === %s", season, ver, st, format(Sys.time(), "%m-%d %H:%M")))
      sp <- specs[[st]]
      if (sp$type == "bin") {
        fit_bin_prod(df_s, st, sp, prod_path(season, ver, st))
      } else {
        d <- ball_ystar(df_s, st, sp$formula,
                        file.path(fit_dir, sprintf("prod_%d_ystar_%s.rds", season, st)))
        fit_ball_prod(d, df_s, sp, prod_path(season, ver, st))
      }
    }
  }
  invisible(NULL)
}

# 診断の要約
prod_report <- function(season, ver) {
  map(names(stage_def), \(st) {
    f <- readRDS(prod_path(season, ver, st))
    s <- f$summary(c("alpha", "sigma"))
    tibble(stage = st, minutes = round(f$time()$total / 60),
           divergent = sum(f$diagnostic_summary()$num_divergent),
           sigma_batter = s$mean[s$variable == "sigma[1]"],
           max_rhat = max(s$rhat), min_ess = min(s$ess_bulk))
  }) |>
    list_rbind()
}
