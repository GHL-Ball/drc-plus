# =============================================================================
# 03_postprocess.R
# 後処理：フィットの結果から最終的な指標を作る
# =============================================================================

# 所要時間は約2時間（1シーズン・1つの版あたり10分前後）。済んだ分は保存され、再実行時は飛ばす。
if (!(exists("df_pa") && exists("rv_all"))) stop("先に scripts/01_prepare.R を実行してください")

data_dir <- Sys.getenv("DRC_DATA_DIR")

# ---- (1) 段階間の相関（四球・三振・本塁打・長打） ----------------------------------
# 指標に使うのは、全シーズン共通の値（pooled）
omega_path <- file.path(data_dir, "omega_by_ver.rds")
omega_by_ver <- if (file.exists(omega_path)) {
  readRDS(omega_path)
} else {
  map(c(pac = "pac", nopa = "nopa"), estimate_omega)
}
saveRDS(omega_by_ver, omega_path)

# シーズンごとの相関と、共通の値
low <- lower.tri(diag(4))
pair_names <- outer(corr_stages, corr_stages, paste, sep = "-")[low]
map(omega_by_ver$nopa$each, \(om) set_names(round(om[low], 2), pair_names)) |>
  bind_rows(.id = "season") |>
  bind_rows(tibble(season = "pooled", !!!set_names(round(omega_by_ver$nopa$pooled[low], 2), pair_names)))
round(omega_by_ver$pac$pooled, 2)

# ---- (2) 得点価値の補正（アウト・単打・二塁打） ---------------------------------------
val_path <- file.path(data_dir, "value_layer_v1.rds")
val_all <- if (file.exists(val_path)) {
  readRDS(val_path)
} else {
  map(seasons, \(s) { message("value: ", s); value_season(s) }) |> list_rbind()
}
saveRDS(val_all, val_path)
# 打者間の差の標準偏差（アウトは毎年0.005前後）
val_all |> distinct(season, event, tau) |> pivot_wider(names_from = event, values_from = tau)

# ---- (3) 指標 ----------------------------------------------------------------
grid <- expand_grid(season = seasons, ver = c("pac", "nopa"))
post_v3 <- pmap(grid, \(season, ver) {
  message(sprintf("post v3: %d / %s  %s", season, ver, format(Sys.time(), "%H:%M")))
  post_season(season, ver, omega_by_ver[[ver]]$pooled)
})

# 主指標は打席数を入れない版（nopa）。打席数を入れた版（pac）は参考列 plus_pa。
#   plus          : 主指標（得点価値の補正を含む）
#   plus_sd、plus_q05、plus_q95 : 事後分布の標準偏差と90%区間
#   plus_no_value : 得点価値の補正を除いた値
#   plus_raw      : 縮小も環境の調整もしない値（実際の頻度 × 得点価値）
res_v3 <- map2(post_v3, grid$ver, \(o, v) mutate(o$res, ver = v)) |> list_rbind()
drc_final <- res_v3 |>
  filter(ver == "nopa") |>
  select(-ver, -batter_idx) |>
  left_join(res_v3 |> filter(ver == "pac") |> select(season, batter_id, plus_pa = plus),
            by = c("season", "batter_id"))
saveRDS(drc_final, file.path(data_dir, "drc_final_v3.rds"))

# ---- 確認 --------------------------------------------------------------------
# シーズンごとの散らばり（300打席以上の野手）
drc_final |>
  filter(n_pa >= 300, !is_pitcher) |>
  summarise(
    n = n(),
    sd_raw = sd(plus_raw), sd_own = sd(plus), sd_own_pa = sd(plus_pa),
    mean_sd = mean(plus_sd),
    .by = season
  )

# Statcast データを使う段階の補正係数（傾き）と、得点価値
map2(post_v3[grid$ver == "nopa"], seasons, \(o, s) {
  tibble(season = s, hr_slope = o$cal$hr[2], hit_slope = o$cal$hit[2], xbh_slope = o$cal$xbh[2],
         w_hr = o$w["hr"], w_x1b = o$w["x1b"], w_bb = o$w["bb"], w_so = o$w["so"], w_out = o$w["out"])
}) |>
  list_rbind()

# 各年の上位5人
drc_final |>
  filter(n_pa >= 300, !is_pitcher) |>
  slice_max(plus, n = 5, by = season) |>
  select(season, batter_name, n_pa, plus_raw, plus, plus_sd, plus_q05, plus_q95,
         plus_pa) |>
  print(n = Inf)

# 得点価値の補正の大きさ（自作 DRC+ の単位）
drc_final |>
  filter(n_pa >= 300, !is_pitcher) |>
  mutate(val_adj = plus - plus_no_value) |>
  summarise(sd_adj = sd(val_adj), min_adj = min(val_adj), max_adj = max(val_adj), .by = season)
