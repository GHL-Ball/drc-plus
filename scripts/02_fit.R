# =============================================================================
# 02_fit.R
# 段階ごとのフィット（6シーズン × 2つの版 × 8段階 = 96個）
# =============================================================================

# 所要時間は1シーズンあたり約4時間半（4コア）。
# 段階ごとに保存するので、途中で止めても、もう一度実行すれば続きから始まる。
if (!(exists("df_pa") && exists("rv_all"))) stop("先に scripts/01_prepare.R を実行してください")

mod_bin <- cmdstan_model(file.path(drc_root, "stan", "stage_binary_v3.stan"))
mod_frac <- cmdstan_model(file.path(drc_root, "stan", "stage_frac_v1.stan"))

# 再開する前に：書き込みの途中で止めた場合、読めないファイルが残ることがある。
# 次の4行で見つけて消してから、下のループを実行する。
# files <- list.files(fit_dir, pattern = "^prod_.*\\.rds$", full.names = TRUE)
# ok <- map_lgl(files, \(f) !inherits(try(readRDS(f), silent = TRUE), "try-error"))
# basename(files[!ok])
# file.remove(files[!ok])

for (season in c(2024, 2025, 2023, 2022, 2026, 2021)) {
  tryCatch(run_season(season),
           error = \(e) message("!!! ", season, " で停止: ", conditionMessage(e)))
  gc()
}

# ---- 進み具合と診断 ------------------------------------------------------------
done <- expand_grid(season = seasons, ver = c("pac", "nopa"), stage = names(stage_def)) |>
  mutate(done = file.exists(prod_path(season, ver, stage)))
done |> summarise(n_done = sum(done), .by = c(season, ver)) |>
  pivot_wider(names_from = ver, values_from = n_done)

# 目安：発散が全体の1%未満、R-hat が1.05以下
diag_all <- expand_grid(season = seasons, ver = c("nopa", "pac")) |>
  pmap(\(season, ver) mutate(prod_report(season, ver), season = season, ver = ver, .before = 1)) |>
  list_rbind()
diag_all |> summarise(minutes = sum(minutes), divergent = sum(divergent),
                      max_rhat = max(max_rhat), min_ess = min(min_ess), .by = c(season, ver))
