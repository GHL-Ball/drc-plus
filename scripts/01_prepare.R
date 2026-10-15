# =============================================================================
# 01_prepare.R
# 関数とデータの読み込み（他のスクリプトの前に必ず実行する）
# =============================================================================

library(here)
drc_root <- here()
r_files <- sort(list.files(file.path(drc_root, "R"), pattern = "\\.R$", full.names = TRUE))
if (length(r_files) == 0) stop("R/ が見つかりません。drc-plus.Rproj からプロジェクトとして開いてください: ", drc_root)
for (f in r_files) {
  source(f, encoding = "UTF-8")
}

# 打席データ（2021〜2026年）と、イベントの9分類
df_pa <- load_pa() |>
  mutate(event = classify_event(events))

# 投手の打席、野手登板
role_tbl <- make_role_tbl(df_pa)
ppp_tbl <- make_ppp_tbl(role_tbl)

# 打席ごとの得点価値
rv_all <- load_rv_all()

# ---- 確認 --------------------------------------------------------------------
# イベントの件数（NA は対象外の打席）
df_pa |>
  count(year, event) |>
  pivot_wider(names_from = year, values_from = n) |>
  print(n = Inf)

# 投手と判定した打者の数
role_tbl |> count(year, is_pitcher) |>
  pivot_wider(names_from = is_pitcher, values_from = n)

# R/PA
rv_all |> summarise(n = n(), r_pa = first(r_pa), .by = year)

# 打席数の上限（各年500前後、上限以上の打者が130人前後）
map(seasons, \(s) {
  d <- make_season(df_pa, role_tbl, s)
  info <- pa_cap_info(d)
  n_pa <- count(d, batter_idx)$n
  tibble(season = s, n_teams = n_distinct(d$home_team), n_slots = info$n_slots,
         pa_cap = round(info$cap, 1), pa_of_last_slot = info$min_top,
         n_batters = length(n_pa), n_at_cap = sum(n_pa >= info$cap))
}) |>
  list_rbind()
