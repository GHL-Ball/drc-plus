# =============================================================================
# 01_data.R
# データの読み込みと、シーズン単位・段階単位のデータ作成
# =============================================================================

# ---- 一球データ（打席の最終投球） ------------------------------------------------
# sabRmetrics で取得した一球データに、独自に列を加えたもの。使う列は次のとおり。
# load_rv_all() では、ほかに run_value_actual と R/PA を使う（README を参照）。
pa_cols <- c(
  # 試合と打席
  "game_id", "year", "at_bat_number", "half_inning", "home_team", "away_team",
  "venue_name",
  # 打者と投手
  "batter_id", "batter_name", "bat_side", "pitcher_id", "pitch_hand",
  "handedness_type",
  # 結果と塁アウト状況
  "events", "outs", "bases",
  # 打球
  "launch_speed", "launch_angle", "hit_coord_x", "hit_coord_y",
  # 環境と走力
  "temp_f_spin", "avg_hp_to_1b"
)

# 打席の最後の投球だけを読む
load_pa <- function() {
  sel <- paste0('"', pa_cols, '"', collapse = ", ")
  con <- dbConnect(SQLite(), db_path, flags = SQLITE_RO)
  on.exit(dbDisconnect(con), add = TRUE)
  dbGetQuery(con, paste("SELECT", sel, "FROM", pbp_table, "WHERE is_last_pitch = 1")) |>
    as_tibble()
}

# 打席ごとの得点価値（rv_pa）と、そのシーズンの R/PA（r_pa）。
# rv_pa は、その打席の全投球の run_value_actual の合計。
load_rv_all <- function() {
  path <- file.path(Sys.getenv("DRC_DATA_DIR"), "rv_all.rds")
  if (file.exists(path)) return(readRDS(path))
  rv_all <- local({
    con <- dbConnect(SQLite(), db_path, flags = SQLITE_RO)
    on.exit(dbDisconnect(con), add = TRUE)
    dbGetQuery(con, paste(
      'SELECT year, game_id, at_bat_number,',
      'SUM(run_value_actual) AS rv_pa, MAX("R/PA") AS r_pa',
      "FROM", pbp_table, "GROUP BY year, game_id, at_bat_number"
    ))
  }) |>
    as_tibble()
  saveRDS(rv_all, path)
  rv_all
}

# ---- イベントの分類 ------------------------------------------------------------
# 敬遠、打撃妨害、犠打、打席途中の交代などは NA（対象外）になる
classify_event <- function(events) {
  case_when(
    events %in% c("strikeout", "strikeout_double_play") ~ "so",
    events == "walk" ~ "bb",
    events == "hit_by_pitch" ~ "hbp",
    events == "single" ~ "x1b",
    events == "double" ~ "x2b",
    events == "triple" ~ "x3b",
    events == "home_run" ~ "hr",
    events == "field_error" ~ "roe",
    events %in% c("field_out", "force_out", "grounded_into_double_play", "double_play",
                  "fielders_choice", "fielders_choice_out", "sac_fly",
                  "sac_fly_double_play", "triple_play") ~ "out",
    .default = NA_character_
  )
}

# ---- 投手の打席と、野手登板 ----------------------------------------------------
# 年度ごとに、打席数より対戦打者数が多い選手を投手とみなす
make_role_tbl <- function(df_pa) {
  df_valid <- df_pa |> filter(!is.na(event))
  df_valid |>
    count(year, batter_id, name = "n_pa_all") |>
    left_join(
      count(df_valid, year, pitcher_id, name = "n_bf"),
      by = join_by(year, batter_id == pitcher_id)
    ) |>
    mutate(
      n_bf = coalesce(n_bf, 0L),
      is_pitcher = n_bf > n_pa_all & !batter_id %in% two_way_ids
    )
}

# 野手登板（投手ではない選手）
make_ppp_tbl <- function(role_tbl) {
  role_tbl |>
    filter(!is_pitcher, n_bf > 0, !batter_id %in% two_way_ids) |>
    transmute(year, pitcher_id = batter_id, is_pos_pitcher = TRUE)
}

# ---- シーズン単位のデータ --------------------------------------------------------
# 番号（*_idx）はシーズン内で振る。2022年以降は投手の打席を除く。
make_season <- function(df_pa, role_tbl, season) {
  df_pa |>
    filter(year == .env$season, !is.na(event)) |>
    left_join(select(role_tbl, year, batter_id, is_pitcher),
              by = c("year", "batter_id")) |>
    filter(.env$season == 2021 | !is_pitcher) |>
    mutate(
      bat_team = if_else(half_inning == "bottom", home_team, away_team),
      fld_team = if_else(half_inning == "bottom", away_team, home_team),
      batter_idx = as.integer(factor(batter_id)),
      pitcher_idx = as.integer(factor(pitcher_id)),
      park_idx = as.integer(factor(venue_name)),
      park_side_idx = as.integer(factor(paste(venue_name, bat_side))),
      fld_team_idx = as.integer(factor(fld_team)),
      x_pitcher_bat = as.integer(is_pitcher),
      x_bat_left = as.integer(bat_side == "L"),
      x_same_hand = as.integer(handedness_type == "Same"),
      x_temp = (temp_f_spin - 70) / 10
    )
}

# 打席数の上限：打順の枠（球団数 × 9。投手が打つ年は野手の打席の割合を掛ける）の数だけ、
# 打席数の多い順に野手を取り、その平均打席数を上限とする
pa_cap_info <- function(d) {
  pa_tbl <- d |> count(batter_idx, is_pitcher, name = "n_pa")
  n_slots <- round(n_distinct(d$home_team) * 9 *
                     sum(pa_tbl$n_pa[!pa_tbl$is_pitcher]) / sum(pa_tbl$n_pa))
  top <- sort(pa_tbl$n_pa[!pa_tbl$is_pitcher], decreasing = TRUE)[seq_len(n_slots)]
  list(n_slots = n_slots, cap = mean(top), min_top = min(top))
}

# 本番用：野手登板のフラグと、上限つきの打席数の対数（上限以上の打者は0）を加える。
# df_pa、role_tbl、ppp_tbl はグローバルのものを使う。
make_season_prod <- function(season) {
  d <- make_season(df_pa, role_tbl, season)
  cap <- pa_cap_info(d)$cap
  d |>
    left_join(ppp_tbl, by = c("year", "pitcher_id")) |>
    mutate(x_pos_pitcher = as.integer(coalesce(is_pos_pitcher, FALSE))) |>
    mutate(n_pa_season = n(), .by = batter_idx) |>
    mutate(pa_cap = cap,
           x_logpa = log1p(pmin(n_pa_season, cap)) - log1p(cap))
}

# 段階のデータ：その段階に進んだ打席だけを残し、当たりを y とする
make_stage <- function(df_season, stage) {
  def <- stage_def[[stage]]
  df_season |>
    filter(event %in% def$pool) |>
    mutate(y = as.integer(event %in% def$success))
}
