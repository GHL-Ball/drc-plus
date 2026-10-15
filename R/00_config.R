# =============================================================================
# 00_config.R
# パッケージ、パス、指標の定義（段階、イベント、モデル式）
# =============================================================================

library(tidyverse)
library(DBI)
library(RSQLite)
library(mgcv)
library(cmdstanr)
library(posterior)
library(here)

# ---- パス --------------------------------------------------------------------
# 個人の環境に依存するパスは .Renviron に書く（.Renviron.example を参照）。
#   DRC_DATA_DIR : フィットや中間結果の保存先（数十GBになる）
#   DRC_DB_PATH  : 一球データの SQLite ファイル
stopifnot(nzchar(Sys.getenv("DRC_DATA_DIR")), nzchar(Sys.getenv("DRC_DB_PATH")))
fit_dir <- file.path(Sys.getenv("DRC_DATA_DIR"), "fits")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
db_path <- Sys.getenv("DRC_DB_PATH")
pbp_table <- "pbp"   # 一球データのテーブル名（1 行が 1 球。2021〜2026 年を 1 つのテーブルに入れる）

# このリポジトリの場所（drc-plus.Rproj から開いた場合の作業ディレクトリ）
drc_root <- here()

seasons <- 2021:2026

# ---- イベントと段階 ------------------------------------------------------------
# 9カテゴリ。この並びが、確率の配列と得点価値のベクトルの並びになる。
cats <- c("hbp", "bb", "so", "hr", "roe", "x1b", "x2b", "x3b", "out")

# 逐次的な二項モデル。pool = その段階に進む打席、success = 「当たり」とするカテゴリ
stage_def <- list(
  hbp = list(pool = c("hbp", "bb", "so", "hr", "roe", "x1b", "x2b", "x3b", "out"),
             success = "hbp"),
  bb  = list(pool = c("bb", "so", "hr", "roe", "x1b", "x2b", "x3b", "out"),
             success = "bb"),
  so  = list(pool = c("so", "hr", "roe", "x1b", "x2b", "x3b", "out"),
             success = "so"),
  hr  = list(pool = c("hr", "roe", "x1b", "x2b", "x3b", "out"), success = "hr"),
  roe = list(pool = c("roe", "x1b", "x2b", "x3b", "out"), success = "roe"),
  hit = list(pool = c("x1b", "x2b", "x3b", "out"), success = c("x1b", "x2b", "x3b")),
  xbh = list(pool = c("x1b", "x2b", "x3b"), success = c("x2b", "x3b")),
  x3b = list(pool = c("x2b", "x3b"), success = "x3b")
)

# 相関を使って縮小し直す段階
corr_stages <- c("bb", "so", "hr", "xbh")

# 常に野手水準の打者として扱う選手（大谷翔平）
two_way_ids <- c(660271)

# ---- 打球単位のモデル式 --------------------------------------------------------
# 本塁打：EV・LA・SA（打者基準 & 球場基準）、気温
gam_hr_formula <- y ~ te(launch_speed, launch_angle, spray_pull, k = c(8, 8, 6)) +
  s(spray_field, park_f, bs = "fs", k = 6, m = 1) + x_temp

# 安打（BIP）、長打：上に加えて一塁到達タイムと守備側チーム。
# 走力の変数は交差検証で選んだ（一塁到達タイムのみの式）。
gam_hit_formula <- y ~ te(launch_speed, launch_angle, spray_pull, k = c(8, 8, 6)) +
  s(avg_hp_to_1b, k = 5) + ti(launch_angle, avg_hp_to_1b, k = c(6, 4)) +
  s(fld_f, bs = "re") + s(spray_field, park_f, bs = "fs", k = 4, m = 1) + x_temp

# 得点価値の補正：同じイベント・同じ塁アウト状況の平均との差を、打球の性質と打者の属性で説明する。
# sg は塁アウト状況の5つのカテゴリ。
val_formula <- resid ~ sg + s(launch_angle, by = sg, k = 6) + s(launch_speed, by = sg, k = 5) +
  s(spray_field, by = sg, k = 5) + sg:hp1b_c + sg:x_bat_left
