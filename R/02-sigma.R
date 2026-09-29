# 02-sigma.R
# Аналитические цели по биологической вариации, сигма-метрика и обратная задача.
#
# Источник формул (открытые публикации):
#   Fraser CG, Petersen PH. Desirable standards for laboratory tests if they are to
#     fulfill medical needs. Clin Chem. 1993;39(7):1447-53.
#   Ricos C, Alvarez V, Cava F, et al. Scand J Clin Lab Invest. 1999;59(7):491-500.
#   Реализация-ориентир: пакет valytics (CRAN), функция ate_from_bv().
#
# Допущение, существенное для всего сопоставления: допустимая суммарная ошибка TEa
# ОДНА И ТА ЖЕ для обеих рамок, то есть сигма-метрика считается по TEa, выведенной
# из биологической вариации. Если лаборатория берёт TEa для сигмы из другого
# источника (CLIA, RiliBAK, EQA), сопоставление относится уже к другой задаче;
# этот случай обсуждается в статье, но не является предметом расчёта.
#
# Формулы:
#   CV_A  = imp_mult  * CV_I
#   Bias  = bias_mult * sqrt(CV_I^2 + CV_G^2)
#   TEa   = k * CV_A + Bias,  k = 1.65
#   sigma = (TEa - |Bias_набл|) / CV_набл
#
# Обратная задача: CV_требуемая(S) = (TEa - |Bias_набл|) / S.
# Два сценария по наблюдаемому смещению:
#   bias = 0      -> CV_треб = TEa / S;
#   bias = допуск -> CV_треб = k * CV_A / S.
#
# Замкнутые формы отношения к допустимой неточности по биологической вариации:
#   bias = допуск: CV_треб / CV_A = k / S                  (не зависит ни от аналита,
#                                                           ни от уровня строгости)
#   bias = 0     : CV_треб / CV_A = (k + 0.5 * R) / S,     R = sqrt(1 + (CV_G/CV_I)^2)
#                  (зависит только от отношения CV_G/CV_I, но не от уровня)
# Уровень строгости сокращается именно в отношении к CV_A; отношение к CV_I
# (CV_треб/CV_I = k * imp_mult / S в сценарии bias = допуск) от уровня ЗАВИСИТ.
#
# Уровни: optimal 0.25/0.125, desirable 0.50/0.25, minimum 0.75/0.375.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# --- запись и чтение CSV ------------------------------------------------------
# Разделитель — точка с запятой. Десятичный разделитель в числах остаётся
# точкой, но файл с разделителем-запятой некорректно разбирается на колонки в
# локали, где запятая служит десятичным разделителем.
write_csv_sc <- function(x, path) {
  write.table(x, path, sep = ";", row.names = FALSE, col.names = TRUE,
              quote = TRUE, qmethod = "double", na = "NA", fileEncoding = "UTF-8")
}
read_csv_sc <- function(path) {
  read.csv(path, sep = ";", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
args_all <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
root <- normalizePath(file.path(here, ".."), mustWork = TRUE)
out_dir <- file.path(root, "out")

bv <- read_csv_sc(file.path(out_dir, "bv_meta.csv"))

# --- из длинной таблицы в широкую: одна строка на аналит --------------------
w <- bv %>%
  select(analyte_id, analyte_name, var_type, median, lower, upper, number_used,
         matrix, source) %>%
  pivot_wider(
    id_cols = c(analyte_id, analyte_name, matrix, source),
    names_from = var_type,
    values_from = c(median, lower, upper, number_used),
    names_sep = "."
  )

# аналиты без CVg: для слагаемого смещения берём CV_G = 0 и помечаем это
w <- w %>%
  mutate(
    cvg_imputed = is.na(median.cvg),
    median.cvg_used = ifelse(is.na(median.cvg), 0, median.cvg),
    # R = sqrt(1 + (CV_G/CV_I)^2) — относительный вклад межиндивидуальной вариации
    rel_cvg = sqrt(1 + (median.cvg_used / median.cvi)^2),
    # границы интервалов входных оценок (нужны для анализа устойчивости, блок S7);
    # при отсутствии границы берётся сама медиана
    cvi_lo = ifelse(is.na(lower.cvi), median.cvi, lower.cvi),
    cvi_hi = ifelse(is.na(upper.cvi), median.cvi, upper.cvi),
    cvg_lo = ifelse(is.na(lower.cvg), median.cvg_used, lower.cvg),
    cvg_hi = ifelse(is.na(upper.cvg), median.cvg_used, upper.cvg)
  )

# --- исключение аналитов, для которых нет входных оценок --------------------
# Ферритин исключён из расчёта: база не публикует для него мета-анализ CV_G
# (в матрице сыворотки CV_G приводит только одно из двух исследований, вошедших
# в мета-анализ CV_I), а подстановка CV_G = 0 создавала бы крайний случай
# «требование сигма-метрики строже», работающий в пользу основного вывода, и
# делала бы этот случай нижней границей приводимых диапазонов. Состав панели
# при этом сохраняется полностью: см. data/panel_analyte_ids.csv, столбец
# used_in_analysis.
EXCLUDED_ANALYTES <- "Ferritin"

w_all <- w
w <- w %>% filter(!analyte_name %in% EXCLUDED_ANALYTES)

# --- уровни аналитических целей ---------------------------------------------
levels_tbl <- tibble::tribble(
  ~level,      ~imp_mult, ~bias_mult,
  "optimal",        0.25,      0.125,
  "desirable",      0.50,      0.250,
  "minimum",        0.75,      0.375
)
K_COVERAGE <- 1.65
# Набор значений коэффициента охвата для анализа чувствительности:
# 1.0 (~84 % одностороннего охвата), 1.65 (~95 %), 1.96 (~97.5 %), 2.0 (~97.7 %),
# 2.58 (~99.5 %).
# Основной расчёт во всех остальных выходных файлах выполняется при K_COVERAGE.
K_SET <- c(1.0, K_COVERAGE, 1.96, 2.0, 2.58)
SIGMA_TARGETS <- c(4, 5, 6)

spec <- w %>%
  crossing(levels_tbl) %>%
  mutate(
    cv_a_allow   = imp_mult * median.cvi,
    bias_allow   = bias_mult * sqrt(median.cvi^2 + median.cvg_used^2),
    tea          = K_COVERAGE * cv_a_allow + bias_allow,
    # отношение смещения к допустимой неточности: bias_allow / CV_A = 0.5 * R
    bias_over_cv_a = bias_allow / cv_a_allow
  )

# --- требование к неточности метода для достижения целевой сигмы ------------
# sigma(c) = (TEa - |bias|) / c  =>  c_req(S) = (TEa - |bias|) / S
req <- spec %>%
  crossing(sigma_target = SIGMA_TARGETS) %>%
  mutate(
    bias_scenario  = "bias=0",
    c_req          = tea / sigma_target,
    c_req_over_cv_a = c_req / cv_a_allow,
    c_req_over_cvi  = c_req / median.cvi,
    # ожидаемые значения замкнутых форм (контроль ниже)
    ratio_cv_a_expected = (K_COVERAGE + 0.5 * rel_cvg) / sigma_target,
    ratio_cvi_expected  = imp_mult * (K_COVERAGE + 0.5 * rel_cvg) / sigma_target
  )

req_real <- spec %>%
  crossing(sigma_target = SIGMA_TARGETS) %>%
  mutate(
    bias_scenario   = "bias=allow",
    # TEa - bias_allow = k * cv_a_allow
    c_req           = (tea - bias_allow) / sigma_target,
    c_req_over_cv_a = c_req / cv_a_allow,
    c_req_over_cvi  = c_req / median.cvi,
    ratio_cv_a_expected = K_COVERAGE / sigma_target,
    ratio_cvi_expected  = K_COVERAGE * imp_mult / sigma_target
  )

req <- bind_rows(req, req_real) %>%
  mutate(
    # "строже целей по биологической вариации" = требование меньше допустимой неточности
    stricter_than_bv = c_req_over_cv_a < 1,
    # порог по CV_G/CV_I, при котором аналит перестаёт быть затронутым (bias = 0):
    # (k + 0.5*R)/S < 1  <=>  CV_G/CV_I < sqrt(4*(S-k)^2 - 1)
    threshold_cvg_over_cvi_bias0 = sqrt(pmax(4 * (sigma_target - K_COVERAGE)^2 - 1, 0)),
    cvg_over_cvi = median.cvg_used / median.cvi
  ) %>%
  arrange(analyte_name, level, bias_scenario, sigma_target)

write_csv_sc(spec, file.path(out_dir, "sigma_specs.csv"))
write_csv_sc(req,  file.path(out_dir, "sigma_requirements.csv"))

# --- контроль замкнутых форм -------------------------------------------------
chk_real <- req %>%
  filter(bias_scenario == "bias=allow") %>%
  mutate(d_cv_a = abs(c_req_over_cv_a - ratio_cv_a_expected),
         d_cvi  = abs(c_req_over_cvi  - ratio_cvi_expected))

chk_zero <- req %>%
  filter(bias_scenario == "bias=0") %>%
  mutate(d_cv_a = abs(c_req_over_cv_a - ratio_cv_a_expected),
         d_cvi  = abs(c_req_over_cvi  - ratio_cvi_expected))

chk_max <- max(c(chk_real$d_cv_a, chk_real$d_cvi, chk_zero$d_cv_a, chk_zero$d_cvi))

# --- чувствительность к коэффициенту охвата k --------------------------------
# CV_A и Bias от k не зависят, поэтому для любого k пересчитывается только TEa.
#   bias = допуск: CV_треб = k * CV_A / S  ->  отношение к CV_A = k / S;
#   bias = 0     : CV_треб = (k * CV_A + Bias) / S -> отношение = (k + 0.5*R) / S.
sens_base <- spec %>% filter(level == "desirable") %>%
  select(analyte_name, cv_a_allow, bias_allow, rel_cvg)

k_sens <- bind_rows(lapply(K_SET, function(k) {
  bind_rows(lapply(SIGMA_TARGETS, function(S) {
    a <- sens_base %>%
      transmute(k = k, sigma_target = S, bias_scenario = "bias=allow",
                c_req = k * cv_a_allow / S,
                ratio = k / S,
                threshold = NA_real_)
    z <- sens_base %>%
      transmute(k = k, sigma_target = S, bias_scenario = "bias=0",
                c_req = (k * cv_a_allow + bias_allow) / S,
                ratio = (k * cv_a_allow + bias_allow) / S / cv_a_allow,
                threshold = sqrt(max(4 * (S - k)^2 - 1, 0)))
    bind_rows(a, z)
  }))
}))

k_sens_summary <- k_sens %>%
  group_by(k, bias_scenario, sigma_target) %>%
  summarise(n_analytes = n(),
            n_stricter = sum(ratio < 1),
            pct_stricter = 100 * sum(ratio < 1) / n(),
            median_ratio = median(ratio),
            q1_ratio = quantile(ratio, 0.25),
            q3_ratio = quantile(ratio, 0.75),
            threshold_cvg_over_cvi = threshold[1],
            median_c_req = median(c_req),
            .groups = "drop") %>%
  arrange(k, bias_scenario, sigma_target) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write_csv_sc(k_sens_summary, file.path(out_dir, "k_sensitivity.csv"))

# --- аналиты, не затронутые расхождением (пороговая логика) -----------------
# Публикуются только производные величины: отношение CV_G/CV_I и отношение
# требования к допустимой неточности; сырые значения CV_I и CV_G не выводятся.
s2 <- bind_rows(lapply(SIGMA_TARGETS, function(S) {
  thr <- sqrt(max(4 * (S - K_COVERAGE)^2 - 1, 0))
  spec %>% filter(level == "desirable") %>%
    mutate(sigma_target  = S,
           threshold     = thr,
           cvg_over_cvi  = median.cvg_used / median.cvi,
           ratio_cv_a    = (K_COVERAGE + 0.5 * rel_cvg) / S,
           above_threshold = cvg_over_cvi > thr) %>%
    filter(above_threshold) %>%
    transmute(sigma_target, threshold = round(thr, 2), analyte_name,
              cvg_over_cvi = round(cvg_over_cvi, 2),
              ratio_cv_a = round(ratio_cv_a, 3))
}))

write_csv_sc(s2, file.path(out_dir, "table_s2_threshold_analytes.csv"))

# --- S7: устойчивость пороговых классификаций к интервалам входных оценок ----
# База EFLM приводит для каждой оценки не только медиану, но и интервал (lower/upper).
# Замкнутые формы от медиан не зависят, а СЧЁТНЫЕ классификации зависят: аналит
# считается не затронутым, когда CV_G/CV_I превышает порог sqrt(4(S-k)^2 - 1).
# Пересчёт: отношение (k + 0,5R)/S вычисляется при CV на границах интервалов,
# отдельно по CV_G, отдельно по CV_I и совместно по обеим оценкам.
R_of <- function(cvi, cvg) sqrt(1 + (cvg / cvi)^2)

s7 <- bind_rows(lapply(SIGMA_TARGETS, function(S) {
  # ratio_max — наиболее благоприятное для гипотезы «не затронут» сочетание границ
  # (CV_I на нижней, CV_G на верхней); ratio_min — противоположное.
  ratio_max <- (K_COVERAGE + 0.5 * R_of(w$cvi_lo, w$cvg_hi)) / S
  ratio_min <- (K_COVERAGE + 0.5 * R_of(w$cvi_hi, w$cvg_lo)) / S
  base      <- (K_COVERAGE + 0.5 * R_of(w$median.cvi, w$median.cvg_used)) / S
  cvg_up    <- (K_COVERAGE + 0.5 * R_of(w$median.cvi, w$cvg_hi)) / S
  cvg_lo_   <- (K_COVERAGE + 0.5 * R_of(w$median.cvi, w$cvg_lo)) / S
  cvi_up    <- (K_COVERAGE + 0.5 * R_of(w$cvi_hi, w$median.cvg_used)) / S
  cvi_lo_   <- (K_COVERAGE + 0.5 * R_of(w$cvi_lo, w$median.cvg_used)) / S
  robust_strict <- ratio_max < 1   # строже допустимого при любом сочетании границ
  robust_unaff  <- ratio_min >= 1  # не затронут при любом сочетании границ
  # запас до порога у аналитов, названных не затронутыми по медианам
  un <- which(base >= 1)
  margin <- if (length(un)) min(100 * (base[un] - 1)) else NA_real_
  tibble(
    sigma_target = S,
    threshold = round(sqrt(max(4 * (S - K_COVERAGE)^2 - 1, 0)), 2),
    n_analytes = nrow(w),
    n_stricter_base = sum(base < 1),
    n_stricter_cvg_at_upper = sum(cvg_up < 1),
    n_stricter_cvg_at_lower = sum(cvg_lo_ < 1),
    n_stricter_cvi_at_upper = sum(cvi_up < 1),
    n_stricter_cvi_at_lower = sum(cvi_lo_ < 1),
    n_stricter_joint_least = sum(ratio_min < 1),
    n_stricter_joint_most = sum(ratio_max < 1),
    n_robust_stricter = sum(robust_strict),
    n_robust_unaffected = sum(robust_unaff),
    n_classification_unstable = sum(!robust_strict & !robust_unaff),
    n_unaffected_base = length(un),
    min_margin_unaffected_pct = round(margin, 2)
  )
}))

write_csv_sc(s7, file.path(out_dir, "table_s7_threshold_sensitivity.csv"))

# --- публикуемый производный перечень по всем аналитам ------------------------
# Сырые CV_I и CV_G не публикуются (условия использования базы EFLM); приводятся
# только производные величины и метаданные о числе первичных исследований.
derived <- req %>%
  filter(level == "desirable", bias_scenario == "bias=0") %>%
  select(analyte_name, sigma_target, ratio_bias0 = c_req_over_cv_a,
         threshold = threshold_cvg_over_cvi_bias0) %>%
  pivot_wider(names_from = sigma_target,
              values_from = c(ratio_bias0, threshold),
              names_prefix = "sigma_") %>%
  left_join(
    req %>% filter(level == "desirable", bias_scenario == "bias=allow",
                   sigma_target == 4) %>%
      transmute(analyte_name, ratio_biasallow_sigma4 = c_req_over_cv_a),
    by = "analyte_name") %>%
  left_join(w %>% transmute(analyte_name,
                            cvg_over_cvi = round(median.cvg_used / median.cvi, 3),
                            number_used_cvi = number_used.cvi,
                            number_used_cvg = number_used.cvg,
                            cvg_imputed),
            by = "analyte_name") %>%
  mutate(across(starts_with("ratio_") | starts_with("threshold_"),
                ~ round(.x, 4))) %>%
  select(analyte_name, cvg_over_cvi,
         ratio_bias0_sigma_4 = ratio_bias0_sigma_4,
         ratio_bias0_sigma_5 = ratio_bias0_sigma_5,
         ratio_bias0_sigma_6 = ratio_bias0_sigma_6,
         threshold_sigma_4 = threshold_sigma_4,
         threshold_sigma_5 = threshold_sigma_5,
         threshold_sigma_6 = threshold_sigma_6,
         ratio_biasallow_sigma4,
         number_used_cvi, number_used_cvg, cvg_imputed) %>%
  arrange(analyte_name)

write_csv_sc(derived, file.path(out_dir, "table1_analytes_derived.csv"))

# --- перечень идентификаторов панели (для воспроизведения состава) ------------
# Публикуется полный состав панели, включая аналиты, не вошедшие в расчёт.
panel_ids <- bv %>%
  select(analyte_id, analyte_name, var_type, number_used, matrix, source) %>%
  pivot_wider(id_cols = c(analyte_id, analyte_name, matrix, source),
              names_from = var_type,
              values_from = number_used, names_prefix = "number_used_") %>%
  transmute(analyte_id, analyte_name, matrix, provenance = source,
            number_used_cvi = number_used_cvi,
            number_used_cvg = number_used_cvg,
            used_in_analysis = !(analyte_name %in% EXCLUDED_ANALYTES)) %>%
  arrange(provenance, analyte_name)

write_csv_sc(panel_ids, file.path(root, "data", "panel_analyte_ids.csv"))

# --- рабочий пример численной проверки (для текста статьи) -------------------
worked <- spec %>%
  filter(level == "desirable", analyte_name %in% c("Sodium", "C-reactive protein (CRP)")) %>%
  transmute(analyte_name, cv_i = median.cvi, cv_g = median.cvg_used,
            cv_a_allow, bias_allow, rel_cvg,
            ratio_bias0_sigma6 = (K_COVERAGE + 0.5 * rel_cvg) / 6,
            ratio_biasallow_sigma6 = K_COVERAGE / 6)

write_csv_sc(worked, file.path(out_dir, "worked_example.csv"))

# --- контроль регрессии основного расчёта (k = 1.65) -------------------------
reg <- k_sens_summary %>% filter(k == K_COVERAGE, bias_scenario == "bias=0")
main <- req %>% filter(level == "desirable", bias_scenario == "bias=0") %>%
  group_by(sigma_target) %>%
  summarise(n_stricter = sum(stricter_than_bv),
            median_ratio = median(c_req_over_cv_a), .groups = "drop")

# --- отчёт ------------------------------------------------------------------
cat("Аналитов собрано в панель  :", nrow(w_all), "\n")
cat("Исключено из расчёта       :", nrow(w_all) - nrow(w),
    "(", paste(setdiff(w_all$analyte_name, w$analyte_name), collapse = ", "), ")\n")
cat("Аналитов в расчёте         :", nrow(w), "\n")
cat("  с CVi и CVg              :", sum(!w$cvg_imputed), "\n")
cat("  с CVi, без CVg (CV_G=0)  :", sum(w$cvg_imputed), "\n")
cat("Матрицы:\n"); print(table(w$matrix))

cat("\nКонтроль замкнутых форм: max |ratio - expected| =",
    format(chk_max, scientific = TRUE), "\n")
cat("  (bias=allow: CV_треб/CV_A = k/S =", paste(round(K_COVERAGE / SIGMA_TARGETS, 4), collapse = " / "), ")\n")
cat("  (bias=0    : CV_треб/CV_A = (k + 0.5*R)/S)\n")

cat("\nПример расчёта (Glucose, desirable):\n")
print(spec %>% filter(analyte_name == "Glucose", level == "desirable") %>%
        select(median.cvi, median.cvg, cv_a_allow, bias_allow, tea, rel_cvg), row.names = FALSE)

cat("\nСводка по сценариям (уровень desirable, все сигмы):\n")
print(req %>%
        filter(level == "desirable") %>%
        group_by(bias_scenario, sigma_target) %>%
        summarise(n = n(),
                  n_stricter = sum(stricter_than_bv),
                  median_ratio = round(median(c_req_over_cv_a), 4),
                  q1 = round(quantile(c_req_over_cv_a, .25), 4),
                  q3 = round(quantile(c_req_over_cv_a, .75), 4),
                  .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

cat("\nСценарий bias=0: аналиты, НЕ требующие неточности строже CV_A (desirable):\n")
print(req %>%
        filter(bias_scenario == "bias=0", level == "desirable", !stricter_than_bv) %>%
        transmute(analyte_name, sigma_target,
                  cvg_over_cvi = round(cvg_over_cvi, 2),
                  threshold = round(threshold_cvg_over_cvi_bias0, 2),
                  ratio = round(c_req_over_cv_a, 3)) %>%
        arrange(sigma_target, desc(ratio)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nИнвариантность к уровню строгости (отношение к CV_A, bias=allow):\n")
print(req %>% filter(bias_scenario == "bias=allow") %>%
        group_by(sigma_target) %>%
        summarise(min_ratio = min(c_req_over_cv_a), max_ratio = max(c_req_over_cv_a),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)

cat("\nЗависимость от уровня строгости (отношение к CV_I, bias=allow) — для контраста:\n")
print(req %>% filter(bias_scenario == "bias=allow", sigma_target == 6) %>%
        group_by(level) %>%
        reframe(ratio_cvi_min = min(c_req_over_cvi), ratio_cvi_max = max(c_req_over_cvi)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nЧувствительность к коэффициенту охвата k:\n")
print(k_sens_summary %>%
        transmute(k, scenario = bias_scenario, sigma = sigma_target,
                  ratio = round(median_ratio, 3),
                  stricter = paste0(n_stricter, "/", n_analytes),
                  threshold = round(threshold_cvg_over_cvi, 2)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nКонтроль регрессии при k = 1.65 (bias = 0, уровень desirable):\n")
main_chk <- main %>% rename(n_stricter_main = n_stricter, median_main = median_ratio)
print(reg %>% select(sigma_target, n_stricter_k = n_stricter, median_k = median_ratio) %>%
        left_join(main_chk, by = "sigma_target") %>%
        mutate(match_n = n_stricter_k == n_stricter_main,
               match_median = abs(median_k - median_main) < 1e-4) %>%
        as.data.frame(), row.names = FALSE)

cat("\nАналиты выше порога CV_G/CV_I (не затронуты расхождением):\n")
print(s2 %>% as.data.frame(), row.names = FALSE)

cat("\nS7. Устойчивость пороговых классификаций к интервалам входных оценок:\n")
print(s7 %>% as.data.frame(), row.names = FALSE)

cat("\nЧисло первичных исследований (CV_I): не более 3 у",
    sum(w$number_used.cvi <= 3), "из", nrow(w), "; равно 1 у", sum(w$number_used.cvi == 1), "\n")

cat("\nРабочий пример для текста статьи:\n")
print(worked %>% as.data.frame(), row.names = FALSE)

cat("\nСформировано: out/k_sensitivity.csv, out/table_s2_threshold_analytes.csv,",
    "out/table_s7_threshold_sensitivity.csv, out/table1_analytes_derived.csv,",
    "data/panel_analyte_ids.csv, out/worked_example.csv\n")
