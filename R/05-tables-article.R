# 05-tables-article.R
# Формирует таблицы в готовом к вставке виде (markdown, русские заголовки,
# десятичная запятая): основные таблицы статьи и дополнительные материалы,
# а также сводку ключевых чисел для текста.
#
# Результат: out/article_tables.md, out/supplementary_tables.md,
#            out/article_numbers.txt

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

args_all <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
root <- normalizePath(file.path(here, ".."), mustWork = TRUE)
out_dir <- file.path(root, "out")

spec  <- read.csv(file.path(out_dir, "sigma_specs.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
req   <- read.csv(file.path(out_dir, "sigma_requirements.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
t2    <- read.csv(file.path(out_dir, "table2_requirement_vs_bv.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
s1    <- read.csv(file.path(out_dir, "table_s1_qc_operating_characteristics.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
s2    <- read.csv(file.path(out_dir, "table_s2_threshold_analytes.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
ks    <- read.csv(file.path(out_dir, "k_sensitivity.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
s5    <- read.csv(file.path(out_dir, "table_s5_level_strictness.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
s6    <- read.csv(file.path(out_dir, "table_s6_arl_start.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
est   <- read.csv(file.path(out_dir, "estimator_comparison.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
work  <- read.csv(file.path(out_dir, "worked_example.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
det   <- read.csv(file.path(out_dir, "qc_power.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
fr    <- read.csv(file.path(out_dir, "qc_false_rejection.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")

K_COVERAGE <- 1.65

# --- утилиты форматирования --------------------------------------------------
ru <- function(x, digits = 2) {
  formatC(round(as.numeric(x), digits), format = "f", digits = digits, decimal.mark = ",")
}
md_table <- function(df) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(hdr, sep, rows)
}
write_utf8 <- function(txt, path) {
  con <- file(path, open = "w", encoding = "UTF-8"); writeLines(txt, con); close(con)
}

# --- Таблица 1: представительные аналиты ------------------------------------
# Состав должен совпадать с таблицей 1 в manuscript.md; полный набор производных
# величин по всем 101 аналиту — в out/table1_analytes_desirable.csv.
rep_analytes <- c(
  "Sodium", "Potassium", "Chloride", "Calcium (Ca)", "Magnesium",
  "Glucose", "Creatinine", "Urea", "Albumin",
  "Triglycerides", "Cholesterol", "HDL cholesterol", "LDL Cholesterol",
  "Alanine transaminase (ALT)", "Aspartate transaminase (AST)",
  "g-glutamyl transferase activity", "C-reactive protein (CRP)",
  "Thyroid stimulating hormone (TSH)", "Thyroxine - free (FT4)",
  "Haemoglobin (Hb)", "Leukocytes", "Haemoglobin A1c (IFCC)"
)

r0 <- req %>% filter(bias_scenario == "bias=0", level == "desirable") %>%
  select(analyte_name, sigma_target, c_req) %>%
  pivot_wider(names_from = sigma_target, values_from = c_req, names_prefix = "cv")

t1 <- spec %>%
  filter(level == "desirable", analyte_name %in% rep_analytes) %>%
  left_join(r0, by = "analyte_name") %>%
  transmute(
    Аналит          = analyte_name,
    `CV_I, %`       = ru(median.cvi, 2),
    `CV_G, %`       = ifelse(is.na(median.cvg), "н/д", ru(median.cvg, 2)),
    `CV_A, %`       = ru(cv_a_allow, 2),
    `Смещение, %`   = ru(bias_allow, 2),
    `TEa, %`        = ru(tea, 2),
    `CV для 4σ, %`  = ru(cv4, 3),
    `CV для 5σ, %`  = ru(cv5, 3),
    `CV для 6σ, %`  = ru(cv6, 3),
    `CV_треб/CV_A при 6σ` = ru((K_COVERAGE + 0.5 * rel_cvg) / 6, 3)
  ) %>%
  arrange(match(Аналит, rep_analytes))

# --- Таблица 2А/2Б/2В --------------------------------------------------------
t2a <- t2 %>% filter(bias_scenario == "bias=allow", level == "desirable") %>%
  arrange(sigma_target) %>%
  transmute(
    `Целевая сигма`              = sigma_target,
    `Аналитов в панели`          = n_analytes,
    `Требуют CV строже целей BV` = n_stricter,
    `Доля, %`                    = ru(pct_stricter, 1),
    `Отношение CV_треб/CV_A`     = paste0(ru(analytic_ratio, 4), " = 1,65/σ")
  ) %>% mutate(across(everything(), as.character))

t2b <- t2 %>% filter(bias_scenario == "bias=0", level == "desirable") %>%
  arrange(sigma_target) %>%
  transmute(
    `Целевая сигма`              = sigma_target,
    `Аналитов в панели`          = n_analytes,
    `Требуют CV строже целей BV` = n_stricter,
    `Доля, %`                    = ru(pct_stricter, 1),
    `Медиана отношения`          = ru(median_ratio, 3),
    `Межквартильный размах`      = paste0(ru(q1_ratio, 3), "–", ru(q3_ratio, 3))
  ) %>% mutate(across(everything(), as.character))

t2c <- NULL  # таблица 2В перенесена в текст статьи (в прозе), см. STATYA.md

# --- Таблица S1: операционные характеристики --------------------------------
s1md <- s1 %>%
  transmute(
    Схема  = scheme,
    `Pfr (95 % ДИ)` = sprintf("%s (%s–%s)", ru(p_per_run, 5), ru(p_lo, 5), ru(p_hi, 5)),
    `ARL, серий (95 % ДИ)` = sprintf("%s (%s–%s)", ru(arl_true, 1), ru(arl_true_lo, 1), ru(arl_true_hi, 1)),
    `1/Pfr` = ru(arl_1_over_pfr, 1),
    `Цензурировано` = ru(censored_frac, 4),
    `Ped 0,5 SD` = ru(p_detect_k_s0.5, 4),
    `Ped 1 SD`   = ru(p_detect_k_s1, 5),
    `Ped 1,5 SD` = ru(p_detect_k_s1.5, 5),
    `Ped 2 SD`   = ru(p_detect_k_s2, 5),
    `Ped 3 SD`   = ru(p_detect_k_s3, 5)
  ) %>% mutate(across(everything(), as.character))

# --- Таблица S2: пороговая логика --------------------------------------------
s2md <- bind_rows(lapply(sort(unique(s2$sigma_target)), function(S) {
  d <- s2 %>% filter(sigma_target == S) %>% arrange(desc(cvg_over_cvi))
  tibble(
    `Целевая сигма` = as.character(S),
    `Порог CV_G/CV_I` = ru(d$threshold[1], 2),
    `Число аналитов выше порога` = as.character(nrow(d)),
    `Аналиты (CV_G/CV_I; CV_треб/CV_A)` = paste(
      sprintf("%s (%s; %s)", d$analyte_name, ru(d$cvg_over_cvi, 2), ru(d$ratio_cv_a, 3)),
      collapse = "; ")
  )
}))

# --- Таблица S3: чувствительность к k ---------------------------------------
s3a <- ks %>% filter(bias_scenario == "bias=allow") %>%
  select(k, sigma_target, median_ratio) %>%
  pivot_wider(names_from = sigma_target, values_from = median_ratio,
              names_prefix = "sigma_") %>%
  transmute(`Коэффициент охвата k` = ru(k, 2),
            `Отношение k/σ при σ = 4` = ru(sigma_4, 4),
            `σ = 5` = ru(sigma_5, 4),
            `σ = 6` = ru(sigma_6, 4),
            `Основной расчёт` = ifelse(abs(k - K_COVERAGE) < 1e-9, "да", ""))

s3b <- ks %>% filter(bias_scenario == "bias=0") %>%
  select(k, sigma_target, n_stricter, n_analytes, median_ratio, threshold_cvg_over_cvi) %>%
  pivot_wider(names_from = sigma_target,
              values_from = c(n_stricter, n_analytes, median_ratio, threshold_cvg_over_cvi),
              names_sep = "_s") %>%
  transmute(`Коэффициент охвата k` = ru(k, 2),
            `Затронуто при σ = 4` = paste0(n_stricter_s4, " из ", n_analytes_s4),
            `σ = 5` = paste0(n_stricter_s5, " из ", n_analytes_s5),
            `σ = 6` = paste0(n_stricter_s6, " из ", n_analytes_s6),
            `Медиана отношения при σ = 4` = ru(median_ratio_s4, 3),
            `Порог CV_G/CV_I при σ = 4` = ru(threshold_cvg_over_cvi_s4, 2),
            `Основной расчёт` = ifelse(abs(k - K_COVERAGE) < 1e-9, "да", ""))
s3a <- s3a %>% mutate(across(everything(), as.character))
s3b <- s3b %>% mutate(across(everything(), as.character))

# --- Таблица S4: сравнение оценивателей Pfr ---------------------------------
s4md <- est %>%
  transmute(
    Схема = scheme,
    `Стационарная частота Pfr` = ru(p_stationary, 5),
    `Оцениватель версии 1.0`   = ru(legacy_p, 5),
    `Отношение (1.0 / стационарная)` = ru(legacy_over_stationary, 3),
    `Точное значение (1-3s)`   = ifelse(is.na(exact_p), "—", ru(exact_p, 5)),
    `Отношение (1.0 / точное)` = ifelse(is.na(legacy_over_exact), "—", ru(legacy_over_exact, 3))
  ) %>% mutate(across(everything(), as.character))

# --- Таблица S5: зависимость от уровня строгости ------------------------------
s5md <- s5 %>%
  arrange(match(level, c("optimal", "desirable", "minimum")), sigma_target) %>%
  transmute(
    `Уровень строгости`          = level,
    `Целевая сигма`              = sigma_target,
    `Требуют CV строже целей BV` = paste0(n_stricter, " из ", n_analytes),
    `Медиана CV_треб/CV_A`       = ru(median_ratio_cv_a, 3),
    `Медиана CV_треб/CV_I`       = ru(median_ratio_cvi, 3),
    `Медиана CV_треб, %`         = ru(median_c_req, 2)
  ) %>% mutate(across(everything(), as.character))

# --- Таблица S6: влияние начала наблюдения на ARL ----------------------------
s6md <- s6 %>%
  transmute(
    Схема = scheme,
    `ARL без предыстории (95 % ДИ)` =
      sprintf("%s (%s–%s)", ru(arl_zero, 1), ru(arl_zero_lo, 1), ru(arl_zero_hi, 1)),
    `ARL при заполненном буфере (95 % ДИ)` =
      sprintf("%s (%s–%s)", ru(arl_warm, 1), ru(arl_warm_lo, 1), ru(arl_warm_hi, 1)),
    `Отношение` = ru(arl_warm_over_zero, 3),
    `1/Pfr`     = ru(one_over_pfr, 1),
    `Цепочек с прогретым буфером` = n_chains_warm
  ) %>% mutate(across(everything(), as.character))

# --- ключевые числа для текста ----------------------------------------------
d0  <- t2 %>% filter(bias_scenario == "bias=0", level == "desirable")
dA  <- t2 %>% filter(bias_scenario == "bias=allow", level == "desirable")
n_an <- d0$n_analytes[1]
g <- function(df, s, col) df[[col]][df$sigma_target == s]

s2_4 <- s2 %>% filter(sigma_target == 4) %>% arrange(desc(cvg_over_cvi))
s2_5 <- s2 %>% filter(sigma_target == 5)
s2_6 <- s2 %>% filter(sigma_target == 6)

k_main <- ks %>% filter(bias_scenario == "bias=0", sigma_target == 4) %>% arrange(k)
k_allow <- ks %>% filter(bias_scenario == "bias=allow", sigma_target == 4) %>% arrange(k)

lvl <- req %>%
  filter(bias_scenario == "bias=0", sigma_target == 6) %>%
  group_by(level) %>%
  summarise(ratio_cv_a = median(c_req_over_cv_a),
            ratio_cvi  = median(c_req_over_cvi),
            c_req      = median(c_req), .groups = "drop") %>%
  arrange(match(level, c("optimal", "desirable", "minimum")))

ped2 <- det %>% filter(shift_sd == 2)
ped3 <- det %>% filter(shift_sd == 3)

nums <- c(
  sprintf("Аналитов в панели: %d (из них с CV_G: %d)",
          n_an, sum(!spec$cvg_imputed[spec$level == "desirable"])),
  sprintf("Контроль замкнутых форм: max |отношение - ожидаемое| = %s",
          format(max(abs(req$c_req_over_cv_a - req$ratio_cv_a_expected)), scientific = TRUE)),
  "",
  "Сценарий bias = допуск (аналитический результат):",
  sprintf("  CV_треб/CV_A = k/sigma = %s (4), %s (5), %s (6); строже целей BV %d из %d при любой сигме",
          ru(1.65/4, 4), ru(1.65/5, 4), ru(1.65/6, 4), dA$n_stricter[1], n_an),
  "  условие расхождения: sigma > k, то есть при sigma = 4-6 выполняется при любом k < 4",
  "",
  "Сценарий bias = 0 (идеализированный), уровень desirable:",
  sprintf("  сигма 4: строже целей BV %d из %d (%s %%), медиана %s (МКР %s-%s)",
          g(d0,4,"n_stricter"), n_an, ru(g(d0,4,"pct_stricter"),1),
          ru(g(d0,4,"median_ratio"),3), ru(g(d0,4,"q1_ratio"),3), ru(g(d0,4,"q3_ratio"),3)),
  sprintf("  сигма 5: строже целей BV %d из %d (%s %%), медиана %s (МКР %s-%s)",
          g(d0,5,"n_stricter"), n_an, ru(g(d0,5,"pct_stricter"),1),
          ru(g(d0,5,"median_ratio"),3), ru(g(d0,5,"q1_ratio"),3), ru(g(d0,5,"q3_ratio"),3)),
  sprintf("  сигма 6: строже целей BV %d из %d (%s %%), медиана %s (МКР %s-%s)",
          g(d0,6,"n_stricter"), n_an, ru(g(d0,6,"pct_stricter"),1),
          ru(g(d0,6,"median_ratio"),3), ru(g(d0,6,"q1_ratio"),3), ru(g(d0,6,"q3_ratio"),3)),
  sprintf("  пороги CV_G/CV_I: %s (4), %s (5), %s (6)",
          ru(g(d0,4,"threshold_cvg_over_cvi"),2), ru(g(d0,5,"threshold_cvg_over_cvi"),2),
          ru(g(d0,6,"threshold_cvg_over_cvi"),2)),
  sprintf("  не затронуты при 4 sigma (%d, исчерпывающе): %s",
          nrow(s2_4), paste(sprintf("%s (%.2f)", s2_4$analyte_name, s2_4$cvg_over_cvi), collapse = "; ")),
  sprintf("  не затронуты при 5 sigma (%d): %s", nrow(s2_5), paste(s2_5$analyte_name, collapse = "; ")),
  sprintf("  не затронуты при 6 sigma (%d): %s", nrow(s2_6), paste(s2_6$analyte_name, collapse = "; ")),
  "",
  "Роль уровня строгости (сигма 6, bias = 0):",
  paste0("  медиана CV_треб/CV_A = ", paste(ru(lvl$ratio_cv_a, 3), collapse = " / "),
         " (optimal / desirable / minimum)"),
  paste0("  медиана CV_треб/CV_I = ", paste(ru(lvl$ratio_cvi, 3), collapse = " / ")),
  paste0("  медиана CV_треб = ", paste(ru(lvl$c_req, 2), collapse = " / "), " %"),
  "",
  "Чувствительность к коэффициенту охвата k (sigma = 4, bias = 0):",
  paste0("  затронуто аналитов: ",
         paste(sprintf("k = %s -> %d из %d", ru(k_main$k, 2), k_main$n_stricter, k_main$n_analytes),
               collapse = "; ")),
  paste0("  медиана отношения: ",
         paste(sprintf("k = %s -> %s", ru(k_main$k, 2), ru(k_main$median_ratio, 3)),
               collapse = "; ")),
  paste0("  при bias = допуск отношение k/sigma при sigma = 4: ",
         paste(sprintf("k = %s -> %s", ru(k_allow$k, 2), ru(k_allow$median_ratio, 4)),
               collapse = "; "), " (все < 1)"),
  "",
  "Численная проверка (сигма 6):",
  paste0("  ", work$analyte_name, ": CV_A = ", ru(work$cv_a_allow, 3),
         "; R = ", ru(work$rel_cvg, 3),
         "; (1,65 + 0,5R)/6 = ", ru(work$ratio_bias0_sigma6, 4),
         "; при смещении на границе допуска k/6 = ", ru(work$ratio_biasallow_sigma6, 4)),
  "",
  "Операционные характеристики:",
  paste0("  точные значения 1-3s: Pfr(N=1) = ", ru(fr$exact_p_per_run[fr$scheme=="1-3s (N=1)"],5),
         " (ARL ", ru(fr$exact_arl[fr$scheme=="1-3s (N=1)"],1), "), Pfr(N=2) = ",
         ru(fr$exact_p_per_run[fr$scheme=="1-3s (N=2)"],5), " (ARL ",
         ru(fr$exact_arl[fr$scheme=="1-3s (N=2)"],1), ")"),
  "  истинная ARL (zero-state, время до первого сигнала от начала контроля): ",
  paste0("    ",
         paste(sprintf("%s = %s (%s-%s)", fr$scheme, ru(fr$arl_true,1), ru(fr$arl_true_lo,1),
                       ru(fr$arl_true_hi,1)), collapse = "; ")),
  paste0("  ARL при заполненном буфере (остаточное время после прогрева 20 серий): ",
         paste(sprintf("%s = %s", fr$scheme, ru(fr$arl_warm,1)), collapse = "; ")),
  paste0("  максимальное расхождение ARL при прогреве и zero-state ARL: ",
         ru(100 * max(abs(fr$arl_warm_over_zero - 1)), 1), " %"),
  paste0("  отношение истинной ARL к 1/Pfr: от ",
         ru(min(fr$arl_true_over_geom), 3), " до ", ru(max(fr$arl_true_over_geom), 3),
         " (расхождение не более ", ru(100 * max(abs(fr$arl_true_over_geom - 1)), 1), " %)"),
  paste0("  доля цензурированных наблюдений: не более ", ru(max(fr$censored_frac), 4)),
  paste0("  оцениватель версии 1.0 против стационарной частоты: ",
         paste(sprintf("%s = %s", est$scheme, ru(est$legacy_over_stationary, 3)), collapse = "; ")),
  paste0("  Ped при сдвиге 1 SD: ",
         paste(sprintf("%s = %s", det$scheme[det$shift_sd==1], ru(det$p_detect_k[det$shift_sd==1],3)),
               collapse = "; ")),
  paste0("  Ped при сдвиге 2 SD: минимум по схемам ", ru(min(ped2$p_detect_k), 5),
         " (нижняя граница 95 % ДИ не ниже ", ru(min(ped2$pd_lo), 5), ")"),
  paste0("  Ped при сдвиге 3 SD: нижняя граница 95 % ДИ не ниже ", ru(min(ped3$pd_lo), 5)),
  "",
  "Различие сценариев объясняется слагаемым смещения, пропорциональным sqrt(CV_I^2 + CV_G^2)."
)

# --- сборка файлов -----------------------------------------------------------
md <- c(
  "# Таблицы для статьи (готовы к вставке)",
  "",
  "Источник данных: EFLM Biological Variation Database, https://biologicalvariation.eu/ (дата доступа 25.09.2026).",
  "Все величины, кроме CV_I и CV_G, рассчитаны автором. Десятичный разделитель — запятая.",
  "",
  "## Таблица 1. Аналитические цели и требования к неточности метода (уровень desirable, сценарий bias = 0)",
  "",
  md_table(t1),
  "",
  "Примечание. CV_I — внутрииндивидуальная, CV_G — межиндивидуальная биологическая вариация;",
  "CV_A — допустимая аналитическая неточность; TEa — допустимая суммарная ошибка;",
  "CV для 4σ/5σ/6σ — максимальная неточность метода, при которой достигается соответствующая сигма при нулевом смещении;",
  "CV_треб/CV_A при 6σ — отношение требования сигма-метрики к допустимой неточности по биологической вариации",
  "(аналитически равно (1,65 + 0,5R)/6, где R = sqrt(1 + (CV_G/CV_I)^2)). н/д — значение CV_G в базе отсутствует.",
  "Полный набор производных величин по всем 101 аналиту — в файле out/table1_analytes_desirable.csv;",
  "исчерпывающие перечни аналитов, не затронутых расхождением, — в таблице S2 дополнительных материалов.",
  "",
  "## Таблица 2А. Соотношение требований при смещении на границе допуска (аналитический результат)",
  "",
  md_table(t2a),
  "",
  "Примечание. В этом сценарии TEa − |смещение| = k × CV_A, поэтому отношение требований не зависит",
  "ни от аналита, ни от уровня строгости и равно k/σ. Условие «строже цели по биологической вариации»",
  "выполняется для всех аналитов при любой σ > k, то есть при σ = 4–6 — для любого k < 4.",
  "",
  "## Таблица 2Б. Соотношение требований при нулевом смещении (идеализированный сценарий, уровень desirable)",
  "",
  md_table(t2b),
  "",
  "Примечание. «Медиана отношения» — медиана величины CV_требуемая / CV_A по панели из 101 аналита;",
  "значение меньше единицы означает, что достижение данной сигмы требует неточности строже допустимой",
  "по биологической вариации. Отношение (k + 0,5R)/σ зависит от отношения CV_G/CV_I, но не от уровня",
  "строгости; аналит перестаёт быть затронутым при CV_G/CV_I > sqrt(4(σ − k)^2 − 1), то есть 4,59 при",
  "σ = 4, 6,62 при σ = 5 и 8,64 при σ = 6. Исчерпывающие перечни таких аналитов — в таблице S2.",
  ""
)
write_utf8(md, file.path(out_dir, "article_tables.md"))

smd <- c(
  "# Дополнительные материалы",
  "",
  "Источник данных: EFLM Biological Variation Database, https://biologicalvariation.eu/ (дата доступа 25.09.2026).",
  "Все величины, кроме CV_I и CV_G, рассчитаны автором. Десятичный разделитель — запятая.",
  "",
  "## Таблица S1. Операционные характеристики схем внутрилабораторного контроля",
  "",
  md_table(s1md),
  "",
  "Примечание. Pfr — стационарная вероятность сигнала на одну серию (доля серий со срабатыванием на",
  "цепочках, достигших стационарного режима: 100 000 цепочек по 50 серий, прогрев 20 серий исключён);",
  "95 % ДИ — по разбросу долей между цепочками. ARL — истинное среднее число серий от начала контроля",
  "без предшествующей истории до первого сигнала (10 000 цепочек, длина подобрана так, чтобы доля",
  "цензурированных наблюдений была порядка 0,1–0,4 %); для правил без памяти (1-3s) ARL = 1/Pfr.",
  "Столбец «1/Pfr» приведён для сопоставления: это обратная стационарная частота сигналов, а не",
  "определение ARL; влияние начала наблюдения на ARL разобрано в таблице S6. Ped — вероятность обнаружения сдвига в пределах 100 серий после его появления",
  "(20 000 цепочек на клетку, сдвиг вводится после 20 серий контрольного состояния), 95 % ДИ Уилсона.",
  "Для схем из одного правила 1-3s точные значения: Pfr = 0,002700 при ARL 370,4 (N = 1) и",
  "Pfr = 0,005392 при ARL 185,4 (N = 2); оценки Монте-Карло попадают в 95 % ДИ этих значений.",
  "Рисунок S1 — out/figS1_qc_power.png.",
  "",
  "## Таблица S2. Аналиты, не затронутые расхождением (пороговая логика)",
  "",
  md_table(s2md),
  "",
  "Примечание. Приведены все аналиты панели, для которых отношение CV_G/CV_I превышает порог",
  "sqrt(4(σ − 1,65)^2 − 1); перечни исчерпывающие. Для этих аналитов требование сигма-метрики",
  "не строже допустимой неточности по биологической вариации (отношение CV_треб/CV_A ≥ 1) при",
  "нулевом смещении. Публикуются только производные величины: сырые значения CV_I и CV_G в",
  "дополнительных материалах не приводятся.",
  "",
  "## Таблица S3А. Чувствительность к коэффициенту охвата k: смещение на границе допуска",
  "",
  md_table(s3a),
  "",
  "## Таблица S3Б. Чувствительность к коэффициенту охвата k: нулевое смещение (σ = 4, 5, 6)",
  "",
  md_table(s3b),
  "",
  "Примечание. k — коэффициент охвата в выражении TEa = k × CV_A + смещение (1,0 ≈ 84 %,",
  "1,65 ≈ 95 %, 1,96 ≈ 97,5 %, 2,0 ≈ 97,7 %, 2,58 ≈ 99,5 % одностороннего охвата). При смещении",
  "на границе допуска отношение требований равно k/σ; при нулевом смещении — (k + 0,5R)/σ.",
  "Полная сводка — out/k_sensitivity.csv.",
  "",
  "## Таблица S4. Сравнение оценивателей вероятности ложного сигнала",
  "",
  md_table(s4md),
  "",
  "Примечание. «Оцениватель версии 1.0» — доля репликатов, в которых срабатывание произошло хотя бы",
  "раз за окно из 50 серий, делённая на 50. Он систематически смещён вниз примерно в",
  "[1 − (1 − p)^50] / (50p) раз и при больших p ограничен сверху величиной 1/50; в расчёте",
  "используется стационарная частота сигналов. Для схемы 1-3s (N = 1) смещение составляет ≈ −5 %,",
  "для полного набора правил при четырёх измерениях — ≈ −58 %.",
  "",
  "## Таблица S5. Зависимость требований от уровня строгости (сценарий bias = 0)",
  "",
  md_table(s5md),
  "",
  "Примечание. Для каждого уровня строгости и каждой целевой сигмы приведены медианы по панели",
  "из 101 аналита. Отношение CV_треб/CV_A одинаково на всех уровнях, поскольку допустимое смещение",
  "пропорционально допустимой неточности с коэффициентом 0,5; отношение к CV_I и абсолютная",
  "величина требуемой неточности пропорциональны множителю уровня m_imp.",
  "",
  "## Таблица S6. Влияние начала наблюдения на среднюю длину серии до сигнала",
  "",
  md_table(s6md),
  "",
  "Примечание. «ARL без предыстории» — zero-state ARL: среднее число серий от начала контроля",
  "(буфер правил пуст) до первого сигнала. «ARL при заполненном буфере» — среднее остаточное",
  "время до сигнала среди цепочек, не давших сигнала за первые 20 серий (буфер правил с памятью",
  "к этому моменту заполнен); счёт остаточного времени начинается заново, поэтому величина",
  "сопоставима с zero-state ARL. Расхождение между определениями не превышает 1 %, то есть",
  "выбор начала наблюдения на результат практически не влияет. Оба значения получены на 10 000",
  "цепочках; доля цензурированных наблюдений не превышает 0,5 %.",
  ""
)
write_utf8(smd, file.path(out_dir, "supplementary_tables.md"))

write_utf8(nums, file.path(out_dir, "article_numbers.txt"))

cat("Сформировано: out/article_tables.md, out/supplementary_tables.md, out/article_numbers.txt\n\n")
cat(paste(nums, collapse = "\n"), "\n")
