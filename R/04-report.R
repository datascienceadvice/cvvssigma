# 04-report.R
# Сводные таблицы (CSV) и рисунки (PNG).
#
# Основной текст статьи: рисунки 1 и 2.
# Дополнительные материалы: таблица S1 (операционные характеристики) и рисунок S1.
#
# Результат: out/table1_analytes_desirable.csv, out/table2_requirement_vs_bv.csv,
#            out/table_s1_qc_operating_characteristics.csv,
#            out/fig1_required_imprecision.png, out/fig2_ratio_to_bv.png,
#            out/figS1_qc_power.png

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

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

spec <- read_csv_sc(file.path(out_dir, "sigma_specs.csv"))
req  <- read_csv_sc(file.path(out_dir, "sigma_requirements.csv"))
fr   <- read_csv_sc(file.path(out_dir, "qc_false_rejection.csv"))
det  <- read_csv_sc(file.path(out_dir, "qc_power.csv"))

K_COVERAGE <- 1.65
LEVELS <- c("optimal", "desirable", "minimum")

# --- Таблица 1: сводка по аналитам (уровень desirable) ----------------------
t1 <- spec %>%
  filter(level == "desirable") %>%
  select(analyte_name, matrix, cvi = median.cvi, cvg = median.cvg,
         n_cvi = number_used.cvi, n_cvg = number_used.cvg,
         cv_a_allow, bias_allow, tea, rel_cvg, cvg_imputed) %>%
  arrange(analyte_name)

# требование к неточности для сигмы 4/5/6 при bias = 0 (сценарий из таблицы 1 статьи)
r0 <- req %>% filter(bias_scenario == "bias=0", level == "desirable") %>%
  select(analyte_name, sigma_target, c_req) %>%
  pivot_wider(names_from = sigma_target, values_from = c_req,
              names_prefix = "cv_required_sigma_")

t1 <- t1 %>% left_join(r0, by = "analyte_name") %>%
  mutate(across(starts_with("cv_required"), ~ round(.x, 3)),
         across(c(cvi, cvg, cv_a_allow, bias_allow, tea, rel_cvg), ~ round(.x, 3)))

write_csv_sc(t1, file.path(out_dir, "table1_analytes_desirable.csv"))

# --- Таблица 2: соотношение требований сигмы и целей по биологической вариации
t2 <- req %>%
  group_by(bias_scenario, level, sigma_target) %>%
  summarise(n_analytes   = n(),
            n_stricter   = sum(stricter_than_bv),
            pct_stricter = 100 * sum(stricter_than_bv) / n(),
            median_ratio = median(c_req_over_cv_a),
            q1_ratio     = quantile(c_req_over_cv_a, 0.25),
            q3_ratio     = quantile(c_req_over_cv_a, 0.75),
            # аналитический предел отношения (замкнутая форма)
            analytic_ratio = ifelse(bias_scenario[1] == "bias=allow",
                                    K_COVERAGE / sigma_target[1],
                                    median((K_COVERAGE + 0.5 * rel_cvg) / sigma_target)),
            threshold_cvg_over_cvi = ifelse(bias_scenario[1] == "bias=allow",
                                            NA_real_,
                                            sqrt(max(4 * (sigma_target[1] - K_COVERAGE)^2 - 1, 0))),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write_csv_sc(t2, file.path(out_dir, "table2_requirement_vs_bv.csv"))

# --- Таблица S1: операционные характеристики схем контроля -------------------
# arl_mc — оценка ARL по Монте-Карло, устойчивая к цензурированию; arl_over_pfr —
# отношение к обратной стационарной частоте сигналов (для правил без памяти должно
# быть близко к 1, для правил с памятью сравнимо с ARL заполненного буфера, S6).
t3 <- fr %>%
  select(scheme, n_per_run, p_per_run, p_lo, p_hi,
         arl_mc, arl_mc_lo, arl_mc_hi, arl_1_over_pfr,
         exact_p_per_run, exact_arl, runs_arl, n_censored_arl, censored_frac) %>%
  left_join(det %>% select(scheme, shift_sd, p_detect_k, pd_lo, pd_hi) %>%
              pivot_wider(names_from = shift_sd,
                          values_from = c(p_detect_k, pd_lo, pd_hi),
                          names_sep = "_s"),
            by = "scheme") %>%
  mutate(
    across(c(p_per_run, p_lo, p_hi, exact_p_per_run), ~ round(.x, 6)),
    across(c(arl_mc, arl_mc_lo, arl_mc_hi, arl_1_over_pfr, exact_arl), ~ round(.x, 1)),
    arl_over_pfr = round(arl_mc / arl_1_over_pfr, 4),
    across(starts_with("p_detect_k_s"), ~ round(.x, 6)),
    across(starts_with("pd_lo_s"), ~ round(.x, 6)),
    across(starts_with("pd_hi_s"), ~ round(.x, 6)),
    censored_frac = round(censored_frac, 6)
  )

write_csv_sc(t3, file.path(out_dir, "table_s1_qc_operating_characteristics.csv"))

# --- Таблица S5: зависимость требований от уровня строгости -------------------
t5 <- req %>%
  filter(bias_scenario == "bias=0") %>%
  group_by(level, sigma_target) %>%
  summarise(n_analytes       = n(),
            n_stricter       = sum(stricter_than_bv),
            median_ratio_cv_a = median(c_req_over_cv_a),
            median_ratio_cvi  = median(c_req_over_cvi),
            median_c_req      = median(c_req),
            .groups = "drop") %>%
  arrange(match(level, c("optimal", "desirable", "minimum")), sigma_target) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write_csv_sc(t5, file.path(out_dir, "table_s5_level_strictness.csv"))

# --- Таблица S6: влияние начала наблюдения на ARL ----------------------------
t6 <- fr %>%
  transmute(scheme, n_per_run,
            arl_zero = round(arl_mc, 2),
            arl_zero_lo = round(arl_mc_lo, 2), arl_zero_hi = round(arl_mc_hi, 2),
            arl_warm = round(arl_warm, 2),
            arl_warm_lo = round(arl_warm_lo, 2), arl_warm_hi = round(arl_warm_hi, 2),
            one_over_pfr = round(arl_1_over_pfr, 2),
            arl_warm_over_zero = round(arl_warm_over_zero, 4),
            n_chains_warm = n_chains_warm,
            warm_censored_frac = round(warm_censored_frac, 6))

write_csv_sc(t6, file.path(out_dir, "table_s6_arl_start.csv"))

# устаревшие артефакты (переименованные рисунки предыдущих версий)
for (f in c("table3_qc_operating_characteristics.csv", "fig3_qc_power.png",
            "fig1_required_imprecision.png")) {
  p <- file.path(out_dir, f)
  if (file.exists(p)) unlink(p)
}

# --- Рисунок 1: аналитическая связь требований BV и сигма-метрики ------------
# Отношение CV_треб/CV_A как функция наблюдаемого смещения, выраженного в долях
# допустимого: ratio(beta) = (k + 0.5 * R * (1 - beta)) / sigma, где
# beta = |bias| / Bias_доп и R = sqrt(1 + (CV_G/CV_I)^2). При beta = 1 все
# аналиты сходятся в k/sigma (свойства аналита исчезают), при beta = 0 разброс
# максимален и определяется отношением CV_G/CV_I. Полоса — диапазон R по панели,
# линия — медиана R, пунктир — граница «строже цели» (ratio = 1), точечная
# горизонталь — аналитический предел k/sigma.
R_panel <- req %>% filter(level == "desirable", bias_scenario == "bias=0") %>%
  distinct(analyte_id, rel_cvg) %>% pull(rel_cvg)
R_rng <- range(R_panel); R_med <- median(R_panel)
cat("\nРисунок 1: R по панели — от", round(R_rng[1], 3), "до", round(R_rng[2], 3),
    ", медиана", round(R_med, 3), "\n")

g1 <- expand.grid(beta = seq(0, 1, by = 0.005), sigma_target = c(4, 5, 6)) %>%
  mutate(
    lo  = (K_COVERAGE + 0.5 * R_rng[1] * (1 - beta)) / sigma_target,
    hi  = (K_COVERAGE + 0.5 * R_rng[2] * (1 - beta)) / sigma_target,
    med = (K_COVERAGE + 0.5 * R_med    * (1 - beta)) / sigma_target,
    sigma_target = factor(sigma_target, levels = c(4, 5, 6))
  )
floor1 <- data.frame(sigma_target = factor(c(4, 5, 6), levels = c(4, 5, 6)),
                     y = K_COVERAGE / c(4, 5, 6))

p1 <- ggplot(g1, aes(x = beta)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
  geom_line(aes(y = med), linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_hline(data = floor1, aes(yintercept = y), linetype = "dotted", linewidth = 0.4) +
  facet_wrap(~ sigma_target, nrow = 1,
             labeller = labeller(sigma_target = function(x) paste0("Sigma = ", x))) +
  labs(x = "Observed bias as a fraction of the allowable bias, |bias| / Bias_allow",
       y = "Required imprecision / allowable imprecision  (CV_req / CV_A)") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "fig1_bv_vs_sigma.png"), p1, width = 9, height = 4, dpi = 300)

# --- Рисунок 2: распределение отношения по панели ---------------------------
# Три панели по целевой сигме; вертикаль ratio = 1 — граница «строже цели»,
# сплошная вертикаль — медиана по панели. Подписи дают число аналитов строже
# границы и медиану.
f2 <- req %>%
  filter(bias_scenario == "bias=0", level == "desirable") %>%
  mutate(sigma_target = factor(sigma_target, levels = c(4, 5, 6)))
ann2 <- f2 %>% group_by(sigma_target) %>%
  summarise(n = n(), k = sum(c_req_over_cv_a < 1),
            med = median(c_req_over_cv_a), .groups = "drop") %>%
  mutate(lab = sprintf("stricter: %d of %d\nmedian %.3f", k, n, med))

p2 <- ggplot(f2, aes(x = c_req_over_cv_a)) +
  geom_histogram(bins = 22, fill = "grey75", colour = "grey35", linewidth = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(data = ann2, aes(xintercept = med), linewidth = 0.4) +
  geom_text(data = ann2, aes(x = Inf, y = Inf, label = lab),
            hjust = 1.08, vjust = 1.15, size = 3) +
  facet_wrap(~ sigma_target, nrow = 1, scales = "free_y",
             labeller = labeller(sigma_target = function(x) paste0("Sigma = ", x))) +
  labs(x = "Required imprecision / allowable imprecision  (CV_req / CV_A)",
       y = "Number of analytes") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "fig2_ratio_to_bv.png"), p2, width = 9, height = 4, dpi = 300)

# --- Рисунок S1: операционные характеристики схем контроля ------------------
p3 <- ggplot(det, aes(x = shift_sd, y = p_detect_k, colour = scheme)) +
  geom_errorbar(aes(ymin = pd_lo, ymax = pd_hi), width = 0.06, linewidth = 0.35) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  labs(x = "Shift in method SD units",
       y = "Probability of detection within 100 runs",
       colour = "QC scheme") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "figS1_qc_power.png"), p3, width = 8.5, height = 5, dpi = 300)

# --- Рисунок S2: требуемая неточность по уровням строгости -------------------
fS2 <- req %>%
  filter(bias_scenario == "bias=0") %>%
  mutate(sigma_target = factor(sigma_target, levels = c(4, 5, 6)),
         level = factor(level, levels = LEVELS))

pS2 <- ggplot(fS2, aes(x = level, y = c_req)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey85") +
  facet_wrap(~ sigma_target, nrow = 1,
             labeller = labeller(sigma_target = function(x) paste0("Sigma = ", x))) +
  labs(x = "Analytical performance level (Fraser hierarchy)",
       y = "Required method imprecision, CV (%)") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "figS2_required_imprecision.png"), pS2, width = 8, height = 4.5, dpi = 300)

# --- отчёт ------------------------------------------------------------------
cat("Сформировано:\n")
cat("  table1_analytes_desirable.csv                    —", nrow(t1), "аналитов\n")
cat("  table2_requirement_vs_bv.csv                     —", nrow(t2), "строк\n")
cat("  table_s1_qc_operating_characteristics.csv        —", nrow(t3), "строк\n")
cat("  table_s5_level_strictness.csv                    —", nrow(t5), "строк\n")
cat("  table_s6_arl_start.csv                           —", nrow(t6), "строк\n")
cat("  fig1_bv_vs_sigma.png\n  fig2_ratio_to_bv.png\n  figS1_qc_power.png\n  figS2_required_imprecision.png\n")

cat("\nТаблица 2 (все сценарии и уровни):\n")
print(t2 %>% select(bias_scenario, level, sigma_target, n_stricter, pct_stricter,
                    median_ratio, analytic_ratio, threshold_cvg_over_cvi) %>%
        as.data.frame(), row.names = FALSE)

cat("\nНаиболее строгие требования (сигма 6, bias=0, desirable):\n")
print(t1 %>% arrange(cv_required_sigma_6) %>%
        select(analyte_name, cvi, cvg, cv_a_allow, cv_required_sigma_6) %>% head(5),
      row.names = FALSE)
