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

args_all <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
root <- normalizePath(file.path(here, ".."), mustWork = TRUE)
out_dir <- file.path(root, "out")

spec <- read.csv(file.path(out_dir, "sigma_specs.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
req  <- read.csv(file.path(out_dir, "sigma_requirements.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
fr   <- read.csv(file.path(out_dir, "qc_false_rejection.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
det  <- read.csv(file.path(out_dir, "qc_power.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")

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

write.csv(t1, file.path(out_dir, "table1_analytes_desirable.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

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

write.csv(t2, file.path(out_dir, "table2_requirement_vs_bv.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# --- Таблица S1: операционные характеристики схем контроля -------------------
t3 <- fr %>%
  select(scheme, n_per_run, p_per_run, p_lo, p_hi,
         arl_true, arl_true_lo, arl_true_hi, arl_1_over_pfr,
         exact_p_per_run, exact_arl, runs_arl, censored_frac) %>%
  left_join(det %>% select(scheme, shift_sd, p_detect_k, pd_lo, pd_hi) %>%
              pivot_wider(names_from = shift_sd,
                          values_from = c(p_detect_k, pd_lo, pd_hi),
                          names_sep = "_s"),
            by = "scheme") %>%
  mutate(
    across(c(p_per_run, p_lo, p_hi, exact_p_per_run), ~ round(.x, 6)),
    across(c(arl_true, arl_true_lo, arl_true_hi, arl_1_over_pfr, exact_arl), ~ round(.x, 1)),
    across(starts_with("p_detect_k_s"), ~ round(.x, 6)),
    across(starts_with("pd_lo_s"), ~ round(.x, 6)),
    across(starts_with("pd_hi_s"), ~ round(.x, 6)),
    censored_frac = round(censored_frac, 4)
  )

write.csv(t3, file.path(out_dir, "table_s1_qc_operating_characteristics.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

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

write.csv(t5, file.path(out_dir, "table_s5_level_strictness.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# --- Таблица S6: влияние начала наблюдения на ARL ----------------------------
t6 <- fr %>%
  transmute(scheme, n_per_run,
            arl_zero = round(arl_true, 2),
            arl_zero_lo = round(arl_true_lo, 2), arl_zero_hi = round(arl_true_hi, 2),
            arl_warm = round(arl_warm, 2),
            arl_warm_lo = round(arl_warm_lo, 2), arl_warm_hi = round(arl_warm_hi, 2),
            one_over_pfr = round(arl_1_over_pfr, 2),
            arl_warm_over_zero = round(arl_warm_over_zero, 4),
            n_chains_warm = n_chains_warm,
            warm_censored_frac = round(warm_censored_frac, 4))

write.csv(t6, file.path(out_dir, "table_s6_arl_start.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# устаревшие артефакты версии 2.0 (таблица 3 и рисунок 3 перенесены в дополнительные материалы)
for (f in c("table3_qc_operating_characteristics.csv", "fig3_qc_power.png")) {
  p <- file.path(out_dir, f)
  if (file.exists(p)) unlink(p)
}

# --- Рисунок 1: требуемая неточность против цели по биологической вариации ---
f1 <- req %>%
  filter(bias_scenario == "bias=0") %>%
  mutate(sigma_target = factor(sigma_target, levels = c(4, 5, 6)),
         level = factor(level, levels = LEVELS))

p1 <- ggplot(f1, aes(x = level, y = c_req, fill = sigma_target)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ sigma_target, nrow = 1,
             labeller = labeller(sigma_target = function(x) paste0("Sigma = ", x))) +
  labs(x = "Analytical performance level (Fraser hierarchy)",
       y = "Required method imprecision, CV (%)",
       fill = "Sigma target") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "fig1_required_imprecision.png"), p1, width = 8, height = 4.5, dpi = 300)

# --- Рисунок 2: отношение требуемой неточности к допустимой по BV -----------
# Вертикаль на 1 — граница «строже цели»; вертикаль на k/sigma — аналитический
# предел для сценария, в котором смещение находится на границе допуска.
f2 <- req %>%
  filter(bias_scenario == "bias=0") %>%
  mutate(sigma_target = factor(sigma_target, levels = c(4, 5, 6)),
         level = factor(level, levels = LEVELS))
lim <- data.frame(sigma_target = factor(c(4, 5, 6), levels = c(4, 5, 6)),
                  x = K_COVERAGE / c(4, 5, 6))

p2 <- ggplot(f2, aes(x = c_req_over_cv_a, fill = level)) +
  geom_histogram(bins = 30, alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(data = lim, aes(xintercept = x), linetype = "solid", linewidth = 0.4) +
  facet_wrap(~ sigma_target, nrow = 1,
             labeller = labeller(sigma_target = function(x) paste0("Sigma = ", x))) +
  labs(x = "Required imprecision / allowable imprecision (biological variation)",
       y = "Number of analytes", fill = "Level") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "fig2_ratio_to_bv.png"), p2, width = 9, height = 4.5, dpi = 300)

# --- Рисунок S1: операционные характеристики схем контроля ------------------
p3 <- ggplot(det, aes(x = shift_sd, y = p_detect_k, colour = scheme)) +
  geom_errorbar(aes(ymin = pd_lo, ymax = pd_hi), width = 0.06, linewidth = 0.35) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  labs(x = "Shift in method SD units",
       y = "Probability of detection within 100 runs",
       colour = "QC scheme") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "figS1_qc_power.png"), p3, width = 8.5, height = 5, dpi = 300)

# --- отчёт ------------------------------------------------------------------
cat("Сформировано:\n")
cat("  table1_analytes_desirable.csv                    —", nrow(t1), "аналитов\n")
cat("  table2_requirement_vs_bv.csv                     —", nrow(t2), "строк\n")
cat("  table_s1_qc_operating_characteristics.csv        —", nrow(t3), "строк\n")
cat("  table_s5_level_strictness.csv                    —", nrow(t5), "строк\n")
cat("  table_s6_arl_start.csv                           —", nrow(t6), "строк\n")
cat("  fig1_required_imprecision.png\n  fig2_ratio_to_bv.png\n  figS1_qc_power.png\n")

cat("\nТаблица 2 (все сценарии и уровни):\n")
print(t2 %>% select(bias_scenario, level, sigma_target, n_stricter, pct_stricter,
                    median_ratio, analytic_ratio, threshold_cvg_over_cvi) %>%
        as.data.frame(), row.names = FALSE)

cat("\nНаиболее строгие требования (сигма 6, bias=0, desirable):\n")
print(t1 %>% arrange(cv_required_sigma_6) %>%
        select(analyte_name, cvi, cvg, cv_a_allow, cv_required_sigma_6) %>% head(5),
      row.names = FALSE)
