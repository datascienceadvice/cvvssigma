# 03-power.R
# Операционные характеристики схем внутрилабораторного контроля на основе правил
# Вестгарда: вероятность ложного сигнала на серию (Pfr), средняя длина серии до
# сигнала (ARL) и вероятность обнаружения сдвига в пределах K серий после его
# появления (Ped).
#
# Обозначения: z — результат контрольного измерения в единицах стандартного
# отклонения метода; сдвиг задаётся в тех же единицах (1.0 = одно SD).
#
# Правила: 1-3s, 2-2s, R4s (в пределах серии), 4-1s, 10x.
# Источник правил: Westgard JO, Barry PL, Hunt MR, Groth T. Clin Chem. 1981;27(3):493-501.
#
# Методика (см. также METODIKA.md, раздел 5).
#
# 1) Pfr. Оценивается как доля СЕРИЙ, в которых схема срабатывает, на длинных
#    цепочках серий, достигших стационарного режима (прогрев burn-in серий
#    исключается из подсчёта — это нужно правилам с памятью: 2-2s, 4-1s, 10x).
#    Оцениватель доли несмещён; 95% ДИ строится по разбросу долей между
#    независимыми цепочками (t-распределение), что учитывает корреляцию серий
#    внутри одной цепочки.
#    Прежняя версия оценивала Pfr как (число репликатов со срабатыванием) /
#    (число репликатов x число серий). Такой оцениватель систематически смещён
#    вниз примерно в [1 - (1-p)^K] / (K*p) раз (при p = 0.0027 и K = 50 это
#    около -6.4%), а при больших p упирается в потолок 1/K; он здесь не
#    используется.
# 2) ARL. Оценивается напрямую как среднее число серий от начала контроля до
#    первого сигнала (zero-state ARL) и отдельно в установившемся режиме
#    (после прогрева B_WARM серий). Оцениватель УСТОЙЧИВ К ЦЕНЗУРИРОВАНИЮ:
#    цепочки, не давшие сигнала за L серий, не отбрасываются, а получают
#    подстановочное значение L + 1/p (для прогретого режима (L - B_WARM) + 1/p),
#    то есть ожидаемое остаточное время. Прежняя версия усредняла только
#    цепочки со сигналом и потому систематически занижала ARL: отброшенными
#    оказывались самые длинные наблюдения. Горизонт L = ARL_RUNS_MULT / Pfr
#    выбран так, чтобы доля цензурирования была порядка 1e-5.
#    ДИ строятся по t-распределению между цепочками.
#    ARL = 1 / Pfr приводится рядом как отдельная величина; для правил без
#    памяти она совпадает с zero-state ARL, для правил с памятью сопоставима с
#    ARL заполненного буфера.
# 3) Ped. Моделируется появление сдвига: цепочка начинается в контрольном
#    состоянии (B_RUNS серий), затем сдвиг вводится и сохраняется; фиксируется,
#    сработала ли схема в пределах K_RUNS серий после появления сдвига.
#    Такой расчёт корректен для правил с памятью: перед сдвигом накоплена
#    «внутриконтрольная» история, а не пустая.
# 4) Для правила 1-3s известны точные значения: вероятность срабатывания за одну
#    серию p(s) = 1 - (1 - p1(s))^N, где p1(s) = Phi(s - 3) + Phi(-3 - s) для
#    одного измерения; ARL = 1/p(s); Ped(K) = 1 - (1 - p(s))^K. Они приводятся
#    рядом с Монте-Карло как контроль корректности.
#
# Воспроизводимость: зерно генератора зафиксировано, параметры вынесены в блок
# ниже; повторный запуск даёт идентичные результаты.

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

args_all <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
root <- normalizePath(file.path(here, ".."), mustWork = TRUE)
out_dir <- file.path(root, "out")

# --- параметры моделирования ------------------------------------------------
SEED          <- 20260925
N_REP_FR      <- 100000L   # цепочек для оценки Pfr
RUNS_FR       <- 50L       # серий в цепочке (Pfr)
BURN_IN       <- 20L       # прогрев, исключаемый из подсчёта Pfr
BLOCK_FR      <- 5000L     # размер блока цепочек (Pfr)
N_REP_DET     <- 20000L    # цепочек на клетку (обнаружение сдвига)
K_RUNS        <- 100L      # окно обнаружения сдвига, серий
B_RUNS        <- 20L       # серий контрольного состояния до появления сдвига
BLOCK_DET     <- 2000L     # размер блока цепочек (обнаружение)
N_REP_ARL     <- 10000L    # цепочек для оценки ARL по Монте-Карло (время до 1-го сигнала)
BLOCK_ARL     <- 1000L     # размер блока цепочек (ARL)
ARL_MIN_RUNS  <- 300L      # минимальная длина цепочки при оценке ARL
ARL_MAX_RUNS  <- 6000L     # максимальная длина цепочки при оценке ARL
ARL_RUNS_MULT <- 12        # длина цепочки = ARL_RUNS_MULT / Pfr (цензурирование ~ e^-12)
B_WARM        <- 20L       # прогрев для оценки ARL в установившемся режиме
SHIFTS        <- c(0.5, 1, 1.5, 2, 3)
CONF          <- 0.95

set.seed(SEED)

# Схемы контроля: имя -> (число контрольных измерений в серии, набор правил)
schemes <- tibble::tribble(
  ~scheme,          ~n_per_run, ~rules,
  "1-3s (N=1)",              1L, "1_3s",
  "1-3s (N=2)",              2L, "1_3s",
  "1-3s / 2-2s (N=2)",       2L, "1_3s,2_2s",
  "1-3s / 2-2s / R4s (N=2)", 2L, "1_3s,2_2s,R_4s",
  "multi-rule (N=2)",        2L, "1_3s,2_2s,R_4s,4_1s",
  "multi-rule (N=4)",        4L, "1_3s,2_2s,R_4s,4_1s,10x"
)

# --- индикаторы срабатывания правил по позициям измерений --------------------
# Возвращает логическую матрицу той же размерности, что Z: TRUE в позиции, на
# которой (с учётом предыдущих позиций) срабатывает хотя бы одно из правил.
# Семантика правил:
#   1-3s  — |z| > 3;
#   2-2s  — два последовательных измерения по одну сторону от нуля и оба > 2 SD
#           (как внутри серии, так и на границе серий);
#   R4s   — размах внутри одной серии > 4 SD;
#   4-1s  — четыре последовательных измерения по одну сторону и все > 1 SD;
#   10x   — десять последовательных измерений по одну сторону от нуля.
rule_indicators <- function(Z, n_per_run, rules) {
  nb <- nrow(Z); M <- ncol(Z)
  ind <- matrix(FALSE, nb, M)
  S <- sign(Z)
  nz <- S != 0

  if ("1_3s" %in% rules) ind <- ind | (abs(Z) > 3)

  if ("2_2s" %in% rules && M >= 2L) {
    A <- abs(Z) > 2
    j <- 2:M; i <- j - 1L
    w <- A[, j, drop = FALSE] & A[, i, drop = FALSE] &
         (S[, j, drop = FALSE] == S[, i, drop = FALSE]) & nz[, j, drop = FALSE]
    ind[, j] <- ind[, j] | w
  }

  if ("4_1s" %in% rules && M >= 4L) {
    B <- abs(Z) > 1
    j <- 4:M
    w <- nz[, j, drop = FALSE]
    for (k in 0:3) w <- w & B[, j - k, drop = FALSE] &
                        (S[, j, drop = FALSE] == S[, j - k, drop = FALSE])
    ind[, j] <- ind[, j] | w
  }

  if ("10x" %in% rules && M >= 10L) {
    j <- 10:M
    w <- nz[, j, drop = FALSE]
    for (k in 1:9) w <- w & (S[, j, drop = FALSE] == S[, j - k, drop = FALSE])
    ind[, j] <- ind[, j] | w
  }

  if ("R_4s" %in% rules && n_per_run >= 2L) {
    n_blocks <- M %/% n_per_run
    for (b in seq_len(n_blocks)) {
      idx <- ((b - 1L) * n_per_run + 1L):(b * n_per_run)
      seg <- Z[, idx, drop = FALSE]
      hit <- (do.call(pmax, as.data.frame(seg)) -
              do.call(pmin, as.data.frame(seg))) > 4
      ind[, idx[n_per_run]] <- ind[, idx[n_per_run]] | hit
    }
  }

  ind
}

# --- срабатывание схемы по сериям -------------------------------------------
# T[i, r] = TRUE, если в серии r цепочки i сработало хотя бы одно правило.
run_triggers <- function(ind, n_per_run) {
  nb <- nrow(ind); R <- ncol(ind) %/% n_per_run
  T <- matrix(FALSE, nb, R)
  for (r in seq_len(R)) {
    idx <- ((r - 1L) * n_per_run + 1L):(r * n_per_run)
    T[, r] <- rowSums(ind[, idx, drop = FALSE]) > 0
  }
  T
}

# --- доли по цепочкам и t-интервал ------------------------------------------
mean_ci_t <- function(x, conf = CONF) {
  n <- length(x); m <- mean(x); s <- sd(x)
  se <- s / sqrt(n)
  tc <- qt(1 - (1 - conf) / 2, df = n - 1L)
  c(mean = m, lo = m - tc * se, hi = m + tc * se)
}

wilson_ci <- function(x, n, conf = CONF) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(p = p, lo = max(centre - half, 0), hi = min(centre + half, 1))
}

# --- Pfr: доля серий со срабатыванием в стационарном режиме ------------------
# Дополнительно на тех же цепочках считается «старый» оцениватель версии 1.0
# (доля репликатов, где срабатывание произошло хотя бы раз за окно RUNS_FR
# серий, делённая на число серий) — для численной демонстрации его смещения.
false_rejection <- function(scheme, n_per_run, rules) {
  rules <- strsplit(rules, ",", fixed = TRUE)[[1]]
  counted <- RUNS_FR - BURN_IN
  total_triggers <- 0; total_rates <- numeric(0); n_chains <- 0L
  legacy_hits <- 0L
  done <- 0L
  while (done < N_REP_FR) {
    nb <- min(BLOCK_FR, N_REP_FR - done)
    Z <- matrix(rnorm(nb * RUNS_FR * n_per_run), nrow = nb)
    T <- run_triggers(rule_indicators(Z, n_per_run, rules), n_per_run)
    seg <- T[, (BURN_IN + 1L):RUNS_FR, drop = FALSE]
    total_triggers <- total_triggers + sum(seg)
    total_rates <- c(total_rates, rowSums(seg) / counted)
    legacy_hits <- legacy_hits + sum(rowSums(T) > 0)
    n_chains <- n_chains + nb
    done <- done + nb
  }
  ci <- mean_ci_t(total_rates)
  p  <- ci[["mean"]]
  tibble(
    scheme = scheme, n_per_run = n_per_run,
    n_replicates = n_chains, runs_per_chain = RUNS_FR, burn_in = BURN_IN,
    runs_counted = counted, triggers = total_triggers,
    p_per_run = p, p_lo = max(ci[["lo"]], 0), p_hi = ci[["hi"]],
    # обратная стационарная частота сигналов (геометрическое приближение ARL)
    arl_1_over_pfr = 1 / p,
    arl_1_over_pfr_lo = 1 / ci[["hi"]],
    arl_1_over_pfr_hi = 1 / max(ci[["lo"]], 1e-12),
    # оцениватель версии 1.0 — приводится только для демонстрации смещения
    legacy_p = legacy_hits / (n_chains * RUNS_FR)
  )
}

# --- ARL по Монте-Карло: среднее число серий от начала контроля до сигнала ---
# Основное определение — zero-state ARL: цепочка стартует без предшествующей истории
# (первая серия — первая серия контроля), фиксируется номер первой серии со
# срабатыванием. Это стандартное для SPC определение «среднего числа выборок от
# начала контроля до сигнала»; для правил без памяти (1-3s) оно совпадает с 1/Pfr.
# Дополнительно считается ARL при заполненном буфере: среди цепочек, не давших
# сигнала за первые B_WARM серий, берётся ОСТАТОЧНОЕ время от конца прогрева до
# первого сигнала (счёт начинается заново, поэтому величина сопоставима с zero-state
# ARL). Длина цепочки выбирается так, чтобы доля цензурированных наблюдений была
# порядка e^-12, а сами цензурированные цепочки получают подстановочное значение
# (см. комментарий к оценивателю ARL в шапке файла).
arl_first_signal <- function(scheme, n_per_run, rules, p_stationary) {
  rules <- strsplit(rules, ",", fixed = TRUE)[[1]]
  L <- as.integer(ceiling(min(max(ARL_MIN_RUNS, ARL_RUNS_MULT / p_stationary),
                            ARL_MAX_RUNS)))
  hits <- 0L; done <- 0L; first_all <- numeric(0)
  warm_hits <- 0L; warm_n <- 0L; warm_all <- numeric(0)
  while (done < N_REP_ARL) {
    nb <- min(BLOCK_ARL, N_REP_ARL - done)
    Z <- matrix(rnorm(nb * L * n_per_run), nrow = nb)
    T <- run_triggers(rule_indicators(Z, n_per_run, rules), n_per_run)

    f <- max.col(T, ties.method = "first")
    none <- rowSums(T) == 0L
    f[none] <- NA_integer_
    ok <- !is.na(f)
    hits <- hits + sum(ok)
    # цензурированные справа цепочки (сигнала нет за L серий) не отбрасываются,
    # а получают ожидаемое остаточное время L + 1/p
    first_all <- c(first_all, f[ok], rep(L + 1 / p_stationary, sum(!ok)))

    # установившийся режим: у цепочек без сигнала в прогреве (f > B_WARM или нет
    # сигнала вовсе) первое срабатывание совпадает с первым срабатыванием вообще,
    # поэтому достаточно вычесть длину прогрева
    qualified <- is.na(f) | f > B_WARM
    if (any(qualified)) {
      okw <- qualified & ok
      warm_n <- warm_n + sum(qualified)
      warm_hits <- warm_hits + sum(okw)
      warm_all <- c(warm_all, f[okw] - B_WARM,
                    rep((L - B_WARM) + 1 / p_stationary, sum(qualified & !ok)))
    }

    done <- done + nb
  }
  ci  <- mean_ci_t(first_all)
  ciw <- mean_ci_t(warm_all)
  tibble(
    scheme = scheme, runs_arl = L, n_replicates_arl = N_REP_ARL,
    arl_hits = hits, n_censored_arl = N_REP_ARL - hits,
    censored_frac = 1 - hits / N_REP_ARL,
    arl_mc = ci[["mean"]],
    arl_mc_lo = ci[["lo"]], arl_mc_hi = ci[["hi"]],
    n_chains_warm = warm_n,
    warm_censored_frac = 1 - warm_hits / max(warm_n, 1L),
    arl_warm = ciw[["mean"]],
    arl_warm_lo = ciw[["lo"]], arl_warm_hi = ciw[["hi"]]
  )
}

# --- Ped: срабатывание в пределах K_RUNS серий после появления сдвига --------
detection <- function(scheme, n_per_run, rules, shift) {
  rules <- strsplit(rules, ",", fixed = TRUE)[[1]]
  runs_total <- B_RUNS + K_RUNS
  hits <- 0L; done <- 0L; first_hits <- numeric(0)
  while (done < N_REP_DET) {
    nb <- min(BLOCK_DET, N_REP_DET - done)
    Z <- matrix(rnorm(nb * runs_total * n_per_run), nrow = nb)
    # сдвиг вводится после B_RUNS серий контрольного состояния
    Z[, (B_RUNS * n_per_run + 1L):(runs_total * n_per_run)] <-
      Z[, (B_RUNS * n_per_run + 1L):(runs_total * n_per_run)] + shift
    T <- run_triggers(rule_indicators(Z, n_per_run, rules), n_per_run)
    seg <- T[, (B_RUNS + 1L):runs_total, drop = FALSE]
    any_hit <- rowSums(seg) > 0
    hits <- hits + sum(any_hit)
    if (any(any_hit)) first_hits <- c(first_hits, apply(seg[any_hit, , drop = FALSE], 1, which.max))
    done <- done + nb
  }
  ci <- wilson_ci(hits, N_REP_DET)
  tibble(
    scheme = scheme, n_per_run = n_per_run, shift_sd = shift,
    n_replicates = N_REP_DET, k_runs = K_RUNS, hits = hits,
    p_detect_k = ci[["p"]], pd_lo = ci[["lo"]], pd_hi = ci[["hi"]],
    mean_first_hit = if (length(first_hits)) mean(first_hits) else NA_real_
  )
}

# --- точные значения для правила 1-3s ---------------------------------------
p1_exact <- function(shift) pnorm(shift - 3) + pnorm(-3 - shift)
p_exact  <- function(shift, n_per_run) 1 - (1 - p1_exact(shift))^n_per_run

cat("Точные значения для правила 1-3s:\n")
cat(sprintf("  сдвиг 0: Pfr(N=1) = %.6f (ARL %.1f), Pfr(N=2) = %.6f (ARL %.1f)\n",
            p_exact(0, 1), 1 / p_exact(0, 1), p_exact(0, 2), 1 / p_exact(0, 2)))
cat(sprintf("  Ped за %d серий, сдвиг 1 SD: N=1 %.4f, N=2 %.4f\n",
            K_RUNS, 1 - (1 - p_exact(1, 1))^K_RUNS, 1 - (1 - p_exact(1, 2))^K_RUNS))
cat(sprintf("  Ped за %d серий, сдвиг 2 SD: N=1 %.6f, N=2 %.6f\n\n",
            K_RUNS, 1 - (1 - p_exact(2, 1))^K_RUNS, 1 - (1 - p_exact(2, 2))^K_RUNS))

# --- расчёт -----------------------------------------------------------------
t0 <- Sys.time()

fr <- bind_rows(lapply(seq_len(nrow(schemes)), function(i)
  false_rejection(schemes$scheme[i], schemes$n_per_run[i], schemes$rules[i])))

det <- bind_rows(lapply(seq_len(nrow(schemes)), function(i)
  bind_rows(lapply(SHIFTS, function(sh)
    detection(schemes$scheme[i], schemes$n_per_run[i], schemes$rules[i], sh)))))

cat("Оценка истинной ARL (время до первого сигнала):\n")
arls <- bind_rows(lapply(seq_len(nrow(schemes)), function(i)
  arl_first_signal(schemes$scheme[i], schemes$n_per_run[i], schemes$rules[i],
                   fr$p_per_run[fr$scheme == schemes$scheme[i]])))

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# точные значения — только для схем, состоящих из одного правила 1-3s
fr <- fr %>%
  left_join(schemes %>% select(scheme, n_per_run, rules), by = c("scheme", "n_per_run")) %>%
  left_join(arls, by = "scheme") %>%
  mutate(
    exact_p_per_run = ifelse(rules == "1_3s", p_exact(0, n_per_run), NA_real_),
    exact_arl       = ifelse(rules == "1_3s", 1 / p_exact(0, n_per_run), NA_real_),
    # относительные расхождения: zero-state ARL, установившийся режим и 1/Pfr
    arl_over_pfr = arl_mc / arl_1_over_pfr,
    arl_warm_over_zero = arl_warm / arl_mc,
    # отклонение МК-оценки от точного аналитического значения (только 1-3s), %
    dev_from_exact_pct = ifelse(is.na(exact_arl), NA_real_, 100 * (arl_mc / exact_arl - 1))
  )
det <- det %>%
  left_join(schemes %>% select(scheme, n_per_run, rules), by = c("scheme", "n_per_run")) %>%
  mutate(
    exact_p_detect = ifelse(rules == "1_3s",
                            1 - (1 - p_exact(shift_sd, n_per_run))^K_RUNS, NA_real_)
  )

# демонстрация смещения оценивателя версии 1.0
est_cmp <- fr %>%
  transmute(scheme, p_stationary = p_per_run, exact_p = exact_p_per_run,
            legacy_p,
            legacy_over_stationary = legacy_p / p_per_run,
            legacy_over_exact = legacy_p / exact_p_per_run)

write.csv(fr,  file.path(out_dir, "qc_false_rejection.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(det, file.path(out_dir, "qc_power.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(est_cmp, file.path(out_dir, "estimator_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(false_rejection = fr, detection = det, arl = arls, estimator = est_cmp,
             params = list(seed = SEED, n_rep_fr = N_REP_FR, runs_fr = RUNS_FR,
                           burn_in = BURN_IN, n_rep_det = N_REP_DET,
                           k_runs = K_RUNS, b_runs = B_RUNS, shifts = SHIFTS,
                           n_rep_arl = N_REP_ARL, arl_runs_mult = ARL_RUNS_MULT,
                           b_warm = B_WARM)),
        file.path(out_dir, "qc_power.rds"))

# --- отчёт ------------------------------------------------------------------
cat("Схемы контроля:\n"); print(schemes, n = nrow(schemes))
cat(sprintf("\nВремя расчёта: %.1f с\n", elapsed))

cat("\nPfr (стационарная частота сигналов) и ARL по схемам:\n")
print(fr %>% transmute(
  scheme,
  Pfr = sprintf("%.5f (%.5f-%.5f)", p_per_run, p_lo, p_hi),
  ARL_mc = sprintf("%.1f (%.1f-%.1f)", arl_mc, arl_mc_lo, arl_mc_hi),
  one_over_Pfr = sprintf("%.1f", arl_1_over_pfr),
  ratio = sprintf("%.3f", arl_over_pfr),
  runs = runs_arl, censored = sprintf("%.6f", censored_frac)
) %>% as.data.frame(), row.names = FALSE)

cat("\nКонтроль корректности ARL для схем 1-3s (точное значение должно попадать в 95% ДИ оценки):\n")
print(fr %>% filter(!is.na(exact_arl)) %>%
        transmute(scheme, ARL_zero = round(arl_mc, 1), ARL_warm = round(arl_warm, 1),
                  ARL_exact = round(exact_arl, 1),
                  in_ci_zero = exact_arl >= arl_mc_lo & exact_arl <= arl_mc_hi,
                  dev_zero_pct = round(100 * (arl_mc / exact_arl - 1), 2),
                  dev_warm_pct = round(100 * (arl_warm / exact_arl - 1), 2)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nВлияние начала наблюдения на ARL (zero-state и установившийся режим):\n")
print(fr %>% transmute(
  scheme,
  ARL_zero = sprintf("%.1f", arl_mc),
  ARL_warm = sprintf("%.1f", arl_warm),
  ratio = sprintf("%.3f", arl_warm_over_zero),
  chains_warm = n_chains_warm,
  warm_censored = sprintf("%.4f", warm_censored_frac)
) %>% as.data.frame(), row.names = FALSE)

cat("\nДемонстрация смещения оценивателя версии 1.0:\n")
print(est_cmp %>%
        transmute(scheme,
                  stationary = round(p_stationary, 5),
                  legacy = round(legacy_p, 5),
                  exact = ifelse(is.na(exact_p), "-", round(exact_p, 5)),
                  ratio_to_stationary = round(legacy_over_stationary, 4)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nКонтроль корректности Pfr (1-3s): МК-оценка должна попадать в точное значение.\n")
print(fr %>% filter(!is.na(exact_p_per_run)) %>%
        transmute(scheme, MC = round(p_per_run, 5), exact = round(exact_p_per_run, 5),
                  in_ci = exact_p_per_run >= p_lo & exact_p_per_run <= p_hi) %>%
        as.data.frame(), row.names = FALSE)

cat(sprintf("\nВероятность обнаружения сдвига в пределах %d серий (95%% ДИ):\n", K_RUNS))
print(det %>% filter(shift_sd %in% c(0.5, 1, 1.5, 2)) %>%
        transmute(scheme, shift = shift_sd,
                  Ped = sprintf("%.3f (%.3f-%.3f)", p_detect_k, pd_lo, pd_hi),
                  exact = ifelse(is.na(exact_p_detect), "-", sprintf("%.3f", exact_p_detect))) %>%
        as.data.frame(), row.names = FALSE)

cat("\nСводка Ped по величине сдвига (Монте-Карло):\n")
print(det %>% select(scheme, shift_sd, p_detect_k) %>%
        pivot_wider(names_from = shift_sd, values_from = p_detect_k,
                    names_prefix = "shift_") %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
        as.data.frame(), row.names = FALSE)
