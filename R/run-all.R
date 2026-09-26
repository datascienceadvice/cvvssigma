# run-all.R — полный конвейер темы A1
#
# Запуск:
#   & 'D:\R\R-4.3.1\bin\Rscript.exe' --vanilla R\run-all.R
#
# Скрипты читают друг за другом: 01 -> 02 -> 03 -> 04 -> 05.
# Данные в data/raw уже выгружены (см. data/raw/README.md): сеть для процессов
# недоступна, автоматическая загрузка невозможна.

args_all <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
root <- normalizePath(file.path(here, ".."), mustWork = TRUE)

scripts <- c("01-ingest.R", "02-sigma.R", "03-power.R", "04-report.R", "05-tables-article.R")

t0 <- Sys.time()
for (s in scripts) {
  cat("\n", strrep("=", 62), "\n== ", s, "\n", strrep("=", 62), "\n", sep = "")
  ts <- Sys.time()
  source(file.path(here, s), echo = FALSE, chdir = FALSE)
  cat(sprintf("-- %s выполнено за %.1f с\n", s, as.numeric(difftime(Sys.time(), ts, units = "secs"))))
}
cat(sprintf("\nКонвейер завершён за %.1f с. Результаты: %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            file.path(root, "out")))
