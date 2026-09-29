# 01-ingest.R
# Разбор обрезанных JSON-ответов EFLM Biological Variation Database в аккуратную таблицу.
#
# Особенности источника:
#   * web_fetch кладёт в файл преамбулу "Fetched <url>" — она уже отрезана при переносе;
#   * каждый ответ обрезан на ~100078 символов и может обрываться посреди JSON,
#     поэтому штатный fromJSON на файле целиком падает. Разбираем только целые объекты.

suppressPackageStartupMessages({
  library(jsonlite)
})

# --- пути -------------------------------------------------------------------
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
raw_dir <- file.path(root, "data", "raw")
out_dir <- file.path(root, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- чтение тела ответа -----------------------------------------------------
read_body <- function(path) {
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  i <- regexpr("[\\[{]", txt)
  if (i < 0) stop("в файле нет JSON: ", path)
  substring(txt, i)
}

# --- выделение целых объектов верхнего уровня внутри массива "data" ---------
# Сканируем от первой '[' (начало массива data), считаем глубину фигурных скобок,
# корректно пропускаем строковые литералы и экранирование.
split_top_objects <- function(txt) {
  chars <- strsplit(txt, "", fixed = TRUE)[[1]]
  start <- match(TRUE, chars == "[")
  if (is.na(start)) stop("не найден массив data")
  objs <- character(0)
  depth <- 0L
  in_str <- FALSE
  esc <- FALSE
  buf_start <- NA_integer_
  for (k in seq.int(start, length(chars))) {
    ch <- chars[k]
    if (in_str) {
      if (esc) esc <- FALSE
      else if (ch == "\\") esc <- TRUE
      else if (ch == "\"") in_str <- FALSE
      next
    }
    if (ch == "\"") { in_str <- TRUE; next }
    if (ch == "{") {
      if (depth == 0L) buf_start <- k
      depth <- depth + 1L
      next
    }
    if (ch == "}") {
      depth <- depth - 1L
      if (depth == 0L && !is.na(buf_start)) {
        objs <- c(objs, paste0(chars[buf_start:k], collapse = ""))
        buf_start <- NA_integer_
      }
      next
    }
  }
  objs
}

# --- разбор одного объекта {analyte, metas} ---------------------------------
pick_name <- function(a) {
  cand <- c(a$display_name, a$full_name, a$description, a$alternate_names)
  cand <- cand[!is.na(cand) & nzchar(trimws(cand))]
  if (length(cand)) trimws(cand[1]) else paste0("analyte_", a$id)
}
num_or_na <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)

parse_one <- function(obj) {
  x <- fromJSON(obj, simplifyVector = FALSE)
  a <- x$analyte
  if (is.null(a$id)) return(NULL)
  rows <- list()
  for (m in x$metas) {
    rows[[length(rows) + 1L]] <- data.frame(
      analyte_id   = as.integer(a$id),
      analyte_name = pick_name(a),
      test_code    = if (is.null(a$test_code)) NA_character_ else a$test_code,
      nlmc_id      = if (is.null(a$nlmc_id)) NA_character_ else a$nlmc_id,
      disciplines  = if (is.null(a$disciplines)) NA_character_ else a$disciplines,
      var_type     = sub("^:", "", m$var_type),
      median       = num_or_na(m$median),
      lower        = num_or_na(m$lower),
      upper        = num_or_na(m$upper),
      number_used  = as.integer(num_or_na(m$number_used)),
      matrix       = if (is.null(m$matrix$matrix_expansion)) NA_character_
                     else m$matrix$matrix_expansion,
      meta_id      = as.integer(num_or_na(m$id)),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# --- конвейер ---------------------------------------------------------------
meta_file <- file.path(raw_dir, "meta_calculations.json")
body <- read_body(meta_file)
objs <- split_top_objects(body)

parsed <- lapply(objs, parse_one)
parsed <- Filter(Negate(is.null), parsed)
bv <- do.call(rbind, parsed)
bv$source <- "bulk"

# --- добор: аналиты, не попавшие в обрезанный дамп --------------------------
# Значения получены через web_fetch из /api/meta_calculations/meta_by_analyte/{id}
# и зафиксированы в data/raw/topup_meta.csv (см. README.md в data/raw).
topup_file <- file.path(raw_dir, "topup_meta.csv")
if (file.exists(topup_file)) {
  tp <- read.csv(topup_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  tp$disciplines <- NA_character_
  tp$source <- "topup"
  bv <- rbind(bv[, names(tp)], tp)
}

# Приоритет у полного дампа: если аналит есть и там, и там — оставляем bulk.
bv <- bv[order(bv$analyte_id, bv$var_type, bv$source != "bulk"), ]
bv <- bv[!duplicated(bv[, c("analyte_id", "var_type")]), ]
bv <- bv[order(bv$analyte_name, bv$var_type), ]
rownames(bv) <- NULL

write_csv_sc(bv, file.path(out_dir, "bv_meta.csv"))

# --- отчёт ------------------------------------------------------------------
wide <- reshape(bv[, c("analyte_id", "analyte_name", "var_type", "median", "number_used")],
                idvar = c("analyte_id", "analyte_name"), timevar = "var_type",
                direction = "wide")

cat("Файл-источник      :", basename(meta_file), "\n")
cat("Размер тела, симв. :", nchar(body), "\n")
cat("Целых объектов     :", length(objs), "\n")
cat("Разобрано записей  :", nrow(bv), "\n")
cat("Уникальных аналитов:", length(unique(bv$analyte_id)), "\n")
cat("  с CVi            :", length(unique(bv$analyte_id[bv$var_type == "cvi"])), "\n")
cat("  с CVg            :", length(unique(bv$analyte_id[bv$var_type == "cvg"])), "\n")
cat("Только CVi, без CVg:", length(setdiff(unique(bv$analyte_id[bv$var_type == "cvi"]),
                                          unique(bv$analyte_id[bv$var_type == "cvg"]))), "\n")

cat("\nПроверка обрезки хвоста:\n")
cat("  последний объект завершён корректно:", !is.na(objs[length(objs)]), "\n")
tail_txt <- substr(body, nchar(body) - 60, nchar(body))
cat("  последние 60 символов тела:", gsub("\n", " ", tail_txt), "\n")

cat("\nМатрицы:\n"); print(table(bv$matrix, useNA = "ifany"))
cat("\nЗаписей по источнику:\n"); print(table(bv$source))
cat("Аналитов по источнику:\n")
print(table(tapply(bv$source, bv$analyte_id, function(s) s[1])))

cat("\nКонтрольная выборка:\n")
key <- c("Glucose", "Sodium", "Potassium", "Creatinine", "Albumin",
         "Triglycerides", "Thyroid stimulating hormone (TSH)", "Calcium (Ca)")
print(wide[wide$analyte_name %in% key, ], row.names = FALSE)
