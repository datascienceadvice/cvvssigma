# cvvssigma

**Consistency of biological-variation goals and the sigma metric in internal quality
control planning: an analytical result and its illustration on a panel of 101 analytes.**

This repository contains the R code, derived tables and figures for a study that compares
two frameworks used together when planning internal quality control in clinical
laboratories — analytical goals derived from biological variation (the Fraser and Petersen
hierarchy) and the Six Sigma metric. Both frameworks are expressed through a single
allowable total error (TEa). The study answers the question analytically under a shared TEa
and illustrates the answer quantitatively on a panel of 101 analytes.

Repository: <https://github.com/datascienceadvice/cvvssigma>

---

## What the study shows

With one and the same TEa used by both frameworks:

- **Bias at the allowable limit.** The ratio of the imprecision required for a target sigma
  to the biological-variation-derived allowable imprecision equals `k / σ`, where `k` is the
  imprecision multiplier. This is an algebraic identity: it depends on neither the analyte
  nor the stringency level, and it would be the same for 5, 101 or 1000 analytes. It is
  below 1 for **all 101 analytes** at σ = 4, 5 and 6 for any `k < 4`
  (0.4125 / 0.3300 / 0.2750 at `k` = 1.65).
- **Zero bias (idealised scenario).** The ratio becomes `(k + 0.5·R) / σ`, where
  `R = sqrt(1 + (CV_G/CV_I)²)`. At the desirable level the sigma requirement is stricter
  than the biological-variation goal for **92 of 101 analytes (91.1 %)** at σ = 4 and for
  **98 of 101 (97.0 %)** at σ = 5 and 6; the median ratios are 0.730, 0.584 and 0.486.
  An analyte escapes only if `CV_G/CV_I` exceeds `sqrt(4(σ − k)² − 1)` — 4.59 at σ = 4,
  6.62 at σ = 5, 8.64 at σ = 6 — which selects 9 analytes at σ = 4 and
  3 (CA 19-9, alpha-fetoprotein, CEA) at σ = 5 and 6.
- **The two frameworks are hierarchical, not contradictory.** Because the ratio is below 1,
  a method that reaches the target sigma automatically meets the biological-variation goal;
  the converse fails. A method operating exactly at the desirable goal attains
  `σ = k = 1.65` with bias at the allowable limit and `σ = k + 0.5·R` with zero bias — a
  panel median of **2.9**, below the threshold of 4 for 92 of 101 analytes. The desirable
  goal is therefore not the binding constraint.
- **The stringency level does not change this.** It is the ratio to CV_A that is
  level-invariant (0.486 at σ = 6 for optimal, desirable and minimum alike), while the
  ratio to CV_I (0.122 / 0.243 / 0.365) and the absolute requirement (0.91 / 1.82 / 2.74 %)
  scale with the level multiplier.
- **Robustness to the imprecision multiplier.** At the allowable bias limit the discrepancy
  holds for σ = 4–6 and any `k < 4`; under zero bias the number of affected analytes at
  σ = 4 falls from 98 of 101 (`k` = 1.0) to 64 of 101 (`k` = 2.58).
- **The threshold lists are not fully robust to input uncertainty.** Substituting the
  lower/upper bounds of the `CV_I` and `CV_G` intervals reported by the database moves the
  number of affected analytes to 59–101 at σ = 4, 77–101 at σ = 5 and 87–101 at σ = 6
  (point estimates 92 / 98 / 98). The classification changes for 42, 24 and 14 analytes
  respectively, and no analyte stays unaffected under every combination of bounds
  (`out/table_s7_threshold_sensitivity.csv`). CEA at σ = 6 sits 0.31 % above the threshold.
- **Quality control schemes.** False-rejection probability per run rises from 0.00270
  (1-3s, N = 1) to 0.04088 (full rule set plus 10x, N = 4); the Monte Carlo average run
  length falls from 369.6 to 25.3 runs; the probability of a signal at a 1 SD shift within
  100 runs rises from 0.899 to 1.000. At a 2 SD shift the lower 95 % confidence bound is at
  least 0.9998 for every scheme. `1/Pfr` differs from the zero-state ARL by no more than
  3.8 % (the largest gap is for rules with memory, where the comparable quantity is the
  filled-buffer ARL), and the ARL starting convention changes the estimate by no more than
  1.7 %.

Formulas, assumptions and the full methodological rationale are in
[`METODIKA.md`](METODIKA.md) (in Russian, matching the manuscript).

---

## Repository contents

```
.
├── README.md                    this file
├── METODIKA.md                  calculation methodology: formulas, assumptions, definitions
├── R/
│   ├── run-all.R                run the whole pipeline with one command
│   ├── 01-ingest.R              parse the raw API responses into out/bv_meta.csv
│   ├── 02-sigma.R               analytical goals, sigma metric, closed-form relations,
│   │                            sensitivity to the imprecision multiplier k, threshold
│   │                            robustness to input intervals, derived per-analyte table
│   ├── 03-power.R               Monte Carlo: false rejection, ARL, shift detection
│   ├── 04-report.R              summary tables (CSV) and figures (PNG)
│   └── 05-tables-article.R      manuscript-ready tables and the key-numbers digest
├── data/
│   ├── panel_analyte_ids.csv    panel composition: analyte ids, matrix, provenance
│   └── raw/                     provenance of the source data (no database content)
└── out/                         derived results: tables, figures, key numbers
```

### How the panel was assembled

The panel contains 101 analytes, assembled from the EFLM Biological Variation Database on
25 September 2026:

- **77 analytes** come from the bulk `/api/meta_calculations/` response. That response could
  only be obtained truncated (100 078 characters), and the truncation is caused by the
  web-access tool used to fetch it, not by the server: three requests to the same endpoint
  with different pagination parameters returned byte-identical bodies, so the server ignores
  the parameters and the full body never arrives. The bulk part is therefore an
  alphabetical prefix of the database, not a sample.
- **24 analytes** were fetched individually through `/api/meta_calculations/meta_by_analyte/{id}`.
  Inclusion rule: the basic routine menu (electrolytes, metabolites, lipids, proteins,
  enzymes, hormones, basic haematology) that the truncated bulk response did not contain.

`data/panel_analyte_ids.csv` lists every analyte with its id, matrix, provenance
(`bulk` / `topup`) and the number of primary studies behind its `CV_I` and `CV_G`. No
`CV_I` / `CV_G` values are included.

---

## What is not included, and why

| Excluded | Reason |
|---|---|
| `manuscript.md` | The article is under review and is not published here. |
| `data/raw/*.json`, `data/raw/topup_meta.csv` | Content of the EFLM Biological Variation Database. The licence permits use in non-commercial research but forbids redistribution and storage in any other electronic retrieval system. |
| `out/bv_meta.csv`, `out/sigma_specs.csv`, `out/sigma_requirements.csv`, `out/table1_analytes_desirable.csv`, `out/worked_example.csv`, `out/article_tables.md` | These pipeline artefacts carry the `CV_I` / `CV_G` values in explicit column form and therefore fall under the same restriction. The published alternative is `out/table1_analytes_derived.csv`, which holds derived quantities only. |
| `issues.txt` | Reviewer correspondence. |
| `docs/NOTES.md` | Internal working notes: environment constraints, licensing checklist, correction log, remaining tasks before submission. |
| `review/` | Editorial synthesis and revision roadmap produced by the internal review panel. |

Everything else in `out/` contains **derived** quantities only (ratios, medians, counts,
thresholds, control-rule operating characteristics) and is published. The exact list of
exclusions is in [`.gitignore`](.gitignore).

Because the source data are excluded, a fresh clone cannot run the pipeline until the data
are obtained. See the next section.

---

## Data source and licence

The study uses the open **EFLM Biological Variation Database**
(<https://biologicalvariation.eu/>), retrieved on 25 September 2026. At that date the
database reported 191 meta-analyses, 3366 primary specifications and 608 references
(endpoint `/api/meta`). The database is
copyright of EFLM; the terms of use state, verbatim:

> This website and its content is copyright of EFLM. All rights reserved. You may not,
> except with our express written permission, distribute or commercially exploit the
> content. You may, however, print or download to a local hard disk extracts for your
> personal and non-commercial use only. If you present the content you must acknowledge
> the website as the source of the material. You cannot transmit the content or store it
> in any other website or other form of electronic retrieval system.

**To reproduce the results**, export the meta-analysis values yourself and place them in
`data/raw/` using the file names expected by `R/01-ingest.R`; the expected inventory and
the provenance of the original extraction are documented in
[`data/raw/README.md`](data/raw/README.md). No database content is distributed with this
repository.

---

## Requirements

- R 4.3.1 (the pipeline was developed and verified on this version)
- R packages: `jsonlite`, `dplyr`, `tidyr`, `ggplot2`

On Windows, `Rscript.exe` is often not on `PATH`; in that case use its full path, for
example:

```powershell
& 'C:\Program Files\R\R-4.3.1\bin\Rscript.exe' --vanilla R\run-all.R
```

---

## Running the pipeline

Whole pipeline with one command:

```powershell
Rscript --vanilla R/run-all.R
```

Individual steps, in order:

```powershell
Rscript --vanilla R/01-ingest.R          # parse raw responses      (~1 s)
Rscript --vanilla R/02-sigma.R           # goals, sigma, k study    (~10 s)
Rscript --vanilla R/03-power.R           # Monte Carlo              (~3.5 min)
Rscript --vanilla R/04-report.R          # tables and figures       (~5 s)
Rscript --vanilla R/05-tables-article.R  # manuscript-ready tables  (~1 s)
```

Total runtime is about 4 min on an idle machine. The Monte Carlo step is the only expensive
one and is CPU-bound; wall-clock time depends strongly on machine load.

Step 2 in the list is the analytical core and does not need the Monte Carlo results; steps
4 and 5 do.

### Reproducibility

The random number generator is seeded (`set.seed(20260925)`, parameter block at the top of
`R/03-power.R`), so repeated runs reproduce `out/qc_false_rejection.csv`,
`out/qc_power.csv`, `out/estimator_comparison.csv` and the figures exactly. Only the
wall-clock time varies.

---

## Outputs

Published with the repository:

| File | Contents |
|---|---|
| `out/supplementary_tables.md` | Supplementary tables S1–S7 |
| `out/article_numbers.txt` | Key numbers quoted in the manuscript |
| `out/session_info.txt` | R version and package versions used for the results |
| `out/table1_analytes_derived.csv` | Per-analyte derived quantities for all 101 analytes (ratios, thresholds, study counts; no `CV_I` / `CV_G` values) |
| `data/panel_analyte_ids.csv` | Panel composition and provenance |
| `out/table2_requirement_vs_bv.csv` | Ratio of sigma requirements to biological-variation goals by scenario and level |
| `out/table_s1_qc_operating_characteristics.csv` | Pfr, Monte Carlo ARL, filled-buffer ARL and signal probability for six control schemes |
| `out/table_s2_threshold_analytes.csv` | Lists of analytes above the `CV_G/CV_I` threshold |
| `out/table_s5_level_strictness.csv` | Dependence of the requirements on the stringency level |
| `out/table_s6_arl_start.csv` | ARL under the two starting conventions |
| `out/table_s7_threshold_sensitivity.csv` | Robustness of the threshold classifications to the input intervals |
| `out/k_sensitivity.csv` | Sensitivity of the results to the imprecision multiplier `k` |
| `out/estimator_comparison.csv` | Naive versus stationary false-rejection estimator |
| `out/qc_false_rejection.csv`, `out/qc_power.csv`, `out/qc_power.rds` | Full Monte Carlo results |
| `out/fig1_required_imprecision.png`, `out/fig2_ratio_to_bv.png`, `out/figS1_qc_power.png` | Figures (300 dpi) |

Generated locally only, and therefore absent from the clone (see
[`.gitignore`](.gitignore)): `out/bv_meta.csv`, `out/sigma_specs.csv`,
`out/sigma_requirements.csv`, `out/table1_analytes_desirable.csv`,
`out/worked_example.csv`, `out/article_tables.md`.

---

## Licence

No licence file is provided with this repository, so no licence is granted for the code or
the documentation and default copyright applies. If you intend to reuse the code, please
get in touch.

The database content is **not** covered in any case: it remains the property of EFLM and is
governed by the terms quoted above.

---

## Citation

If you use this code or the derived results, please cite the manuscript (details to be
added on publication) and acknowledge the EFLM Biological Variation Database as the source
of the underlying data.

## Contact

Questions and corrections: open an issue in this repository.
