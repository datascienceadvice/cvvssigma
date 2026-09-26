# cvvssigma

**Consistency of biological-variation goals and the sigma metric in internal quality
control planning: analysis of 101 analytes.**

This repository contains the R code, derived tables and figures for a study that compares
two frameworks used together when planning internal quality control in clinical
laboratories — analytical goals derived from biological variation (the Fraser and Petersen
hierarchy) and the Six Sigma metric. Both frameworks are expressed through a single
allowable total error (TEa); the study asks whether the imprecision requirements they imply
can be satisfied simultaneously.

Repository: <https://github.com/datascienceadvice/cvvssigma>

---

## What the study shows

With one and the same TEa used by both frameworks:

- **Bias at the allowable limit.** The ratio of the imprecision required for a target sigma
  to the biological-variation-derived allowable imprecision equals `k / σ`, where `k` is the
  coverage factor. It depends on neither the analyte nor the stringency level, and it is
  below 1 for **all 101 analytes** at σ = 4, 5 and 6 for any `k < 4`
  (0.4125 / 0.3300 / 0.2750 at `k` = 1.65).
- **Zero bias (idealised scenario).** The ratio becomes `(k + 0.5·R) / σ`, where
  `R = sqrt(1 + (CV_G/CV_I)²)`. At the desirable level the sigma requirement is stricter
  than the biological-variation goal for **92 of 101 analytes (91.1 %)** at σ = 4 and for
  **98 of 101 (97.0 %)** at σ = 5 and 6; the median ratios are 0.730, 0.584 and 0.486.
  An analyte escapes only if `CV_G/CV_I` exceeds `sqrt(4(σ − k)² − 1)` — 4.59 at σ = 4,
  6.62 at σ = 5, 8.64 at σ = 6 — which selects exactly 9 analytes at σ = 4 and exactly
  3 (CA 19-9, alpha-fetoprotein, CEA) at σ = 5 and 6.
- **The stringency level does not remove the discrepancy.** It is the ratio to CV_A that is
  level-invariant (0.486 at σ = 6 for optimal, desirable and minimum alike), while the
  ratio to CV_I (0.122 / 0.243 / 0.365) and the absolute requirement (0.91 / 1.82 / 2.74 %)
  scale with the level multiplier.
- **Robustness to the coverage factor.** At the allowable bias limit the discrepancy holds
  for σ = 4–6 and any `k < 4`; under zero bias the number of affected analytes at σ = 4
  falls from 98 of 101 (`k` = 1.0) to 64 of 101 (`k` = 2.58).
- **Quality control schemes.** False-rejection probability per run rises from 0.00270
  (1-3s, N = 1) to 0.04088 (full rule set plus 10x, N = 4); the zero-state average run
  length falls from 365.4 to 25.6 runs; the probability of detecting a 1 SD shift within
  100 runs rises from 0.899 to 1.000. At a 2 SD shift the lower 95 % confidence bound of
  the detection probability is at least 0.9998 for every scheme. `1/Pfr` differs from the
  true ARL by no more than 4.7 %, and the choice of the ARL starting convention (empty
  buffer versus filled buffer) changes it by no more than 1 %.

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
│   │                            sensitivity to the coverage factor k
│   ├── 03-power.R               Monte Carlo: false rejection, ARL, shift detection
│   ├── 04-report.R              summary tables (CSV) and figures (PNG)
│   └── 05-tables-article.R      manuscript-ready tables and the key-numbers digest
├── data/raw/                    provenance of the source data (no database content)
├── docs/NOTES.md                working notes: environment constraints, licensing,
│                                correction log, remaining tasks before submission
└── out/                         derived results: tables, figures, key numbers
```

---

## What is not included, and why

| Excluded | Reason |
|---|---|
| `manuscript.md` | The article is under review and is not published here. |
| `data/raw/*.json`, `data/raw/topup_meta.csv` | Content of the EFLM Biological Variation Database. The licence permits use in non-commercial research but forbids redistribution and storage in any other electronic retrieval system. |
| `out/bv_meta.csv`, `out/sigma_specs.csv`, `out/sigma_requirements.csv`, `out/table1_analytes_desirable.csv`, `out/worked_example.csv`, `out/article_tables.md` | These pipeline artefacts carry the `CV_I` / `CV_G` values in explicit column form and therefore fall under the same restriction. |
| `issues.txt` | Reviewer correspondence. |

Everything else in `out/` contains **derived** quantities only (ratios, medians, counts,
thresholds, control-rule operating characteristics) and is published. The exact list of
exclusions is in [`.gitignore`](.gitignore).

Because the source data are excluded, a fresh clone cannot run the pipeline until the data
are obtained. See the next section.

---

## Data source and licence

The study uses the open **EFLM Biological Variation Database**
(<https://biologicalvariation.eu/>), retrieved on 25 September 2026. The database is
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
Rscript --vanilla R/02-sigma.R           # goals, sigma, k study    (~2 s)
Rscript --vanilla R/03-power.R           # Monte Carlo              (1.5-2.5 min)
Rscript --vanilla R/04-report.R          # tables and figures       (~5 s)
Rscript --vanilla R/05-tables-article.R  # manuscript-ready tables  (~1 s)
```

Total runtime is 1.5–3 min on an idle machine. The Monte Carlo step is the only expensive
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
| `out/supplementary_tables.md` | Supplementary tables S1–S6 |
| `out/article_numbers.txt` | Key numbers quoted in the manuscript |
| `out/table2_requirement_vs_bv.csv` | Ratio of sigma requirements to biological-variation goals by scenario and level |
| `out/table_s1_qc_operating_characteristics.csv` | Pfr, ARL and detection probability for six control schemes |
| `out/table_s2_threshold_analytes.csv` | Exhaustive lists of analytes above the `CV_G/CV_I` threshold |
| `out/table_s5_level_strictness.csv` | Dependence of the requirements on the stringency level |
| `out/table_s6_arl_start.csv` | ARL under the two starting conventions |
| `out/k_sensitivity.csv` | Sensitivity of the results to the coverage factor `k` |
| `out/estimator_comparison.csv` | Legacy versus current false-rejection estimator |
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
