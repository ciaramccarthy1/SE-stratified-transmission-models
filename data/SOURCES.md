# Data sources

Provenance for every file in `data/`. Files marked **auto-fetched** are
downloaded by `codes/fetch_data.R` and cached locally. Files marked
**generated** are produced by `codes/prepare_model_inputs.R` from those
raw inputs.

---

## Contact matrices

### `base_matrix.csv` — auto-fetched
- **Description**: Long-form age × IMD contact matrix. 5 IMD quintiles × 16 five-year age bands (0–4, 5–9, …, 70–74, 75+), participant × contact. Columns: `Participant_IMD, Contact_IMD, Participant_age_group, Contact_age_group, mean, lower, upper`.
- **Source**: [lucy-gf/imd_matrices](https://github.com/lucy-gf/imd_matrices), built from the LSHTM Reconnect Social Contact Survey ([Zenodo](https://zenodo.org/records/17339866)).
- **URL**: <https://raw.githubusercontent.com/lucy-gf/imd_matrices/main/matrices/base_matrix.csv>
- **Fetched by**: `codes/fetch_data.R`
- **Notes**: Already balanced for reciprocity in the source repo.

### `Mas50.csv` — generated
- **Description**: 50×50 wide contact matrix for the model's 10 age bands × 5 IMD. `cm[i, j]` = mean contacts of contact-group `j` by participant-group `i`, where `ig = is*na + ia` (IMD-major).
- **Inputs**: `base_matrix.csv`
- **Generator**: `codes/prepare_model_inputs.R`
- **Aggregation method**:
  - Contact-side merge of {a₁, a₂} → A: **sum** of rates (a participant in I contacts both sub-bands).
  - Participant-side merge of {p₁, p₂} → P: **uniform mean** of rates (one participant in exactly one sub-band).
  - **TODO**: participant-side should be population-weighted using ONS sub-band populations (already fetched as `ons_lsoa_syoa_2022-2024.xlsx`; just wire into `prepare_model_inputs.R`).

---

## Demographics

### `ons_lsoa_syoa_2022-2024.xlsx` — auto-fetched
- **Description**: ONS LSOA-level population by single year of age (0–90+) and sex, mid-2022 through mid-2024. One sheet per year; `prepare_model_inputs.R` uses the Mid-2024 sheet.
- **Source**: [ONS — Lower super output area mid-year population estimates](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates).
- **Fetched by**: `codes/fetch_data.R`

### `iod2025_lsoa_ranks_deciles.csv` — auto-fetched
- **Description**: English Indices of Deprivation 2025, File 7. Per LSOA: all 7 IMD domain ranks/scores/deciles plus population denominators. Decile 1 = most deprived.
- **Source**: [GOV.UK — English indices of deprivation 2025](https://www.gov.uk/government/statistics/english-indices-of-deprivation-2025), File 7.
- **Fetched by**: `codes/fetch_data.R`

### `demographics_10age.csv` — generated
- **Description**: Population by IMD quintile × 10 model age bands (0–4, 5–14, 15–19, 20–29, 30–39, 40–49, 50–59, 60–64, 65–74, 75+). One row per (IMD, age) cell; columns `Age, IMD, Population, tot_pop, Proportion`.
- **Inputs**: `ons_lsoa_syoa_2022-2024.xlsx` (Mid-2024 sheet), `iod2025_lsoa_ranks_deciles.csv`.
- **Generator**: `codes/prepare_model_inputs.R`
- **Method**: Sum female + male at each single year of age per LSOA; bin into model age bands; join LSOAs to IMD deciles; convert decile → quintile (`ceiling(decile/2)`); aggregate population by (quintile, band). Inner join keeps the ~33.7k English LSOAs present in both sources.

---

## RSV vaccination uptake

### `ukhsa_rsv_uptake_jan2026.html` — auto-fetched
- **Description**: Raw HTML of the UKHSA "Respiratory Syncytial Virus (RSV) older adults vaccination coverage in England" January 2026 report.
- **Source**: [GOV.UK — RSV older adults vaccination coverage (Jan 2026 report)](https://www.gov.uk/government/statistics/respiratory-syncytial-virus-rsv-older-adults-vaccination-coverage-in-england/respiratory-syncytial-virus-rsv-older-adults-vaccination-coverage-in-england-january-2026-report).
- **Fetched by**: `codes/fetch_data.R`
- **Notes**: Report is updated periodically; URL points to the January 2026 snapshot. Includes routine + catch-up cohorts.

### `rsv_uptake_by_imd_decile.csv` — generated
- **Description**: Decile-level RSV uptake percentages parsed from the UKHSA report. Columns: `decile, uptake_pct`. Decile 1 = most deprived.
- **Inputs**: `ukhsa_rsv_uptake_jan2026.html`
- **Generator**: `codes/fetch_data.R` (HTML table parsed inline via regex; no rvest dependency).

---

## Per-age parameters — NOT in `data/`

Per-age values for `u` (susceptibility), `y` (clinical fraction), `m` (mortality given clinical), `h` (hospitalisation fraction given clinical), `mH` (mortality given hospitalised), `rrep` (reporting rate) live in `codes/pars*_.r` files, with source citations inline. They are **not** automatically pullable — the source papers' tables need manual extraction.

| Disease | Pars files | Cited sources |
|---|---|---|
| COVID-19 | `parsC_.r`, `parsCv_.r` | Davies 2020 (Nat Med), Verity 2020 (IFR), Knock 2021 |
| Influenza | `parsF_.r`, `parsFv_.r` | Baguelin 2013, LG Global IFR |
| RSV | `parsR_.r`, `parsRv_.r` | Hodgson 2020 (Lancet ID), Henderson 1979, Waterlow 2021 |

Each per-age vector currently inherits values from the 9-band era; band shifts at 11/12, 17/18, 29/30, 69/70 mean these need re-aggregating from the source papers. See `TODO(10-age scaffold)` and `TODO(H scaffold)` markers in the pars files.

---

## How to refresh inputs

```sh
# Download fresh raw inputs:
Rscript codes/fetch_data.R

# Regenerate derived files from raw inputs:
Rscript codes/prepare_model_inputs.R
```

To force a re-download of any fetched file, delete it from `data/` first (the cache check in `fetch_data.R` is just file existence).
