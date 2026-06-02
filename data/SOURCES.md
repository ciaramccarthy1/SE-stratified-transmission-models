# Data sources

Provenance for every file in `data/`. Files marked **auto-fetched** are
downloaded by `codes/fetch_data.R` and cached locally. Files marked
**generated** are produced by `codes/prepare_model_inputs.R` from
either raw or other inherited inputs.

---

## Contact matrices

### `base_matrix.csv` — auto-fetched
- **Description**: Long-form age × IMD contact matrix. 5 IMD quintiles × 16 five-year age bands (0–4, 5–9, …, 70–74, 75+), participant × contact. Columns: `Participant_IMD, Contact_IMD, Participant_age_group, Contact_age_group, mean, lower, upper`.
- **Source**: [lucy-gf/imd_matrices](https://github.com/lucy-gf/imd_matrices), built from the LSHTM Reconnect Social Contact Survey ([Zenodo](https://zenodo.org/records/17339866)).
- **URL**: <https://raw.githubusercontent.com/lucy-gf/imd_matrices/main/matrices/base_matrix.csv>
- **Fetched by**: `codes/fetch_data.R`
- **Notes**: Already balanced for reciprocity in the source repo.

### `Mas50_urban.csv` — generated
- **Description**: 50×50 wide contact matrix for the model's 10 age bands × 5 IMD, urban. `cm[i, j]` = mean contacts of contact-group `j` by participant-group `i`, where `ig = is*na + ia` (IMD-major).
- **Inputs**: `base_matrix.csv`
- **Generator**: `codes/prepare_model_inputs.R`
- **Aggregation method**:
  - Contact-side merge of {a₁, a₂} → A: **sum** of rates (a participant in I contacts both sub-bands).
  - Participant-side merge of {p₁, p₂} → P: **uniform mean** of rates (one participant in exactly one sub-band).
  - **TODO**: participant-side should be population-weighted using ONS sub-band populations.

### `Mas45_urban.csv` — legacy, inherited
- **Description**: 45×45 wide contact matrix for the old 9-band structure. Age bands: 0–4, 5–11, 12–17, 18–25, 26–34, 35–49, 50–69, 70–79, 80+.
- **Source**: Inherited from upstream `JAN-Filipe/SE-stratified-transmission-models`. **Original provenance unknown** — predates this fork. Likely Polymod-derived but augmented for IMD; the augmentation methodology is undocumented.
- **Status**: No longer used by the 10-band model — kept for reference / diff comparison.

---

## Demographics

### `demographics2021.csv` — legacy, inherited
- **Description**: Population by age × IMD × urban/rural for England, 2021. 9 age bands: 0–4, 5–11, 12–17, 18–29, 30–39, 40–49, 50–59, 60–69, 70+. Columns: `Age, IMD, rural, Population, tot_pop, Proportion`.
- **Source**: Inherited from upstream `JAN-Filipe/SE-stratified-transmission-models`. **Original provenance unknown** — likely derived from ONS 2021 Census small-area data combined with LSOA-IMD lookups and LSOA-urban/rural classification, but the build pipeline is undocumented.
- **Status**: Used only as input to the scaffold script for the 10-band rebanding.
- **TODO**: Rebuild from primary ONS sources for traceability; document the LSOA → IMD → urban/rural lookups used.

### `demographics2021_10age.csv` — generated
- **Description**: 10-band rebanding of `demographics2021.csv` aligned to the model's age bands (0–4, 5–14, 15–19, 20–29, 30–39, 40–49, 50–59, 60–64, 65–74, 75+).
- **Inputs**: `demographics2021.csv`
- **Generator**: `codes/prepare_model_inputs.R`
- **Method**: Linear-uniform splits within bands that don't align between the 9-band source and 10-band target. The 70+ source band is assumed to span 70–89 effectively (width 20) for splitting into 65–74 contributions and 75+.
- **TODO**: Rebuild from ONS 5-year-band data for exact splits; remove the uniform-within-band assumption.

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
# Download fresh raw inputs (Reconnect matrix only, for now):
Rscript codes/fetch_data.R

# Regenerate derived files from raw inputs:
Rscript codes/prepare_model_inputs.R
```

To force a re-download of any fetched file, delete it from `data/` first
(the cache check in `fetch_data.R` is just file existence).
