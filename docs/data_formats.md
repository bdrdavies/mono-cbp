# Data Formats

This document describes the input and output data formats used by mono-cbp.

## Input Data Formats

### Catalogue CSV

The catalogue contains eclipse and orbital parameters for each eclipsing binary system.

#### Standard Format

Required columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `tess_id` | int | TESS Input Catalog ID | 260128333 |
| `period` | float | Orbital period (days) | 5.234 |
| `bjd0` | float | Reference epoch (BJD - 2457000) | 1325.456 |
| `sectors` | str | TESS sectors (comma-separated) | "1,2,3,4" |
| `prim_pos` | float | Primary eclipse phase position [0-1] | 0.0 |
| `prim_width` | float | Primary eclipse phase width [0-1] | 0.05 |
| `sec_pos` | float | Secondary eclipse phase position [0-1] | 0.5 |
| `sec_width` | float | Secondary eclipse phase width [0-1] | 0.03 |

#### TEBC Format

When using `--tebc` flag or `TEBC=True`, `mono-cbp` assumes that the catalogue contains double Gaussian (2g) and polyfit (pf) variants:

| Column | Type | Description |
|--------|------|-------------|
| `prim_pos_2g`, `prim_pos_pf` | float | Primary eclipse position (2g preferred) |
| `prim_width_2g`, `prim_width_pf` | float | Primary eclipse width (2g preferred) |
| `sec_pos_2g`, `sec_pos_pf` | float | Secondary eclipse position (2g preferred) |
| `sec_width_2g`, `sec_width_pf` | float | Secondary eclipse width (2g preferred) |

The code preferentially uses `*_2g` values, falling back to `*_pf` only if 2g values are invalid (null or zero). This is because the `*_2g` widths are generally wider than the `*_pf` ones, so the eclipse removal is more conservative.

### Light Curve Files

Light curves should be in one of two formats with the following naming convention:

**Filename convention:** `TIC_{tic_id}_{sector}.{txt|npz}`
- Sector is always 2 digits with leading zeros (e.g., `06`, `10`)
- Examples: `TIC_260128333_06.txt`, `TIC_146530594_09.npz`

#### Format 1: Text files (`.txt`)

Space or tab-delimited ASCII file with header row.

**Header:** `TIME FLUX FLUX_ERR PHASE ECL_MASK`

**Columns:**
1. `TIME`: BJD - 2457000
2. `FLUX`: Normalized flux (median = 1.0)
3. `FLUX_ERR`: Flux uncertainty
4. `PHASE`: Orbital phase [0-1] (optional, calculated if missing)
5. `ECL_MASK`: Eclipse mask [0/1] where 0=out-of-eclipse, 1=in-eclipse (optional)

**Example:**
```
# TIME FLUX FLUX_ERR PHASE ECL_MASK
1544.451890067498 1.002350687980652 0.001127733965404 0.360831473378058 0.0
1544.472723991579 1.000355958938599 0.001121848821640 0.361600594988388 0.0
1544.493557915660 0.999854624271393 0.001116090640426 0.362369715598718 0.0
...
```

#### Format 2: NumPy arrays (`.npz`)

Binary NumPy archive format containing multiple arrays.

**Required keys:**
- `time`: Time array (BTJD - 2457000)
- `flux`: Normalized flux array
- `flux_err`: Flux uncertainty array

**Optional keys:**
- `phase`: Orbital phase array [0-1]
- `eclipse_mask`: Eclipse mask array (boolean)

**Example:**
```python
import numpy as np

data = {
    'time': np.array([...]),           # BTJD - 2457000
    'flux': np.array([...]),           # Normalized flux (median = 1.0)
    'flux_err': np.array([...]),       # Flux uncertainty
    'phase': np.array([...]),          # Orbital phase (optional)
    'eclipse_mask': np.array([...])    # Boolean eclipse mask (optional)
}
np.savez('TIC_260128333_06.npz', **data)
```


## Output Data Formats

### Transit Finding Output

#### Event File

**Default filename:** `transit_events.txt`

Space-delimited text file with header.

**Columns for 'cb' (cosine+biweight) method:**
```
TIC SECTOR TIME PHASE DEPTH DURATION SNR WIN_LEN_MAX_SNR DET_DEPENDENCE SKYE_FLAG
```

**Columns for 'cp' (cosine+pspline) method:**
```
TIC SECTOR TIME PHASE DEPTH DURATION SNR SKYE_FLAG
```

**Column descriptions:**

| Column | Description | Units/Format |
|--------|-------------|--------------|
| `TIC` | TIC ID | string |
| `SECTOR` | TESS sector number | string (no leading zeros) |
| `TIME` | Event time | BJD - 2457000 |
| `PHASE` | Orbital phase of event | [0-1] |
| `DEPTH` | Transit depth | Fractional flux |
| `DURATION` | Transit duration | Days |
| `SNR` | Signal-to-noise ratio | Rounded to 2 decimals |
| `WIN_LEN_MAX_SNR` | Biweight window length with max SNR | Days (cb method only) |
| `DET_DEPENDENCE` | Detrending dependence flag | 0=robust, 1=dependent (cb method only) |
| `SKYE_FLAG` | Systematic artifact flag | 1=flagged (if sector_times provided) |

**Example:**
```
TIC SECTOR TIME PHASE DEPTH DURATION SNR WIN_LEN_MAX_SNR DET_DEPENDENCE SKYE_FLAG
260128333 6 1471.5741463040006 0.22426096261004425 0.0010733890846215122 0.22916727320694008 3.67 1.7 1 0
260128333 6 1483.9074888870716 0.06851414847620152 0.0031062932199997784 0.2708334206708969 14.62 2.3 0 0
```

#### Event Snippets (Optional)

**Directory:** `{output_dir}/event_snippets/`

**Filename convention:** `TIC_{tic}_{sector}_{event_no}.npz`

NumPy archive containing extracted event data for model comparison.

**Keys:**
- `tic`: TIC ID (int)
- `sector`: Sector number (int)
- `event_no`: Event number (int)
- `time`: Time array for event window
- `flux`: Detrended flux array
- `flux_err`: Flux error array
- `event_time`: Mid-transit time
- `event_width`: Event duration

#### Vetting Plots (Optional)

If `generate_vetting_plots=True`:
- Diagnostic plots saved for each detected event
- Shows detrended light curve with event highlighted
- Displays residuals and statistical metrics

### Model Comparison Output

**Default filename:** `model_comparison_results.csv`

CSV file containing model comparison results for each event.

**Columns:**
```
filename,best_fit,aic_transit,aic_sinusoidal,aic_linear,aic_step,rmse_transit,rmse_sinusoidal,rmse_linear,rmse_step
```

**Column descriptions:**

| Column | Description |
|--------|-------------|
| `filename` | Event snippet filename |
| `best_fit` | Best-fit model classification |
| `aic_transit` | Akaike Information Criterion (AIC) for transit model |
| `aic_sinusoidal` | AIC for sinusoidal model |
| `aic_linear` | AIC for linear model |
| `aic_step` | AIC for step function model |
| `rmse_transit` | Root mean square error (RMSE) for transit model |
| `rmse_sinusoidal` | RMSE for sinusoidal model |
| `rmse_linear` | RMSE for linear model |
| `rmse_step` | RMSE for step function model |

**Classification categories:**

| Code | Meaning | Criteria |
|------|---------|----------|
| `T` | Transit | Transit is best fit, AIC diff ≥2, RMSE ≤ threshold |
| `AT` | Ambiguous Transit | Transit is best fit, AIC diff ≥2, RMSE > threshold |
| `Sin` | Sinusoid | Sinusoidal model is best fit |
| `ASin` | Ambiguous Sinusoid | Sinusoidal is best but RMSE > threshold |
| `L` | Linear | Linear model is best fit |
| `AL` | Ambiguous Linear | Linear is best but RMSE > threshold |
| `St` | Step | Step function is best fit |
| `ASt` | Ambiguous Step | Step is best but RMSE > threshold |
| `A` | Ambiguous | AIC difference <2 between models |

**Example:**
```csv
filename,best_fit,aic_transit,aic_sinusoidal,aic_linear,aic_step,rmse_transit,rmse_sinusoidal,rmse_linear,rmse_step
TIC_260128333_6_1,A,-506.75,-503.21,-499.98,-505.23,1.078,1.195,1.277,1.202
TIC_260128333_6_2,T,-568.30,-524.37,-345.00,-523.33,0.940,1.401,2.412,1.343
```

### Injection-Retrieval Output

#### Results File

**Default filename:** `inj-ret_results.csv`

CSV file containing injection and recovery results.

**Columns:**
```
tics,sectors,injected_times,injected_depths,injected_durations,injected_snrs,recovered,recovered_times,recovered_depths,recovered_durations,recovered_snrs
```

**Column descriptions:**

| Column | Description | Notes |
|--------|-------------|-------|
| `tics` | TIC ID | int |
| `sectors` | Sector number | int |
| `injected_times` | Injection time | BJD - 2457000 |
| `injected_depths` | Injected transit depth | Fractional flux |
| `injected_durations` | Injected transit duration | Days |
| `injected_snrs` | SNR of injected transit |  |
| `recovered` | Recovery flag | True/False |
| `recovered_times` | Detected event time | NaN if not recovered |
| `recovered_depths` | Detected transit depth | NaN if not recovered |
| `recovered_durations` | Detected transit duration | NaN if not recovered |
| `recovered_snrs` | Detected SNR | NaN if not recovered |

**Example:**
```csv
tics,sectors,injected_times,injected_depths,injected_durations,injected_snrs,recovered,recovered_times,recovered_depths,recovered_durations,recovered_snrs
260128333,7,1511.7512196791108,0.001,0.1,,False,,,,
260128333,6,1486.990820220012,0.001,0.1,2.698,False,,,,
260128333,6,1486.990820220012,0.0014678,0.1,3.964,True,1486.991,0.00147,0.099,3.95
```

#### Statistics File

**Default filename:** `inj-ret_results_stats.csv`

CSV file containing summary statistics for each transit model.

**Columns:**
```
model_idx,depth,duration,n_injections,n_recoveries,recovery_rate
```

**Column descriptions:**

| Column | Description |
|--------|-------------|
| `model_idx` | Transit model index |
| `depth` | Model transit depth (fractional) |
| `duration` | Model transit duration (days) |
| `n_injections` | Total number of injections for this model |
| `n_recoveries` | Number of successful recoveries |
| `recovery_rate` | Recovery rate (n_recoveries / n_injections) |

**Example:**
```csv
model_idx,depth,duration,n_injections,n_recoveries,recovery_rate
0,0.001,0.1,2,0,0.0
1,0.001,0.25,2,0,0.0
7,0.0014678,0.1,2,2,1.0
15,0.0014678,0.25,2,2,1.0
```