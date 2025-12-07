# Configuration Guide

Complete guide to configuring mono-cbp.

## Overview

mono-cbp uses a hierarchical configuration system:
1. Default values defined in `mono_cbp/config/defaults.py`
2. User overrides passed as dictionaries
3. Command-line arguments (for CLI usage)

## Configuration Structure

```python
config = {
    'transit_finding': { ... },
    'model_comparison': { ... },
    'injection_retrieval': { ... }
}
```

## Transit Finding Configuration

### Basic Parameters

```python
config = {
    'transit_finding': {
        'edge_cutoff': 0.0,
        'mad_threshold': 3.0,
        'detrending_method': 'cb',
        'generate_vetting_plots': False,
        'generate_skye_plots': False,
        'generate_event_snippets': True,
        'save_event_snippets': True,
        'cadence_minutes': 30,
        'cosine': {
            'win_len_max': 12,
            'win_len_min': 1,
            'fap_threshold': 1e-2,
            'poly_order': 2,
        },
        'biweight': {
            'win_len_max': 3,
            'win_len_min': 1,
        },
        'pspline': {
            'max_splines': 25,
        },
        'filters': {
            'min_snr': 5,
            'max_duration_days': 1,
            'det_dependence_threshold': 18,
        }
    }
}
```

**Parameters:**

- `edge_cutoff` (float, default: 0.0)
  - Timespan of flux values to exclude at the edges of data gaps
  - In units of days 

- `mad_threshold` (float, default: 3.0)
  - Multiplier for Median Absolute Deviation (MAD) threshold for event detection
  - Higher = fewer false positives, lower sensitivity
  - Recommended range: 3.0 - 5.0

- `detrending_method` (str, default: 'cb')
  - `'cb'`: cosine + biweight (provides detrending-dependence metric)
  - `'cp'`: cosine + pspline

- `generate_vetting_plots` (bool, default: False)
  - Whether to create diagnostic plots for each detection
  - Useful for visual inspection but slow for large datasets

- `generate_skye_plots` (bool, default: False)
  - Whether to generate per-sector Skye metric plots
  - Essentially a histogram of event counts across the timespan of a given sector

- `generate_event_snippets` (bool, default: True)
  - Save event snippet arrays in-memory for analysis

- `save_event_snippets` (bool, default: True)
  - Save event data to disk for later model comparison

- `cadence_minutes` (int, default: 30)
  - Cadence of data in minutes

### Detrending Parameters

#### Cosine Filter

```python
config = {
    'transit_finding': {
        'cosine': {
            'win_len_max': 12.0,
            'win_len_min': 1.0,
            'fap_threshold': 0.01,
            'poly_order': 2,
        }
    }
}
```

**Parameters:**
- `win_len_max` (float, default: 12.0): Maximum detrending window length (days)
- `win_len_min` (float, default: 1.0): Minimum detrending window length (days)
- `fap_threshold` (float, default: 0.01): False alarm probability threshold for Lomb-Scargle periodogram peak
- `poly_order` (int, default: 2): Polynomial order for initial detrending before calculating Lomb-Scargle periodogram

The algorithm automatically selects optimal window size between these bounds by calculating the maximum window length at which no significant periodicities are present in the light curve.

#### Biweight Filter (CB method only)

```python
config = {
    'transit_finding': {
        'biweight': {
            'win_len_max': 3.0,
            'win_len_min': 1.0,
        }
    }
}
```

**Parameters:**
- `win_len_max` (float, default: 3.0): Maximum detrending window length for biweight grid (days)
- `win_len_min` (float, default: 1.0): Miniumum detrending window length for biweight grid (days)

#### Penalised Spline (`pspline`) Filter

```python
config = {
    'transit_finding': {
        'pspline': {
            'max_splines': 25,
        }
    }
}
```

**Parameters:**
- `max_splines` (int, default: 25): Maximum number of splines for penalized spline detrending

### Filtering Parameters

```python
config = {
    'transit_finding': {
        'filters': {
            'min_snr': 5.0,
            'max_duration_days': 1.0,
            'det_dependence_threshold': 18,
        }
    }
}
```

**Parameter:**

- `min_snr` (float, default: 5.0): Signal-to-noise ratio threshold

- `max_duration_days` (float, default: 1.0): Maximum TCE duration (days)

- `det_detection_threshold` (int, default: 18): Threshold for number of biweight-detrended light curves that a TCE is detected for flagging event as "detrending-dependent". E.g., for the default value, if an event is detected in < 18 biweight-detrended light curves, then it is flagged as detrending-dependent. Likewise, if an event is detected in > 18 biweight-detrended light curves, then it is flagged as detrending-independent.

## Model Comparison Configuration

```python
config = {
    'model_comparison': {
        'tune': 1000,
        'draws': 1000,
        'chains': 4,
        'cores': 4,
        'target_accept': 0.99,
        'sigma_threshold': 3,
        'aic_threshold': 2,
        'rmse_threshold': 1.2,
        'save_plots': False,
        'plot_dir': None,
    }
}
```

**Parameters:**

- `tune` (int, default: 1000)
  - Number of tuning/burn-in steps for each Markov chain
  - More tuning = better convergence but slower

- `draws` (int, default: 1000)
  - Number of posterior samples per chain
  - More draws = better statistics but slower

- `chains` (int, default: 4)
  - Number of independent Markov chains
  - More chains = better convergence diagnostics

- `cores` (int, default: 4)
  - CPU cores for parallel chain execution
  - Set to number of physical cores

- `target_accept` (float, default: 0.99)
  - Target acceptance rate for NUTS sampler
  - Higher = smaller step size, slower but more accurate

- `sigma_threshold` (float, default: 3.0)
  - Number of standard deviations for outlier removal before fitting models

- `aic_threshold` (float, default: 2.0)
  - AIC units difference threshold for distinguishing between models

- `rmse_threshold` (float, default: 1.2)
  - Threshold for root mean square error (RMSE) for quantifying model fit quality

- `save_plots` (bool, default: False)
  - Whether to save model comparison diagnostic plots

- `plot_dir` (str, default: None)
  - Directory for saving diagnostic plots

## Injection-Retrieval Configuration

```python
config = {
    'injection_retrieval': {
        'n_injections': 1000,
    }
}
```

**Parameters:**

- `n_injections` (int, default: 1000)
  - Total number of synthetic transits to inject for each transit model
  - More injections = better statistics but slower

## Complete Configuration Example

```python
from mono_cbp import MonoCBPPipeline

# Complete configuration
config = {
    'transit_finding': {
        'edge_cutoff': 0.0,
        'mad_threshold': 3.0,
        'detrending_method': 'cb',
        'generate_vetting_plots': True,
        'generate_skye_plots': False,
        'generate_event_snippets': True,
        'save_event_snippets': True,
        'cadence_minutes': 30,
        'cosine': {
            'win_len_max': 12.0,
            'win_len_min': 1.0,
            'fap_threshold': 0.01,
            'poly_order': 2,
        },
        'biweight': {
            'win_len_max': 3.0,
            'win_len_min': 1.0,
        },
        'pspline': {
            'max_splines': 25,
        },
        'filters': {
            'min_snr': 5.0,
            'max_duration_days': 1.0,
            'det_dependence_threshold': 18,
        }
    },
    'model_comparison': {
        'tune': 1000,
        'draws': 1000,
        'chains': 4,
        'cores': 4,
        'target_accept': 0.99,
        'sigma_threshold': 3,
        'aic_threshold': 2,
        'rmse_threshold': 1.2,
        'save_plots': False,
        'plot_dir': None,
    },
    'injection_retrieval': {
        'n_injections': 200,
    }
}

# Use configuration
pipeline = MonoCBPPipeline(
    catalogue_path='./catalouges/TEBC_morph_05_P_7.csv',
    config=config
)
```

### Loading from JSON

```python
import json
from mono_cbp import MonoCBPPipeline

with open('config.json', 'r') as f:
    config = json.load(f)

pipeline = MonoCBPPipeline(
    catalogue_path='./catalouges/TEBC_morph_05_P_7.csv',
    config=config
)
```

## Command-Line Configuration

Override configuration via command-line:

```bash
mono-cbp find-transits \
    --data-dir ./data \
    --threshold 4.0 \
    --method cb \
```
