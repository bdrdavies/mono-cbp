# mono-cbp: Search for Monotransits of Circumbinary Planets

A comprehensive Python pipeline for detecting circumbinary planets in TESS eclipsing binary light curves through the identification of single transit events (monotransits).

## Overview

**mono-cbp** is a research-grade pipeline designed to systematically search for circumbinary planets by detecting individual transit signatures in TESS eclipsing binary systems. It combines automated eclipse masking, advanced transit detection algorithms, and Bayesian model comparison to distinguish genuine transit events from systematic artifacts.

## Key Features

- **Eclipse Masking**: Automatically mask primary and secondary eclipses in eclipsing binary light curves with configurable eclipse phase and width parameters
- **Transit Detection**: Advanced detrending (cosine filtering, biweight, wotan) and single-event detection using the monofind algorithm
- **Bayesian Model Comparison**: Probabilistic event classification using PyMC and Exoplanet (transit vs. systematic artifacts)
- **Injection-Retrieval Testing**: Completeness analysis through synthetic transit injection and recovery statistics
- **Modular Architecture**: Use individual components independently or run the complete integrated pipeline
- **Configuration-Driven**: Easily customize parameters via Python dictionaries without modifying code
- **Command-Line Interface**: Shell scripts and CLI subcommands for batch processing and reproducibility

## Installation

### Requirements

- Python 3.8 or higher
- pip package manager
- ~500 MB disk space (with dependencies)

### From Source

```bash
git clone https://github.com/bdrdavies/mono-cbp.git
cd mono-cbp
pip install -e .
```

### Core Dependencies

The main dependencies include:

**Numerical & Data Processing:**
- `numpy` - Numerical computations
- `pandas` - Data manipulation and CSV handling
- `scipy` - Scientific computing utilities

**Astronomical Data:**
- `astropy` - Astronomical calculations and units
- `lightkurve` - TESS light curve handling

**Detrending & Filtering:**
- `wotan` - Advanced detrending methods
- `matplotlib` - Visualization

**Bayesian Inference:**
- `pymc` (v5.12.0+) - Probabilistic modeling and MCMC sampling
- `exoplanet` - Exoplanet-specific models and utilities

See [requirements.txt](requirements.txt) for the complete dependency list and exact versions.

### Troubleshooting Installation

If you encounter issues with PyMC or Exoplanet installation:
```bash
# Install pre-compiled wheels (recommended)
pip install --upgrade pip setuptools wheel
pip install -e .

# Or install with conda (if using Anaconda)
conda install -c conda-forge pymc exoplanet
pip install -e .
```

## Quick Start

### 1. Complete Pipeline

Run the full pipeline with all steps:

```python
from mono_cbp import MonoCBPPipeline

# Initialize pipeline with your data
pipeline = MonoCBPPipeline(
    catalogue_path='catalogues/tebc.csv',
    data_dir='./lightcurves',
    config={'transit_finding': {'mad_threshold': 4.0}}
)

# Execute complete pipeline
# Eclipse masking → Detrending → Transit finding → Model comparison
results = pipeline.run(
    find_transits=True,
    vet_candidates=True,
    generate_plots=True
)

# Access results
print(f"Transit events found: {len(results['events'])}")
print(f"Classified as transits: {results['summary']['transit_count']}")
```

### 2. Individual Components

#### Eclipse Masking

Automatically mask primary and secondary eclipses based on orbital ephemerides:

```python
from mono_cbp import EclipseMasker
from mono_cbp.utils import load_catalogue

# Load catalogue with eclipse parameters
catalogue = load_catalogue('catalogues/tebc.csv')

# Create masker and process all systems
masker = EclipseMasker(catalogue, data_dir='./lightcurves')
masker.mask_all()

# Visualize results for a specific system
masker.plot_phase_fold(tic_id=319011894, save_fig=True)
```

#### Transit Detection

Find potential transit events in detrended light curves:

```python
from mono_cbp import TransitFinder

finder = TransitFinder(
    eclipse_params='catalogues/tebc.csv',
    config={
        'transit_finding': {
            'mad_threshold': 4.0,              # Detection sensitivity
            'detrending_method': 'cb',         # 'cb', 'cp', or 'wotan'
            'generate_vetting_plots': True
        }
    }
)

# Process directory of light curves
results = finder.process_directory(
    input_dir='./lightcurves',
    output_file='detected_transits.txt',
    plot_output_dir='./diagnostic_plots'
)

print(f"Events detected: {len(results['events'])}")
```

#### Bayesian Model Comparison

Classify detected events using Bayesian inference (distinguishes transits from artifacts):

```python
from mono_cbp import ModelComparator

comparator = ModelComparator(
    config={
        'model_comparison': {
            'tune': 1000,
            'draws': 1000,
            'chains': 4,
            'cores': 4
        }
    }
)

# Classify a single event
result = comparator.compare_event('events/snippet_1.npz')
print(f"Classification: {result['best_fit']}")  # 'T', 'AT', 'A', 'AN', or 'N'
print(f"Transit probability: {result['transit_prob']:.2%}")

# Batch classify multiple events
classifications = comparator.classify_directory(
    event_dir='./events',
    output_file='classifications.csv'
)
```

#### Injection-Retrieval Testing

Test detection completeness by injecting synthetic transits and recovering them:

```python
from mono_cbp import TransitInjector

injector = TransitInjector(
    transit_models_path='catalogues/transit_models.npz',
    eclipse_params='catalogues/tebc.csv',
    config={'transit_finding': {'mad_threshold': 4.0}}
)

# Run injection-retrieval analysis
results = injector.run_injection_retrieval(
    data_dir='./lightcurves',
    n_injections=100,
    output_file='injection_results.csv'
)

# Summary statistics
print(f"Recovery rate: {injector.stats['recovery_rate']:.1%}")
print(f"False positive rate: {injector.stats['false_pos_rate']:.1%}")
print(f"SNR threshold: {injector.stats['snr_threshold']:.2f}")
```

## Configuration

The pipeline uses a hierarchical configuration system that allows fine-grained control over all processing steps. Default settings are defined in [mono_cbp/config/defaults.py](mono_cbp/config/defaults.py).

### Default Configuration

```python
config = {
    'transit_finding': {
        # Detection sensitivity (MAD-based threshold)
        'mad_threshold': 3.0,                 # Sigma threshold for event detection
        'detrending_method': 'cb',            # 'cb' (cosine), 'cp' (biweight), or 'wotan'
        'generate_vetting_plots': False,      # Generate diagnostic plots

        # Cosine detrending parameters
        'cosine': {
            'win_len_max': 12,                # Maximum window length (hours)
            'win_len_min': 1,                 # Minimum window length (hours)
        },

        # Biweight detrending parameters
        'biweight': {
            'win_len_max': 3,                 # Maximum window length (hours)
            'win_len_min': 1,                 # Minimum window length (hours)
        },

        # Event filtering
        'filters': {
            'min_snr': 5,                     # Minimum signal-to-noise ratio
            'max_duration_days': 1,           # Maximum event duration (days)
        }
    },

    'model_comparison': {
        # MCMC sampling parameters (PyMC)
        'tune': 1000,                         # Tuning steps
        'draws': 1000,                        # Posterior samples
        'chains': 4,                          # Number of chains
        'cores': 4,                           # CPU cores to use

        # Classification thresholds
        'aic_threshold': 10,                  # AIC difference threshold for model preference
        'rmse_threshold': 0.005,              # RMSE threshold for artifact detection
    },

    'injection_retrieval': {
        'n_injections_per_file': 10,          # Synthetic transits per light curve
        'injection_depths': [0.001, 0.005, 0.01],  # Transit depths to test (fractional)
    }
}
```

### Custom Configuration

Override defaults by passing a custom configuration dictionary:

```python
from mono_cbp import TransitFinder

# Custom configuration for high sensitivity
custom_config = {
    'transit_finding': {
        'mad_threshold': 2.5,                 # Lower threshold = more sensitive
        'generate_vetting_plots': True,
        'filters': {
            'min_snr': 3,                     # Lower SNR limit
            'max_duration_days': 2,
        }
    },
    'model_comparison': {
        'tune': 2000,                         # More tuning steps
        'draws': 2000,                        # More posterior samples
    }
}

finder = TransitFinder(
    eclipse_params='catalogues/tebc.csv',
    config=custom_config
)
```

### Configuration Hierarchy

Configurations are merged hierarchically:
1. **Hard defaults**: Built-in defaults in `mono_cbp/config/defaults.py`
2. **Custom config**: User-provided configuration (merges with defaults)
3. **Component-level**: Individual component initialization parameters (override previous levels)

This ensures backward compatibility while allowing full customization.

## Input Data Format

### Catalogue CSV

The catalogue should contain eclipse and orbital parameters for each system. Required columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `tess_id` | int | TESS Input Catalog ID | 319011894 |
| `period` | float | Orbital period (days) | 5.234 |
| `bjd0` | float | Reference epoch (BJD - 2457000) | 1325.456 |
| `sectors` | str | TESS sectors (comma-separated) | 1,2,3,4 |
| `prim_pos` | float | Primary eclipse phase position [0-1] | 0.0 |
| `prim_width` | float | Primary eclipse phase width [0-1] | 0.05 |
| `sec_pos` | float | Secondary eclipse phase position [0-1] | 0.5 |
| `sec_width` | float | Secondary eclipse phase width [0-1] | 0.03 |

Example catalogue snippet:
```csv
tess_id,period,bjd0,sectors,prim_pos,prim_width,sec_pos,sec_width
319011894,5.234,1325.456,"1,2,3,4",0.0,0.05,0.5,0.03
327837473,3.891,1255.123,"5,6,7",0.02,0.04,0.48,0.04
```

### Light Curve Files

Light curves should be in one of two formats:

#### Format 1: Text files (`.txt`)

Space or tab-delimited with required columns:
- `TIME`: Observation time (BJD - 2457000)
- `FLUX`: Normalized flux [dimensionless]
- `FLUX_ERR`: Flux uncertainty [dimensionless]
- `PHASE`: Orbital phase [0-1] (optional, calculated if missing)
- `ECL_MASK`: Eclipse mask [0/1] (optional, calculated if missing)

**Filename convention:** `TIC_{tic_id}_S{sector}.txt`

Example:
```
TIME FLUX FLUX_ERR PHASE ECL_MASK
1325.456 0.9998 0.0001 0.0 0
1325.457 0.9997 0.0001 0.01 0
1325.458 0.9996 0.0001 0.02 0
...
```

#### Format 2: NumPy arrays (`.npz`)

Binary format containing arrays accessible by name:

```python
import numpy as np

data = {
    'time': np.array([...]),           # BJD - 2457000
    'flux': np.array([...]),           # Normalized flux
    'flux_err': np.array([...]),       # Flux uncertainty
    'phase': np.array([...]),          # Orbital phase (optional)
    'eclipse_mask': np.array([...])    # Eclipse mask (optional)
}
np.savez('TIC_319011894_S1.npz', **data)
```

**Filename convention:** `TIC_{tic_id}_S{sector}.npz`

### Directory Structure Example

```
project_root/
├── catalogues/
│   ├── tebc.csv              # Main catalogue
│   └── transit_models.npz    # Precomputed transit models
├── lightcurves/
│   ├── TIC_319011894_S1.txt
│   ├── TIC_319011894_S2.npz
│   ├── TIC_327837473_S1.txt
│   └── ...
├── results/
│   ├── detected_transits.txt
│   ├── classifications.csv
│   └── diagnostic_plots/
└── README.md
```

## Output Formats

### Transit Finding Output

**Event file** (`.txt`, space-delimited):
```
TIC_ID SECTOR TIME PHASE DEPTH DURATION_HOURS SNR
319011894 1 1327.234 0.234 0.0045 2.3 7.2
319011894 1 1329.456 0.456 0.0052 2.1 8.1
327837473 1 1328.123 0.145 0.0038 1.9 6.5
```

**Columns:**
- `TIC_ID`: System identifier
- `SECTOR`: TESS sector
- `TIME`: Event time (BJD - 2457000)
- `PHASE`: Orbital phase of event [0-1]
- `DEPTH`: Transit depth [fractional]
- `DURATION_HOURS`: Event duration [hours]
- `SNR`: Signal-to-noise ratio

**Optional vetting plots** (if `generate_vetting_plots=True`):
- Diagnostic plots for each detected event
- Shows detrended light curve with event highlighted
- Displays residuals and statistical metrics

**Optional event snippets** (`.npz`):
- Extracted event data for model comparison
- Contains time, flux, and uncertainty around event

### Model Comparison Output

**Classifications file** (`.csv`):
```
EVENT_ID,TIC_ID,SECTOR,BEST_FIT,TRANSIT_PROB,AIC_TRANSIT,AIC_ARTIFACT,NOTES
1,319011894,1,T,0.95,-245.3,-235.1,High confidence
2,319011894,1,AT,0.68,-230.5,-233.2,Ambiguous
3,327837473,1,N,0.12,-200.1,-215.3,Systematic artifact
```

**Classification categories:**
- `T`: Unambiguous transit (high posterior probability for transit model)
- `AT`: Ambiguous transit (moderate confidence, could be artifact)
- `A`: Ambiguous (equal probability for multiple models)
- `AN`: Ambiguous non-transit (leans toward artifact)
- `N`: Not a transit (high confidence it's an artifact)

**Columns:**
- `EVENT_ID`: Unique event identifier
- `TIC_ID`: System identifier
- `SECTOR`: TESS sector
- `BEST_FIT`: Most likely classification
- `TRANSIT_PROB`: Posterior probability of transit model
- `AIC_TRANSIT`: AIC for transit model fit
- `AIC_ARTIFACT`: AIC for artifact model fit
- `NOTES`: Additional interpretation notes

### Injection-Retrieval Output

**Results file** (`.csv`):
```
DEPTH_INJECTED,DEPTH_RECOVERED,SNR_INJECTED,RECOVERED,TIC_ID,SECTOR
0.001,0.0009,3.2,1,319011894,1
0.001,,,0,319011894,1
0.005,0.0048,7.1,1,319011894,1
...
```

**Summary statistics** (printed or logged):
```
Injection-Retrieval Statistics
================================
Total injections: 100
Successful recoveries: 89
Recovery rate: 89.0%
False positive rate: 2.1%
SNR threshold (50% recovery): 4.2

By depth:
  0.001 frac: 72% recovery
  0.005 frac: 91% recovery
  0.010 frac: 98% recovery
```

**Columns:**
- `DEPTH_INJECTED`: Synthetic transit depth [fractional]
- `DEPTH_RECOVERED`: Detected transit depth (if recovered)
- `SNR_INJECTED`: SNR of injected transit
- `RECOVERED`: Binary flag (1 = recovered, 0 = missed)
- `TIC_ID`: System identifier
- `SECTOR`: TESS sector

## Command-Line Interface

The package provides convenient shell scripts for batch processing:

### Using Shell Scripts

Shell scripts in the `scripts/` directory provide command-line interfaces:

```bash
# Run the complete pipeline
./scripts/run_pipeline.sh \
  --catalogue catalogues/tebc.csv \
  --data-dir ./lightcurves \
  --output-dir ./results \
  --mad-threshold 4.0

# Eclipse masking only
./scripts/mask_eclipses.sh \
  --catalogue catalogues/tebc.csv \
  --data-dir ./lightcurves \
  --output-dir ./results

# Transit finding only
./scripts/find_transits.sh \
  --catalogue catalogues/tebc.csv \
  --data-dir ./lightcurves \
  --output-dir ./results \
  --plots

# Model comparison on detected events
./scripts/compare_models.sh \
  --event-dir ./results/event_snippets \
  --output-file ./results/classifications.csv

# Injection-retrieval testing
./scripts/inject_retrieve.sh \
  --catalogue catalogues/tebc.csv \
  --data-dir ./lightcurves \
  --models catalogues/transit_models.npz \
  --output-file ./results/injection_stats.csv
```

For detailed script options, examine the scripts or run with `--help`:
```bash
./scripts/run_pipeline.sh --help
```

## Project Structure

The repository is organized as follows:

```
mono-cbp/
├── mono_cbp/                      # Main Python package
│   ├── __init__.py                # Package initialization, public API
│   ├── __version__.py             # Version information
│   ├── cli.py                     # Command-line interface
│   ├── pipeline.py                # Main orchestration class
│   │
│   ├── eclipse_masking/           # Eclipse masking module
│   │   ├── __init__.py
│   │   ├── masker.py             # EclipseMasker class
│   │   └── utils.py
│   │
│   ├── transit_finding/           # Transit detection module
│   │   ├── __init__.py
│   │   ├── finder.py             # TransitFinder class
│   │   └── utils.py
│   │
│   ├── model_comparison/          # Bayesian model comparison module
│   │   ├── __init__.py
│   │   ├── comparator.py         # ModelComparator class
│   │   └── models.py             # Transit/artifact models
│   │
│   ├── injection_retrieval/       # Injection-retrieval module
│   │   ├── __init__.py
│   │   ├── injector.py           # TransitInjector class
│   │   └── utils.py
│   │
│   ├── utils/                     # Utility modules
│   │   ├── data.py               # Data loading, catalogue handling
│   │   ├── eclipses.py           # Eclipse calculations
│   │   ├── detrending.py         # Detrending algorithms
│   │   ├── monofind.py           # Monofind algorithm
│   │   ├── transit_models.py     # Transit model management
│   │   ├── plotting.py           # Visualization utilities
│   │   └── common.py             # General utilities
│   │
│   └── config/                    # Configuration
│       └── defaults.py            # Default configuration
│
├── examples/                      # Jupyter notebooks
│   ├── 01_complete_pipeline.ipynb
│   ├── 02_eclipse_masking.ipynb
│   ├── 03_transit_finding.ipynb
│   ├── 04_model_comparison.ipynb
│   └── 05_injection_retrieval.ipynb
│
├── scripts/                       # Shell scripts for CLI
│   ├── run_pipeline.sh
│   ├── mask_eclipses.sh
│   ├── find_transits.sh
│   ├── compare_models.sh
│   └── inject_retrieve.sh
│
├── catalogues/                    # Data files
│   ├── tebc.csv                  # TESS Eclipsing Binary Catalogue
│   └── transit_models.npz        # Precomputed transit models
│
├── data/                          # Sample light curves
│   ├── TIC_319011894_S1.txt
│   ├── TIC_319011894_S2.npz
│   └── ...
│
├── docs/                          # Documentation
│   └── api_reference.md          # Full API documentation
│
├── README.md                      # This file
├── setup.py                       # Package installation
├── requirements.txt               # Dependencies
└── LICENSE                        # MIT License
```

## Examples & Tutorials

Five Jupyter notebooks in the `examples/` directory demonstrate each component:

1. **[01_complete_pipeline.ipynb](examples/01_complete_pipeline.ipynb)** - End-to-end workflow on sample data
2. **[02_eclipse_masking.ipynb](examples/02_eclipse_masking.ipynb)** - Eclipse masking demonstration
3. **[03_transit_finding.ipynb](examples/03_transit_finding.ipynb)** - Transit detection walkthrough
4. **[04_model_comparison.ipynb](examples/04_model_comparison.ipynb)** - Bayesian classification details
5. **[05_injection_retrieval.ipynb](examples/05_injection_retrieval.ipynb)** - Completeness testing

## Full API Documentation

Complete API reference is available in [docs/api_reference.md](docs/api_reference.md), including:
- Detailed docstrings for all public classes and methods
- Parameter descriptions and type information
- Return value documentation
- Code examples

## Citation

If you use mono-cbp in your research, please cite:

```bibtex
@software{mono_cbp,
  author = {Miller, Ben},
  title = {mono-cbp: A Pipeline for Detecting Circumbinary Planets through Monotransit Detection in TESS Eclipsing Binary Light Curves},
  year = {2024},
  url = {https://github.com/bdrdavies/mono-cbp},
  note = {Research software, University of Warwick}
}
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for:
- Code style guidelines
- Testing requirements
- Pull request process
- Development setup

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

This pipeline builds upon the following open-source projects:

- **[monofind](https://github.com/samuelgill/monofind)** - Single transit detection algorithm
- **[wotan](https://github.com/hippke/wotan)** - Advanced light curve detrending
- **[PyMC](https://github.com/pymc-devs/pymc)** - Probabilistic programming and MCMC sampling
- **[Exoplanet](https://github.com/exoplanet-dev/exoplanet)** - Exoplanet-specific PyMC models
- **[lightkurve](https://docs.lightkurve.org/)** - TESS light curve processing
- **[astropy](https://www.astropy.org/)** - Core astronomical calculations

## Support & Contact

For questions, issues, or feature requests:
- Open an issue on [GitHub Issues](https://github.com/bdrdavies/mono-cbp/issues)
- Check existing [Discussions](https://github.com/bdrdavies/mono-cbp/discussions)
- Review the [Full Documentation](docs/api_reference.md)

For direct inquiries, contact the development team at the University of Warwick.

## Roadmap

Planned enhancements:
- [ ] GPU acceleration for MCMC sampling
- [ ] Extended mission (EM) light curve support
- [ ] Machine learning-based vetting (in development)
- [ ] Interactive vetting dashboard
- [ ] Expanded detrending method library
