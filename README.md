# mono-cbp: Search for Monotransits of Circumbinary Planets

A comprehensive Python pipeline for detecting circumbinary planets in TESS eclipsing binary light curves through the identification of single transit events.

## Features

- **Eclipse Masking**: Automatically mask primary and secondary eclipses in eclipsing binary light curves
- **Transit Detection**: Advanced detrending and event detection using monofind algorithm
- **Model Comparison**: Bayesian model comparison for event classification (transit vs. systematic artifacts)
- **Injection-Retrieval**: Completeness testing through synthetic transit injection
- **Modular Design**: Use individual components or run the complete pipeline
- **Configuration-Driven**: Easily customize parameters without modifying code

## Installation

### From Source

```bash
git clone https://github.com/bdrdavies/mono-cbp.git
cd mono-cbp
pip install -e .
```

### Dependencies

The main dependencies are:
- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `astropy`
- `lightkurve`
- `pymc`
- `exoplanet`
- `wotan`

See [requirements.txt](requirements.txt) for the complete list.

## Quick Start

### Complete Pipeline

```python
from mono_cbp import MonoCBPPipeline

# Initialize pipeline
pipeline = MonoCBPPipeline(
    catalogue_path='eclipse_params.csv',
    data_dir='./lightcurves',
    orbit_params_path='orbit_params.csv',
    sector_times_path='sector_times.csv'
)

# Run full pipeline
# Eclipse masking is always performed automatically
results = pipeline.run(
    find_transits=True,
    vet_candidates=True
)
```

### Individual Components

#### Eclipse Masking

```python
from mono_cbp import EclipseMasker

masker = EclipseMasker(catalogue, data_dir='./data')
masker.mask_all()
masker.plot_phase_fold(tic_id=319011894, save_fig=True)
```

#### Transit Finding

```python
from mono_cbp import TransitFinder

finder = TransitFinder(
    eclipse_params='eclipse_params.csv',
    orbit_params='orbit_params.csv',
    sector_times='sector_times.csv',
    config={'transit_finding': {'mad_threshold': 4.0}}
)

results = finder.process_directory(
    'data',
    output_file='transit_events.txt',
    plot_output_dir='./plots'
)
```

#### Model Comparison

```python
from mono_cbp import ModelComparator

comparator = ModelComparator()
result = comparator.compare_event('event_snippet.npz')
print(f"Classification: {result['best_fit']}")
```

#### Injection-Retrieval

```python
from mono_cbp import TransitInjector

injector = TransitInjector(
    transit_models_path='transit_models.npz',
    eclipse_params='eclipse_params.csv',
    orbit_params='orbit_params.csv'
)

results = injector.run_injection_retrieval(
    'data',
    n_files=100,
    output_file='injection_results.csv'
)
print(f"Recovery rate: {injector.stats['recovery_rate']:.2%}")
```

## Configuration

The pipeline uses a hierarchical configuration system. Default settings are in `mono_cbp/config/defaults.py`:

```python
config = {
    'transit_finding': {
        'mad_threshold': 3.0,
        'detrending_method': 'cb',  # 'cb' or 'cp'
        'generate_vetting_plots': False,
        'cosine': {
            'win_len_max': 12,
            'win_len_min': 1,
        },
        'biweight': {
            'win_len_max': 3,
            'win_len_min': 1,
        },
        'filters': {
            'min_snr': 5,
            'max_duration_days': 1,
        }
    },
    'model_comparison': {
        'tune': 1000,
        'draws': 1000,
        'chains': 4,
        'cores': 4,
    }
}
```

Override defaults by passing a config dictionary:

```python
my_config = {
    'transit_finding': {
        'mad_threshold': 4.0,
        'generate_vetting_plots': True
    }
}

pipeline = MonoCBPPipeline(
    catalogue_path='catalogue.csv',
    config=my_config
)
```

## Input Data Format

### Catalogue CSV

Required columns:
- `tess_id`: TIC ID
- `period`: Orbital period (days)
- `bjd0`: Reference epoch (BJD - 2457000)
- `sectors`: TESS sectors (comma-separated)
- `prim_pos`: Primary eclipse phase position
- `prim_width`: Primary eclipse phase width
- `sec_pos`: Secondary eclipse phase position
- `sec_width`: Secondary eclipse phase width

### Light Curve Files

Supported formats:
- `.txt` files with columns: TIME, FLUX, FLUX_ERR, PHASE, ECL_MASK
- `.npz` files with arrays: time, flux, flux_err, phase, eclipse_mask

Filename convention: `TIC_{tic_id}_S{sector}.txt` or `.npz`

## Output

### Transit Finding
- Text file with detected events (TIC, sector, time, phase, depth, duration, SNR)
- Optional vetting plots showing detected events
- Optional event snippets (.npz) for model comparison

### Model Comparison
- CSV file with classifications:
  - `T`: Unambiguous transit
  - `AT`: Ambiguous transit
  - `A`: Ambiguous
  - `AN`: Ambiguous non-transit
  - `N`: Not a transit

### Injection-Retrieval
- CSV file with injection/recovery statistics

## Command-Line Interface

```bash
# Run full pipeline
mono-cbp run --catalogue catalogue.csv --data-dir ./data

# Eclipse masking only
mono-cbp mask-eclipses --catalogue catalogue.csv --data-dir ./data

# Transit finding only
mono-cbp find-transits --eclipse-params eclipse.csv --data-dir ./data

# Model comparison
mono-cbp compare-models --event-dir ./event_snippets

# Injection-retrieval
mono-cbp inject-retrieve --models transit_models.npz --data-dir ./data
```

## Documentation

Full documentation is available at [link to docs] (coming soon).

## Citation

If you use this software in your research, please cite:

```bibtex
@software{mono_cbp,
  author = {Miller, Ben},
  title = {mono-cbp: Circumbinary Planet Detection Pipeline},
  year = {2024},
  url = {https://github.com/bdrdavies/mono-cbp}
}
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Transit detection algorithm based on [monofind](https://github.com/samuelgill/monofind)
- Detrending uses [wotan](https://github.com/hippke/wotan)
- Model fitting with [PyMC](https://github.com/pymc-devs/pymc) and [exoplanet](https://github.com/exoplanet-dev/exoplanet)

## Contact

For questions or issues, please open an issue on GitHub or contact [your email].
