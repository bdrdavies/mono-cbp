# Quick Start Guide

This guide aims to get you up and running with `mono-cbp` ASAP.

## Basic Usage

### 1. Prepare Your Data

Organize your light curve files and create a catalogue CSV:

```
project/
├── data/
│   ├── TIC_123456789_01.txt
│   ├── TIC_123456789_02.txt
│   └── ...
├── catalogues
    ├── catalogue.csv
    └── sector_times.csv
```

#### Catalogue Format

`catalogue.csv` should contain:

```csv
tess_id,period,bjd0,sectors,prim_pos,prim_width,sec_pos,sec_width
123456789,2.5,1234.5,"1,2,3",0.0,0.05,0.5,0.05
```

- `tess_id`: TIC ID of the target
- `period`: Priod of the eclipsing binary
- `bjd0`: Reference epoch (time of mid-eclipse for primary eclipse)
- `sectors`: TESS sectors that the target was observed in. Note that this column is only necessary if downloading the data from MAST with `mono_cbp.utils.catalogue_to_lc_files`.
- `prim_pos`: Position of the primary eclipse in phase space (should be very close to either 0 or 1)
- `prim_width`: Width of primary eclipse in phase space
- `sec_pos`: Position of secondary eclipse in phase space. For a circularised orbit, should be very close to 0.5
- `sec_width`: Width of secondary eclipse in phase space

### 2. Run the Complete Pipeline

```python
from mono_cbp import MonoCBPPipeline

# Initialise
pipeline = MonoCBPPipeline(
    catalogue_path='catalogue.csv',
    data_dir='./data',
    sector_times_path='sector_times.csv'
)

# Run (eclipse masking and transit finding always performed)
results = pipeline.run(
    vet_candidates=True,      # Optional: run model comparison vetting
    injection_retrieval=False  # Optional: run injection-retrieval testing
)

# Check results
print(f"Transit finding results: {results['transit_finding']}")
print(f"Vetting results: {results['vetting']}")
```

### 3. Run Individual Components

#### Eclipse Masking Only

```python
from mono_cbp import EclipseMasker
import pandas as pd

catalogue = pd.read_csv('catalogue.csv')
masker = EclipseMasker(catalogue, data_dir='./data')

# Mask all systems
masker.mask_all()

# View results
masker.plot_bin_phase_fold(tic_id=123456789, save_fig=True)
```

#### Transit Finding Only

```python
from mono_cbp import TransitFinder
import pandas as pd

catalogue = pd.read_csv('catalogue.csv')
finder = TransitFinder(
    catalogue=catalogue,
    sector_times='sector_times.csv'
)

# Process all files
results = finder.process_directory(
    'data',
    output_file='transit_events.txt',
    plot_output_dir='./plots'
)

print(f"Detected {len(results)} events")
```

#### Model Comparison Only

```python
from mono_cbp import ModelComparator

comparator = ModelComparator()

# Compare a single event
result = comparator.compare_event('event_snippet.npz')
print(f"Classification: {result['best_fit']}")

# Batch compare multiple events
results = comparator.compare_events(
    'event_snippets/',
    output_file='classifications.csv'
)

print(f"Processed {len(results)} events")
```

#### Injection-Retrieval Testing

```python
from mono_cbp import TransitInjector
import pandas as pd

catalogue = pd.read_csv('catalogue.csv')
injector = TransitInjector(
    transit_models_path='transit_models.npz',
    catalogue=catalogue
)

# Run injection-retrieval
results = injector.run_injection_retrieval(
    'data',
    n_injections=100,
    output_file='injection_retrieval_results.csv'
)

print(f"Injection-retrieval complete: {len(results)} results")
```

## Customizing Configuration

Override default settings:

```python
config = {
    'transit_finding': {
        'mad_threshold': 4.0,  # More conservative threshold
        'generate_vetting_plots': True,
        'filters': {
            'min_snr': 7,  # Higher SNR requirement
            'max_duration_days': 0.5
        }
    },
    'model_comparison': {
        'tune': 2000,  # More tuning steps
        'draws': 2000,
        'cores': 8
    }
}

pipeline = MonoCBPPipeline(
    catalogue_path='catalogue.csv',
    config=config
)
```

## Command-Line Interface

Run from the terminal (if CLI is available):

```bash
# Full pipeline
mono-cbp run --catalogue catalogue.csv --data-dir ./data

# Individual components
mono-cbp mask-eclipses --catalogue catalogue.csv --data-dir ./data
mono-cbp find-transits --catalogue catalogue.csv --data-dir ./data
mono-cbp compare-models --event-dir ./events
mono-cbp inject-retrieve --models models.npz --catalogue catalogue.csv --data-dir ./data

# With options
mono-cbp find-transits \
    --catalogue catalogue.csv \
    --data-dir ./data \
    --threshold 4.0 \
    --output transit_events.txt
```

**Note:** For detailed configuration options, use Python API or edit configuration files rather than CLI flags.

## Next Steps

- Read the [User Guide](user_guide.md) for detailed explanations
- Explore the Jupyter notebooks in examples/ for common use cases
- Check [Configuration](configuration.md) for all available options
- Review the [API Reference](api_reference.md) for complete documentation