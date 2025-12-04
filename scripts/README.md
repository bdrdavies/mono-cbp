# mono-cbp Scripts

Shell scripts for running mono-cbp from the command line.

## Prerequisites

1. Install mono-cbp:
```bash
pip install -e .
```

2. Make scripts executable (Unix/macOS):
```bash
chmod +x scripts/*.sh
```

## Scripts

### `run_pipeline.sh`
Run the complete mono-cbp pipeline (eclipse masking → transit finding → model comparison).

**Usage:**
```bash
./scripts/run_pipeline.sh [catalogue] [data_dir]
```

**Example:**
```bash
./scripts/run_pipeline.sh data/catalogues/eclipse_params.csv data/lightcurves
```

**Default paths:**
- Catalogue: `data/catalogues/eclipse_params.csv`
- Data directory: `data/lightcurves`
- Output: `results/`

---

### `mask_eclipses.sh`
Mask primary and secondary eclipses in eclipsing binary light curves.

**Usage:**
```bash
./scripts/mask_eclipses.sh [catalogue] [input_dir] [output_dir]
```

**Example:**
```bash
./scripts/mask_eclipses.sh \
    data/catalogues/eclipse_params.csv \
    data/lightcurves/raw \
    data/lightcurves/masked
```

**Default paths:**
- Catalogue: `data/catalogues/eclipse_params.csv`
- Input: `data/lightcurves/raw`
- Output: `data/lightcurves/masked`

---

### `find_transits.sh`
Detect transit events using the monofind algorithm.

**Usage:**
```bash
./scripts/find_transits.sh [data_dir] [output_file]
```

**Example:**
```bash
./scripts/find_transits.sh data/lightcurves/masked results/transit_events.txt
```

**Configuration:**
Edit the script to change:
- `MAD_THRESHOLD`: Detection threshold (default: 3.5)
- `METHOD`: Detrending method, `cb` or `cp` (default: `cb`)

**Default paths:**
- Data directory: `data/lightcurves/masked`
- Output: `results/transit_events.txt`

---

### `compare_models.sh`
Run Bayesian model comparison to classify detected events.

**Usage:**
```bash
./scripts/compare_models.sh [event_dir] [output_file]
```

**Example:**
```bash
./scripts/compare_models.sh results/event_snippets results/classifications.csv
```

**Default paths:**
- Event directory: `results/event_snippets`
- Output: `results/classifications.csv`

**Note:** Requires event snippets from `find_transits.sh`

---

### `inject_retrieve.sh`
Run injection-retrieval testing to characterize completeness.

**Usage:**
```bash
./scripts/inject_retrieve.sh [data_dir] [output_file]
```

**Example:**
```bash
./scripts/inject_retrieve.sh data/lightcurves/masked results/injection_results.csv
```

**Configuration:**
Edit the script to change:
- `N_INJECTIONS`: Number of injections (default: 100)
- `MODELS_FILE`: Transit models file (default: `data/transit_models.npz`)

**Default paths:**
- Data directory: `data/lightcurves/masked`
- Output: `results/injection_results.csv`

---

## Example Workflows

### Complete pipeline from scratch

```bash
# 1. Mask eclipses
./scripts/mask_eclipses.sh

# 2. Find transits
./scripts/find_transits.sh

# 3. Compare models
./scripts/compare_models.sh

# Or run all at once:
./scripts/run_pipeline.sh
```

### Reprocess with different parameters

Edit `find_transits.sh` to change threshold:
```bash
MAD_THRESHOLD=4.0  # More conservative
```

Then run:
```bash
./scripts/find_transits.sh data/lightcurves/masked results/transits_conservative.txt
```

### Test detection completeness

```bash
# Create transit models first (Python)
python -c "
from mono_cbp.utils import create_transit_models
import numpy as np

models = create_transit_models(
    radii=np.linspace(0.5, 2.5, 15),
    periods=np.logspace(0, 2, 25),
    stellar_params={'R_star': 1.0, 'M_star': 1.5}
)

np.savez('data/transit_models.npz', **models)
print('Transit models created')
"

# Run injection-retrieval
./scripts/inject_retrieve.sh
```

## Using Python API Directly

For more control, use the Python API instead:

```python
from mono_cbp import MonoCBPPipeline

pipeline = MonoCBPPipeline(
    catalogue_path='data/catalogues/eclipse_params.csv',
    data_dir='data/lightcurves',
    config={'transit_finding': {'mad_threshold': 4.0}}
)

results = pipeline.run()
```

See the [examples/](../examples/) directory for Jupyter notebooks.

## Using mono-cbp CLI

All scripts use the `mono-cbp` command-line tool. You can also call it directly:

```bash
# Run complete pipeline
mono-cbp run --catalogue catalogue.csv --data-dir data/

# Individual steps
mono-cbp mask-eclipses --catalogue catalogue.csv --data-dir data/raw
mono-cbp find-transits --eclipse-params eclipse.csv --data-dir data/masked
mono-cbp compare-models --event-dir event_snippets/
mono-cbp inject-retrieve --models models.npz --data-dir data/
```

For help:
```bash
mono-cbp --help
mono-cbp run --help
```

## Configuration Files

Create a `config.json` file for custom parameters:

```json
{
  "transit_finding": {
    "mad_threshold": 4.0,
    "detrending_method": "cb",
    "filters": {
      "min_snr": 7.0,
      "max_duration_days": 0.8
    }
  },
  "model_comparison": {
    "tune": 1500,
    "draws": 1500,
    "chains": 4
  }
}
```

Then run:
```bash
./scripts/run_pipeline.sh
```

The scripts will automatically use `config.json` if it exists.

## Troubleshooting

### Permission denied

On Unix/macOS:
```bash
chmod +x scripts/*.sh
```

### Command not found: mono-cbp

Install mono-cbp:
```bash
pip install -e .
```

Or use full path:
```bash
python -m mono_cbp.cli run --catalogue catalogue.csv --data-dir data/
```

### No events detected

- Check MAD threshold (try lower value like 3.0)
- Verify eclipse masking worked correctly
- Check input data quality

### Model comparison crashes

- Reduce MCMC samples in configuration
- Check available memory
- Process fewer events at once

## Windows Users

On Windows, use:
```bash
bash scripts/run_pipeline.sh
```

Or call the Python CLI directly:
```bash
python -m mono_cbp.cli run --catalogue catalogue.csv --data-dir data/
```

## Getting Help

- See [documentation](../docs/)
- Check [examples](../examples/)
- Open an issue on GitHub
- Run scripts with no arguments for usage information
