# API Reference

Complete reference for all mono-cbp modules and classes.

## Core Modules

### Pipeline

#### `MonoCBPPipeline`

Main pipeline class coordinating all components.

```python
from mono_cbp import MonoCBPPipeline

pipeline = MonoCBPPipeline(
    catalogue_path: str,
    data_dir: str = './data',
    output_dir: str = './results',
    sector_times_path: str = '../../catalogues/sector_times.csv',
    TEBC: bool = False,
    transit_models_path: str = None,
    config: dict = None
)
```

**Parameters:**
- `catalogue_path` (str): Path to CSV containing catalogue of EBs
- `data_dir` (str, optional): Directory containing light curve files (default: './data')
- `output_dir` (str, optional): Directory to produce output files (default: './results')
- `sector_times_path` (str, optional): Path to file containing the start and end times of TESS sectors (for Skye metric calculation, default: '../../catalogues/sector_times.csv')
- `TEBC` (bool, optional): Flag for TESS Eclipsing Binary Catalogue format (default: False)
- `transit_models_path` (str, optional): Path to transit models for injection-retrieval (default: None)
- `config` (dict, optional): Configuration dictionary (default: None)

**Methods:**

##### `run()`
```python
results = pipeline.run(
    find_transits: bool = True,
    vet_candidates: bool = True,
    injection_retrieval: bool = False,
    **kwargs
) -> dict
```

Run the complete pipeline. Eclipse masking is always performed first.

**Parameters:**
- `find_transits` (bool, optional): Run transit finding (default: True)
- `vet_candidates` (bool, optional): Run model comparison vetting (default: True)
- `injection_retrieval` (bool, optional): Run injection-retrieval test (default: False)
- `**kwargs`: Additional arguments passed to pipeline steps:
  - `mask_eclipses_kwargs`: Arguments for eclipse masking
  - `find_transits_kwargs`: Arguments for transit finding
  - `vet_candidates_kwargs`: Arguments for vetting
  - `injection_retrieval_kwargs`: Arguments for injection-retrieval

**Returns:** Dictionary with results from each pipeline step

##### `mask_eclipses()`
```python
pipeline.mask_eclipses(**kwargs) -> None
```

Execute eclipse masking step (appends eclipse mask as a column to light curve files).

##### `find_transits()`
```python
results = pipeline.find_transits(
    output_file: str = 'transit_events.txt',
    output_dir: str = None,
    plot_output_dir: str = None
) -> pd.DataFrame
```

Execute transit finding step.

**Parameters:**
- `output_file` (str, optional): Filename for transit events output (default: 'transit_events.txt')
- `output_dir` (str, optional): Directory to produce `output_file` (default: None)
- `plot_output_dir` (str, optional): Directory to save diagnostic plots (default: None)

**Returns:** DataFrame with transit detection results

##### `vet_candidates()`
```python
results = pipeline.vet_candidates(
    event_snippets: list = None,
    event_snippets_dir: str = None,
    output_file: str = 'vetting_results.csv',
    output_dir: str = None
) -> pd.DataFrame
```

Execute model comparison vetting.

**Parameters**
- `event_snippets` (list, optional): List of event snippet dictionaries to process in-memory. If provided, takes precedence over `event_snippets_dir` (default: None)
- `event_snippets_dir` (str, optional): Directory with event snippet files (if event snippets stored on disk, default: None)
- `output_file` (str, optional): Output filename (default: 'vetting_results.csv')
- `output_dir` (str, optional): Output directory (default: None)

**Returns:** DataFrame with classification results

##### `run_injection_retrieval()`
```python
results = pipeline.run_injection_retrieval(
    n_injections: int = 100,
    output_file: str = 'injection_results.csv',
    output_dir: str = None
) -> pd.DataFrame
```

Run injection-retrieval testing for all transit models in `transit_models.npz`. Tests each transit model by injecting it into `n_injections` randomly selected light curves. The total number of tests will be `n_injections × number_of_models`.

**Parameters:**
- `n_injections` (int, optional): Number of injection-retrieval tests to perform per transit model. If there are fewer files than requested injections, files will be randomly sampled with replacement (default: 100)
- `output_file` (str, optional): Output filename for results (default: 'injection_results.csv')
- `output_dir` (str, optional): Output directory. If None, defaults to data_dir (default: None)

**Returns:** DataFrame with injection-retrieval results (one row per injection test)

**Raises:**
- `ValueError`: If transit injector not initialized (requires `transit_models_path` in constructor)

##### `plot_bin_phase_fold()`
```python
pipeline.plot_bin_phase_fold(
    tic_id: int,
    save_fig: bool = False,
    save_path: str = '.'
) -> None
```

Plot phase-folded light curve showing eclipse mask for visual assessment.

**Parameters:**
- `tic_id` (int): TIC identifier to plot
- `save_fig` (bool, optional): Whether to save figure to disk (default: False)
- `save_path` (str, optional): Directory to save figure (default: '.')

**Returns:** None

##### `plot_events()`
```python
pipeline.plot_events(
    tic_id: int,
    event_number: int = None,
    save_fig: bool = False,
    save_path: str = '.',
    figsize: tuple = (12, 4)
) -> None
```

Plot detected transit events for a given TIC ID.

**Parameters:**
- `tic_id` (int): TIC ID to plot events for
- `event_number` (int, optional): Specific event number to plot (1-indexed). If None, plots all events for this TIC (default: None)
- `save_fig` (bool, optional): Whether to save the figure (default: False)
- `save_path` (str, optional): Directory to save figure (default: '.')
- `figsize` (tuple, optional): Figure size as (width, height) per event (default: (12, 4))

**Returns:** None

**Raises:**
- `ValueError`: If no events found for the given TIC ID or event number
- `RuntimeError`: If transit finding has not been run yet

**Note:** Requires `generate_event_snippets: True` in config during transit finding.

**Example:**
```python
# Plot all events for a TIC
pipeline.plot_events(260128333)

# Plot only the first event
pipeline.plot_events(260128333, event_number=1, save_fig=True)
```

### Eclipse Masking

#### `EclipseMasker`

Masks eclipses in EB light curves. Calculates eclipse masks from binary ephemeris and eclipse parameters in the catalogue. Supports input files with or without pre-calculated phase data.

```python
from mono_cbp import EclipseMasker

masker = EclipseMasker(
    catalogue: str or pd.DataFrame,
    data_dir: str = './data',
    TEBC: bool = False
)
```

**Parameters:**
- `catalogue` (str or DataFrame): Path to or DataFrame of catalogue with eclipse parameters (prim_pos, prim_width, sec_pos, sec_width) and ephemerides (period, bjd0). If a DataFrame is passed that has already been processed (contains standard eclipse parameter columns), it will be used as-is.
- `data_dir` (str, optional): Directory containing light curve files (default: './data')
- `TEBC` (bool, optional): If True, processes TEBC catalogue format with *_2g and *_pf columns and converts to standard eclipse parameter columns. If a DataFrame is passed that already has standard columns, TEBC processing is skipped (default: False)

**Raises:**
- `FileNotFoundError`: If data_dir does not exist
- `NotADirectoryError`: If data_dir is not a directory
- `ValueError`: If catalogue file cannot be loaded or processed

**Methods:**

##### `mask_all()`
```python
masker.mask_all(
    force: bool = False
) -> None
```

Mask eclipses for all systems in catalogue. Skips files that cannot be processed (e.g., missing TIC ID in catalogue, unsupported format).

**Parameters:**
- `force` (bool, optional): If True, recalculate masks even if they already exist in input files (default: False)

##### `mask_file()`
```python
masker.mask_file(
    file: str,
    force: bool = False
) -> None
```

Mask eclipses in a single light curve file. Loads light curve data (time, flux, flux_err) and optionally phase. If phase is not provided, it is calculated from the catalogue ephemeris. Computes eclipse masks for primary and secondary eclipses and saves the results with an appended eclipse_mask column.

**Input file format:**
- **.npz**: Must contain 'time', 'flux', 'flux_err' keys. Optionally 'phase' and 'eclipse_mask'.
- **.txt**: Must have columns for time, flux, flux_err (and optionally phase, eclipse_mask). Phase can be calculated from ephemeris if not provided.

**Parameters:**
- `file` (str): Name of file to process (including file extension). Filename should be in format 'TIC_<TICID>_<sector>.<ext>'
- `force` (bool, optional): If True, recalculate mask even if it already exists (default: False)

**Returns:**
- `None`. Returns early without processing if an eclipse mask already exists in the input file and force=False.

**Raises:**
- `FileNotFoundError`: If file does not exist
- `ValueError`: If file format is unsupported, TIC ID not found in catalogue, or required data is missing

##### `plot_bin_phase_fold()`
```python
masker.plot_bin_phase_fold(
    tic_id: int,
    save_fig: bool = False,
    save_path: str = '.'
) -> None
```

Plot phase-folded light curve with eclipse masks highlighted. Loads and combines all data files for the given TIC ID. If phase is not provided in the input files, it is calculated from the catalogue ephemeris. If eclipse masks exist, applies polynomial detrending to out-of-eclipse data before plotting for visibility.

**Parameters:**
- `tic_id` (int): TIC identifier
- `save_fig` (bool, optional): If True, save the figure to disk (default: False)
- `save_path` (str, optional): Directory to save figure (default: '.')

**Raises:**
- `FileNotFoundError`: If no files found for the given TIC ID
- `ValueError`: If no valid data can be loaded for plotting

**Note:**
- If the TIC ID is not found in the catalogue, the period will be displayed as 'unknown'.
- Data is automatically sorted by time before plotting.
- Polynomial detrending is applied to out-of-eclipse data if eclipse masks are available.
- Plots with eclipse mask distinction if mask data is present, otherwise plots without distinction.
- If `save_fig` is False, displays the plot to screen. If True, saves to file and closes the plot.

### Transit Finding

#### `TransitFinder`

Detrends input data and identifies threshold crossing events (TCEs) using the `monofind` algorithm. Supports flexible catalogue input (file paths or DataFrames) with optional TEBC format processing.

```python
from mono_cbp import TransitFinder

finder = TransitFinder(
    catalogue: str or pd.DataFrame = None,
    sector_times: str or pd.DataFrame = None,
    config: dict = None,
    TEBC: bool = False
)
```

**Parameters:**
- `catalogue` (str or DataFrame, optional): Path to or DataFrame of catalogue with eclipse parameters (prim_pos, prim_width, sec_pos, sec_width) and binary ephemerides (period, bjd0). If a DataFrame is passed that has already been processed (contains standard eclipse parameter columns), it will be used as-is (default: None)
- `sector_times` (str or DataFrame, optional): Path to or DataFrame of sector times CSV for Skye metric calculation (default: None)
- `config` (dict, optional): Configuration dictionary. Uses defaults if None (default: None)
- `TEBC` (bool, optional): If True, processes TEBC catalogue format with *_2g and *_pf columns and converts to standard eclipse parameter columns. If a DataFrame is passed that already has standard columns, TEBC processing is skipped (default: False)

**Raises:**
- `ValueError`: If catalogue path cannot be loaded or catalogue data is invalid.
- `FileNotFoundError`: If catalogue or sector_times file path does not exist.
- `KeyError`: If configuration dictionary is missing required keys.

**Methods:**

##### `process_directory()`
```python
results = finder.process_directory(
    data_dir: str,
    output_file: str = 'output.txt',
    output_dir: str = None,
    plot_output_dir: str = None
) -> pd.DataFrame
```

Process all light curves in a directory.

**Parameters:**
- `data_dir` (str): Directory containing light curve files
- `output_file` (str, optional): Output filename (default: 'output.txt')
- `output_dir` (str, optional): Directory to save output file (default: None)
- `plot_output_dir` (str, optional): Directory to save plots (default: None)

**Returns:** DataFrame of detected events. Returns empty DataFrame if no events are detected.

**Raises:**
- `FileNotFoundError`: If data_dir does not exist.
- `PermissionError`: If no read permission for data_dir or write permission for output directories.

##### `process_file()`
```python
events = finder.process_file(
    file_path: str,
    plot_output_dir: str = None
) -> list
```

Process a single light curve file.

**Parameters:**
- `file_path` (str): Path to light curve file
- `plot_output_dir` (str, optional): Directory to save plots (default: None)

**Returns:** List of detected event dictionaries. Each event dictionary contains keys: 'time', 'phase', 'depth', 'width', 'duration', 'start_time', 'end_time', 'snr', 'tic', 'sector'. Returns empty list if file cannot be processed or no events are detected.

**Raises:**
- `ValueError`: If filename cannot be parsed.
- `FileNotFoundError`: If file_path does not exist.

##### `save_results()`
```python
results_df = finder.save_results(
    output_file: str = 'output.txt',
    output_dir: str = None
) -> pd.DataFrame
```

Save detected events to a text file and return as DataFrame.

**Parameters:**
- `output_file` (str, optional): Output filename (default: 'output.txt')
- `output_dir` (str, optional): Directory to save output file (defaults to current directory)

**Returns:** DataFrame with all detected events

**Output Format:**
- Text file with space-separated columns
- Columns depend on detrending method:
  - `cb` method: TIC, SECTOR, TIME, PHASE, DEPTH, DURATION, SNR, WIN_LEN_MAX_SNR, DET_DEPENDENCE, SKYE_FLAG (if available)
  - `cp` method: TIC, SECTOR, TIME, PHASE, DEPTH, DURATION, SNR, SKYE_FLAG (if available)
- Event snippets saved to `event_snippets/` subdirectory if `generate_event_snippets` and `save_event_snippets` are enabled in config

**Note:** All detected events are saved without filtering. Use `filter_events()` to apply quality cuts.

##### `filter_events()` (static method)
```python
filtered_df = TransitFinder.filter_events(
    events_df: pd.DataFrame,
    min_snr: float = None,
    max_duration_days: float = None,
    det_dependence_flag: int = None,
    skye_flag: int = None
) -> pd.DataFrame
```

Filter detected events based on quality criteria. This is a static method that can be called without creating a TransitFinder instance.

**Parameters:**
- `events_df` (DataFrame): DataFrame of detected events (from transit finding)
- `min_snr` (float, optional): Minimum signal-to-noise ratio threshold (default: None)
- `max_duration_days` (float, optional): Maximum transit duration in days (default: None)
- `det_dependence_flag` (int, optional): Filter by detrending dependence flag (0 = keep robust events, 1 = keep only flagged events, default: None)
- `skye_flag` (int, optional): Filter by Skye metric flag (0 = keep non-systematic events, 1 = keep systematic events, default: None)

**Returns:** Filtered DataFrame containing only events that pass all specified criteria

**Examples:**
```python
# Apply standard quality filters
filtered = TransitFinder.filter_events(events_df, min_snr=5.0, max_duration_days=1.0)

# Filter for high-quality events (not flagged by any metric)
filtered = TransitFinder.filter_events(
    events_df,
    min_snr=5.0,
    det_dependence_flag=0,  # Robust to detrending
    skye_flag=0             # Not systematic artifact
)

# Filter by SNR only
filtered = TransitFinder.filter_events(events_df, min_snr=7.0)

# Can also call as instance method
finder = TransitFinder(...)
filtered = finder.filter_events(events_df, min_snr=5.0)
```

**Note:** This method does not modify the original DataFrame. Flags are only available if corresponding analysis was performed (e.g., `det_dependence_flag` requires `cb` detrending method, `skye_flag` requires sector times).

**Implementation Notes:**

The TransitFinder class includes several internal constants for configuration:
- `LONG_CADENCE_THRESHOLD_DAYS` (0.0138889): Threshold for detecting long cadence data (~2 minutes)
- `CADENCE_MINUTES_TO_DAYS` (1440): Conversion factor from minutes to days
- `VAR_MAD_WINDOW` (100): Window size for variable median absolute deviation calculation
- `PROGRESS_INTERVAL` (10): Progress logging interval for file processing
- `EVENT_WINDOW_HALF_WIDTH` (0.5): Half-width of event window in days for short-duration events
- `EVENT_GROUPING_TOLERANCE` (10): Cadence-point tolerance for grouping events across detrending window lengths
- `SKYE_FIGURE_SIZE` ((12, 10)): Figure size for Skye metric histograms

Internal helper methods:
- `_load_npz()`: Loads light curve data from NPZ format files
- `_load_txt()`: Loads light curve data from TXT format files
- `_parse_filename()`: Extracts TIC ID and sector from filename
- `_get_eclipse_params()`: Retrieves eclipse parameters from catalogue
- `_extract_event_data()`: Calculates SNR and event metadata
- `_create_event_snippet()`: Creates time-windowed event data for model comparison
- `_calculate_skye_metric()`: Calculates systematic artifact flags across sectors

### Model Comparison

#### `ModelComparator`

Bayesian model comparison for event classification.

```python
from mono_cbp import ModelComparator

comparator = ModelComparator(
    config: dict = None
)
```

**Parameters:**
- `config` (dict, optional): Configuration for MCMC sampling, thresholds, and plotting (default: None)

**Methods:**

##### `compare_event()`
```python
result = comparator.compare_event(
    event_input: str or dict,
    save_plot: bool = False,
    plot_dir: str = None
) -> dict
```

Compare models for a single event using Bayesian model fitting and AIC comparison.

**Parameters:**
- `event_input` (str or dict): Event data input, either:
  - **File path** (str): Path to event snippet (.npz file)
  - **Dictionary**: Event data with required keys: `'time'`, `'flux'`, `'flux_err'`, `'event_time'`, `'event_width'`. Optional keys: `'tic'`, `'sector'`, `'event_no'` (for logging/plotting)
- `save_plot` (bool, optional): Whether to save diagnostic plot showing model fits (default: False)
- `plot_dir` (str, optional): Directory to save plot. If None, defaults to `'./plots'` (default: None)

**Returns:** Dictionary with classification and model comparison results:
- `filename` (str): Event filename (for file input) or constructed filename (for dict input)
- `best_fit` (str): Classification category:
  - `'T'`: Unambiguous transit (transit model is best fit, RMSE ≈ 1)
  - `'AT'`: Ambiguous transit (transit model is best fit, RMSE > 1)
  - `'A'`: Ambiguous (AIC difference < 2 between best models)
  - `'AN'`: Ambiguous non-transit (transit not best fit, RMSE ≈ 1)
  - `'N'`: Not a transit (transit not best fit, RMSE > 1)
- `aic_transit` (float): AIC for transit model
- `aic_sinusoidal` (float): AIC for sinusoidal model
- `aic_linear` (float): AIC for linear model
- `aic_step` (float): AIC for step function model
- `rmse_transit` (float): Root mean square error for transit model
- `rmse_sinusoidal` (float): RMSE for sinusoidal model
- `rmse_linear` (float): RMSE for linear model
- `rmse_step` (float): RMSE for step function model

**Model Fitting:**
This method fits four models to the event data:
1. **Transit model**: Exoplanet transit using PyMC3/Exoplanet
2. **Sinusoidal model**: Sinusoidal variation (for eclipsing binary remnants)
3. **Linear model**: Linear trend
4. **Step function model**: Step function (for systematic artifacts)

**Example:**
```python
# From file
result = comparator.compare_event('event_snippets/TIC_12345_S10_1.npz')

# From memory
event_data = {
    'time': time_array,
    'flux': flux_array,
    'flux_err': flux_err_array,
    'event_time': 1234.56,
    'event_width': 0.1,
    'tic': 12345,
    'sector': 10,
    'event_no': 1
}
result = comparator.compare_event(event_data, save_plot=True, plot_dir='./plots')

print(f"Classification: {result['best_fit']}")
print(f"Transit RMSE: {result['rmse_transit']:.2f}")
```

##### `compare_events()`
```python
results = comparator.compare_events(
    events_input: str or list,
    output_file: str = 'model_comparison_results.csv',
    output_dir: str = None,
    save_plots: bool = None,
    plot_dir: str = None
) -> pd.DataFrame
```

Compare models for multiple events (batch processing).

**Parameters:**
- `events_input` (str or list): Event data input, either:
  - **Directory path** (str): Path to directory containing event snippet files (.npz)
  - **List of dictionaries**: List of event data dictionaries with required keys: `'time'`, `'flux'`, `'flux_err'`, `'event_time'`, `'event_width'`. Optional keys: `'tic'`, `'sector'`, `'event_no'`
- `output_file` (str, optional): Output CSV filename (default: 'model_comparison_results.csv')
- `output_dir` (str, optional): Directory to save output file. If None and `events_input` is a directory, defaults to `events_input` directory (default: None)
- `save_plots` (bool, optional): Whether to save diagnostic plots for each event. If None, uses value from config (default: None)
- `plot_dir` (str, optional): Directory to save plots. If None, uses value from config or falls back to `output_dir` (default: None)

**Returns:** DataFrame with classification results for all events

**DataFrame Columns:**
- `filename` (str): Event filename or constructed name
- `best_fit` (str): Classification (T, AT, A, AN, N)
- `aic_transit` (float): AIC for transit model
- `aic_sinusoidal` (float): AIC for sinusoidal model
- `aic_linear` (float): AIC for linear model
- `aic_step` (float): AIC for step function model
- `rmse_transit` (float): RMSE for transit model
- `rmse_sinusoidal` (float): RMSE for sinusoidal model
- `rmse_linear` (float): RMSE for linear model
- `rmse_step` (float): RMSE for step function model

**Behavior:**
- Processes events in batches with progress logging (every 10 events)
- Handles errors gracefully - failed events are logged and skipped
- Saves results to CSV if `output_dir` is specified
- Logs classification summary statistics at completion

**Examples:**
```python
# Process directory of event snippet files
results_df = comparator.compare_events(
    events_input='results/event_snippets',
    output_file='vetting_results.csv',
    output_dir='results',
    save_plots=True,
    plot_dir='results/model_plots'
)

# Process in-memory event data (from TransitFinder)
finder = TransitFinder(...)
finder.process_directory(...)
event_snippets = finder.results['event_snippets']

results_df = comparator.compare_events(
    events_input=event_snippets,
    output_file='vetting_results.csv',
    output_dir='results'
)

# View classification summary
print(results_df['best_fit'].value_counts())
```

**Note:** Private methods `_fit_transit_model()`, `_fit_sinusoidal_model()`, `_fit_linear_model()`, and `_fit_step_model()` are used internally for model fitting.

### Injection-Retrieval

#### `TransitInjector`

Synthetic transit injection and recovery testing.

```python
from mono_cbp import TransitInjector

injector = TransitInjector(
    transit_models_path: str,
    catalogue: pd.DataFrame or str = None,
    config: dict = None,
    TEBC: bool = False
)
```

**Parameters:**
- `transit_models_path` (str): Path to transit models (.npz file)
- `catalogue` (DataFrame or str, optional): Catalogue with eclipse parameters. Can be either a DataFrame or path to CSV file. If TEBC format, set `TEBC=True` (default: None)
- `config` (dict, optional): Configuration dictionary (default: None)
- `TEBC` (bool, optional): If True, use TEBC catalogue format with eclipse parameter selection logic (default: False)

**Methods:**

##### `run_injection_retrieval()`
```python
results = injector.run_injection_retrieval(
    data_dir: str,
    n_injections: int = 100,
    output_file: str = 'inj-ret_results.csv',
    output_dir: str = None
) -> pd.DataFrame
```

Run injection-retrieval tests for all transit models in `transit_models.npz`. Tests each transit model by injecting it into `n_injections` randomly selected light curves. The total number of tests will be `n_injections × number_of_models`.

**Parameters:**
- `data_dir` (str): Directory containing light curve files
- `n_injections` (int, optional): Number of injection-retrieval tests to perform per transit model. If there are fewer files than requested injections, files will be randomly sampled with replacement (default: 100)
- `output_file` (str, optional): Output filename (default: 'inj-ret_results.csv')
- `output_dir` (str, optional): Output directory. If None, defaults to data_dir (default: None)

**Returns:** DataFrame with injection-retrieval results containing (one row per injection test):
- TIC ID and sector
- Injected transit parameters (time, depth, duration, SNR)
- Recovery status (boolean)
- Recovered parameters (time, depth, duration, SNR) if detected, NaN otherwise

**Behavior:**
- Light curves are inverted to ensure detected events are injected transits, not real signals
- Gaps in light curves are temporarily filled to allow injections across gaps
- Each test uses a randomly selected injection time from the available light curve
- Per-model recovery statistics are automatically saved to `{output_file}_stats.csv`
- Per-model statistics are also accessible via `injector.stats` after completion

**Output Files:**
- `{output_file}`: Individual injection-retrieval test results (one row per test)
- `{output_file}_stats.csv`: Per-model recovery statistics with columns:
  - `model_idx` (int): Model index
  - `depth` (float): Transit depth
  - `duration` (float): Transit duration in days
  - `n_injections` (int): Number of injections for this model
  - `n_recoveries` (int): Number of successful recoveries
  - `recovery_rate` (float): Recovery rate (n_recoveries / n_injections)

**Note:** The `_inject_transit()` method is used internally for synthetic transit injection.

## Utility Modules

### `mono_cbp.utils.data`

Catalogue and data file utilities.

#### Functions

##### `load_catalogue()`
```python
from mono_cbp.utils.data import load_catalogue

catalogue = load_catalogue(
    catalogue_path: str,
    TEBC: bool = False
) -> pd.DataFrame
```

Load eclipsing binary catalogue from CSV file.

**Parameters:**
- `catalogue_path` (str): Path to catalogue CSV
- `TEBC` (bool): Whether file is in TESS Eclipsing Binary Catalog format (default: False)

**Returns:** DataFrame with catalogue data

##### `get_row()`
```python
from mono_cbp.utils.data import get_row

row = get_row(
    catalogue: pd.DataFrame,
    tic_id: int
) -> dict or None
```

Get catalogue entry for a specific TIC ID.

**Returns:** Dictionary with row data, or None if not found

##### `bin_to_long_cadence()`
```python
from mono_cbp.utils.data import bin_to_long_cadence

time_binned, flux_binned, flux_err_binned = bin_to_long_cadence(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]
```

Bin short-cadence TESS data to 30-minute long cadence using lightkurve.

**Parameters:**
- `time` (np.ndarray): Time values to bin
- `flux` (np.ndarray): Flux values to bin
- `flux_err` (np.ndarray): Flux error values to bin

**Returns:** Tuple of (binned_time, binned_flux, binned_flux_err)

**Note:** Bins to 30 minutes (hardcoded), not a configurable parameter.

### `mono_cbp.utils.eclipses`

Eclipse calculation and masking utilities.

#### Functions

##### `time_to_phase()`
```python
from mono_cbp.utils.eclipses import time_to_phase

phase = time_to_phase(
    time: np.ndarray,
    period: float,
    bjd0: float,
    centre: float = 0.5
) -> np.ndarray
```

Calculate orbital phase from time.

**Parameters:**
- `time` (np.ndarray): Time values (BJD)
- `period` (float): Orbital period in days
- `bjd0` (float): Reference epoch (BJD)
- `centre` (float, optional): Centre of the phase fold (default: 0.5)

**Returns:** Orbital phase values (0-1)

**Note:** When `centre=0.5`, phase 0.5 corresponds to the reference epoch. When `centre=0.0`, phase 0.0 corresponds to the reference epoch.

##### `get_eclipse_mask()`
```python
from mono_cbp.utils.eclipses import get_eclipse_mask

mask = get_eclipse_mask(
    phase: np.ndarray,
    eclipse_pos: float,
    eclipse_width: float
) -> np.ndarray
```

Create boolean mask for eclipse events.

**Parameters:**
- `phase` (np.ndarray): Orbital phase values
- `eclipse_pos` (float): Eclipse center position (0-1)
- `eclipse_width` (float): Eclipse width in phase units

**Returns:** Boolean array (True = in eclipse)

##### `get_eclipse_indices()`
```python
from mono_cbp.utils.eclipses import get_eclipse_indices

indices = get_eclipse_indices(
    phase: np.ndarray,
    eclipse_pos: float,
    eclipse_width: float
) -> np.ndarray
```

Get indices of points within eclipse.

**Returns:** Array of indices for points in eclipse

### `mono_cbp.utils.detrending`

Light curve detrending functions.

#### Functions

##### `detrend()`
```python
from mono_cbp.utils.detrending import detrend

detrended, trend = detrend(
    time: np.ndarray,
    flux: np.ndarray,
    eclipse_mask: np.ndarray = None,
    method: str = 'cb',
    config: dict = None
) -> tuple[np.ndarray, np.ndarray]
```

Detrend light curve using specified method (cosine-biweight or cosine-plus).

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `eclipse_mask` (np.ndarray, optional): Boolean mask for eclipses
- `method` (str): Detrending method ('cb' or 'cp', default: 'cb')
- `config` (dict, optional): Configuration parameters

**Returns:** Tuple of (detrended_flux, trend_flux)

##### `cosine_detrend()`
```python
from mono_cbp.utils.detrending import cosine_detrend

detrended = cosine_detrend(
    flux: np.ndarray,
    win_len: float = 5.0
) -> np.ndarray
```

Apply cosine filter detrending to flux.

**Parameters:**
- `flux` (np.ndarray): Flux values to detrend
- `win_len` (float): Window length in days

**Returns:** Detrended flux array

##### `run_multi_biweight()`
```python
from mono_cbp.utils.detrending import run_multi_biweight

detrended = run_multi_biweight(
    flux: np.ndarray,
    win_len_max: float = 3.0,
    win_len_min: float = 1.0
) -> np.ndarray
```

Apply biweight detrending with optimal window selection.

**Returns:** Detrended flux array

### `mono_cbp.utils.monofind`

Core monofind algorithm implementation for detecting single transit events.

#### Functions

##### `monofind()`
```python
from mono_cbp.utils.monofind import monofind

peaks, meta = monofind(
    time: np.ndarray,
    flux: np.ndarray,
    mad: float = 3.0,
    var_mad: np.ndarray = None
) -> tuple
```

Detect single transit-like events in detrended light curve using MAD threshold.

**Parameters:**
- `time` (np.ndarray): Time values of the light curve (masked)
- `flux` (np.ndarray): Detrended flux values
- `mad` (float, optional): MAD multiplier for detection threshold (default: 3.0)
- `var_mad` (np.ndarray, optional): Variable MAD values calculated using rolling window for adaptive thresholding. If None, uses constant MAD (default: None)

**Returns:** Tuple of (event_indices, metadata_dict) where:
- `event_indices` (list): Indices of detected event peaks in the masked time/flux arrays
- `metadata_dict` (dict): Dictionary containing:
  - `'depths'`: Event depths (flux dip magnitudes)
  - `'widths'`: Event widths (durations in days)
  - `'start_times'`: Event start times
  - `'end_times'`: Event end times

**Algorithm:**
Identifies dips in the light curve that exceed the MAD-based threshold, groups nearby points into events, and calculates event properties.

### `mono_cbp.utils.plotting`

Plotting utilities for light curves and events.

#### Functions

##### `plot_event()`
```python
from mono_cbp.utils.plotting import plot_event

fig = plot_event(
    time: np.ndarray,
    flux: np.ndarray,
    event_time: float,
    window: float = 1.0,
    title: str = None
) -> matplotlib.figure.Figure
```

Plot zoomed view of detected event.

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `event_time` (float): Time of event center
- `window` (float): Window size in days (default: 1.0)
- `title` (str, optional): Plot title

**Returns:** matplotlib Figure object

### `mono_cbp.utils.common`

General utility functions.

#### Functions

##### `setup_logging()`
```python
from mono_cbp.utils.common import setup_logging

logger = setup_logging(
    level: int = logging.INFO,
    log_file: str = None
) -> logging.Logger
```

Configure logging for mono-cbp package.

**Parameters:**
- `level` (int, optional): Logging level (e.g., logging.INFO, logging.DEBUG) (default: logging.INFO)
- `log_file` (str, optional): Path to log file. If None, logs to console only (default: None)

**Returns:** Configured logger instance

##### `get_snr()`
```python
from mono_cbp.utils.common import get_snr

snr = get_snr(
    depth: float,
    error: float,
    duration: float,
    cadence: int = 30
) -> float
```

Calculate signal-to-noise ratio for a transit event.

**Parameters:**
- `depth` (float): Transit depth (fractional flux decrease)
- `error` (float): Combined measurement uncertainty per cadence point
- `duration` (float): Transit duration in days
- `cadence` (int, optional): Observation cadence in minutes (default: 30)

**Returns:** Signal-to-noise ratio value

**Formula:**
SNR = (depth / error) × sqrt(N_cadences) where N_cadences = duration × (1440 / cadence)

### `mono_cbp.utils.transit_models`

Transit model creation and handling.

#### Functions

##### `create_transit_models()`
```python
from mono_cbp.utils.transit_models import create_transit_models

models = create_transit_models(
    radii: np.ndarray,
    periods: np.ndarray,
    stellar_params: dict,
    config: dict = None
) -> dict
```

Create synthetic transit models for injection-retrieval testing.

**Returns:** Dictionary with model parameters and synthetic transits

##### `save_transit_models()`
```python
from mono_cbp.utils.transit_models import save_transit_models

save_transit_models(
    models: dict,
    output_path: str
) -> None
```

Save transit models to .npz file.

##### `load_transit_models()`
```python
from mono_cbp.utils.transit_models import load_transit_models

models = load_transit_models(
    models_path: str
) -> dict
```

Load transit models from .npz file.

**Returns:** Dictionary with transit models

## Configuration

### Default Configuration

```python
from mono_cbp.config import get_default_config

config = get_default_config()
```

Returns the default configuration dictionary.

### Configuration Structure

```python
config = {
    'transit_finding': {
        'edge_cutoff': float,
        'mad_threshold': float,
        'detrending_method': str,  # 'cb' or 'cp'
        'generate_vetting_plots': bool,
        'generate_skye_plots': bool,
        'generate_event_snippets': bool,
        'save_event_snippets': bool,
        'cadence_minutes': int,
        'cosine': {
            'win_len_max': float,
            'win_len_min': float,
            'fap_threshold': float,
            'poly_order': int,
        },
        'biweight': {
            'win_len_max': float,
            'win_len_min': float,
        },
        'pspline': {
            'max_splines': int,
        },
        'filters': {
            'min_snr': 5,
            'max_duration_days': 1,
            'det_dependence_threshold': 18,
        }
    },
    'model_comparison': {
        'tune': int,
        'draws': int,
        'chains': int,
        'cores': int,
        'target_accept': float,
        'sigma_threshold': float,
        'aic_threshold': float,
        'save_plots': bool,
        'plot_dir': str,
    },
    'injection_retrieval': {
        'n_injections': int,
    }
}
```

## Type Hints

The package uses type hints throughout. Key types:

```python
from typing import Union, Optional, List, Dict, Tuple
import numpy as np
import pandas as pd

PathLike = Union[str, Path]
ArrayLike = Union[np.ndarray, List, Tuple]
```
