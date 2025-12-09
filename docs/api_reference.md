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
    vet_candidates: bool = True,
    injection_retrieval: bool = False,
    **kwargs
) -> dict
```

Run the complete pipeline. Eclipse masking and transit finding are always performed.

**Parameters:**
- `vet_candidates` (bool, optional): Run model comparison vetting (default: True)
- `injection_retrieval` (bool, optional): Run injection-retrieval test (default: False)
- `**kwargs`: Additional arguments passed to pipeline steps:
  - `mask_eclipses_kwargs`: Arguments for eclipse masking
  - `find_transits_kwargs`: Arguments for transit finding
  - `vet_candidates_kwargs`: Arguments for vetting
  - `injection_retrieval_kwargs`: Arguments for injection-retrieval

**Returns:** Dictionary with results from each pipeline step

**Note:** Eclipse masking and transit finding are mandatory steps that always run. Only vetting and injection-retrieval are optional.

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
    n_injections: int = None,
    output_file: str = 'inj-ret_results.csv',
    output_dir: str = None,
    plot_completeness: bool = False,
    completeness_kwargs: dict = None
) -> pd.DataFrame
```

Run injection-retrieval testing for all transit models in `transit_models.npz`. Tests each transit model by injecting it into `n_injections` randomly selected light curves. The total number of tests will be `n_injections × number_of_models`.

**Parameters:**
- `n_injections` (int, optional): Number of injection-retrieval tests to perform per transit model. If there are fewer files than requested injections, files will be randomly sampled with replacement. If None, uses config value (default: None)
- `output_file` (str, optional): Output filename for results (default: 'inj-ret_results.csv')
- `output_dir` (str, optional): Output directory. If None, defaults to pipeline's output_dir (default: None)
- `plot_completeness` (bool, optional): Whether to generate completeness plot (default: False)
- `completeness_kwargs` (dict, optional): Keyword arguments to pass to `TransitInjector.plot_completeness()`. Common options:
  - `figsize` (tuple): Figure size (default: (5, 4))
  - `cmap` (str): Colormap name (default: 'viridis')
  - `save_fig` (bool): Whether to save figure (default: False)
  - `output_path` (str): Path to save figure (default: 'completeness.png')
  - `dpi` (int): DPI for saved figure (default: 300)

**Returns:** DataFrame with injection-retrieval results containing (one row per injection test):
- TIC ID and sector
- Injected transit parameters (time, depth, duration, SNR)
- Recovery status (boolean)
- Recovered parameters (time, depth, duration, SNR) if detected, NaN otherwise

**Additional Outputs:**
- `{output_file}_stats.csv`: Per-model recovery statistics automatically saved

**Raises:**
- `ValueError`: If transit injector not initialised (requires `transit_models_path` in constructor)

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
- **.npz**: Must have columns for 'time', 'flux', 'flux_err'. Optionally 'phase' and 'eclipse_mask'.
- **.txt**: Must have columns for 'time', 'flux', 'flux_err'. Optionally 'phase' and 'eclipse_mask'.

In both cases, 'phase' and 'eclipse_mask' can be calculated from the input catalogue.

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

**Note:** This method does not modify the original DataFrame. Flags are only available if corresponding analysis was performed (e.g., `det_dependence_flag` requires `cb` detrending method, `skye_flag` requires sector times).

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
- `best_fit` (str): Classification category (based on best AIC model and RMSE):
  - `'T'`: Transit (transit model best fit, AIC diff ≥ 2, RMSE ≤ threshold)
  - `'AT'`: Ambiguous transit (transit model best fit, AIC diff ≥ 2, RMSE > threshold)
  - `'Sin'`: Sinusoid (sinusoid model best fit, AIC diff ≥ 2, RMSE ≤ threshold)
  - `'ASin'`: Ambiguous sinusoid (sinusoid model best fit, AIC diff ≥ 2, RMSE > threshold)
  - `'L'`: Linear (linear model best fit, AIC diff ≥ 2, RMSE ≤ threshold)
  - `'AL'`: Ambiguous linear (linear model best fit, AIC diff ≥ 2, RMSE > threshold)
  - `'St'`: Step (step model best fit, AIC diff ≥ 2, RMSE ≤ threshold)
  - `'ASt'`: Ambiguous step (step model best fit, AIC diff ≥ 2, RMSE > threshold)
  - `'A'`: Ambiguous (AIC difference < 2 from best model to all others)
- `aic_transit` (float): AIC for transit model
- `aic_sinusoidal` (float): AIC for sinusoidal model
- `aic_linear` (float): AIC for linear model
- `aic_step` (float): AIC for step function model
- `rmse_transit` (float): RMSE for transit model
- `rmse_sinusoidal` (float): RMSE for sinusoidal model
- `rmse_linear` (float): RMSE for linear model
- `rmse_step` (float): RMSE for step function model

**Model Fitting:**
This method fits four models to the event data:
1. **Transit model**: Exoplanet transit using [exoplanet](https://docs.exoplanet.codes/en/latest/)
2. **Sinusoidal model**: Sinusoidal variation
3. **Linear model**: Linear trend
4. **Step function model**: Step function. Performs a 2nd order polynomial fit across the largest flux jump, if the jump is > 3 sigma outlier

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
- `best_fit` (str): Classification (T, AT, Sin, ASin, L, AL, St, ASt, A)
- `aic_transit` (float): AIC for transit model
- `aic_sinusoidal` (float): AIC for sinusoidal model
- `aic_linear` (float): AIC for linear model
- `aic_step` (float): AIC for step function model
- `rmse_transit` (float): RMSE for transit model
- `rmse_sinusoidal` (float): RMSE for sinusoidal model
- `rmse_linear` (float): RMSE for linear model
- `rmse_step` (float): RMSE for step function model

**Behaviour:**
- Processes events in batches with progress logging (every 10 events)
- Failed events are logged and skipped
- Saves results to CSV if `output_dir` is specified
- Logs classification summary statistics at completion

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
    n_injections: int = None,
    output_file: str = 'inj-ret_results.csv',
    output_dir: str = None
) -> pd.DataFrame
```

Run injection-retrieval tests for all transit models in `transit_models.npz`. Tests each transit model by injecting it into `n_injections` randomly selected light curves. The total number of tests will be `n_injections × number_of_models`.

**Parameters:**
- `data_dir` (str): Directory containing light curve files
- `n_injections` (int, optional): Number of injection-retrieval tests to perform per transit model. If there are fewer files than requested injections, files will be randomly sampled with replacement (default: None)
- `output_file` (str, optional): Output filename (default: 'inj-ret_results.csv')
- `output_dir` (str, optional): Output directory. If None, defaults to data_dir (default: None)

**Returns:** DataFrame with injection-retrieval results containing (one row per injection test):
- TIC ID and sector
- Injected transit parameters (time, depth, duration, SNR)
- Recovery status (boolean)
- Recovered parameters (time, depth, duration, SNR) if detected, NaN otherwise

**Behaviour:**
- Light curves are inverted to ensure detected events are injected transits, not real signals
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

##### `plot_completeness()`
```python
fig, ax = injector.plot_completeness(
    stats_file: str = None,
    figsize: tuple = (5, 4),
    cmap: str = 'viridis',
    save_fig: bool = False,
    output_path: str = 'completeness.png',
    dpi: int = 300,
    font_family: str = 'sans-serif',
    font_size: int = 8
) -> tuple
```

Plot completeness matrix showing recovery rate as a function of transit depth and duration.

**Parameters:**
- `stats_file` (str, optional): Path to existing stats CSV file. If None, uses `self.stats` (default: None)
- `figsize` (tuple, optional): Figure size (width, height) in inches (default: (5, 4))
- `cmap` (str, optional): Colormap name for heatmap (default: 'viridis')
- `save_fig` (bool, optional): Whether to save the figure (default: False)
- `output_path` (str, optional): Path to save figure if `save_fig=True` (default: 'completeness.png')
- `dpi` (int, optional): DPI for saved figure (default: 300)
- `font_family` (str, optional): Font family for plot text (default: 'sans-serif')
- `font_size` (int, optional): Font size in points (default: 8)

**Returns:** Tuple of (fig, ax) matplotlib figure and axes objects, or (None, None) if no statistics available

**Note:** This method should be called after `run_injection_retrieval()` or with a valid `stats_file` path.

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

Bin short-cadence TESS data to 30-minute long cadence using [lightkurve](https://lightkurve.github.io/lightkurve/).

**Parameters:**
- `time` (np.ndarray): Time values to bin
- `flux` (np.ndarray): Flux values to bin
- `flux_err` (np.ndarray): Flux error values to bin

**Returns:** Tuple of (binned_time, binned_flux, binned_flux_err)

### `mono_cbp.utils.eclipses`

Eclipse calculation and masking utilities.

#### Functions

##### `time_to_phase()`
```python
from mono_cbp.utils.eclipses import time_to_phase

phase = time_to_phase(
    time: np.ndarray,
    period: float,
    t0: float,
    centre: float = 0.5
) -> np.ndarray
```

Convert time measurements to orbital phase of the eclipsing binary.

**Parameters:**
- `time` (np.ndarray): Time values (BJD)
- `period` (float): Orbital period in days
- `bjd0` (float): Reference epoch (BJD)
- `centre` (float, optional): Centre of the phase fold (default: 0.5)

**Returns:** Orbital phase values (0-1)

**Note:** When `centre=0.5`, phase 0.0 corresponds to the reference epoch. When `centre=0.0`, phase 0.5 corresponds to the reference epoch.

##### `get_eclipse_mask()`
```python
from mono_cbp.utils.eclipses import get_eclipse_mask

mask = get_eclipse_mask(
    phases: np.ndarray,
    pos: float,
    width: float
) -> np.ndarray
```

Create boolean mask for eclipse events.

**Parameters:**
- `phase` (np.ndarray): Orbital phase values
- `pos` (float): Mid-clipse position (0-1)
- `width` (float): Eclipse width in phase units

**Returns:** Boolean array (True = in-eclipse)

##### `get_eclipse_indices()`
```python
from mono_cbp.utils.eclipses import get_eclipse_indices

indices = get_eclipse_indices(
    phase: np.ndarray,
    pos: float,
    width: float
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

flat_lcs, trend_lcs, bi_win_lens, cos_success = detrend(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
    method: str,
    fname: str,
    cos_win_len_max: int = 12,
    cos_win_len_min: int = 1,
    fap_threshold: float = 1e-2,
    poly_order: int = 2,
    mask: list = [],
    edge_cutoff: int = 0,
    max_splines: int = 25,
    bi_win_len_max: int = 3,
    bi_win_len_min: int = 1
) -> tuple
```

Detrend a TESS EB light curve using different methods.

**Method Options:**
- `'cb'`: Iterative cosine + multi-biweight detrending
- `'cp'`: Iterative cosine + penalised spline detrending

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `flux_err` (np.ndarray): Flux error values
- `method` (str): Detrending method to use ('cb' or 'cp')
- `fname` (str): Filename to print to user output
- `cos_win_len_max` (int, optional): Maximum window length for cosine detrending in days (default: 12)
- `cos_win_len_min` (int, optional): Minimum window length for cosine detrending in days (default: 1)
- `fap_threshold` (float, optional): False alarm probability threshold (default: 1e-2)
- `poly_order` (int, optional): Polynomial order for trend fitting (default: 2)
- `mask` (list, optional): Boolean mask for data points to exclude from fitting (default: [])
- `edge_cutoff` (int, optional): Amount of data at edges to exclude in days (default: 0)
- `max_splines` (int, optional): Maximum number of splines for penalised spline fitting (default: 25)
- `bi_win_len_max` (int, optional): Maximum window length for biweight detrending in days (default: 3)
- `bi_win_len_min` (int, optional): Minimum window length for biweight detrending in days (default: 1)

**Returns:** Tuple of (detrended_flux, trend_flux, biweight_window_lengths, cosine_success_count)

##### `cosine_detrend()`
```python
from mono_cbp.utils.detrending import cosine_detrend

result = cosine_detrend(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
    win_len_max: int = 12,
    win_len_min: int = 1,
    threshold: float = 1e-2,
    poly_order: int = 2,
    mask: np.ndarray = None,
    edge_cutoff: int = 0
) -> tuple or np.ndarray
```

Performs iterative cosine detrending on the input light curve.

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `flux_err` (np.ndarray): Flux error values
- `win_len_max` (int, optional): Maximum window length for cosine fitting in days (default: 12)
- `win_len_min` (int, optional): Minimum window length for cosine fitting in days (default: 1)
- `threshold` (float, optional): Threshold for false alarm probability (default: 1e-2)
- `poly_order` (int, optional): Order of polynomial for initial detrending (default: 2)
- `mask` (np.ndarray, optional): Boolean mask for data points to exclude from fitting (default: None)
- `edge_cutoff` (int, optional): Amount of data at edges to exclude in time units (default: 0)

**Returns:** Tuple of (detrended_flux, fitted_trend, window_length) if successful, just the input flux if unsuccessful

**Raises:**
- `ValueError`: If maximum window length is smaller than minimum window length

##### `run_multi_biweight()`
```python
from mono_cbp.utils.detrending import run_multi_biweight

biweight_lcs, biweight_trends, win_len_grid = run_multi_biweight(
    time: np.ndarray,
    flux: np.ndarray,
    max_win_len: int = 3,
    min_win_len: int = 1,
    edge_cutoff: int = 0
) -> tuple
```

Run biweight detrending over a range of window lengths.

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `max_win_len` (int, optional): Maximum window length for biweight detrending in days (default: 3)
- `min_win_len` (int, optional): Minimum window length for biweight detrending in days (default: 1)
- `edge_cutoff` (int, optional): Amount of data at edges to exclude in days (default: 0)

**Returns:** Tuple of (detrended_light_curves, fitted_trends, window_length_grid)

**Raises:**
- `ValueError`: If window length inputs are invalid

##### `slider_detrend()`
```python
from mono_cbp.utils.detrending import slider_detrend

flatten_lc, trend_lc = slider_detrend(
    time: np.ndarray,
    flux: np.ndarray,
    win_len: float,
    mask: np.ndarray = None,
    edge_cutoff: int = 0
) -> tuple
```

Apply sliding window detrending to the light curve (biweight).

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `win_len` (float): Length of the sliding window in days
- `mask` (np.ndarray, optional): Boolean mask for data points to exclude from fitting (default: None)
- `edge_cutoff` (int, optional): Amount of data at edges to exclude in time units (default: 0)

**Returns:** Tuple of (detrended_light_curve, fitted_trend)

##### `poly_normalise()`
```python
from mono_cbp.utils.detrending import poly_normalise

normalised_flux = poly_normalise(
    time: np.ndarray,
    flux: np.ndarray,
    order: int = 2
) -> np.ndarray
```

Normalises the flux using a polynomial fit.

Used for reducing the half-sector periodicity from TESS data before calculating the [Lomb-Scargle periodogram](https://docs.astropy.org/en/stable/timeseries/lombscargle.html).

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `order` (int, optional): Order of the polynomial (default: 2)

**Returns:** Normalised flux array

##### `get_period_max_power()`
```python
from mono_cbp.utils.detrending import get_period_max_power

period = get_period_max_power(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray
) -> float
```

Calculates the period with maximum power from the [Lomb-Scargle periodogram](https://docs.astropy.org/en/stable/timeseries/lombscargle.html).

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `flux_err` (np.ndarray): Flux error values

**Returns:** Period with maximum power

##### `get_fap()`
```python
from mono_cbp.utils.detrending import get_fap

fap = get_fap(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray
) -> float
```

Calculates the false alarm probability of the peak of the [Lomb-Scargle periodogram](https://docs.astropy.org/en/stable/timeseries/lombscargle.html).

**Parameters:**
- `time` (np.ndarray): Time values
- `flux` (np.ndarray): Flux values
- `flux_err` (np.ndarray): Flux error values

**Returns:** False alarm probability at the peak of the periodogram

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

Find individual threshold crossing events (TCEs) in timeseries data based on the Median Absolute Deviation (MAD) of the light curve.

**Parameters:**
- `time` (np.ndarray): Time values of the light curve
- `flux` (np.ndarray): Flux values of the light curve
- `mad` (float, optional): MAD multiplier for event detection (default: 3.0)
- `var_mad` (np.ndarray, optional): MAD values calculated using a rolling window for adaptive thresholding. If None, uses constant MAD (default: None)

**Returns:** Tuple of (event_indices, metadata_dict) where:
- `event_indices` (np.ndarray): Indices of detected event peaks
- `metadata_dict` (dict): Dictionary containing:
  - `'threshold'`: Detection threshold value
  - `'depths'`: Event depths (flux dip magnitudes)
  - `'widths'`: Event widths (durations in days)
  - `'start_times'`: Event start times
  - `'end_times'`: Event end times

**Algorithm:**
Identifies dips in the light curve that exceed the MAD-based threshold, groups nearby points into events, and calculates event properties.

##### `get_var_mad()`
```python
from mono_cbp.utils.monofind import get_var_mad

var_mad = get_var_mad(
    flux: np.ndarray,
    npoints: int
) -> np.ndarray
```

Calculates the MAD of a light curve over a running window.

**Parameters:**
- `flux` (np.ndarray): Flux values of a light curve
- `npoints` (int): Number of data points in a given window

**Returns:** MAD of the light curve over a running window

##### `get_gaps_indices()`
```python
from mono_cbp.utils.monofind import get_gaps_indices

gaps_indices = get_gaps_indices(
    time: np.ndarray,
    break_tolerance: float
) -> np.ndarray
```

Array indices where 'time' has gaps longer than 'break_tolerance'.

**Parameters:**
- `time` (np.ndarray): Array of time values
- `break_tolerance` (float): Threshold for the gap distance in days

**Returns:** Indices on the time axis where gaps occur

##### `consecutive()`
```python
from mono_cbp.utils.monofind import consecutive

chunks = consecutive(
    data: np.ndarray,
    stepsize: int = 1
) -> list
```

Split an array into consecutive chunks if gap between values is above a given step size.

**Parameters:**
- `data` (np.ndarray): Data to find gaps in
- `stepsize` (int, optional): The size of the step to cluster consecutive data (default: 1)

**Returns:** List of arrays of clustered data

### `mono_cbp.utils.plotting`

Plotting utilities for light curves and events.

#### Functions

##### `plot_event()`
```python
from mono_cbp.utils.plotting import plot_event

fig = plot_event(
    time: np.ndarray,
    event_time: float,
    flat_flux: np.ndarray,
    raw_flux: np.ndarray,
    flux_err: np.ndarray,
    trend: np.ndarray,
    fname: str,
    mad: float,
    var_mad: np.ndarray,
    depth: float,
    width: float,
    phase: float,
    SNR: float,
    peaks: np.ndarray,
    event_no: int,
    ecl_mask: np.ndarray = None,
    output_dir: str = None,
    mask: list = [],
    figsize: tuple = (20/3, 8),
    save: bool = True,
    return_fig: bool = False
) -> matplotlib.figure.Figure or None
```

Plot the light curve for a detected event.

**Parameters:**
- `time` (np.ndarray): Array of time values
- `event_time` (float): Time of the event in days
- `flat_flux` (np.ndarray): Flattened flux array
- `raw_flux` (np.ndarray): Raw flux array
- `flux_err` (np.ndarray): Flux error array
- `trend` (np.ndarray): Trend array
- `fname` (str): Filename
- `mad` (float): Threshold multiplier of Median Absolute Deviation
- `var_mad` (np.ndarray): Variable Median Absolute Deviation
- `depth` (float): Depth of the event
- `width` (float): Width of the event
- `phase` (float): Binary phase of the event
- `SNR` (float): Signal-to-noise ratio of the event
- `peaks` (np.ndarray): Array of detected peak indices
- `event_no` (int): Event number
- `ecl_mask` (np.ndarray, optional): Eclipse mask array (default: None)
- `output_dir` (str, optional): Output directory (default: None)
- `mask` (list, optional): List of boolean masks for the data (default: [])
- `figsize` (tuple, optional): Figure size (default: (20/3, 8))
- `save` (bool, optional): Whether to save the plot (default: True)
- `return_fig` (bool, optional): Whether to return the figure object (default: False)

**Returns:** matplotlib Figure object if return_fig is True, otherwise None

##### `plot_no_events()`
```python
from mono_cbp.utils.plotting import plot_no_events

fig = plot_no_events(
    time: np.ndarray,
    flat_flux: np.ndarray,
    raw_flux: np.ndarray,
    flux_err: np.ndarray,
    trend: np.ndarray,
    fname: str,
    mad: float,
    var_mad: np.ndarray,
    ecl_mask: np.ndarray = None,
    output_dir: str = None,
    mask: list = [],
    figsize: tuple = (20/3, 8),
    save: bool = True,
    return_fig: bool = False
) -> matplotlib.figure.Figure or None
```

Plot light curves with no events detected.

**Parameters:**
- `time` (np.ndarray): Time array
- `flat_flux` (np.ndarray): Flattened flux array
- `raw_flux` (np.ndarray): Raw flux array
- `flux_err` (np.ndarray): Flux error array
- `trend` (np.ndarray): Trend array
- `fname` (str): Filename
- `mad` (float): Threshold multiplier of Median Absolute Deviation
- `var_mad` (np.ndarray): Variable Median Absolute Deviation
- `ecl_mask` (np.ndarray, optional): Eclipse mask array (default: None)
- `output_dir` (str, optional): Output directory (default: None)
- `mask` (list, optional): Mask for the data (default: [])
- `figsize` (tuple, optional): Figure size (default: (20/3, 8))
- `save` (bool, optional): Whether to save the plot to disk (default: True)
- `return_fig` (bool, optional): Whether to return the figure object (default: False)

**Returns:** matplotlib Figure object if return_fig is True, otherwise None

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

Set up logging for mono-cbp package.

**Parameters:**
- `level` (int, optional): Logging level (e.g., logging.INFO, logging.DEBUG) (default: logging.INFO)
- `log_file` (str, optional): Path to log file. If None, logs to console only (default: None)

**Returns:** Configured logger instance

**Note:** Clears existing handlers and sets up both console and file handlers (if log_file is specified).

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

Calculate signal-to-noise ratio for a threshold crossing event (TCE).

**Parameters:**
- `depth` (float): TCE depth (fractional flux decrease)
- `error` (float): Combined measurement uncertainty
- `duration` (float): TCE duration in days
- `cadence` (int, optional): Observation cadence in minutes (default: 30)

**Returns:** Signal-to-noise ratio value

**Formula:**
SNR = (depth / error) × sqrt(duration / cadence_days) where cadence_days = cadence_minutes / 1440

### `mono_cbp.utils.transit_models`

Transit model creation and handling for injection-retrieval testing.

#### Functions

##### `create_transit_models()`
```python
from mono_cbp.utils.transit_models import create_transit_models

models_dict = create_transit_models(
    depth_range: tuple = (1e-3, 1e-2),
    duration_range: tuple = (0.1, 1.0),
    num_depths: int = 7,
    num_durations: int = 7,
    time_range: tuple = (-1, 1),
    cadence_minutes: float = 30,
    impact_parameter: float = 0.0,
    period: float = 10.0,
    limb_dark_coeffs: tuple = (0.3, 0.2)
) -> dict
```

Create a grid of synthetic transit light curves for injection-retrieval testing.

This function generates transit models across a grid of transit depths and durations, which can be used to characterize pipeline completeness and detection efficiency.

**Parameters:**
- `depth_range` (tuple, optional): (min, max) transit depth in fractional flux units (default: (1e-3, 1e-2) = 0.1% to 1%)
- `duration_range` (tuple, optional): (min, max) transit duration in days (default: (0.1, 1.0))
- `num_depths` (int, optional): Number of depth values to sample (default: 7)
- `num_durations` (int, optional): Number of duration values to sample (default: 7)
- `time_range` (tuple, optional): (start, end) time range in days centered on transit (default: (-1, 1))
- `cadence_minutes` (float, optional): Observation cadence in minutes (default: 30)
- `impact_parameter` (float, optional): Impact parameter (0 = center of limb) (default: 0.0)
- `period` (float, optional): Orbital period in days (arbitrary, >2×duration) (default: 10.0)
- `limb_dark_coeffs` (tuple, optional): Quadratic limb darkening coefficients (u1, u2) (default: (0.3, 0.2))

**Returns:** Dictionary containing:
- `'time'`: Time array (same for all models)
- `'models'`: List of transit model dictionaries, each containing:
  - `'flux'`: Normalized flux array
  - `'depth'`: Transit depth
  - `'duration'`: Transit duration in days
  - `'impact_parameter'`: Impact parameter
  - `'ror'`: Radius ratio (planet radius / star radius)
- `'num_depths'`: Number of depth values
- `'num_durations'`: Number of duration values
- `'depth_range'`: Depth range tuple
- `'duration_range'`: Duration range tuple
- `'cadence_minutes'`: Cadence in minutes

##### `save_transit_models()`
```python
from mono_cbp.utils.transit_models import save_transit_models

save_transit_models(
    models_dict: dict,
    filepath: str
) -> None
```

Save transit models to an .npz file.

**Parameters:**
- `models_dict` (dict): Dictionary from `create_transit_models()`
- `filepath` (str): Path to save .npz file

##### `load_transit_models()`
```python
from mono_cbp.utils.transit_models import load_transit_models

models_dict = load_transit_models(
    filepath: str
) -> dict
```

Load transit models from an .npz file.

**Parameters:**
- `filepath` (str): Path to .npz file created by `save_transit_models()`

**Returns:** Dictionary with same structure as `create_transit_models()` output

## Configuration

The mono-cbp package uses a hierarchical configuration system to control pipeline behaviour. Configuration can be provided when initialising pipeline components or by using the configuration utilities.

### Configuration Functions

#### `get_default_config()`

```python
from mono_cbp.config import get_default_config

config = get_default_config()
```

Get a copy of the default configuration.

**Returns:** Deep copy of the default configuration dictionary

**Note:** Returns a deep copy to prevent accidental modification of the default configuration.

#### `merge_config()`

```python
from mono_cbp.config import merge_config

merged = merge_config(
    user_config: dict,
    default_config: dict = None
) -> dict
```

Merge user configuration with default configuration.

**Parameters:**
- `user_config` (dict): User-provided configuration
- `default_config` (dict, optional): Base configuration. Uses `DEFAULT_CONFIG` if None (default: None)

**Returns:** Merged configuration dictionary

**Behaviour:**
- Recursively merges nested dictionaries
- User values override default values
- Creates a deep copy to prevent mutation of inputs

**Example:**
```python
from mono_cbp.config import get_default_config, merge_config

# Start with defaults
config = get_default_config()

# Override specific values
user_config = {
    'transit_finding': {
        'mad_threshold': 4.0,
        'cosine': {
            'win_len_max': 15
        }
    }
}

# Merge configurations
final_config = merge_config(user_config, config)
```

### Configuration Structure

The default configuration dictionary with all available options:

```python
config = {
    'transit_finding': {
        'edge_cutoff': 0.0,                     # Edge cutoff in days
        'mad_threshold': 3.0,                   # MAD threshold multiplier
        'detrending_method': 'cb',              # 'cb' (cosine+biweight) or 'cp' (cosine+pspline)
        'generate_vetting_plots': False,        # Generate diagnostic plots for events
        'generate_skye_plots': False,           # Generate Skye metric histograms
        'generate_event_snippets': True,        # Generate event data snippets
        'save_event_snippets': True,            # Save event snippets to disk
        'cadence_minutes': 30,                  # Cadence in minutes
        'cosine': {
            'win_len_max': 12,                  # Maximum window length (days)
            'win_len_min': 1,                   # Minimum window length (days)
            'fap_threshold': 1e-2,              # False alarm probability threshold
            'poly_order': 2,                    # Polynomial order
        },
        'biweight': {
            'win_len_max': 3,                   # Maximum window length (days)
            'win_len_min': 1,                   # Minimum window length (days)
        },
        'pspline': {
            'max_splines': 25,                  # Maximum number of splines
        },
        'filters': {
            'min_snr': 5,                       # Minimum SNR for filtering
            'max_duration_days': 1,             # Maximum duration in days
            'det_dependence_threshold': 18,     # Detrending dependence threshold
        }
    },
    'model_comparison': {
        'tune': 1000,                           # MCMC tuning steps
        'draws': 1000,                          # MCMC draws per chain
        'chains': 4,                            # Number of MCMC chains
        'cores': 4,                             # Number of CPU cores
        'target_accept': 0.99,                  # MCMC target acceptance rate
        'sigma_threshold': 3,                   # Sigma threshold for outlier removal
        'aic_threshold': 2,                     # AIC difference threshold
        'rmse_threshold': 1.2,                  # RMSE threshold for classification
        'save_plots': False,                    # Save model comparison plots
        'plot_dir': None,                       # Directory for plots (None = use output_dir)
    },
    'injection_retrieval': {
        'n_injections': 1000,                   # Number of injections per model
    }
}
```

### Configuration Sections

#### `transit_finding`

Controls the transit finding process including detrending and event detection.

**Top-level parameters:**
- `edge_cutoff` (float): Amount of data at light curve edges to exclude (in days)
- `mad_threshold` (float): MAD multiplier for event detection threshold
- `detrending_method` (str): Detrending method - `'cb'` (cosine+biweight) or `'cp'` (cosine+pspline)
- `generate_vetting_plots` (bool): Generate diagnostic plots showing detected events
- `generate_skye_plots` (bool): Generate histograms of Skye metric distributions
- `generate_event_snippets` (bool): Extract event data snippets for vetting
- `save_event_snippets` (bool): Save event snippets to disk (requires `generate_event_snippets=True`)
- `cadence_minutes` (int): Observation cadence in minutes

**`cosine` subsection:**
Parameters for iterative cosine detrending step.
- `win_len_max` (float): Maximum window length in days
- `win_len_min` (float): Minimum window length in days
- `fap_threshold` (float): False alarm probability threshold for periodogram
- `poly_order` (int): Polynomial order for initial detrending

**`biweight` subsection:**
Parameters for biweight detrending (used with `detrending_method='cb'`).
- `win_len_max` (float): Maximum window length in days
- `win_len_min` (float): Minimum window length in days

**`pspline` subsection:**
Parameters for penalised spline detrending (used with `detrending_method='cp'`).
- `max_splines` (int): Maximum number of splines for fitting

**`filters` subsection:**
Quality filter thresholds for detected events (applied via `TransitFinder.filter_events()`).
- `min_snr` (float): Minimum signal-to-noise ratio
- `max_duration_days` (float): Maximum transit duration in days
- `det_dependence_threshold` (int): Detrending dependence threshold (number of window lengths)

#### `model_comparison`

Controls Bayesian model comparison for candidate vetting.

- `tune` (int): Number of MCMC tuning steps
- `draws` (int): Number of MCMC draws per chain
- `chains` (int): Number of MCMC chains
- `cores` (int): Number of CPU cores for parallel sampling
- `target_accept` (float): MCMC target acceptance rate
- `sigma_threshold` (float): Sigma threshold for outlier removal before fitting
- `aic_threshold` (float): AIC difference threshold for model selection (models differing by less are ambiguous)
- `rmse_threshold` (float): RMSE threshold for classification quality
- `save_plots` (bool): Save model comparison diagnostic plots
- `plot_dir` (str or None): Directory for plots (None uses output_dir)

#### `injection_retrieval`

Controls injection-retrieval testing.

- `n_injections` (int): Number of injection tests per transit model

### Using Configuration

Configuration can be passed when initialising pipeline components:

```python
from mono_cbp import MonoCBPPipeline
from mono_cbp.config import get_default_config, merge_config

# Get default config
config = get_default_config()

# Customize specific values
user_config = {
    'transit_finding': {
        'mad_threshold': 4.0,
    }
}
config = merge_config(user_config, config)

# Initialise pipeline with custom config
pipeline = MonoCBPPipeline(
    catalogue_path='catalogue.csv',
    data_dir='./data',
    config=config
)
```