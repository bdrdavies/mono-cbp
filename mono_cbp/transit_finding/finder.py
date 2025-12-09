"""Finding additional transit events in eclipsing binary light curves."""

import os
import numpy as np
import pandas as pd
import logging
import warnings
import matplotlib.pyplot as plt
from ..utils import (
    bin_to_long_cadence, get_var_mad, monofind, split_tol,
    time_to_phase, get_eclipse_mask, get_snr
)
from ..utils.detrending import detrend
from ..utils.plotting import plot_no_events, plot_event
from ..config import get_default_config, merge_config
from ..utils.data import process_tebc_catalogue
from ..utils import load_catalogue

logger = logging.getLogger('mono_cbp.transit_finding')
warnings.simplefilter("ignore")

# Constants
MINUTES_TO_DAYS = 1 / 1440
VAR_MAD_WINDOW = 100           # data points
PROGRESS_INTERVAL = 10
EVENT_WINDOW_HALF_WIDTH = 0.5  # days (used for event windows when duration < 1 day)
EVENT_GROUPING_TOLERANCE = 10  # data points
SKYE_FIGURE_SIZE = (12, 10)

class TransitFinder:
    """Class for finding additional transit events in eclipsing binary light curves.

    This class handles the Threshold Crossing Event (TCE) detection pipeline including:
    - Light curve detrending
    - Event detection with monofind
    - SNR estimation and filtering
    - Skye metric calculation for systematic artifact identification
    - Generation of event snippets for model comparison

    Attributes:
        catalogue (pd.DataFrame or None): Catalogue with binary ephemerides and eclipse parameters
        sector_times (pd.DataFrame or None): Sector times for Skye metric calculation
        config (dict): Full configuration dictionary
        transit_config (dict): Transit finding specific configuration parameters
        results (dict): Detected events and metadata including event times, depths, durations, etc.
        stats (dict): Statistics on processing including total files, events, and cosine successes
    """

    def __init__(self, catalogue=None, sector_times=None, config=None, TEBC=False):
        """Initialise TransitFinder.

        Args:
            catalogue (str or pd.DataFrame, optional): Path to or DataFrame of catalogue with
                binary ephemerides and eclipse parameters. Must contain columns: tess_id, period, bjd0,
                prim_pos, prim_width, sec_pos, sec_width (or *_pf/*_2g if TEBC=True for eclipse parameters).
            sector_times (str or pd.DataFrame, optional): Path to or DataFrame of sector times for Skye metric
            config (dict, optional): Configuration dictionary. Uses default_config if None.
            TEBC (bool, optional): If True, processes TEBC catalogue format with *_2g and *_pf columns
                and converts to standard eclipse parameter columns. If a DataFrame is passed that already
                has standard columns, TEBC processing is skipped. Defaults to False.
        """
        # Load catalogue
        if catalogue is not None:
            if isinstance(catalogue, str):
                self.catalogue = load_catalogue(catalogue, TEBC=TEBC)
                logger.info(f"Loaded catalogue from {catalogue}")
            else:
                # If DataFrame is passed directly and has TEBC columns, process it only if not already processed
                if TEBC:
                    # Check if catalogue has already been processed (has standard eclipse parameter columns)
                    required_cols = {'prim_pos', 'prim_width', 'sec_pos', 'sec_width'}
                    if not required_cols.issubset(catalogue.columns):
                        catalogue = process_tebc_catalogue(catalogue)
                    else:
                        logger.info("Catalogue already has standard eclipse parameter columns, skipping TEBC processing")
                self.catalogue = catalogue
        else:
            self.catalogue = None
            logger.warning("No catalogue provided")

        if sector_times is not None:
            if isinstance(sector_times, str):
                self.sector_times = pd.read_csv(sector_times, comment='#')
                logger.info(f"Loaded sector times from {sector_times}")
            else:
                self.sector_times = sector_times
        else:
            self.sector_times = None
            logger.warning("No sector times provided - Skye metric will not be available")

        # Configuration
        self.config = merge_config(config, get_default_config()) if config else get_default_config()
        self.transit_config = self.config['transit_finding']

        # Results storage
        self.results = self._initialise_results()
        self.stats = self._initialise_stats()

    @staticmethod
    def _initialise_results():
        """Initialise empty results dictionary for storing detected transit events.

        Returns:
            dict: Empty results dictionary with keys:
                - tics/sectors: TIC ID and sector of each detected event
                - event_times/event_phases: Time and orbital phase of detected events
                - event_depths/event_durations: Transit depth (fraction) and duration (days)
                - event_snrs: Signal-to-noise ratio for each event
                - win_len_max_SNR: Biweight window length with maximum SNR (cb method only)
                - det_dependence: Flag indicating detrending method dependence (0=robust, 1=dependent)
                - start_times/end_times: Event window boundaries
                - skye_flags: Systematic artifact flags from Skye metric
                - event_snippets: Event data windows for model comparison
        """
        return {
            'tics': [],
            'sectors': [],
            'event_times': [],
            'event_phases': [],
            'event_depths': [],
            'event_durations': [],
            'event_snrs': [],
            'win_len_max_SNR': [],
            'det_dependence': [],
            'start_times': [],
            'end_times': [],
            'skye_flags': [],
            'event_snippets': []
        }

    @staticmethod
    def _initialise_stats():
        """Initialise empty statistics dictionary for tracking processing metrics.

        Returns:
            dict: Empty statistics dictionary that is updated throughout processing with keys:
                - total_files: Number of light curve files processed
                - total_events: Total number of transit events detected across all files
                - cosine_successes: Number of files where cosine detrending succeeded
        """
        return {
            'total_files': 0,
            'total_events': 0,
            'cosine_successes': 0
        }

    def _clear_results(self):
        """Clear all stored results and statistics.

        Called automatically at the start of process_directory to prevent accumulation.
        """
        self.results = self._initialise_results()
        self.stats = self._initialise_stats()

    def process_directory(self, data_dir, output_file='output.txt', output_dir=None, plot_output_dir=None):
        """Process all light curve files in a directory.

        Args:
            data_dir (str): Directory containing light curve files (.npz or .txt)
            output_file (str, optional): Name of output file. Defaults to 'output.txt'.
            output_dir (str, optional): Directory to save output file. Defaults to current working directory.
            plot_output_dir (str, optional): Directory to save plots. Defaults to None.

        Returns:
            pd.DataFrame: DataFrame of detected events. Returns empty DataFrame if no events
            are detected.
        """
        # Clear previous results
        self._clear_results()

        logger.info(f"Processing files in {data_dir}")

        # Get list of files
        files = [f for f in os.listdir(data_dir) if f.endswith('.npz') or f.endswith('.txt')]
        self.stats['total_files'] = len(files)
        logger.info(f"Found {len(files)} files to process")

        # Process each file
        for i, file in enumerate(files):
            if i % PROGRESS_INTERVAL == 0:
                logger.info(f"Progress: {i}/{len(files)}")

            file_path = os.path.join(data_dir, file)
            self.process_file(file_path, plot_output_dir=plot_output_dir)

        # Calculate Skye metric
        if self.sector_times is not None:
            self._calculate_skye_metric(plot_output_dir=plot_output_dir)

        # Save results
        results_df = self._save_results(output_file, output_dir)

        # Log statistics
        logger.info(f"Processing complete:")
        logger.info(f"  Total files: {self.stats['total_files']}")
        logger.info(f"  Total events: {self.stats['total_events']}")
        logger.info(f"  Cosine detrending successes: {self.stats['cosine_successes']}/{self.stats['total_files']}")

        return results_df

    def _load_npz(self, file_path):
        """Load light curve data from npz file.

        Args:
            file_path (str): Path to npz file

        Returns:
            tuple: (time, flux, flux_err, phase, ecl_mask) where:
                - time: np.ndarray of time values
                - flux: np.ndarray of flux values
                - flux_err: np.ndarray of flux error values
                - phase: np.ndarray of binary phase values or None if not present in file
                - ecl_mask: np.ndarray (bool) of eclipse mask or None if not present in file
        """
        data = np.load(file_path)
        time = data['time']
        flux = data['flux']
        flux_err = data['flux_err']
        phase = data.get('phase', None)
        ecl_mask_raw = data.get('eclipse_mask', None)
        ecl_mask = ecl_mask_raw.astype(bool) if ecl_mask_raw is not None else None
        return time, flux, flux_err, phase, ecl_mask

    def _load_txt(self, file_path):
        """Load light curve data from txt file.

        Args:
            file_path (str): Path to txt file

        Returns:
            tuple: (time, flux, flux_err, phase, ecl_mask) where:
                - time: np.ndarray of time values
                - flux: np.ndarray of flux values
                - flux_err: np.ndarray of flux error values
                - phase: np.ndarray of binary phase values
                - ecl_mask: np.ndarray (bool) of eclipse mask or None if not present in file
        """
        data = np.loadtxt(file_path, skiprows=1)
        time = data[:, 0]
        flux = data[:, 1]
        flux_err = data[:, 2]
        phase = data[:, 3]
        ecl_mask = data[:, 4].astype(bool) if data.shape[1] > 4 else None
        return time, flux, flux_err, phase, ecl_mask

    def process_file(self, file_path, plot_output_dir=None):
        """Process a single light curve file.

        Args:
            file_path (str): Path to light curve file
            plot_output_dir (str, optional): Directory to save plots. Defaults to None.

        Returns:
            list: List of detected event dictionaries. Each event dictionary contains keys:
                'time', 'phase', 'depth', 'width', 'duration', 'start_time', 'end_time', 'snr',
                'tic', 'sector'. Returns empty list if file cannot be processed or no events are detected.
        """
        # Parse filename
        file = os.path.basename(file_path)
        split_file = os.path.splitext(file)

        # Load file data
        if split_file[1] == '.npz':
            time, flux, flux_err, phase, ecl_mask = self._load_npz(file_path)
        elif split_file[1] == '.txt':
            time, flux, flux_err, phase, ecl_mask = self._load_txt(file_path)
        else:
            logger.warning(f"Unsupported file format: {file}")
            return []

        # Extract TIC and sector
        tic, sector = self._parse_filename(file)
        if tic is None:
            logger.warning(f"Could not parse TIC/sector from {file}, skipping")
            return []

        # Get eclipse parameters
        prim_pos, prim_width, sec_pos, sec_width = self._get_eclipse_params(tic)

        # Bin to long cadence if needed
        if np.median(np.gradient(time[~np.isnan(flux)])) < self.transit_config['cadence_minutes'] * MINUTES_TO_DAYS:
            time, flux, flux_err = bin_to_long_cadence(time, flux, flux_err)
            if self.catalogue is not None:
                row = self.catalogue[self.catalogue['tess_id'] == tic]
                if not row.empty:
                    phase = time_to_phase(time, row['period'].values[0], row['bjd0'].values[0])
                    if prim_pos is not None:
                        eclipse_mask_p = get_eclipse_mask(phase, prim_pos, prim_width)
                        eclipse_mask_s = get_eclipse_mask(phase, sec_pos, sec_width)
                        ecl_mask = np.logical_or(eclipse_mask_p, eclipse_mask_s)

        # Create mask
        nan_mask = ~np.isnan(flux * time * flux_err)
        if ecl_mask is not None:
            mask = nan_mask & ~ecl_mask
        else:
            mask = nan_mask

        # Detrend
        method = self.transit_config['detrending_method']
        edge_cutoff = self.transit_config['edge_cutoff']

        detrend_result = detrend(
            time, flux, flux_err, method=method, fname=split_file[0],
            mask=mask, edge_cutoff=edge_cutoff,
            cos_win_len_max=self.transit_config['cosine']['win_len_max'],
            cos_win_len_min=self.transit_config['cosine']['win_len_min'],
            fap_threshold=self.transit_config['cosine']['fap_threshold'],
            poly_order=self.transit_config['cosine']['poly_order'],
            max_splines=self.transit_config['pspline']['max_splines'],
            bi_win_len_max=self.transit_config['biweight']['win_len_max'],
            bi_win_len_min=self.transit_config['biweight']['win_len_min']
        )

        # Update stats
        self.stats['cosine_successes'] += detrend_result[3]

        # Find events based on method
        if method == 'cp':
            events = self._process_cp_events(
                time, flux, flux_err, phase, ecl_mask, mask,
                detrend_result, tic, sector, split_file[0],
                plot_output_dir
            )
        elif method == 'cb':
            events = self._process_cb_events(
                time, flux, flux_err, phase, ecl_mask, mask,
                detrend_result, tic, sector, split_file[0],
                plot_output_dir
            )
        else:
            logger.error(f"Unknown detrending method: {method}")
            events = []

        self.stats['total_events'] += len(events)
        return events

    def _process_cp_events(self, time, flux, flux_err, phase, ecl_mask, mask,
                          detrend_result, tic, sector, fname, plot_output_dir):
        """Process and extract TCEs using cosine + pspline detrending.

        Args:
            time (np.ndarray): Time array.
            flux (np.ndarray): Flux array.
            flux_err (np.ndarray): Flux error array.
            phase (np.ndarray): Binary phase array or None.
            ecl_mask (np.ndarray): Eclipse mask array or None.
            mask (np.ndarray): Combined mask array indicating valid out-of-eclipse data points.
            detrend_result (tuple): Result from the detrending function.
            tic (int): TIC ID.
            sector (int): Sector number.
            fname (str): Filename for the light curve file.
            plot_output_dir (str): Directory to save plots.

        Returns:
            list: Detected event dictionaries with keys:
                'time', 'phase', 'depth', 'width', 'duration', 'start_time', 'end_time', 'snr', 'tic', 'sector'
        """
        flatten_lc, trend_lc, _, _ = detrend_result
        mad = self.transit_config['mad_threshold']
        var_mad = get_var_mad(flatten_lc, VAR_MAD_WINDOW)

        # Run monofind
        peaks, meta = monofind(time[mask], flatten_lc, mad=mad, var_mad=var_mad)

        if len(peaks) == 0:
            # Plot no events if requested
            if self.transit_config['generate_vetting_plots'] and plot_output_dir:
                plot_no_events(time, flatten_lc, flux, flux_err, trend_lc,
                             fname, mad=mad, var_mad=var_mad, ecl_mask=ecl_mask,
                             mask=mask, output_dir=plot_output_dir)
            return []

        # Process each event
        # NOTE: All events are saved, filtering can be done later
        events = []
        for j in range(len(peaks)):
            event_data = self._extract_event_data(
                time, flux_err, phase, flatten_lc, mask,
                peaks[j], meta, j, tic, sector
            )

            # Save ALL events (no filtering)
            events.append(event_data)
            self.results['tics'].append(str(tic))
            self.results['sectors'].append(str(sector))
            self.results['event_times'].append(event_data['time'])
            self.results['event_phases'].append(event_data['phase'])
            self.results['event_depths'].append(event_data['depth'])
            self.results['event_durations'].append(event_data['duration'])
            self.results['event_snrs'].append(event_data['snr'])

            if self.transit_config['generate_vetting_plots'] and plot_output_dir:
                plot_event(time, event_data['time'], flatten_lc, flux, flux_err,
                         trend_lc, fname, mad, var_mad, event_data['depth'],
                         event_data['width'], event_data['phase'],
                         event_data['snr'], peaks, j+1, ecl_mask=ecl_mask,
                         mask=mask, output_dir=plot_output_dir)

        return events

    def _process_cb_events(self, time, flux, flux_err, phase, ecl_mask, mask,
                          detrend_result, tic, sector, fname, plot_output_dir):
        """Process and extract TCEs using cosine + biweight detrending.

        Args:
            time (np.ndarray): Time array.
            flux (np.ndarray): Flux array.
            flux_err (np.ndarray): Flux error array.
            phase (np.ndarray): Binary phase array or None.
            ecl_mask (np.ndarray): Eclipse mask array or None.
            mask (np.ndarray): Combined mask array indicating valid out-of-eclipse data points.
            detrend_result (tuple): Result from the detrending function.
            tic (int): TIC ID.
            sector (int): Sector number.
            fname (str): Filename for the light curve file.
            plot_output_dir (str): Directory to save plots.

        Returns:
            list: Detected event dictionaries (one per group) with keys:
                'time', 'phase', 'depth', 'width', 'duration', 'start_time', 'end_time', 'snr', 'tic', 'sector',
                plus: 'win_len' (biweight window of max SNR)
        """
        flatten_lcs, trend_lcs, bi_win_lens, _ = detrend_result
        mad = self.transit_config['mad_threshold']

        # Storage for all detected events across all biweight windows
        all_peaks = []
        event_data_all = []

        # Loop over biweight windows
        for index, lc in enumerate(flatten_lcs):
            win_len = bi_win_lens[index]
            var_mad = get_var_mad(lc, VAR_MAD_WINDOW)
            peaks, meta = monofind(time[mask], lc, mad=mad, var_mad=var_mad)
            all_peaks.append(peaks)

            # Extract event data for each peak
            for j in range(len(peaks)):
                event_data = self._extract_event_data(
                    time, flux_err, phase, lc, mask,
                    peaks[j], meta, j, tic, sector
                )
                event_data['win_len'] = win_len
                event_data['flat_lc_idx'] = index
                event_data['var_mad'] = var_mad
                event_data_all.append(event_data)

        # Flatten peaks list
        all_peaks_flat = [p for ps in all_peaks for p in ps]

        if len(all_peaks_flat) == 0:
            # Plot no events if requested
            if self.transit_config['generate_vetting_plots'] and plot_output_dir:
                plot_no_events(time, flatten_lcs[-1], flux, flux_err, trend_lcs[-1],
                             fname, mad=mad, var_mad=get_var_mad(flatten_lcs[-1], VAR_MAD_WINDOW),
                             ecl_mask=ecl_mask, mask=mask, output_dir=plot_output_dir)
            return []

        # Group events detected at similar times across different window lengths
        all_peaks_flat_sorted = np.sort(all_peaks_flat)
        time_sorted_idx = np.argsort([e['time'] for e in event_data_all])
        event_data_sorted = [event_data_all[i] for i in time_sorted_idx]

        events_grouped = split_tol(all_peaks_flat_sorted, EVENT_GROUPING_TOLERANCE)

        # Select highest SNR event from each group
        # NOTE: All events are saved, filtering can be done later with filter_events()
        events = []
        start_idx = 0
        for group in events_grouped:
            group_events = event_data_sorted[start_idx:start_idx+len(group)]
            snrs = [e['snr'] for e in group_events]
            max_snr_idx = np.argmax(snrs)
            best_event = group_events[max_snr_idx]

            # Calculate detrending dependence flag
            det_dep = 1 if len(group) <= self.transit_config['filters']['det_dependence_threshold'] else 0

            # Save ALL events (no filtering)
            events.append(best_event)
            self.results['tics'].append(str(tic))
            self.results['sectors'].append(str(sector))
            self.results['event_times'].append(best_event['time'])
            self.results['event_phases'].append(best_event['phase'])
            self.results['event_depths'].append(best_event['depth'])
            self.results['event_durations'].append(best_event['duration'])
            self.results['event_snrs'].append(best_event['snr'])
            self.results['win_len_max_SNR'].append(round(best_event['win_len'], 1))
            self.results['det_dependence'].append(det_dep)

            # Save event snippet if requested
            if self.transit_config['generate_event_snippets']:
                snippet = self._create_event_snippet(
                    time, flux_err, flatten_lcs[best_event['flat_lc_idx']],
                    mask, best_event, tic, sector, len(events)
                )
                self.results['event_snippets'].append(snippet)

            # Plot ONLY if event passes filters
            if self.transit_config['generate_vetting_plots'] and plot_output_dir:
                flat_lc_idx = best_event['flat_lc_idx']
                peaks = all_peaks[flat_lc_idx]
                plot_event(time, best_event['time'], flatten_lcs[flat_lc_idx],
                         flux, flux_err, trend_lcs[flat_lc_idx], fname, mad,
                         best_event['var_mad'], best_event['depth'],
                         best_event['duration'], best_event['phase'],
                         best_event['snr'], peaks, len(events),
                         ecl_mask=ecl_mask, mask=mask, output_dir=plot_output_dir)

            start_idx += len(group)

        return events

    def _extract_event_data(self, time, flux_err, phase, flatten_lc, mask,
                           peak_idx, meta, event_idx, tic, sector):
        """Extract event metadata and calculate SNR for a single detected event.

        Args:
            time (np.ndarray): Time array.
            flux_err (np.ndarray): Flux error array.
            phase (np.ndarray): Binary phase array or None.
            flatten_lc (np.ndarray): Detrended/flattened light curve.
            mask (np.ndarray): Boolean mask of valid points.
            peak_idx (int): Index of detected events in masked arrays.
            meta (dict): Metadata dictionary from monofind with keys:
                - start_times: Event start times.
                - end_times: Event end times.
                - widths: Event widths/durations.
                - depths: Transit depths (fraction).
            event_idx (int): Index into meta arrays.
            tic (int): TIC ID.
            sector (int): Sector number.

        Returns:
            dict: Event properties dictionary with keys:
                'time', 'phase', 'depth', 'width', 'duration', 'start_time', 'end_time', 'snr', 'tic', 'sector'
        """
        event_time = time[mask][peak_idx]
        event_phase = phase[mask][peak_idx] if phase is not None else None

        # Calculate SNR
        in_transit = (time[mask] >= meta['start_times'][event_idx]) & \
                     (time[mask] <= meta['end_times'][event_idx])
        event_flux_err = np.median(flux_err[mask][in_transit])

        # Calculate scatter
        event_duration = meta['widths'][event_idx]
        if (meta['end_times'][event_idx] + event_duration) - \
           (meta['start_times'][event_idx] - event_duration) >= 1:
            window = (time[mask][~in_transit] >= meta['start_times'][event_idx] - event_duration) & \
                     (time[mask][~in_transit] <= meta['end_times'][event_idx] + event_duration)
        else:
            mid_time = meta['start_times'][event_idx] + (meta['end_times'][event_idx] - meta['start_times'][event_idx]) / 2
            window = (time[mask][~in_transit] >= mid_time - EVENT_WINDOW_HALF_WIDTH) & \
                     (time[mask][~in_transit] <= mid_time + EVENT_WINDOW_HALF_WIDTH)

        event_scatter = np.std(flatten_lc[~in_transit][window])
        event_err = np.sqrt(event_flux_err**2 + event_scatter**2)

        cadence_minutes = self.transit_config['cadence_minutes']
        cadence_days = cadence_minutes * MINUTES_TO_DAYS
        event_duration_days = meta['widths'][event_idx]  # Time-span (includes gaps)
        event_duration_cadence = len(time[mask][in_transit]) * cadence_days  # For SNR only
        snr = round(get_snr(meta['depths'][event_idx], event_err, event_duration_cadence, cadence=cadence_minutes), 2)

        return {
            'time': event_time,
            'phase': event_phase,
            'depth': meta['depths'][event_idx],
            'width': meta['widths'][event_idx],
            'duration': event_duration_days,
            'start_time': meta['start_times'][event_idx],
            'end_time': meta['end_times'][event_idx],
            'snr': snr,
            'tic': tic,
            'sector': sector
        }

    def _create_event_snippet(self, time, flux_err, flatten_lc, mask, event_data, tic, sector, event_no):
        """Create a data snippet around transit event for downstream model comparison.

        Args:
            time (np.ndarray): Time array.
            flux_err (np.ndarray): Flux error array.
            flatten_lc (np.ndarray): Detrended light curve.
            mask (np.ndarray): Boolean mask of valid points.
            event_data (dict): Event dictionary containing:
                - start_time, end_time: Event window boundaries.
                - duration: Event duration (days).
                - time: Event time.
                - width: Event width from monofind.
            tic (int): TIC ID.
            sector (int): Sector number.
            event_no (int): Event number (for file naming)

        Returns:
            dict: Event snippet with keys:
                'tic', 'sector', 'event_no': Identifiers
                'time', 'flux', 'flux_err': Data arrays in snippet window
                'event_time': Time of mid-transit
                'event_width': Width from monofind
        """
        duration = event_data['duration']
        start_time = event_data['start_time']
        end_time = event_data['end_time']

        if (end_time + duration) - (start_time - duration) >= 1:
            event_window = (time[mask] >= start_time - duration) & \
                          (time[mask] <= end_time + duration)
        else:
            mid_time = start_time + (end_time - start_time) / 2
            event_window = (time[mask] >= mid_time - EVENT_WINDOW_HALF_WIDTH) & \
                          (time[mask] <= mid_time + EVENT_WINDOW_HALF_WIDTH)

        return {
            'tic': tic,
            'sector': sector,
            'event_no': event_no,
            'time': time[mask][event_window],
            'flux': flatten_lc[event_window],
            'flux_err': flux_err[mask][event_window],
            'event_time': event_data['time'],
            'event_width': duration
        }

    def _get_eclipse_params(self, tic):
        """Get eclipse parameters for a given TIC ID.
        
        Args:
            tic (int): TIC ID
        
        Returns:
            tuple: (prim_pos, prim_width, sec_pos, sec_width) or
                   (None, None, None, None) if TIC not found or catalogue is None
        """
        if self.catalogue is None:
            return None, None, None, None

        row = self.catalogue[self.catalogue['tess_id'] == tic]
        if row.empty:
            return None, None, None, None

        return (row['prim_pos'].values[0], row['prim_width'].values[0],
                row['sec_pos'].values[0], row['sec_width'].values[0])

    def _parse_filename(self, filename):
        """Parse TIC ID and sector number from filename.

        Expected filename format: TIC_<TICID>_<SECTOR>.<EXT>
        Sector should be 2 digits (e.g., 02 or 10) but returned without leading zeros.

        Args:
            filename (str): Filename to parse

        Returns:
            tuple: (tic_id, sector) where tic_id is int and sector is str without leading zeros

        Raises:
            ValueError: If filename format is invalid
        """
        try:
            parts = filename.split('_')
            tic_id = int(parts[1])
            sector_part = parts[2].split('.')[0]
            sector_num = int(sector_part)  # Convert to int to remove leading zero
            sector = str(sector_num)
            return tic_id, sector
        except (ValueError, IndexError) as e:
            raise ValueError(f"Failed to parse TIC ID and sector from filename '{filename}': {e}")

    def _calculate_skye_metric(self, plot_output_dir=None):
        """Calculate Skye metric to flag systematic artifacts.

        The Skye metric identifies events clustered in time across different targets within a sector.
        Inspired by RoboVetter: Thompson et al. (2018) (https://iopscience.iop.org/article/10.3847/1538-4365/aab4f9)
        and pterodactyls: Fernandes et al. (2022) (https://iopscience.iop.org/article/10.3847/1538-3881/ac7b29)

        Args:
            plot_output_dir (str, optional): Directory to save Skye plots. Defaults to None.
        """
        if self.sector_times is None:
            logger.warning("Sector times not provided, skipping Skye metric calculation")
            return

        if len(self.results['event_times']) == 0:
            logger.warning("No events detected, skipping Skye metric calculation")
            return

        event_times = np.array(self.results['event_times'])
        sectors = np.array(self.results['sectors'])
        unique_sectors = np.unique(sectors)

        skye_flags = np.zeros(len(event_times))

        for sector in unique_sectors:
            sector_idx = sectors == sector
            sector_event_times = event_times[sector_idx]

            if len(sector_event_times) == 0:
                continue

            sector_num = int(sector)
            sector_start = self.sector_times['start_time'].values[sector_num - 1]
            sector_end = self.sector_times['end_time'].values[sector_num - 1]

            # Calculate Skye threshold
            event_times_sorted = np.sort(sector_event_times)
            event_times_grouped = split_tol(event_times_sorted, 0.1)
            counts = np.array([len(g) for g in event_times_grouped])

            # Pad counts with zeros to match expected bin count for standard deviation calculation
            num_expected_bins = int(np.ceil((sector_end - sector_start) / 0.1))
            counts_zeros = np.concatenate((counts, np.zeros(max(0, num_expected_bins - len(counts)))))
            skye_threshold = round(3 * np.std(counts_zeros))

            # Minimum threshold
            if skye_threshold <= 2:
                skye_threshold = 3

            # Flag events above threshold
            sector_skye_flags = np.zeros(len(sector_event_times))
            for i, event_time in enumerate(sector_event_times):
                for j, group in enumerate(event_times_grouped):
                    if event_time in group and counts[j] >= skye_threshold:
                        sector_skye_flags[i] = 1

            # Update overall flags in original order
            original_sector_indices = np.where(sector_idx)[0]
            for i, idx in enumerate(original_sector_indices):
                skye_flags[idx] = sector_skye_flags[i]

            logger.info(f"Sector {sector}: Skye threshold = {skye_threshold}, "
                       f"flagged {int(np.sum(sector_skye_flags))} of {len(sector_event_times)} events")

            # Generate Skye plot if requested
            if self.transit_config['generate_skye_plots'] and plot_output_dir:
                bin_widths_plot = np.arange(start=sector_start, stop=sector_end, step=0.15)
                counts_plot, bins_plot = np.histogram(event_times_sorted, bins=bin_widths_plot)

                fig = plt.figure(figsize=SKYE_FIGURE_SIZE)
                plt.axhline(y=skye_threshold, color='red', linestyle='--', linewidth=5)
                plt.stairs(counts_plot, bins_plot, fill=True, color='black', alpha=0.5)
                plt.title(f'Sector {sector}: {len(sector_event_times)} events',
                         fontsize=30, pad=20)
                plt.ylabel('Number of events', fontsize=30)
                plt.xlabel('Time - 2457000 (BJD)', fontsize=30)
                plt.xticks(fontsize=30)
                plt.yticks(fontsize=30)

                # Save plot
                hist_dir = os.path.join(plot_output_dir, 'skye_histograms')
                os.makedirs(hist_dir, exist_ok=True)
                plot_path = os.path.join(hist_dir, f'Skye_S{sector}_hist.png')
                fig.savefig(plot_path)
                plt.close(fig)
                logger.info(f"Saved Skye histogram to {plot_path}")

        self.results['skye_flags'] = [int(flag) for flag in skye_flags]

    def _save_results(self, output_file='output.txt', output_dir=None):
        """Save results to disk.

        The output format and columns depend on the detrending method:
        - 'cb': includes 'win_len' and 'det_dependence' columns
        - 'cp': basic event information only
        - All methods: optionally include 'skye_flag' if Skye metric was calculated

        Args:
            output_file (str, optional): Output filename. Defaults to 'output.txt'.
            output_dir (str, optional): Output directory. Defaults to current working directory if None.

        Returns:
            pd.DataFrame: Results dataframe with columns depending on detrending method.
            Columns include: tic, sector, time, phase, depth, duration, snr, and optionally
            win_len, det_dependence (for 'cb' method), and skye_flag (if calculated).
        """
        if output_dir is None:
            output_dir = os.getcwd()

        output_path = os.path.join(output_dir, output_file)

        # Build output array
        method = self.transit_config['detrending_method']
        output_data = np.column_stack([
            self.results['tics'],
            self.results['sectors'],
            self.results['event_times'],
            self.results['event_phases'],
            self.results['event_depths'],
            self.results['event_durations'],
            self.results['event_snrs']
        ])

        # Define column headers
        if method == 'cb':
            output_data = np.column_stack([
                output_data,
                self.results['win_len_max_SNR'],
                self.results['det_dependence']
            ])
            # Add skye flags if available
            if len(self.results['skye_flags']) > 0:
                output_data = np.column_stack([output_data, self.results['skye_flags']])
                header = "TIC SECTOR TIME PHASE DEPTH DURATION SNR WIN_LEN_MAX_SNR DET_DEPENDENCE SKYE_FLAG"
            else:
                header = "TIC SECTOR TIME PHASE DEPTH DURATION SNR WIN_LEN_MAX_SNR DET_DEPENDENCE"
        elif method == 'cp':
            # Add skye flags if available
            if len(self.results['skye_flags']) > 0:
                output_data = np.column_stack([output_data, self.results['skye_flags']])
                header = "TIC SECTOR TIME PHASE DEPTH DURATION SNR SKYE_FLAG"
            else:
                header = "TIC SECTOR TIME PHASE DEPTH DURATION SNR"
        else:
            header = "TIC SECTOR TIME PHASE DEPTH DURATION SNR"

        # Sort by TIC then by time
        sort_idx = np.lexsort((output_data[:, 2].astype(float), output_data[:, 0]))
        output_data = output_data[sort_idx]

        # Save to file with header
        np.savetxt(output_path, output_data, fmt="%s", header=header, comments='')
        logger.info(f"Saved results to {output_path}")

        # Save event snippets if generated and save_event_snippets is enabled
        if (self.transit_config['generate_event_snippets'] and
            self.transit_config.get('save_event_snippets', True) and
            len(self.results['event_snippets']) > 0):
            snippet_dir = os.path.join(output_dir, 'event_snippets')
            os.makedirs(snippet_dir, exist_ok=True)
            for snippet in self.results['event_snippets']:
                snippet_file = f"TIC_{snippet['tic']}_{snippet['sector']}_{snippet['event_no']}.npz"
                snippet_path = os.path.join(snippet_dir, snippet_file)
                np.savez(snippet_path, **snippet)
            logger.info(f"Saved {len(self.results['event_snippets'])} event snippets to {snippet_dir}")
        elif (self.transit_config['generate_event_snippets'] and
              not self.transit_config.get('save_event_snippets', True) and
              len(self.results['event_snippets']) > 0):
            logger.info(f"Generated {len(self.results['event_snippets'])} event snippets in memory (not saved to disk)")

        # Return as DataFrame - columns must match output_data shape
        if method == 'cb':
            if len(self.results['skye_flags']) > 0:
                columns = ['tic', 'sector', 'time', 'phase', 'depth', 'duration', 'snr',
                          'win_len', 'det_dependence', 'skye_flag']
            else:
                columns = ['tic', 'sector', 'time', 'phase', 'depth', 'duration', 'snr',
                          'win_len', 'det_dependence']
        elif method == 'cp':
            if len(self.results['skye_flags']) > 0:
                columns = ['tic', 'sector', 'time', 'phase', 'depth', 'duration', 'snr', 'skye_flag']
            else:
                columns = ['tic', 'sector', 'time', 'phase', 'depth', 'duration', 'snr']
        else:
            columns = ['tic', 'sector', 'time', 'phase', 'depth', 'duration', 'snr']

        df = pd.DataFrame(output_data, columns=columns)
        return df

    @staticmethod
    def filter_events(events_df, min_snr=None, max_duration_days=None,
                      det_dependence_flag=None, skye_flag=None):
        """Filter detected events based on quality criteria.

        Args:
            events_df (pd.DataFrame): DataFrame of detected events from transit finding
            min_snr (float, optional): Minimum SNR threshold. Events with SNR < min_snr are filtered out.
            max_duration_days (float, optional): Maximum duration in days. Events with duration > max_duration_days are filtered out.
            det_dependence_flag (int, optional): Filter by detrending dependence flag.
                If 0, keep only events detected consistently across detrending windows.
                If 1, keep only events flagged as detrending-dependent.
            skye_flag (int, optional): Filter by Skye metric flag.
                If 0, keep only events not flagged by Skye metric.
                If 1, keep only events flagged by Skye metric.

        Returns:
            pd.DataFrame: Filtered events DataFrame
        """
        filtered_df = events_df.copy()

        # Apply SNR filter
        if min_snr is not None:
            filtered_df = filtered_df[filtered_df['snr'].astype(float) >= min_snr]

        # Apply duration filter
        if max_duration_days is not None:
            filtered_df = filtered_df[filtered_df['duration'].astype(float) <= max_duration_days]

        # Apply detrending dependence filter
        if det_dependence_flag is not None and 'det_dependence' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['det_dependence'].astype(int) == det_dependence_flag]

        # Apply Skye metric filter
        if skye_flag is not None and 'skye_flag' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['skye_flag'].astype(int) == skye_flag]

        return filtered_df
