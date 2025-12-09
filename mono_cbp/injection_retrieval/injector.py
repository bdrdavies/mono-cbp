"""Transit injection and retrieval for completeness testing.

This module injects synthetic transit signals into light curves and attempts
to recover them to assess detection completeness as a function of SNR, depth, duration, etc.
"""

import os
import numpy as np
import pandas as pd
import logging
import warnings
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
from ..utils import (
    bin_to_long_cadence, get_var_mad, monofind,
    time_to_phase, get_eclipse_mask, get_snr,
    load_transit_models
)
from ..utils.detrending import detrend
from ..config import get_default_config, merge_config
from ..utils.data import choose_eclipse_params

logger = logging.getLogger('mono_cbp.injection_retrieval')
warnings.simplefilter("ignore")


class TransitInjector:
    """Class for injection-retrieval testing of transit detection.

    This class:
    1. Loads light curves
    2. Inverts them (to ensure detected events aren't real signals)
    3. Injects synthetic transits
    4. Attempts to recover them
    5. Records recovery statistics

    Attributes:
        transit_models (dict): Precomputed transit models for injection
        catalogue (pd.DataFrame): Catalogue with eclipse and orbital parameters
        config (dict): Configuration parameters
        transit_config (dict): Transit finding configuration
        results (dict): Recovery results
        stats (list): Recovery statistics
    """

    def __init__(self, transit_models_path, catalogue=None, config=None, TEBC=False):
        """Initialise the TransitInjector.

        Args:
            transit_models_path (str): Path to .npz file with transit models
            catalogue (str or pd.DataFrame, optional): Catalogue with eclipse and orbital parameters.
                Must contain columns: tess_id, period, bjd0, prim_pos, prim_width, sec_pos, sec_width.
            config (dict, optional): Configuration dictionary
            TEBC (bool, optional): If True, use TEBC catalogue format with *_2g and *_pf columns.
                Defaults to False.
        """
        # Load transit models
        self.transit_models = load_transit_models(transit_models_path)
        logger.info(f"Loaded {len(self.transit_models['models'])} transit models from {transit_models_path}")

        # Load catalogue
        if catalogue is not None:
            if isinstance(catalogue, str):
                from ..utils import load_catalogue
                self.catalogue = load_catalogue(catalogue, TEBC=TEBC)
            else:
                # If DataFrame is passed directly, apply TEBC processing if needed
                if TEBC:
                    catalogue = catalogue.copy()
                    eclipse_params = catalogue.apply(choose_eclipse_params, axis=1)
                    catalogue['prim_pos'] = eclipse_params.apply(lambda x: x['prim_pos'])
                    catalogue['prim_width'] = eclipse_params.apply(lambda x: x['prim_width'])
                    catalogue['sec_pos'] = eclipse_params.apply(lambda x: x['sec_pos'])
                    catalogue['sec_width'] = eclipse_params.apply(lambda x: x['sec_width'])
                    logger.info("Applied TEBC eclipse parameter mapping to DataFrame")
                self.catalogue = catalogue
        else:
            self.catalogue = None

        # Configuration
        self.config = merge_config(config, get_default_config()) if config else get_default_config()
        self.transit_config = self.config['transit_finding']

        # Results storage
        self.results = {
            'tics': [],
            'sectors': [],
            'injected_times': [],
            'injected_depths': [],
            'injected_durations': [],
            'injected_snrs': [],
            'recovered': [],
            'recovered_times': [],
            'recovered_depths': [],
            'recovered_durations': [],
            'recovered_snrs': []
        }

        self.stats = []  # List of dicts with statistics for each transit model

    def run_injection_retrieval(self, data_dir, n_injections=None, output_file='inj-ret_results.csv',
                                output_dir=None):
        """Run injection-retrieval tests for all transit models.

        Tests each transit model in transit_models.npz by injecting it into n_injections
        randomly selected light curves. The total number of tests will be
        n_injections * number_of_models.

        Args:
            data_dir (str): Directory containing light curve files
            n_injections (int, optional): Number of injection-retrieval tests to perform
                per transit model. If there are fewer files than requested injections,
                files will be randomly sampled with replacement. Defaults to config value.
            output_file (str, optional): Output filename. Defaults to 'inj-ret_results.csv'.
            output_dir (str, optional): Output directory. Defaults to data_dir.

        Returns:
            pd.DataFrame: Results dataframe with one row per injection test
        """
        n_models = len(self.transit_models['models'])
        n_injections = self.config['injection_retrieval']['n_injections'] if n_injections is None else n_injections
        total_tests = n_injections * n_models
        logger.info(f"Starting injection-retrieval tests: {n_injections} injections x {n_models} models = {total_tests} total tests")

        # Get all available files (both .npz and .txt)
        all_files = [f for f in os.listdir(data_dir) if f.endswith('.npz') or f.endswith('.txt')]

        if len(all_files) == 0:
            raise ValueError(f"No light curve files found in {data_dir}")

        # Loop through all transit models
        test_count = 0
        for model_idx, model in enumerate(self.transit_models['models']):
            # Extract model parameters
            flux_model = model['flux']
            depth_model = model['depth']
            duration_model = model['duration']

            logger.info(f"Model {model_idx+1}/{n_models}: depth={depth_model:.4f}, duration={duration_model:.4f}")

            # Track the starting index for this model
            model_start_idx = len(self.results['recovered'])

            # Sample files with replacement if needed
            if len(all_files) < n_injections:
                logger.info(f"  Sampling {n_injections} files with replacement from {len(all_files)} available")
                selected_files = np.random.choice(all_files, size=n_injections, replace=True)
            else:
                logger.info(f"  Sampling {n_injections} files without replacement from {len(all_files)} available")
                selected_files = np.random.choice(all_files, size=n_injections, replace=False)

            # Process each file with this model
            for file in selected_files:
                test_count += 1
                if test_count % 50 == 0:
                    logger.info(f"  Progress: {test_count}/{total_tests} total tests")

                file_path = os.path.join(data_dir, file)
                self.process_file(file_path, flux_model, depth_model, duration_model)

            # Calculate and store statistics
            model_end_idx = len(self.results['recovered'])
            model_recoveries = sum(self.results['recovered'][model_start_idx:model_end_idx])
            model_injections = model_end_idx - model_start_idx
            model_recovery_rate = model_recoveries / model_injections if model_injections > 0 else 0.0

            self.stats.append({
                'model_idx': model_idx,
                'depth': depth_model,
                'duration': duration_model,
                'n_injections': model_injections,
                'n_recoveries': model_recoveries,
                'recovery_rate': model_recovery_rate
            })

            logger.info(f"  Model recovery rate: {model_recovery_rate:.2%} ({model_recoveries}/{model_injections})")

        logger.info(f"Injection-retrieval complete: {total_tests} total tests")

        # Save results
        results_df = self.save_results(output_file, output_dir or data_dir)

        # Save statistics
        stats_file = output_file.replace('.csv', '_stats.csv')
        self.save_stats(stats_file, output_dir or data_dir)

        return results_df

    def process_file(self, file_path, flux_model, depth_model, duration_model):
        """Process a single file through injection-retrieval test.

        Mirrors the structure of finder.py's process_file method, with the addition of
        lightcurve inversion and transit injection.

        Args:
            file_path (str): Path to light curve file
            flux_model (array): Transit model flux
            depth_model (float): Injected transit depth
            duration_model (float): Injected transit duration
        """
        # Parse filename
        file = os.path.basename(file_path)
        split_file = os.path.splitext(file)

        if split_file[1] == '.npz':
            # Load from npz
            data = np.load(file_path, allow_pickle=True)
            time, flux, flux_err = data['time'], data['flux'], data['flux_err']
            phase = data.get('phase', None)
            ecl_mask_raw = data.get('eclipse_mask', None)
            ecl_mask = ecl_mask_raw.astype(bool) if ecl_mask_raw is not None else None
        elif split_file[1] == '.txt':
            # Load from txt
            data = np.loadtxt(file_path, skiprows=1)
            time, flux, flux_err, phase = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
            ecl_mask = data[:, 4].astype(bool) if data.shape[1] > 4 else None
        else:
            logger.warning(f"Unsupported file format: {file}")
            return

        # Extract TIC and sector
        tic, sector = self._parse_filename(file)
        if tic is None:
            logger.warning(f"Could not parse TIC/sector from {file}, skipping")
            return

        # Bin to 30-minute cadence if needed
        if np.median(np.gradient(time[~np.isnan(flux)])) < self.transit_config['cadence_minutes'] / (60 * 24):
            time, flux, flux_err = bin_to_long_cadence(time, flux, flux_err)
            if self.catalogue is not None:
                row = self.catalogue[self.catalogue['tess_id'] == tic]
                if not row.empty:
                    phase = time_to_phase(time, row['period'].values[0], row['bjd0'].values[0])
                    prim_mask = get_eclipse_mask(phase, row['prim_pos'].values[0], row['prim_width'].values[0])
                    sec_mask = get_eclipse_mask(phase, row['sec_pos'].values[0], row['sec_width'].values[0])
                    ecl_mask = np.logical_or(prim_mask, sec_mask)

        # Create mask
        # Exclude NaN values and eclipse regions
        nan_mask = ~np.isnan(flux * time * flux_err)
        if ecl_mask is not None:
            mask = nan_mask & ~ecl_mask  # Out-of-eclipse regions
        else:
            mask = nan_mask

        # Invert lightcurve
        # Ensures recovered events are injected, not real
        flux = ((flux - 1) * -1) + 1

        # Inject transit at random time
        mask_indices = np.where(mask)[0]
        if len(mask_indices) == 0:
            logger.warning(f"{file}: No out-of-eclipse data available for injection")
            return
        inj_idx = np.random.choice(mask_indices)
        inj_time = time[inj_idx]
        flux_inj, inj_time = self._inject_transit(time, flux, flux_model, inj_time)

        logger.debug(f"{file}: Injected at t={inj_time:.2f}, depth={depth_model:.4f}, dur={duration_model:.4f}")

        # Detrend
        method = self.transit_config['detrending_method']
        detrend_result = detrend(
            time, flux_inj, flux_err, method=method, fname=split_file[0],
            mask=mask, edge_cutoff=self.transit_config['edge_cutoff'],
            cos_win_len_max=self.transit_config['cosine']['win_len_max'],
            cos_win_len_min=self.transit_config['cosine']['win_len_min'],
            fap_threshold=self.transit_config['cosine']['fap_threshold'],
            poly_order=self.transit_config['cosine']['poly_order'],
            max_splines=self.transit_config['pspline']['max_splines'],
            bi_win_len_max=self.transit_config['biweight']['win_len_max'],
            bi_win_len_min=self.transit_config['biweight']['win_len_min']
        )

        # Process events and check for recovery
        if method == 'cp':
            self._process_cp_injection(
                time, flux_err, mask,
                detrend_result, tic, sector,
                depth_model, duration_model, inj_time
            )
        elif method == 'cb':
            self._process_cb_injection(
                time, flux_err, mask,
                detrend_result, tic, sector,
                depth_model, duration_model, inj_time
            )
        else:
            logger.error(f"Unknown detrending method: {method}")

    def _process_cp_injection(self, time, flux_err, mask,
                              detrend_result, tic, sector,
                              depth_model, duration_model, inj_time):
        """Process injection-retrieval using cosine + pspline detrending."""
        flatten_lc, _, _, _ = detrend_result
        mad = self.transit_config['mad_threshold']
        var_mad = get_var_mad(flatten_lc, 100)

        # Run monofind on detrended lightcurve
        peaks, meta = monofind(time[mask], flatten_lc, mad=mad, var_mad=var_mad)

        if len(peaks) == 0:
            # No events found - injection not recovered
            self._record_injection_result(
                time, tic, sector, inj_time, depth_model, duration_model,
                flux_err, flatten_lc, mask, None, None, None, None
            )
            return

        # Check if any detected event matches the injection
        event_times = time[mask][peaks]
        time_diffs = np.abs(event_times - inj_time)

        if np.min(time_diffs) >= duration_model / 2:
            # Closest event is too far from injection - not recovered
            self._record_injection_result(
                time, tic, sector, inj_time, depth_model, duration_model,
                flux_err, flatten_lc, mask, None, None, None, None
            )
            return

        # Event found within recovery window
        closest_idx = np.argmin(time_diffs)
        recovered_time = event_times[closest_idx]

        # Extract event measurements
        in_transit = (time >= meta['start_times'][closest_idx]) & \
                    (time <= meta['end_times'][closest_idx])
        recovered_duration = np.sum(in_transit) * (30 / 1440)

        # Calculate depth and SNR
        recovered_depth, recovered_snr = self._measure_event_depth_snr(
            flux_err, flatten_lc, mask, in_transit
        )

        # Record result
        self._record_injection_result(
            time, tic, sector, inj_time, depth_model, duration_model,
            flux_err, flatten_lc, mask, recovered_time, recovered_depth,
            recovered_duration, recovered_snr
        )

    def _process_cb_injection(self, time, flux_err, mask,
                              detrend_result, tic, sector,
                              depth_model, duration_model, inj_time):
        """Process injection-retrieval using cosine + biweight detrending."""
        flatten_lcs, _, _, _ = detrend_result
        mad = self.transit_config['mad_threshold']

        # Storage for all detected events across biweight windows
        all_peaks = []
        all_meta = []
        all_flatten_lcs = []

        # Loop over biweight windows
        for lc in flatten_lcs:
            var_mad = get_var_mad(lc, 100)
            peaks, meta = monofind(time[mask], lc, mad=mad, var_mad=var_mad)
            all_peaks.append(peaks)
            all_meta.append(meta)
            all_flatten_lcs.append(lc)

        # Flatten peaks
        all_peaks_flat = [p for ps in all_peaks for p in ps]

        if len(all_peaks_flat) == 0:
            # No events found - injection not recovered
            self._record_injection_result(
                time, tic, sector, inj_time, depth_model, duration_model,
                flux_err, flatten_lcs[0], mask, None, None, None, None
            )
            return

        # Find closest event across all windows
        # Track: (window_idx, peak_idx_in_window)
        all_event_times = []
        all_event_info = []
        for window_idx, peaks in enumerate(all_peaks):
            for peak_idx_in_window, peak_idx in enumerate(peaks):
                event_time = time[mask][peak_idx]
                all_event_times.append(event_time)
                all_event_info.append((window_idx, peak_idx_in_window))

        all_event_times = np.array(all_event_times)
        time_diffs = np.abs(all_event_times - inj_time)

        if np.min(time_diffs) >= duration_model / 2:
            # Closest event is too far from injection - not recovered
            self._record_injection_result(
                time, tic, sector, inj_time, depth_model, duration_model,
                flux_err, flatten_lcs[0], mask, None, None, None, None
            )
            return

        # Event found within recovery window
        closest_event_idx = np.argmin(time_diffs)
        window_idx, peak_idx_in_window = all_event_info[closest_event_idx]
        recovered_time = all_event_times[closest_event_idx]

        # Extract event measurements from the best biweight window length
        meta = all_meta[window_idx]
        flatten_lc = all_flatten_lcs[window_idx]

        # Use the peak index within this window length's metadata
        in_transit = (time >= meta['start_times'][peak_idx_in_window]) & \
                    (time <= meta['end_times'][peak_idx_in_window])
        recovered_duration = np.sum(in_transit) * (30 / 1440)

        # Calculate depth and SNR
        recovered_depth, recovered_snr = self._measure_event_depth_snr(
            flux_err, flatten_lc, mask, in_transit
        )

        # Record result
        self._record_injection_result(
            time, tic, sector, inj_time, depth_model, duration_model,
            flux_err, flatten_lc, mask, recovered_time, recovered_depth,
            recovered_duration, recovered_snr
        )

    def _measure_event_depth_snr(self, flux_err, flatten_lc, mask, in_transit):
        """Measure event depth and SNR from detrended lightcurve."""
        if np.sum(in_transit) == 0:
            return np.nan, np.nan

        # Create out-of-transit mask
        out_transit = ~in_transit & mask

        # Map to masked time array
        in_transit_on_masked = in_transit[mask]
        out_transit_on_masked = out_transit[mask]

        if np.sum(in_transit_on_masked) == 0 or np.sum(out_transit_on_masked) == 0:
            return np.nan, np.nan

        # Measure depth from detrended lightcurve
        detrended_flux_in = flatten_lc[in_transit_on_masked]
        detrended_flux_out = flatten_lc[out_transit_on_masked]

        # Remove NaN values
        detrended_flux_in = detrended_flux_in[~np.isnan(detrended_flux_in)]
        detrended_flux_out = detrended_flux_out[~np.isnan(detrended_flux_out)]

        if len(detrended_flux_in) == 0 or len(detrended_flux_out) == 0:
            return np.nan, np.nan

        median_in = np.median(detrended_flux_in)
        median_out = np.median(detrended_flux_out)
        depth = np.abs(median_out - median_in)

        # Calculate SNR
        event_flux_err = np.median(flux_err[in_transit])
        event_scatter = np.std(flatten_lc[out_transit_on_masked])
        event_err = np.sqrt(event_flux_err**2 + event_scatter**2)
        event_duration = np.sum(in_transit) * (30 / 1440)
        snr = get_snr(depth, event_err, event_duration)

        return depth, snr

    def _record_injection_result(self, time, tic, sector, inj_time, depth_model, duration_model,
                                 flux_err, flatten_lc, mask, recovered_time, recovered_depth,
                                 recovered_duration, recovered_snr):
        """Record injection-retrieval results."""
        # Calculate injected SNR
        # nearby = points within 1 day of injection (on full time array)
        nearby = np.abs(time - inj_time) < 1.0

        # in_transit = points within the injected transit (EXCLUDE from scatter calculation!)
        in_transit = np.abs(time - inj_time) < (duration_model / 2)

        # nearby_and_masked = points that are nearby, out-of-eclipse, AND not in the injected transit
        nearby_and_masked = nearby & mask & ~in_transit

        if np.any(nearby_and_masked):
            # Get flux error from nearby points
            inj_flux_err = np.median(flux_err[nearby_and_masked])

            # Get scatter from detrended lightcurve
            # flatten_lc corresponds to time[mask], so we need to map nearby_and_masked to the masked indices
            # Create a mapping: which elements of the masked array correspond to nearby_and_masked?
            masked_indices = np.where(mask)[0]
            nearby_masked_indices = np.where(nearby_and_masked)[0]

            # Find positions in the flattened array that correspond to nearby points
            flatten_lc_positions = np.isin(masked_indices, nearby_masked_indices)

            if np.any(flatten_lc_positions):
                inj_scatter = np.std(flatten_lc[flatten_lc_positions])
            else:
                # Fallback: use flux_err as scatter estimate
                inj_scatter = np.std(flux_err[nearby_and_masked])

            inj_err = np.sqrt(inj_flux_err**2 + inj_scatter**2)
            injected_snr = get_snr(depth_model, inj_err, duration_model)
        else:
            # If no nearby out-of-eclipse out-of-transit points, expand search to 2 days
            nearby_wider = np.abs(time - inj_time) < 2.0
            nearby_wider_and_masked = nearby_wider & mask & ~in_transit

            if np.any(nearby_wider_and_masked):
                inj_flux_err = np.median(flux_err[nearby_wider_and_masked])

                # Map to flatten_lc
                masked_indices = np.where(mask)[0]
                nearby_wider_masked_indices = np.where(nearby_wider_and_masked)[0]
                flatten_lc_positions = np.isin(masked_indices, nearby_wider_masked_indices)

                if np.any(flatten_lc_positions):
                    inj_scatter = np.std(flatten_lc[flatten_lc_positions])
                else:
                    inj_scatter = np.std(flux_err[nearby_wider_and_masked])

                inj_err = np.sqrt(inj_flux_err**2 + inj_scatter**2)
                injected_snr = get_snr(depth_model, inj_err, duration_model)
            else:
                # Last resort: use all out-of-eclipse, out-of-transit points
                out_of_transit_masked = mask & ~in_transit
                if np.any(out_of_transit_masked):
                    inj_flux_err = np.median(flux_err[out_of_transit_masked])
                    inj_scatter = np.std(flatten_lc)  # Use all detrended data
                    inj_err = np.sqrt(inj_flux_err**2 + inj_scatter**2)
                    injected_snr = get_snr(depth_model, inj_err, duration_model)
                else:
                    injected_snr = np.nan

        # Record results
        self.results['tics'].append(tic)
        self.results['sectors'].append(sector)
        self.results['injected_times'].append(inj_time)
        self.results['injected_depths'].append(depth_model)
        self.results['injected_durations'].append(duration_model)
        self.results['injected_snrs'].append(injected_snr)
        self.results['recovered'].append(recovered_time is not None)
        self.results['recovered_times'].append(recovered_time if recovered_time is not None else np.nan)
        self.results['recovered_depths'].append(recovered_depth if recovered_depth is not None else np.nan)
        self.results['recovered_durations'].append(recovered_duration if recovered_duration is not None else np.nan)
        self.results['recovered_snrs'].append(recovered_snr if recovered_snr is not None else np.nan)

        if recovered_time is not None:
            logger.debug(f"  -> Recovered! dt={abs(recovered_time - inj_time):.3f}d")
        else:
            logger.debug(f"  -> Not recovered")

    def _inject_transit(self, time, flux, flux_model, inj_time):
        """Inject a transit model into the light curve.

        Args:
            time (array): Time array
            flux (array): Flux array
            flux_model (array): Transit model flux (centered at 0)
            inj_time (float): Injection time

        Returns:
            tuple: (flux with injection, actual injection time)
        """
        # Find closest cadence to injection time
        injection_mid = np.argmin(np.abs(time - inj_time))
        actual_inj_time = time[injection_mid]

        # Calculate injection window
        half_model_len = len(flux_model) // 2
        injection_start = max(0, injection_mid - half_model_len)
        injection_end = injection_start + len(flux_model)

        # Handle boundary conditions
        model = flux_model.copy()
        if injection_end > len(time):
            injection_end = len(time)
            model = model[:injection_end - injection_start]

        if injection_start == 0:
            model = model[half_model_len - injection_mid:]

        # Inject into the lightcurve
        flux_inj = flux.copy()
        min_length = min(len(flux_inj[injection_start:injection_end]), len(model))
        flux_inj[injection_start:injection_start + min_length] += model[:min_length]

        return flux_inj, actual_inj_time

    def _parse_filename(self, filename):
        """Parse TIC and sector from filename."""
        try:
            parts = filename.split('_')
            tic = int(parts[1])
            sector_str = parts[2].replace('.npz', '').replace('.txt', '')
            sector = int(sector_str.lstrip('S').lstrip('0') or '0')
            return tic, sector
        except (IndexError, ValueError):
            return None, None

    def save_results(self, output_file='inj-ret_results.csv', output_dir='.'):
        """Save injection-retrieval results.

        Args:
            output_file (str, optional): Output filename. Defaults to 'inj-ret_results.csv'.
            output_dir (str, optional): Output directory. Defaults to current working directory ('.').

        Returns:
            pd.DataFrame: Results DataFrame
        """
        df = pd.DataFrame(self.results)

        # Handle case where output_file contains a directory path
        output_file_dir = os.path.dirname(output_file)
        if output_file_dir:
            # output_file contains a directory, use it
            output_dir = output_file_dir
            output_file = os.path.basename(output_file)

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_file)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved results to {output_path}")
        return df

    def save_stats(self, output_file='inj-ret_stats.csv', output_dir='.'):
        """Save recovery statistics.

        Args:
            output_file (str, optional): Output filename. Defaults to 'inj-ret_stats.csv'.
            output_dir (str, optional): Output directory. Defaults to current working directory ('.').

        Returns:
            pd.DataFrame: Injection-recovery statistics DataFrame
        """
        df = pd.DataFrame(self.stats)

        # Handle case where output_file contains a directory path
        output_file_dir = os.path.dirname(output_file)
        if output_file_dir:
            # output_file contains a directory, use it
            output_dir = output_file_dir
            output_file = os.path.basename(output_file)

        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_file)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved recovery statistics to {output_path}")
        return df

    def plot_completeness(self, stats_file=None, figsize=(5, 4), cmap='viridis', save_fig=False,
                          output_path='completeness.png', dpi=300, font_family='sans-serif',
                          font_size=8):
        """Plot completeness matrix from injection-retrieval statistics.

        Creates a heatmap showing recovery rate as a function of transit depth and duration.

        Args:
            stats_file (str, optional): Path to existing stats CSV file. If None, uses self.stats. Defaults to None.
            figsize (tuple, optional): Figure size (width, height) in inches. Defaults to (10, 10).
            cmap (str, optional): Colormap name. Defaults to 'viridis'.
            save_fig (bool, optional): Whether to save the figure. Defaults to False.
            output_path (str, optional): Path to save figure. Defaults to 'completeness.png'.
            dpi (int, optional): DPI for saved figure. Defaults to 300.
            font_family (str, optional): Font family for plot. Defaults to 'sans-serif'.
            font_size (int, optional): Font size in points. Defaults to 8.

        Returns:
            tuple: (fig, ax) matplotlib figure and axes objects
        """

        # Load stats from file if provided
        if stats_file:
            stats_df = pd.read_csv(stats_file)

        elif self.stats:
            stats_df = pd.DataFrame(self.stats)
        
        else:
            logger.warning("No statistics available. Run injection-retrieval first.")
            return None, None

        # Configure matplotlib
        plt.rc('font', family=font_family, size=font_size)
        plt.rc('axes', titlesize=font_size, labelsize=font_size)
        plt.rc('xtick', labelsize=font_size)
        plt.rc('ytick', labelsize=font_size)
        plt.rc('legend', fontsize=font_size - 1)
        plt.rc('figure', titlesize=font_size)

        # Get and sort unique depths and durations
        depths_unique = np.sort(stats_df['depth'].unique())
        durations_unique = np.sort(stats_df['duration'].unique())

        # Create recovery rate matrix
        recovery_matrix = np.zeros((len(depths_unique), len(durations_unique)))
        for i, depth in enumerate(depths_unique):
            for j, duration in enumerate(durations_unique):
                matching = stats_df[(stats_df['depth'] == depth) & (stats_df['duration'] == duration)]
                if not matching.empty:
                    recovery_matrix[i, j] = matching['recovery_rate'].values[0]

        # Create figure and plot
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Create heatmap
        im = ax.imshow(recovery_matrix, aspect='auto', cmap=cmap, origin='lower', vmin=0, vmax=1)

        # Add text annotations
        for i in range(len(depths_unique)):
            for j in range(len(durations_unique)):
                text = ax.text(j, i, f"{recovery_matrix[i, j]:.2f}",
                        ha='center', va='center', color='white', fontsize=font_size)
                text.set_path_effects([withStroke(linewidth=1, foreground='black')])

        # Set ticks and labels
        ax.set_xticks(range(len(durations_unique)))
        ax.set_yticks(range(len(depths_unique)))
        ax.set_xticklabels([f"{d:.2f}" for d in durations_unique])
        ax.set_yticklabels([f"{d*100:.2f}" for d in depths_unique])

        ax.set_xlabel("Duration (days)")
        ax.set_ylabel("Depth (%)")

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, label="Recovery Rate")
        cbar.ax.set_yticklabels([f"{x:.2f}" for x in cbar.get_ticks()])
        cbar.ax.tick_params(left=False, right=True)

        fig.tight_layout()

        # Save or show
        if save_fig:
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
            logger.info(f"Saved completeness plot to {output_path}")
        else:
            plt.show()

        return fig, ax
