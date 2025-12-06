"""Bayesian model comparison for vetting and classifying TCEs detected with TransitFinder."""

import numpy as np
import matplotlib.pyplot as plt
import pymc as pm
import exoplanet as xo
from pymc.model import Model
import pytensor.tensor as pyt
import os
import scipy.stats as stats
import warnings
import pandas as pd
import logging
from ..config import get_default_config, merge_config

logger = logging.getLogger('mono_cbp.model_comparison')
warnings.filterwarnings("ignore", category=UserWarning)

class ModelComparator:
    """Class for comparing models to classify TCEs.

    This class fits multiple models (transit, sinusoid, linear, step) to
    detected events and uses the Akaike Information Criterion (AIC) and RMSE
    of the residual light curve fits to classify them.

    Classifications (based on AIC difference and RMSE):
        - 'T': Transit (transit model best fit, AIC difference >= 2, RMSE <= rmse_threshold)
        - 'AT': Ambiguous transit (transit model best fit, AIC difference >= 2, RMSE > rmse_threshold)
        - 'Sin': Sinusoid (sinusoid model best fit, AIC difference >= 2, RMSE <= rmse_threshold)
        - 'ASin': Ambiguous sinusoid (sinusoid model best fit, AIC difference >= 2, RMSE > rmse_threshold)
        - 'L': Linear (linear model best fit, AIC difference >= 2, RMSE <= rmse_threshold)
        - 'AL': Ambiguous linear (linear model best fit, AIC difference >= 2, RMSE > rmse_threshold)
        - 'St': Step (step model best fit, AIC difference >= 2, RMSE <= rmse_threshold)
        - 'ASt': Ambiguous step (step model best fit, AIC difference >= 2, RMSE > rmse_threshold)
        - 'A': Ambiguous (AIC difference < 2, no single clear best fit)

    Attributes:
        config (dict): Configuration parameters
        model_config (dict): Specific configuration for model comparison
    """

    def __init__(self, config=None):
        """Initialise ModelComparator.

        Args:
            config (dict, optional): Configuration dictionary. Uses defaults if None.
        """
        self.config = merge_config(config, get_default_config()) if config else get_default_config()
        self.model_config = self.config['model_comparison']
        logger.info("Initialised ModelComparator")

    def compare_event(self, event_input, save_plot=False, plot_dir=None):
        """Compare models for a single event.

        Args:
            event_input (str or dict): Either a path to event data file (.npz) or a dictionary
                                       with keys: 'time', 'flux', 'flux_err', 'event_time', 'event_width'.
                                       Optional dict keys: 'tic', 'sector', 'event_no'
            save_plot (bool, optional): Whether to save diagnostic plot. Defaults to False.
            plot_dir (str, optional): Directory to save plot. Defaults to None.

        Returns:
            dict: Dictionary with classification and model comparison results

        Raises:
            TypeError: If event_input is neither a string nor a dictionary.
        """
        # Determine if input is file path or in-memory data
        if isinstance(event_input, str):
            # Load from file
            logger.info(f"Processing {os.path.basename(event_input)}")
            data = np.load(event_input)
            time, flux, flux_err = data['time'], data['flux'], data['flux_err']
            time = np.asarray(time, dtype=np.float64)
            flux = np.asarray(flux, dtype=np.float64)
            flux_err = np.asarray(flux_err, dtype=np.float64)
            event_time = data['event_time']
            event_width = data['event_width']
            filename = os.path.basename(event_input)
        elif isinstance(event_input, dict):
            # Load from in-memory dictionary
            time = np.asarray(event_input['time'], dtype=np.float64)
            flux = np.asarray(event_input['flux'], dtype=np.float64)
            flux_err = np.asarray(event_input['flux_err'], dtype=np.float64)
            event_time = event_input['event_time']
            event_width = event_input['event_width']

            # Create filename for logging/plotting
            tic = event_input.get('tic', 'unknown')
            sector = event_input.get('sector', 'unknown')
            event_no = event_input.get('event_no', 0)
            filename = f"TIC_{tic}_{sector}_{event_no}"
            logger.info(f"Processing TIC {tic} Sector {sector} Event {event_no}")
        else:
            raise TypeError("event_input must be either a file path (str) or event data (dict)")

        # Mask NaN values
        nan_mask = ~np.isnan(flux * time * flux_err)
        time, flux, flux_err = time[nan_mask], flux[nan_mask], flux_err[nan_mask]

        # Remove outliers
        flux_median = np.median(flux)
        flux_std = np.std(flux)
        sigma_threshold = self.model_config['sigma_threshold']
        outliers = (flux - flux_median) > sigma_threshold * flux_std
        time, flux, flux_err = time[~outliers], flux[~outliers], flux_err[~outliers]

        # Estimate event depth
        in_transit = (time > (event_time - event_width / 2)) & (time < (event_time + event_width / 2))
        event_depth = abs(np.median(flux[~in_transit]) - np.median(flux[in_transit]))

        # Fit models
        logger.debug("Fitting transit model...")
        trace_transit, n_params_transit = self._fit_transit_model(
            time, flux, flux_err, event_time, event_width, event_depth
        )

        logger.debug("Fitting sinusoidal model...")
        trace_sinusoidal, n_params_sinusoidal = self._fit_sinusoidal_model(
            time, flux, flux_err
        )

        logger.debug("Fitting linear model...")
        trace_linear, n_params_linear = self._fit_linear_model(
            time, flux, flux_err
        )

        logger.debug("Fitting step model...")
        trace_step, n_params_step = self._fit_step_model(
            time, flux, flux_err
        )

        # Extract posterior samples
        q16_transit, q50_transit, q84_transit = np.percentile(
            trace_transit.posterior["__light_curve"].values, [16, 50, 84], axis=(0, 1)
        )
        q16_sinusoidal, q50_sinusoidal, q84_sinusoidal = np.percentile(
            trace_sinusoidal.posterior["__sinusoid"].values, [16, 50, 84], axis=(0, 1)
        )
        q16_linear, q50_linear, q84_linear = np.percentile(
            trace_linear.posterior["__linear"].values, [16, 50, 84], axis=(0, 1)
        )
        q16_step, q50_step, q84_step = np.percentile(
            trace_step.posterior["__step"].values, [16, 50, 84], axis=(0, 1)
        )

        # Define model uncertainties
        transit_model_err = abs(q84_transit - q16_transit) / 2
        sinusoidal_model_err = abs(q84_sinusoidal - q16_sinusoidal) / 2
        linear_model_err = abs(q84_linear - q16_linear) / 2
        step_model_err = abs(q84_step - q16_step) / 2

        # Calculate total errors
        transit_total_err = np.sqrt(flux_err**2 + transit_model_err**2)
        sinusoidal_total_err = np.sqrt(flux_err**2 + sinusoidal_model_err**2)
        linear_total_err = np.sqrt(flux_err**2 + linear_model_err**2)
        step_total_err = np.sqrt(flux_err**2 + step_model_err**2)

        # Calculate residuals
        residuals_transit = (flux - q50_transit) / transit_total_err
        residuals_sinusoidal = (flux - q50_sinusoidal) / sinusoidal_total_err
        residuals_linear = (flux - q50_linear) / linear_total_err
        residuals_step = (flux - q50_step) / step_total_err

        # Calculate RMSE
        rmse_transit = np.sqrt(np.mean(residuals_transit**2))
        rmse_sinusoidal = np.sqrt(np.mean(residuals_sinusoidal**2))
        rmse_linear = np.sqrt(np.mean(residuals_linear**2))
        rmse_step = np.sqrt(np.mean(residuals_step**2))

        # Calculate log-likelihoods
        log_likelihoods = np.array([
            np.sum(stats.norm.logpdf(flux, loc=q50_transit, scale=transit_total_err)),
            np.sum(stats.norm.logpdf(flux, loc=q50_sinusoidal, scale=sinusoidal_total_err)),
            np.sum(stats.norm.logpdf(flux, loc=q50_linear, scale=linear_total_err)),
            np.sum(stats.norm.logpdf(flux, loc=q50_step, scale=step_total_err))
        ])

        n_params = np.array([n_params_transit, n_params_sinusoidal, n_params_linear, n_params_step])

        # Calculate AIC
        aic_arr = 2 * n_params - 2 * log_likelihoods
        aic_order = np.argsort(aic_arr)

        # Classify event
        rmse_threshold = self.model_config['rmse_threshold']
        best_fit = self._classify_event(aic_arr, aic_order, rmse_transit, rmse_sinusoidal, rmse_linear, rmse_step, rmse_threshold)

        logger.info(f"Classification: {best_fit}, RMSE_transit: {rmse_transit:.2f}")

        # Prepare results
        results = {
            'filename': filename,
            'best_fit': best_fit,
            'aic_transit': aic_arr[0],
            'aic_sinusoidal': aic_arr[1],
            'aic_linear': aic_arr[2],
            'aic_step': aic_arr[3],
            'rmse_transit': rmse_transit,
            'rmse_sinusoidal': rmse_sinusoidal,
            'rmse_linear': rmse_linear,
            'rmse_step': rmse_step,
        }

        # Save plot if requested
        if save_plot:
            self._save_comparison_plot(
                time, flux, flux_err,
                q50_transit, q16_transit, q84_transit,
                q50_sinusoidal, q16_sinusoidal, q84_sinusoidal,
                q50_linear, q16_linear, q84_linear,
                q50_step, q16_step, q84_step,
                filename, plot_dir
            )

        return results

    def compare_events(self, events_input, output_file='model_comparison_results.csv', output_dir=None,
                      save_plots=None, plot_dir=None):
        """Compare models for multiple events.

        Args:
            events_input (str or list): Either a directory path containing event .npz files,
                                       or a list of event data dictionaries with keys:
                                       'time', 'flux', 'flux_err', 'event_time', 'event_width'
            output_file (str, optional): Output CSV filename. Defaults to 'model_comparison_results.csv'.
            output_dir (str, optional): Output directory. Defaults to events_input if it's a directory.
            save_plots (bool, optional): Whether to save comparison plots. Defaults to config value.
            plot_dir (str, optional): Directory for plots. Defaults to config value or output_dir.

        Returns:
            pd.DataFrame: DataFrame of all results

        Raises:
            TypeError: If events_input is neither a string nor a list.
        """
        # Use config values if not explicitly provided
        if save_plots is None:
            save_plots = self.model_config.get('save_plots', False)
        if plot_dir is None:
            plot_dir = self.model_config.get('plot_dir', None)
            # Fall back to output_dir if plot_dir not specified
            if plot_dir is None and save_plots:
                plot_dir = output_dir
        # Determine if input is directory or list of events
        if isinstance(events_input, str):
            # Process directory of .npz files
            data_dir = events_input
            # If output_dir not specified, check if output_file contains a directory path
            if output_dir is None:
                output_file_dir = os.path.dirname(output_file)
                if output_file_dir:
                    # output_file contains a directory, use it
                    output_dir = output_file_dir
                    output_file = os.path.basename(output_file)
                else:
                    # output_file is just a filename, default to data_dir
                    output_dir = data_dir

            files = [f for f in os.listdir(data_dir) if f.endswith('.npz')]
            logger.info(f"Processing {len(files)} files from {data_dir}")

            all_results = []
            for i, file in enumerate(files):
                if i % 10 == 0:
                    logger.info(f"Progress: {i}/{len(files)}")

                file_path = os.path.join(data_dir, file)
                try:
                    results = self.compare_event(file_path, save_plot=save_plots, plot_dir=plot_dir)
                    all_results.append(results)
                except Exception as e:
                    logger.error(f"Failed to process {file}: {e}")

        elif isinstance(events_input, list):
            # Process list of in-memory event dictionaries
            logger.info(f"Processing {len(events_input)} events from memory")

            # Handle output_file with directory path
            if output_dir is None:
                output_file_dir = os.path.dirname(output_file)
                if output_file_dir:
                    output_dir = output_file_dir
                    output_file = os.path.basename(output_file)
                else:
                    output_dir = '.'  # Current directory as default

            all_results = []
            for i, event_data in enumerate(events_input):
                if i % 10 == 0:
                    logger.info(f"Progress: {i}/{len(events_input)}")

                try:
                    results = self.compare_event(event_data, save_plot=save_plots, plot_dir=plot_dir)
                    all_results.append(results)
                except Exception as e:
                    logger.error(f"Failed to process event {i}: {e}")

        else:
            raise TypeError("events_input must be either a directory path (str) or list of event dictionaries")

        # Save results
        df = pd.DataFrame(all_results)
        if output_dir:
            # Ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, output_file)
            df.to_csv(output_path, index=False)
            logger.info(f"Saved results to {output_path}")

        # Log summary statistics
        if len(df) > 0:
            classifications = df['best_fit'].value_counts()
            logger.info(f"Classification summary:\n{classifications}")

        return df

    def _fit_transit_model(self, time, flux, flux_err, event_time, event_width, event_depth):
        """Fit transit model to data."""
        with Model():
            mean = pm.Normal("mean", mu=1.0, sigma=0.01)
            # Transit mid-time prior
            t0 = pm.Normal("t0", mu=event_time, sigma=0.05)

            # Limb-darkening coefficients prior
            u = xo.distributions.quad_limb_dark("u", initval=np.array([0.3, 0.2]))

            # Event depth prior
            depth = pm.Normal("event_depth", mu=event_depth, sigma=event_depth/2)

            # Impact parameter prior
            b = pm.Uniform("b", lower=0, upper=1)

            # Approximate radius ratio from transit depth
            ror = xo.LimbDarkLightCurve(u).get_ror_from_approx_transit_depth(depth, b)

            # Duration prior
            duration = pm.Normal("duration", mu=event_width, sigma=0.1)

            # Orbit (arbitrary period longer than snippet)
            period = 10
            orbit = xo.orbits.SimpleTransitOrbit(period=period, duration=duration, t0=t0, ror=ror, b=b)

            # Light curve
            light_curve = (xo.LimbDarkLightCurve(u).get_light_curve(
                orbit=orbit, r=ror, t=time)[:, 0] + mean)

            pm.Deterministic("__light_curve", light_curve)
            pm.Normal("obs_transit", mu=light_curve, sigma=flux_err, observed=flux)

            trace = pm.sample(
                tune=self.model_config['tune'],
                draws=self.model_config['draws'],
                cores=self.model_config['cores'],
                chains=self.model_config['chains'],
                init="adapt_full",
                target_accept=self.model_config['target_accept'],
                progressbar=False
            )

        return trace, 8  # 8 parameters

    def _fit_sinusoidal_model(self, time, flux, flux_err):
        """Fit sinusoidal model to data."""
        with Model():
            mean = pm.Normal("mean", mu=1.0, sigma=0.01)
            amplitude = pm.Normal("amplitude", mu=(np.max(flux) - np.min(flux)), sigma=0.1)
            phase = pm.Uniform("phase", lower=0, upper=2 * np.pi)
            frequency = pm.Uniform("frequency", lower=1/(time[-1]-time[0]), upper=1/((time[-1]-time[0])*0.5))

            sinusoid = amplitude * pm.math.sin(2 * np.pi * frequency * time + phase) + mean

            pm.Deterministic("__sinusoid", sinusoid)
            pm.Normal("obs_sin", mu=sinusoid, sigma=flux_err, observed=flux)

            trace = pm.sample(
                tune=self.model_config['tune'],
                draws=self.model_config['draws'],
                cores=self.model_config['cores'],
                chains=self.model_config['chains'],
                init="adapt_full",
                target_accept=self.model_config['target_accept'],
                progressbar=False
            )

        return trace, 4  # 4 parameters

    def _fit_linear_model(self, time, flux, flux_err):
        """Fit linear model to data."""
        with Model():
            coeff_guess = np.polyfit(time, flux, 1)
            coeff_0 = pm.Normal("coeff_1", mu=coeff_guess[1], sigma=0.1)
            coeff_1 = pm.Normal("coeff_2", mu=coeff_guess[0], sigma=0.1)
            linear = coeff_0 + coeff_1 * time

            pm.Deterministic("__linear", linear)
            pm.Normal("obs_linear", mu=linear, sigma=flux_err, observed=flux)

            trace = pm.sample(
                tune=self.model_config['tune'],
                draws=self.model_config['draws'],
                cores=self.model_config['cores'],
                chains=self.model_config['chains'],
                init="adapt_full",
                target_accept=self.model_config['target_accept'],
                progressbar=False
            )

        return trace, 2  # 2 parameters

    def _fit_step_model(self, time, flux, flux_err):
        """Fit step/polynomial model to data."""
        with Model():
            # Determine if there's a significant step
            flux_diff = np.diff(flux)
            max_diff_idx = np.argmax(abs(flux_diff))
            step_time = time[max_diff_idx]
            diff_std = np.std(flux_diff)

            if abs(flux_diff[max_diff_idx]) > 3 * diff_std:
                # Fit piecewise quadratic
                coeff_guess_1 = np.polyfit(time[time <= step_time], flux[time <= step_time], 2)
                coeff_guess_2 = np.polyfit(time[time > step_time], flux[time > step_time], 2)

                coeff_01 = pm.Normal("coeff_01", mu=coeff_guess_1[2], sigma=0.1)
                coeff_11 = pm.Normal("coeff_11", mu=coeff_guess_1[1], sigma=0.1)
                coeff_21 = pm.Normal("coeff_21", mu=coeff_guess_1[0], sigma=0.1)
                polynomial_1 = coeff_01 + coeff_11 * time[time <= step_time] + coeff_21 * time[time <= step_time]**2

                coeff_02 = pm.Normal("coeff_02", mu=coeff_guess_2[2], sigma=0.1)
                coeff_12 = pm.Normal("coeff_12", mu=coeff_guess_2[1], sigma=0.1)
                coeff_22 = pm.Normal("coeff_22", mu=coeff_guess_2[0], sigma=0.1)
                polynomial_2 = coeff_02 + coeff_12 * time[time > step_time] + coeff_22 * time[time > step_time]**2

                step = pyt.concatenate([polynomial_1, polynomial_2])
                n_params = 6
            else:
                # Fit simple quadratic
                coeff_guess = np.polyfit(time, flux, 2)
                coeff_0 = pm.Normal("coeff_0", mu=coeff_guess[2], sigma=0.1)
                coeff_1 = pm.Normal("coeff_1", mu=coeff_guess[1], sigma=0.1)
                coeff_2 = pm.Normal("coeff_2", mu=coeff_guess[0], sigma=0.1)
                step = coeff_0 + coeff_1 * time + coeff_2 * time**2
                n_params = 3

            pm.Deterministic("__step", step)
            pm.Normal("obs_step", mu=step, sigma=flux_err, observed=flux)

            trace = pm.sample(
                tune=self.model_config['tune'],
                draws=self.model_config['draws'],
                cores=self.model_config['cores'],
                chains=self.model_config['chains'],
                init="adapt_full",
                target_accept=self.model_config['target_accept'],
                progressbar=False
            )

        return trace, n_params

    def _classify_event(self, aic_arr, aic_order, rmse_transit, rmse_sinusoidal, rmse_linear, rmse_step, rmse_threshold):
        """Classify event based on AIC and RMSE.

        Args:
            aic_arr (array): AIC values for [transit, sinusoid, linear, step]
            aic_order (array): Indices sorted by AIC
            rmse_transit (float): RMSE for transit model
            rmse_sinusoidal (float): RMSE for sinusoidal model
            rmse_linear (float): RMSE for linear model
            rmse_step (float): RMSE for step model
            rmse_threshold (float): RMSE threshold for ambiguity classification

        Returns:
            str: Classification result
        """
        # Model names mapping (index in aic_arr)
        model_names = ["Transit", "Sinusoidal", "Linear", "Step"]
        rmse_values = [rmse_transit, rmse_sinusoidal, rmse_linear, rmse_step]

        # Get the best fit model (lowest AIC)
        best_fit_idx = aic_order[0]
        best_fit_model = model_names[best_fit_idx]
        best_fit_aic = aic_arr[best_fit_idx]

        # Get all other models' AICs
        other_models_aic = np.array([aic_arr[i] for i in range(len(aic_arr)) if i != best_fit_idx])

        # Check if best fit satisfies AIC criterion: min_aic - all_others <= -2
        is_unambiguous = np.all(best_fit_aic - other_models_aic <= -2)

        if is_unambiguous:
            # Classify based on model type and RMSE
            rmse_best = rmse_values[best_fit_idx]

            if best_fit_model == "Transit":
                return "T" if rmse_best <= rmse_threshold else "AT"
            elif best_fit_model == "Sinusoidal":
                return "Sin" if rmse_best <= rmse_threshold else "ASin"
            elif best_fit_model == "Linear":
                return "L" if rmse_best <= rmse_threshold else "AL"
            elif best_fit_model == "Step":
                return "St" if rmse_best <= rmse_threshold else "ASt"
        else:
            # Ambiguous classification
            return "A"

    def _save_comparison_plot(self, time, flux, flux_err,
                             q50_transit, q16_transit, q84_transit,
                             q50_sinusoidal, q16_sinusoidal, q84_sinusoidal,
                             q50_linear, q16_linear, q84_linear,
                             q50_step, q16_step, q84_step,
                             filename, plot_dir):
        """Save diagnostic plot comparing all model fits."""
        if plot_dir is None:
            plot_dir = './plots'

        os.makedirs(plot_dir, exist_ok=True)

        # Find gaps to avoid connecting lines
        gap_indices = np.where(np.diff(time) > (np.median(np.diff(time)) * 10))[0]

        plt.figure(figsize=(10, 5))
        plt.errorbar(time, flux, yerr=flux_err, fmt='k.', label="Observed Flux")

        # Plot model fits without connecting over gaps
        for i in range(len(gap_indices) + 1):
            start_idx = 0 if i == 0 else gap_indices[i-1] + 1
            end_idx = len(time) if i == len(gap_indices) else gap_indices[i] + 1

            plt.plot(time[start_idx:end_idx], q50_transit[start_idx:end_idx], 'r-',
                    label="Transit" if i == 0 else "")
            plt.fill_between(time[start_idx:end_idx], q16_transit[start_idx:end_idx],
                           q84_transit[start_idx:end_idx], alpha=0.3, color='r')
            plt.plot(time[start_idx:end_idx], q50_sinusoidal[start_idx:end_idx], 'b-',
                    label="Sinusoid" if i == 0 else "")
            plt.fill_between(time[start_idx:end_idx], q16_sinusoidal[start_idx:end_idx],
                           q84_sinusoidal[start_idx:end_idx], alpha=0.3, color='b')
            plt.plot(time[start_idx:end_idx], q50_linear[start_idx:end_idx], 'g-',
                    label="Linear" if i == 0 else "")
            plt.fill_between(time[start_idx:end_idx], q16_linear[start_idx:end_idx],
                           q84_linear[start_idx:end_idx], alpha=0.3, color='g')
            plt.plot(time[start_idx:end_idx], q50_step[start_idx:end_idx], '-', c='orange',
                    label="Step" if i == 0 else "")
            plt.fill_between(time[start_idx:end_idx], q16_step[start_idx:end_idx],
                           q84_step[start_idx:end_idx], alpha=0.3, color='orange')

        plt.legend()
        plt.xlabel("Time - 2457000 (BTJD days)")
        plt.ylabel("Normalised Flux")

        plot_path = os.path.join(plot_dir, filename.replace(".npz", ".png"))
        plt.savefig(plot_path)
        plt.close()
        logger.info(f"Saved comparison plot to {plot_path}")