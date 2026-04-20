"""Bayesian model comparison for vetting and classifying TCEs detected with TransitFinder."""

import numpy as np
import matplotlib.pyplot as plt
import batman
import lmfit
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
        q50_transit, n_params_transit, transit_params = self._fit_transit_model(
            time, flux, flux_err, event_time, event_width, event_depth
        )

        logger.debug("Fitting sinusoidal model...")
        q50_sinusoidal, n_params_sinusoidal = self._fit_sinusoidal_model(
            time, flux, flux_err
        )

        logger.debug("Fitting linear model...")
        q50_linear, n_params_linear = self._fit_linear_model(
            time, flux, flux_err
        )

        logger.debug("Fitting step model...")
        q50_step, n_params_step = self._fit_step_model(
            time, flux, flux_err
        )

        # Calculate residuals normalised by flux_err only.
        # Using total_err (which includes eval_uncertainty from lmfit's covariance
        # matrix) is inconsistent: for ill-constrained parameters the covariance
        # can be huge, artificially reducing the RMSE and simultaneously inflating
        # the AIC via the -log(sigma) term in the Gaussian likelihood. flux_err is the
        # only noise term that was used during fitting, so it is the correct
        # normaliser for both RMSE and log-likelihood.
        residuals_transit = (flux - q50_transit) / flux_err
        residuals_sinusoidal = (flux - q50_sinusoidal) / flux_err
        residuals_linear = (flux - q50_linear) / flux_err
        residuals_step = (flux - q50_step) / flux_err

        # Calculate RMSE
        rmse_transit = np.sqrt(np.mean(residuals_transit**2))
        rmse_sinusoidal = np.sqrt(np.mean(residuals_sinusoidal**2))
        rmse_linear = np.sqrt(np.mean(residuals_linear**2))
        rmse_step = np.sqrt(np.mean(residuals_step**2))

        # Calculate log-likelihoods using flux_err (consistent with fitting)
        log_likelihoods = np.array([
            np.sum(stats.norm.logpdf(flux, loc=q50_transit, scale=flux_err)),
            np.sum(stats.norm.logpdf(flux, loc=q50_sinusoidal, scale=flux_err)),
            np.sum(stats.norm.logpdf(flux, loc=q50_linear, scale=flux_err)),
            np.sum(stats.norm.logpdf(flux, loc=q50_step, scale=flux_err))
        ])

        n_params = np.array([n_params_transit, n_params_sinusoidal, n_params_linear, n_params_step])

        # Calculate AIC
        aic_arr = 2 * n_params - 2 * log_likelihoods
        aic_order = np.argsort(aic_arr)

        # Classify event
        aic_threshold = self.model_config['aic_threshold']
        rmse_threshold = self.model_config['rmse_threshold']
        best_fit = self._classify_event(aic_arr, aic_order, rmse_transit, rmse_sinusoidal, rmse_linear, rmse_step, rmse_threshold, aic_threshold)

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

        # Add transit model parameters when transit is the preferred classification
        if best_fit in ('T', 'AT'):
            results.update(transit_params)

        # Save plot if requested
        if save_plot:
            self._save_comparison_plot(
                time, flux, flux_err,
                q50_transit, q50_sinusoidal,
                q50_linear, q50_step,
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
        """Fit transit model using batman + lmfit."""
        PERIOD = 100.0  # fixed, arbitrary - must exceed snippet duration

        rp_init = np.sqrt(max(event_depth, 1e-6))
        # For b=0: T ~ per*(1+rp)/(pi*a) => a ~ per*(1+rp)/(pi*T)
        a_init = max(PERIOD * (1.0 + rp_init) / (np.pi * event_width), 3.0)

        # Create batman model once; update params in-place during fitting
        bp = batman.TransitParams()
        bp.t0 = event_time
        bp.per = PERIOD
        bp.rp = rp_init
        bp.a = a_init
        bp.inc = 90.0
        bp.ecc = 0.0
        bp.w = 90.0
        bp.u = [0.3, 0.2]
        bp.limb_dark = "quadratic"
        bat_model = batman.TransitModel(bp, time)

        # `t` is the independent variable passed by lmfit; bat_model already
        # holds the same time array so it is not used directly in the body.
        #
        # Limb darkening uses the Kipping (2013) triangular reparametrisation:
        #   q1 = (u1+u2)^2,  q2 = u1 / (2*(u1+u2))
        # with q1, q2 within [0, 1].  This guarantees the physical constraints
        # u1 >= 0 and u1+u2 <= 1 by construction, avoiding negative stellar
        # intensity at the limb which causes batman to produce flux > 1.
        def transit_func(t, t0, rp, a, b, q1, q2, mean):
            bp.t0 = t0
            bp.rp = abs(rp)   # rp must be > 0
            bp.a = abs(a)     # a must be > 1
            # Convert impact parameter b to inclination: b = a*cos(inc)
            # Clamp b < a to keep the argument of arccos in (-1, 1)
            bp.inc = np.degrees(np.arccos(min(abs(b), abs(a) - 1e-6) / abs(a)))
            sq_q1 = np.sqrt(abs(q1))
            bp.u = [2.0 * sq_q1 * q2, sq_q1 * (1.0 - 2.0 * q2)]
            # batman returns 1 out-of-transit; mean shifts the baseline
            return bat_model.light_curve(bp) + (mean - 1.0)

        out_of_transit = ~((time > event_time - event_width / 2) &
                           (time < event_time + event_width / 2))
        baseline = np.median(flux[out_of_transit]) if out_of_transit.any() else np.median(flux)

        # Initial q1, q2 from u1=0.3, u2=0.2: q1=(0.5)^2=0.25, q2=0.3/(2*0.5)=0.3
        lm = lmfit.Model(transit_func, independent_vars=['t'])
        params = lm.make_params(
            t0=dict(value=event_time,
                    min=event_time - event_width, max=event_time + event_width),
            rp=dict(value=rp_init, min=0.001, max=0.5),
            a=dict(value=a_init, min=2.0, max=1000.0),
            b=dict(value=0.3, min=0.0, max=1.0),
            q1=dict(value=0.25, min=0.0, max=1.0),
            q2=dict(value=0.3, min=0.0, max=1.0),
            mean=dict(value=baseline, min=0.8, max=1.2),
        )

        result = lm.fit(flux, params, t=time, weights=1.0/flux_err,
                        method='least_squares')

        best_fit = result.best_fit

        # Extract best-fit parameter values and 1-sigma uncertainties
        def _stderr(name):
            s = result.params[name].stderr
            return float(s) if s is not None else np.nan

        rp_val = float(abs(result.params['rp'].value))
        a_val  = float(abs(result.params['a'].value))
        b_val  = float(abs(result.params['b'].value))

        # Build the 3x3 covariance submatrix for (rp, a, b) if available
        cov_sub = None
        if result.covar is not None:
            try:
                idx = {name: i for i, name in enumerate(result.var_names)}
                ii = [idx['rp'], idx['a'], idx['b']]
                cov_sub = result.covar[np.ix_(ii, ii)]
                if not np.isfinite(cov_sub).all():
                    cov_sub = None
            except (KeyError, IndexError):
                cov_sub = None

        duration, duration_err = self._compute_transit_duration(
            PERIOD, rp_val, a_val, b_val,
            _stderr('rp'), _stderr('a'), _stderr('b'),
            cov=cov_sub,
        )

        transit_params = {
            'transit_t0':           float(result.params['t0'].value),
            'transit_t0_err':       _stderr('t0'),
            'transit_rp':           rp_val,
            'transit_rp_err':       _stderr('rp'),
            'transit_a':            a_val,
            'transit_a_err':        _stderr('a'),
            'transit_b':            b_val,
            'transit_b_err':        _stderr('b'),
            'transit_q1':           float(result.params['q1'].value),
            'transit_q1_err':       _stderr('q1'),
            'transit_q2':           float(result.params['q2'].value),
            'transit_q2_err':       _stderr('q2'),
            'transit_mean':         float(result.params['mean'].value),
            'transit_mean_err':     _stderr('mean'),
            'transit_duration':     duration,
            'transit_duration_err': duration_err,
        }

        return best_fit, 7, transit_params  # 7 free parameters (t0, rp, a, b, u1, u2, mean)

    def _fit_sinusoidal_model(self, time, flux, flux_err):
        """Fit sinusoidal model to data."""
        T = time[-1] - time[0]

        def sinusoid_func(t, amplitude, frequency, phase, mean):
            return amplitude * np.sin(2 * np.pi * frequency * t + phase) + mean

        lm = lmfit.Model(sinusoid_func)
        params = lm.make_params(
            amplitude=dict(value=0.5 * (np.max(flux) - np.min(flux)), min=0.0),
            # Period constrained between 0.5*T and T (one or two cycles in snippet)
            frequency=dict(value=1.5 / T, min=1.0 / T, max=2.0 / T),
            # Phase unconstrained: bounding a circular parameter causes boundary artefacts
            phase=dict(value=0.0),
            mean=dict(value=np.median(flux)),
        )

        result = lm.fit(flux, params, t=time, weights=1.0 / flux_err,
                        method='least_squares')

        best_fit = result.best_fit

        return best_fit, 4  # amplitude, frequency, phase, mean

    def _fit_linear_model(self, time, flux, flux_err):
        """Fit linear model to data."""
        def linear_func(t, slope, intercept):
            return slope * t + intercept

        lm = lmfit.Model(linear_func)
        params = lm.make_params(
            slope=dict(value=0.0),
            intercept=dict(value=np.median(flux)),
        )

        result = lm.fit(flux, params, t=time, weights=1.0 / flux_err,
                        method='least_squares')

        best_fit = result.best_fit

        return best_fit, 2  # slope, intercept

    def _fit_step_model(self, time, flux, flux_err):
        """Fit step/polynomial model to data."""
        flux_diff = np.diff(flux)
        max_diff_idx = np.argmax(abs(flux_diff))
        step_time = time[max_diff_idx]
        diff_std = np.std(flux_diff)

        flux_mean = np.median(flux)

        if abs(flux_diff[max_diff_idx]) > 3 * diff_std:
            # Piecewise quadratic: different polynomial each side of step_time
            def step_func(t, c01, c11, c21, c02, c12, c22):
                return np.where(
                    t <= step_time,
                    c01 + c11 * t + c21 * t**2,
                    c02 + c12 * t + c22 * t**2,
                )

            lm = lmfit.Model(step_func)
            params = lm.make_params(
                c01=dict(value=flux_mean), c11=dict(value=0.0), c21=dict(value=0.0),
                c02=dict(value=flux_mean), c12=dict(value=0.0), c22=dict(value=0.0),
            )
            n_params = 6
        else:
            # Simple quadratic
            def step_func(t, c0, c1, c2):
                return c0 + c1 * t + c2 * t**2

            lm = lmfit.Model(step_func)
            params = lm.make_params(
                c0=dict(value=flux_mean),
                c1=dict(value=0.0),
                c2=dict(value=0.0),
            )
            n_params = 3

        result = lm.fit(flux, params, t=time, weights=1.0 / flux_err,
                        method='least_squares')

        best_fit = result.best_fit

        return best_fit, n_params

    @staticmethod
    def _compute_transit_duration(period, rp, a, b, rp_err, a_err, b_err, cov=None):
        """Compute transit duration and propagated 1-sigma uncertainty.

        Uses the circular-orbit formula:
            T14 = (P/π) * arcsin( sqrt((1+rp)² - b²) / sqrt(a² - b²) )

        Args:
            period (float): Orbital period (days)
            rp (float): Radius ratio Rp/R★
            a (float): Semi-major axis in stellar radii
            b (float): Impact parameter
            rp_err, a_err, b_err (float): 1-sigma uncertainties (NaN if unavailable)
            cov (ndarray, optional): 3x3 covariance matrix for (rp, a, b). When
                provided, the full correlated propagation σ² = gᵀCg is used,
                which avoids over-estimating the uncertainty due to anti-correlated
                parameters (e.g. a and b). Falls back to diagonal propagation if None.

        Returns:
            tuple: (duration, duration_err) in days
        """
        num2 = (1.0 + rp)**2 - b**2
        den2 = a**2 - b**2

        if num2 <= 0.0 or den2 <= 0.0:
            return np.nan, np.nan

        num = np.sqrt(num2)
        den = np.sqrt(den2)
        arg = num / den

        if arg >= 1.0:
            return np.nan, np.nan

        duration = (period / np.pi) * np.arcsin(arg)

        # Gradient of T w.r.t. (rp, a, b)
        dT_darg = (period / np.pi) / np.sqrt(1.0 - arg**2)
        g = np.array([
            dT_darg * (1.0 + rp) / (num * den),           # dT/drp
            dT_darg * (-a * num / den**3),                 # dT/da
            dT_darg * (b * ((1.0 + rp)**2 - a**2) / (num * den**3)),  # dT/db
        ])

        if cov is not None:
            var = float(g @ cov @ g)
        else:
            # Diagonal (independent) fallback
            errs = np.array([rp_err, a_err, b_err])
            finite = np.isfinite(errs)
            var = float(np.sum((g[finite] * errs[finite])**2))

        duration_err = np.sqrt(var) if var > 0.0 else np.nan
        return duration, duration_err

    def _classify_event(self, aic_arr, aic_order, rmse_transit, rmse_sinusoidal, rmse_linear, rmse_step, rmse_threshold, aic_threshold=2):
        """Classify event based on AIC and RMSE.

        Args:
            aic_arr (array): AIC values for [transit, sinusoid, linear, step]
            aic_order (array): Indices sorted by AIC
            rmse_transit (float): RMSE for transit model
            rmse_sinusoidal (float): RMSE for sinusoidal model
            rmse_linear (float): RMSE for linear model
            rmse_step (float): RMSE for step model
            rmse_threshold (float): RMSE threshold for ambiguity classification
            aic_threshold (float): Minimum AIC difference for an unambiguous classification

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

        # Check if best fit satisfies AIC criterion: best_aic is at least aic_threshold lower than all others
        is_unambiguous = np.all(best_fit_aic - other_models_aic <= -aic_threshold)

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
                             q50_transit, q50_sinusoidal,
                             q50_linear, q50_step,
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
            plt.plot(time[start_idx:end_idx], q50_sinusoidal[start_idx:end_idx], 'b-',
                    label="Sinusoid" if i == 0 else "")
            plt.plot(time[start_idx:end_idx], q50_linear[start_idx:end_idx], 'g-',
                    label="Linear" if i == 0 else "")
            plt.plot(time[start_idx:end_idx], q50_step[start_idx:end_idx], '-', c='orange',
                    label="Step" if i == 0 else "")

        plt.legend()
        plt.xlabel("Time - 2457000 (BTJD days)")
        plt.ylabel("Normalised Flux")

        plot_path = os.path.join(plot_dir, filename.replace(".npz", ".png"))
        plt.savefig(plot_path)
        plt.close()
        logger.info(f"Saved comparison plot to {plot_path}")