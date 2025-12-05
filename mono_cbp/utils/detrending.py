"""Collection of functions for detrending TESS light curves."""

import numpy as np
import logging
import sys
import os
from contextlib import contextmanager
from wotan import flatten
from astropy.timeseries import LombScargle
from .monofind import get_gaps_indices

logger = logging.getLogger('mono_cbp.utils.detrending')


@contextmanager
def suppress_stdout():
    """Context manager to suppress stdout temporarily.
       This is used to suppress output from the iterative cosine detrending.
    """
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def poly_normalise(time, flux, order=2):
    """Normalises the flux using a polynomial fit.

    Used for reducing the half-sector periodicity from
    TESS data before calculating the Lomb-Scargle periodogram.

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        order (int, optional): The order of the polynomial. Defaults to 2.

    Returns:
        np.ndarray: The normalised flux
    """
    nan_mask = ~np.isnan(flux)
    polyfit = np.polyfit(time[nan_mask], flux[nan_mask], order)

    if order == 2:
        flux_poly = flux / (polyfit[0]*time**2+polyfit[1]*time+polyfit[2])

    if order == 1:
        flux_poly = flux / (polyfit[0]*time +polyfit[1])

    return flux_poly


def poly_normalise_gaps(time, flux, gap_size, order=2):
    """Performs a polynomial fit, but splits the data across a specified gap length.

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        gap_size (float): The size of the gaps to split the data across (in days)
        order (int, optional): The order of the polynomial. Defaults to 2.

    Returns:
        np.ndarray: The normalised flux
    """
    gaps_idx = get_gaps_indices(time, gap_size)

    norm_flux = []
    for i in range(len(gaps_idx)-1):
        flux_chunk = flux[gaps_idx[i]:gaps_idx[i+1]]
        time_chunk = time[gaps_idx[i]:gaps_idx[i+1]]
        norm_flux.append(poly_normalise(time_chunk, flux_chunk, order=order))

    return np.concatenate(norm_flux)


def calculate_LSP(time, flux, flux_err):
    """Calculates the Lomb-Scargle periodogram.

    The frequency range is set to search for periods between 2 and 27 days.

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        flux_err (array_like): Array of flux error values

    Returns:
        tuple: A tuple containing the LombScargle object, frequency array, and power spectrum
    """
    nan_mask = ~np.isnan(flux * time * flux_err)

    ls = LombScargle(time[nan_mask], flux[nan_mask], dy=flux_err[nan_mask])
    freq = np.linspace(1/27, 1/2, 1000)
    power = ls.power(frequency=freq)

    return ls, freq, power


def get_period_max_power(time, flux, flux_err):
    """Calculates the period with maximum power from the Lomb-Scargle periodogram.

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        flux_err (array_like): Array of flux error values

    Returns:
        float: The period with maximum power
    """
    _, freq, power = calculate_LSP(time, flux, flux_err)
    max_power_period = 1/freq[np.argmax(power)]

    return max_power_period


def get_fap(time, flux, flux_err):
    """Calculates the false alarm probability of the peak of the Lomb-Scargle periodogram.

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        flux_err (array_like): Array of flux error values

    Returns:
        float: The false alarm probability at the peak of the periodogram
    """
    ls, _, power = calculate_LSP(time, flux, flux_err)

    return ls.false_alarm_probability(power.max())


def cosine_detrend(time, flux, flux_err, win_len_max=12, win_len_min=1, threshold=1e-2, poly_order=2, mask=None, edge_cutoff=0):
    """Performs iterative cosine detrending on the input light curve.

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        flux_err (array_like): Array of flux error values
        win_len_max (int, optional): Maximum window length for cosine fitting (in days). Defaults to 12.
        win_len_min (int, optional): Minimum window length for cosine fitting (in days). Defaults to 1.
        threshold (float, optional): Threshold for false alarm probability. Defaults to 1e-2.
        poly_order (int, optional): Order of polynomial for initial detrending. Defaults to 2.
        mask (array_like, optional): Boolean mask for data points to exclude from the wotan fitting functions. Defaults to None.
        edge_cutoff (int, optional): Defines the amount of data at the edges to exclude (in units of time). Defaults to 0.

    Returns:
        tuple or np.ndarray: Detrended light curve and fitted trend (if successful), or just detrended flux

    Raises:
        ValueError: If maximum window length is smaller than minimum window length
    """
    # Check that the input window lengths are sensible
    if win_len_max <= win_len_min:
        raise ValueError("Maximum window length must be greater than minimum window length.")

    # If the periodicity is longer than a TESS sector or if there is no significant periodicity
    # Fit polynomial (accounting for gaps of size > 1 day) and re-calculate LSP
    if get_period_max_power(time, flux, flux_err) > 27 or get_fap(time, flux, flux_err) > threshold:

        logger.debug('Fitting polynomial to remove trend before L-S periodogram...')

        flux_poly = poly_normalise_gaps(time, flux, 1, order=poly_order)

        # If still no significant periodicity, return the original flux
        if get_period_max_power(time, flux_poly, flux_err) > 27 or get_fap(time, flux_poly, flux_err) > threshold:
            return flux

        else:
            with suppress_stdout():
                cosine_lc, cosine_fit = flatten(time, flux, window_length=win_len_max, return_trend=True, method='cosine', mask=mask, edge_cutoff=edge_cutoff, robust=True)
            win_len = win_len_max

            # Iteratively reduce window length
            while get_fap(time, cosine_lc, flux_err) < threshold:
                win_len -= 0.1

                # Check that the window length does not go below the minimum
                if win_len < win_len_min:
                    win_len += 0.1
                    break

                with suppress_stdout():
                    cosine_lc, cosine_fit = flatten(time, flux, window_length=win_len, return_trend=True, method='cosine', mask=mask, edge_cutoff=edge_cutoff, robust=True)

            return cosine_lc, cosine_fit, win_len

    else:
        with suppress_stdout():
            cosine_lc, cosine_fit = flatten(time, flux, window_length=win_len_max, return_trend=True, method='cosine', mask=mask, edge_cutoff=edge_cutoff, robust=True)
        win_len = win_len_max

        # Iteratively reduce window length
        while get_fap(time, cosine_lc, flux_err) < threshold:
            win_len -= 0.1

            # Check that the window length does not go below the minimum
            if win_len < win_len_min:
                win_len += 0.1
                break

            with suppress_stdout():
                cosine_lc, cosine_fit = flatten(time, flux, window_length=win_len, return_trend=True, method='cosine', mask=mask, edge_cutoff=edge_cutoff, robust=True)

        return cosine_lc, cosine_fit, win_len


def slider_detrend(time, flux, win_len, mask=None, edge_cutoff=0):
    """Apply sliding window detrending to the light curve (biweight).

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        win_len (float): Length of the sliding window (in days)
        mask (array_like, optional): Boolean mask for data points to exclude from the fitting. Defaults to None.
        edge_cutoff (int, optional): Defines the amount of data at the edges to exclude (in units of time). Defaults to 0.

    Returns:
        tuple: Detrended light curve and fitted trend
    """
    flatten_lc, trend_lc = flatten(time, flux, window_length=win_len, return_trend=True, method='biweight', mask=mask, edge_cutoff=edge_cutoff)

    return flatten_lc, trend_lc


def run_multi_biweight(time, flux, max_win_len=3, min_win_len=1, edge_cutoff=0):
    """Run biweight detrending over a range of window lengths.

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        max_win_len (int, optional): Maximum window length for biweight detrending (in days). Defaults to 3.
        min_win_len (int, optional): Minimum window length for biweight detrending (in days). Defaults to 1.
        edge_cutoff (int, optional): Defines the amount of data at the edges to exclude (in days). Defaults to 0.

    Returns:
        tuple: Detrended light curves, fitted trends, and window length grid

    Raises:
        ValueError: If window length inputs are invalid
    """
    # Check that the input window lengths are sensible
    if ((max_win_len * 10) % 1 != 0) or ((min_win_len * 10) % 1 != 0):
        raise ValueError("Invalid window length range. Please input limits that are divisible by 0.1.")

    if max_win_len == 0 or min_win_len == 0:
        raise ValueError("Window length grid cannot include 0.")

    if max_win_len <= min_win_len:
        raise ValueError("Maximum window length must be greater than minimum window length.")

    # Define the window length grid to loop over
    win_len_grid = np.arange(start=min_win_len, stop=max_win_len+0.1, step=0.1)

    # Initialise storage of the detrending output
    biweight_lcs = []
    biweight_trends = []

    # Loop over window lengths
    for wlen in win_len_grid:
        flatten_lc, trend_lc = flatten(time, flux, window_length=wlen, return_trend=True, method='biweight', edge_cutoff=edge_cutoff)
        biweight_lcs.append(flatten_lc)
        biweight_trends.append(trend_lc)

    return biweight_lcs, biweight_trends, win_len_grid


def detrend(
        time,
        flux,
        flux_err,
        method,
        fname,
        cos_win_len_max=12,
        cos_win_len_min=1,
        fap_threshold=1e-2,
        poly_order=2,
        mask=[],
        edge_cutoff=0,
        max_splines=25,
        bi_win_len_max=3,
        bi_win_len_min=1
):
    """Detrend a TESS EB light curve using different methods.

    Method 1: Iterative cosine + multi-biweight detrending
    Method 2: Iterative cosine + penalised spline detrending

    Args:
        time (array_like): Array of time values
        flux (array_like): Array of flux values
        flux_err (array_like): Array of flux error values
        method (str): Detrending method to use ('cp' or 'cb')
        fname (str): Filename to print to user output
        cos_win_len_max (int, optional): Maximum window length for cosine detrending in days. Defaults to 12.
        cos_win_len_min (int, optional): Minimum window length for cosine detrending in days. Defaults to 1.
        fap_threshold (float, optional): False alarm probability threshold. Defaults to 1e-2.
        poly_order (int, optional): Polynomial order for trend fitting. Defaults to 2.
        mask (list, optional): Boolean mask for data points to exclude from the fitting. Defaults to [].
        edge_cutoff (int, optional): Defines the amount of data at the edges to exclude in days. Defaults to 0.
        max_splines (int, optional): Maximum number of splines for penalised spline fitting. Defaults to 25.
        bi_win_len_max (int, optional): Maximum window length for biweight detrending in days. Defaults to 3.
        bi_win_len_min (int, optional): Minimum window length for biweight detrending in days. Defaults to 1.

    Returns:
        tuple: Detrended light curve, fitted trend, biweight window lengths (if applicable), and number of cosine detrending successes
    """

    cos_success = 0

    if len(mask) == 0:
        mask = ~np.isnan(flux * time * flux_err)

    # Iterative cosine detrending
    cos_detrend = cosine_detrend(time[mask], flux[mask],
                                 flux_err[mask], win_len_max=cos_win_len_max,
                                 win_len_min=cos_win_len_min, threshold=fap_threshold,
                                 poly_order=poly_order, mask=None, edge_cutoff=edge_cutoff)

    # Cosine + pspline
    if method == 'cp':
        if isinstance(cos_detrend, tuple):
            cos_flux, cos_trend, cos_wlen = cos_detrend
            flatten_lc, pspline_trend = flatten(time[mask], cos_flux, return_trend=True, method='pspline', max_splines=max_splines)
            trend_lc = cos_trend * pspline_trend
            logger.info(f"Cosine window length: {round(cos_wlen,1)} days")
            cos_success += 1
            logger.debug(f"{fname} cos + pspline")
        else:
            flatten_lc, trend_lc = flatten(time[mask], flux[mask], return_trend=True, method='pspline', max_splines=max_splines)
            logger.debug(f'{fname} pspline only')

        return flatten_lc, trend_lc, [], cos_success

    # Cosine + biweight
    if method == 'cb':
        if isinstance(cos_detrend, tuple):
            cos_flux, cos_trend, cos_wlen = cos_detrend
            flatten_lcs, biweight_trends, biweight_win_lens = run_multi_biweight(time[mask], cos_flux,
                                                                                 max_win_len=bi_win_len_max, min_win_len=bi_win_len_min)
            trend_lcs = cos_trend * biweight_trends
            logger.info(f"{fname} cos + biweight")
            logger.info(f"Cosine window length: {round(cos_wlen,1)} days")
            cos_success += 1
        else:
            cos_flux = cos_detrend
            flatten_lcs, trend_lcs, biweight_win_lens = run_multi_biweight(time[mask], cos_flux,
                                                                           max_win_len=bi_win_len_max, min_win_len=bi_win_len_min)
            logger.debug(f'{fname} biweight only')

        return flatten_lcs, trend_lcs, biweight_win_lens, cos_success
