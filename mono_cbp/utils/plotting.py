"""Functions for generating plots."""

import matplotlib.pyplot as plt
import numpy as np
import os
import logging

logger = logging.getLogger('mono_cbp.utils.plotting')


def plot_no_events(
        time,
        flat_flux,
        raw_flux,
        flux_err,
        trend,
        fname,
        mad,
        var_mad,
        ecl_mask=None,
        output_dir=None,
        mask=[],
        figsize=(20/3, 8),
        save=True,
        return_fig=False
):
    """Plot light curves with no events detected.

    Args:
        time (array_like): Time array
        flat_flux (array_like): Flattened flux array
        raw_flux (array_like): Raw flux array
        flux_err (array_like): Flux error array
        trend (array_like): Trend array
        fname (str): Filename
        mad (float): Threshold multiplier of Median Absolute Deviation
        var_mad (array_like): Variable Median Absolute Deviation
        ecl_mask (array_like, optional): Eclipse mask array. Defaults to None.
        output_dir (str, optional): Output directory. Defaults to None.
        mask (list, optional): Mask for the data. Defaults to [].
        figsize (tuple, optional): Figure size. Defaults to (20/3, 8).
        save (bool, optional): Whether to save the plot to disk. Defaults to True.
        return_fig (bool, optional): Whether to return the figure object. Defaults to False.

    Returns:
        matplotlib.figure.Figure or None: Figure object if return_fig is True
    """
    # If no mask is provided, use all non-NaN data points
    if len(mask) == 0:
        mask = ~np.isnan(raw_flux * time * flux_err)

    # Initialise plot
    fig, axs = plt.subplots(2, 1, figsize=figsize)

    # Detrended light curve
    axs[0].errorbar(time[mask], flat_flux, yerr=flux_err[mask], fmt='none', ecolor='black', elinewidth=1, zorder=1)
    axs[0].plot(time[mask], flat_flux, 'o', markerfacecolor='grey', markeredgecolor='black', markersize=5, zorder=2, label='Data')
    if ecl_mask is not None:
        axs[0].fill_between(time, 0, 1, where=ecl_mask, color='grey', alpha=0.7, transform=axs[0].get_xaxis_transform())

    # Light curve with detrending
    axs[1].errorbar(time[mask], raw_flux[mask], yerr=flux_err[mask], fmt='none', ecolor='black', elinewidth=1, zorder=1)
    axs[1].plot(time[mask], raw_flux[mask], 'o', markerfacecolor='grey', markeredgecolor='black', markersize=5, zorder=2)
    if ecl_mask is not None:
        axs[1].fill_between(time, 0, 1, where=ecl_mask, color='grey', alpha=0.7, transform=axs[1].get_xaxis_transform())

    # Plot the detrending without plotting it over any gaps in the data
    threshold = 0.5
    dt = np.diff(time[mask])
    split_indices = np.where(dt > threshold)[0] + 1
    segments = np.split(np.arange(len(time[mask])), split_indices)

    for seg in segments:
        axs[1].plot(time[mask][seg], trend[seg], linewidth=2, c='red')
        axs[0].plot(time[mask][seg], ((1-mad*var_mad)[seg]), 'r--')

    axs[0].set_ylabel('Detrended Flux')
    axs[1].set_ylabel('Normalised Raw Flux')
    axs[1].set_xlabel('Time - 2457000 (BTJD days)')
    plt.tight_layout()

    # Only save if requested
    if save:
        if output_dir:
            save_dir = output_dir
            os.makedirs(save_dir, exist_ok=True)
        else:
            save_dir = os.path.join(os.getcwd(), 'output_plots', 'no_events')
            os.makedirs(save_dir, exist_ok=True)
        plot_path = os.path.join(save_dir, fname[:-3]+'_'+fname[-2:]+'.png')
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved no-events plot to {plot_path}")
        plt.close(fig)
        if return_fig:
            return fig
    else:
        if return_fig:
            return fig
        else:
            plt.close(fig)


def plot_event(
        time,
        event_time,
        flat_flux,
        raw_flux,
        flux_err,
        trend,
        fname,
        mad,
        var_mad,
        depth,
        width,
        phase,
        SNR,
        peaks,
        event_no,
        ecl_mask=None,
        output_dir=None,
        mask=[],
        figsize=(20/3, 8),
        save=True,
        return_fig=False
):
    """Plot the light curve for a detected event.

    Args:
        time (array_like): Array of time values
        event_time (float): Time of the event (in days)
        flat_flux (array_like): Flattened flux array
        raw_flux (array_like): Raw flux array
        flux_err (array_like): Flux error array
        trend (array_like): Trend array
        fname (str): Filename
        mad (float): Threshold multiplier of Median Absolute Deviation
        var_mad (array_like): Variable Median Absolute Deviation
        depth (float): Depth of the event
        width (float): Width of the event
        start_time (float): Start time of the event
        end_time (float): End time of the event
        phase (float): Binary phase of the event
        SNR (float): Signal-to-noise ratio of the event
        peaks (array_like): Array of detected peak indices
        event_no (int): Event number
        ecl_mask (array_like, optional): Eclipse mask array. Defaults to None.
        output_dir (str, optional): Output directory. Defaults to None.
        mask (list, optional): List of boolean masks for the data. Defaults to [].
        figsize (tuple, optional): Figure size. Defaults to (20/3, 8).
        save (bool, optional): Whether to save the plot. Defaults to True.
        return_fig (bool, optional): Whether to return the figure object. Defaults to False.

    Returns:
        matplotlib.figure.Figure or None: Figure object if return_fig is True
    """
    # Initialise the plot
    fig, axs = plt.subplots(3, 1, figsize=figsize)

    # Event picked up in flattened light curve
    axs[0].errorbar(time[mask], flat_flux, yerr=flux_err[mask], fmt='none', ecolor='black', elinewidth=1, zorder=1)
    axs[0].plot(time[mask], flat_flux, 'o', markerfacecolor='grey', markeredgecolor='black', markersize=5, zorder=2, label='Data')
    if ecl_mask is not None:
        axs[0].fill_between(time, 0, 1, where=ecl_mask, color='grey', alpha=0.7, transform=axs[0].get_xaxis_transform())
    axs[0].set_xlim(event_time - 1, event_time + 1)

    # Rescale the y axis
    xbound = axs[0].get_xbound()
    idx = np.where(np.logical_and(time[mask]>=xbound[0], time[mask]<=xbound[1]))
    f_new = flat_flux[idx]
    # Filter out NaN and Inf values
    f_new_valid = f_new[np.isfinite(f_new)]
    if len(f_new_valid) > 0:
        f_max, f_min = np.max(f_new_valid), np.min(f_new_valid)
        axs[0].set_ylim(f_min-0.001, f_max+0.001)
    else:
        logger.warning(f"No valid flux values found in zoom range for event at {event_time}")
    axs[0].scatter(event_time, axs[0].get_ylim()[0]+0.05*(axs[0].get_ylim()[1]-axs[0].get_ylim()[0]), s=50, c='blue', marker=6)
    axs[0].text(
        0.84, 0.07,
        r'$\delta = ' + str(round(depth*100,2)) + r'\%$' + '\n' +
        r'$t_{\mathrm{dur}} = ' + str(round(width,2)) + r'\,\mathrm{d}$' + '\n' +
        r'$\phi_b = ' + str(round(phase,2)) + r'$' + '\n' +
        r'$\mathrm{SNR} = ' + str(round(SNR,1)) + r'$',
        fontsize=10,
        ha='left',
        transform=axs[0].transAxes,
        bbox=dict(facecolor='white', edgecolor='white', pad=5.0, alpha=0.1, zorder=0)
    )

    # Plot the whole detrended light curve
    axs[1].errorbar(time[mask], flat_flux, yerr=flux_err[mask], fmt='none', ecolor='black', elinewidth=1, zorder=1)
    axs[1].plot(time[mask], flat_flux, 'o', markerfacecolor='grey', markeredgecolor='black', markersize=5, zorder=2)
    if ecl_mask is not None:
        axs[1].fill_between(time, 0, 1, where=ecl_mask, color='grey', alpha=0.7, transform=axs[1].get_xaxis_transform())

    # Rescale y axis
    xbound = axs[1].get_xbound()
    idx = np.where(np.logical_and(time[mask]>=xbound[0], time[mask]<=xbound[1]))
    f_new = flat_flux[idx]
    # Filter out NaN and Inf values
    f_new_valid = f_new[np.isfinite(f_new)]
    if len(f_new_valid) > 0:
        f_max, f_min = np.max(f_new_valid), np.min(f_new_valid)
        axs[1].set_ylim(f_min-0.001, f_max+0.001)
    else:
        logger.warning(f"No valid flux values found in full detrended range for event at {event_time}")
    axs[1].scatter(time[mask][peaks], (axs[1].get_ylim()[0]+0.05*(axs[1].get_ylim()[1]-axs[1].get_ylim()[0]))*np.ones(len(time[mask][peaks])), s=50, c='gray', marker=6)
    axs[1].scatter(event_time, axs[1].get_ylim()[0]+0.05*(axs[1].get_ylim()[1]-axs[1].get_ylim()[0]), s=50, c='blue', marker=6)

    # Plot the whole raw light curve with removed trend
    axs[2].errorbar(time[mask], raw_flux[mask], yerr=flux_err[mask], fmt='none', ecolor='black', elinewidth=1, zorder=1)
    axs[2].plot(time[mask], raw_flux[mask], 'o', markerfacecolor='grey', markeredgecolor='black', markersize=5, zorder=2)
    if ecl_mask is not None:
        axs[2].fill_between(time, 0, 1, where=ecl_mask, color='grey', alpha=0.7, transform=axs[2].get_xaxis_transform())

    # Plot the detrending without plotting it over any gaps in the data
    threshold = 0.5
    dt = np.diff(time[mask])
    split_indices = np.where(dt > threshold)[0] + 1
    segments = np.split(np.arange(len(time[mask])), split_indices)

    for seg in segments:
        axs[2].plot(time[mask][seg], trend[seg], linewidth=2, c='red')
        axs[1].plot(time[mask][seg], ((1-mad*var_mad)[seg]), 'r--')

    axs[0].set_ylabel('Detrended Flux')
    axs[1].set_ylabel('Detrended Flux')
    axs[2].set_ylabel('Normalised Raw Flux')
    axs[2].set_xlabel('Time - 2457000 (BTJD days)')
    plt.tight_layout()
    axs[2].scatter(time[mask][peaks], (axs[2].get_ylim()[0]+0.05*(axs[2].get_ylim()[1]-axs[2].get_ylim()[0]))*np.ones(len(time[mask][peaks])), s=50, c='gray', marker=6)
    axs[2].scatter(event_time, axs[2].get_ylim()[0]+0.05*(axs[2].get_ylim()[1]-axs[2].get_ylim()[0]), s=50, c='blue', marker=6)

    # Only save if requested
    if save:
        if output_dir:
            save_dir = output_dir
            os.makedirs(save_dir, exist_ok=True)
        else:
            save_dir = os.path.join(os.getcwd(), 'output_plots', 'events')
            os.makedirs(save_dir, exist_ok=True)
        plot_path = os.path.join(save_dir, fname[:-3]+'_'+fname[-2:]+'_'+str(event_no)+'.png')
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved event plot to {plot_path}")
        plt.close(fig)
        if return_fig:
            return fig
    else:
        if return_fig:
            return fig
        else:
            plt.close(fig)