"""Collection of utility functions for converting time to phase and obtaining eclipse masks."""

import numpy as np


def time_to_phase(times, period, t0, centre=0.5):
    """Converts an array of time values to phase.

    Args:
        times (array_like): Array of time values
        period (float): Period of the binary
        t0 (float): Time of mid-eclipse
        centre (float, optional): Centre of the phase fold. Defaults to 0.5.

    Returns:
        np.ndarray: Array of phase values
    """
    start_phase = centre - 0.5
    t0 += start_phase * period
    phases = ((times - t0) / period) % 1.0 + start_phase

    return phases


def get_eclipse_indices(phases, pos, width):
    """Obtain the indices of the positions in phase space where an eclipse is.

    Values that are in-eclipse are labelled as True.

    Args:
        phases (array_like): Phase values to calculate indices
        pos (float): The position in phase space of mid-eclipse
        width (float): The width of the eclipse in phase space

    Returns:
        np.ndarray: The indices of in-eclipse data
    """
    # Check that the position and width values are not NaN, and width is not zero
    # Note: pos==0 is valid (eclipse at phase 0)
    if np.isnan(pos) or np.isnan(width) or width == 0:
        return np.array([], dtype=int)

    if pos > 0.95:
        # For eclipses near phase 1, wrap around by subtracting 1 from the upper bound
        idx = np.where(np.logical_or(phases >= pos - width/2, phases <= pos + width/2 - 1))
    elif pos < 0.05:
        # For eclipses near phase 0, wrap around by adding 1 to the lower bound
        idx = np.where(np.logical_or(phases >= pos - width/2 + 1, phases <= pos + width/2))
    else:
        idx = np.where(np.logical_and(phases >= pos - width/2, phases <= pos + width/2))

    return idx[0]


def get_eclipse_mask(phases, pos, width):
    """Obtain a mask of the phases where an eclipse occurs.

    Args:
        phases (array_like): Phase values to calculate mask
        pos (float): The position in phase space of mid-eclipse
        width (float): The width of the eclipse in phase space

    Returns:
        np.ndarray: A boolean mask where True indicates an eclipse
    """
    idx = get_eclipse_indices(phases, pos, width)
    mask = np.zeros_like(phases, dtype=bool)
    mask[idx] = True

    return mask
