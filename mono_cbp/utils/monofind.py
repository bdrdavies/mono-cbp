"""Collection of functions to find threshold crossing events (TCEs).

Contains monofind and its associated helper functions.
"""

from scipy.stats import median_abs_deviation
import numpy as np


def get_gaps_indices(time, break_tolerance):
    """Array indices where 'time' has gaps longer than 'break_tolerance'.

    Args:
        time (array_like): Array of time values
        break_tolerance (float): Threshold for the gap distance in days

    Returns:
        np.ndarray: Indices on the time axis where gaps occur
    """
    gaps = np.diff(time)
    gaps_indices = np.where(gaps > break_tolerance)
    gaps_indices = np.add(gaps_indices, 1)
    gaps_indices = np.concatenate(gaps_indices).ravel()
    gaps_indices = np.append(np.array([0]), gaps_indices)
    gaps_indices = np.append(gaps_indices, np.array([len(time)+1]))
    return gaps_indices


def split_tol(test_list, tol):
    """Take a list and split it into sub-lists with similar values defined by a tolerance.

    Args:
        test_list (list): List to split
        tol (int): Threshold tolerance

    Returns:
        list: List containing split lists
    """
    result = []
    res = []
    start = test_list[0]
    for ele in test_list:
        if ele - start > tol:
            result.append(res)
            res = []
            start = ele
        res.append(ele)
    result.append(res)
    return result


def create_mad_indices(time, npoints):
    """Create a running window to calculate the MAD over a light curve.

    Args:
        time (array): Time values of a light curve
        npoints (int): Number of data points to define your window

    Returns:
        list: The start and stop indices to calculate the MAD over a running window
    """
    output = []

    for i in range(len(time)):
        if i < npoints:
            start = 0
            stop = i + npoints
        elif i >= len(time) - npoints:
            start = i - npoints
            stop = len(time) - 1
        else:
            start = i - npoints
            stop = i + npoints

        output.append([start, stop])

    return output


def get_var_mad(flux, npoints):
    """Calculates the MAD of a light curve over a running window.

    Args:
        flux (array): Flux values of a light curve
        npoints (int): Number of data points in a given window

    Returns:
        np.ndarray: The MAD of the light curve over a running window
    """
    mad_indices = create_mad_indices(flux, npoints)

    var_mad = []
    for i in range(len(mad_indices)):
        val = median_abs_deviation(flux[mad_indices[i][0]:mad_indices[i][1]])
        var_mad.append(val)

    return np.asarray(var_mad)


def consecutive(data, stepsize=1):
    """Split an array into consecutive chunks if gap between values is above a given step size.

    Args:
        data (array): Data to find gaps in
        stepsize (int, optional): The size of the step to cluster consecutive data. Defaults to 1.

    Returns:
        list: A list of arrays of clustered data
    """
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)


def calculate_threshold(flux, mad, var_mad=None):
    """Calculate event detection threshold.

    Args:
        flux (array): Flux values of a light curve
        mad (float): Median absolute deviation multiplier
        var_mad (array, optional): Rolling window variable MAD. Defaults to None.

    Returns:
        tuple: The calculated threshold array and meta threshold value
    """
    if var_mad is not None:
        threshold = mad * var_mad
        meta_threshold = np.ones(len(threshold)) - threshold
    else:
        threshold = np.ones(len(flux)) * mad * median_abs_deviation(flux[~np.isnan(flux)])
        meta_threshold = 1 - threshold[0]

    return threshold, meta_threshold


def detect_events(flux, threshold):
    """Detect events as three consecutive points below threshold.

    Args:
        flux (array): Flux values of a light curve
        threshold (array): Detection threshold

    Returns:
        np.ndarray: Indices of detected events
    """
    mask = (flux < 1 - threshold)
    # Find indices where three consecutive points are below threshold
    events = np.where(mask[:-2] & mask[1:-1] & mask[2:])[0]
    return events


def group_events(events, stepsize=4, min_separation=30):
    """Group events into consecutive chunks and filter close events.

    Args:
        events (array): Indices of detected events
        stepsize (int, optional): Step size for grouping events. Defaults to 4.
        min_separation (int, optional): Minimum cadences between events. Defaults to 30.

    Returns:
        np.ndarray: Indices of grouped events
    """
    if len(events) == 0:
        return np.array([])
    events_split = consecutive(events, stepsize=stepsize)
    mask = np.full(len(events_split), True)
    for i in range(len(events_split)):
        # Fill in missing integers
        events_split[i] = np.arange(events_split[i][0], events_split[i][-1] + 1)
        # Add next two cadences
        events_split[i] = np.append(events_split[i], [events_split[i][-1] + 1, events_split[i][-1] + 2])
        # Find mid-point
        events_split[i] = int(np.median(events_split[i]))
        if i > 0 and events_split[i] - events_split[i-1] < min_separation:
            mask[i-1] = False
    peaks = np.asarray(events_split)[mask]
    return peaks


def extract_event_metadata(time, flux, peaks):
    """Extract metadata for each event: depth, width, start/end times.

    Args:
        time (array): Time values of the light curve
        flux (array): Flux values of the light curve
        peaks (array): Indices of detected events

    Returns:
        tuple: Arrays of depths, widths, start times, and end times
    """
    under_trend = flux < 1
    depths, widths, start_times, end_times = [], [], [], []
    for t_mid_idx in peaks:
        t_max, t_min = [], []
        # Find end time
        for k in range(len(time[t_mid_idx:])):
            if len(t_max) == 1:
                end_idx = t_mid_idx + k - 1
                break
            elif under_trend[t_mid_idx + k]:
                continue
            else:
                t_max.append(time[t_mid_idx + k])
        # Find start time
        for j in range(t_mid_idx)[::-1]:
            if len(t_min) == 1:
                start_idx = j + 1
                break
            elif not under_trend[j]:
                t_min.append(time[j])
        try:
            widths.append(t_max[0] - t_min[0])
            start_times.append(t_min[0])
            end_times.append(t_max[0])
        except IndexError:
            if len(t_max) == 0:
                widths.append(time[-1] - t_min[0])
                start_times.append(t_min[0])
                end_times.append(time[-1])
            if len(t_min) == 0:
                widths.append(t_max[0] - time[0])
                start_times.append(time[0])
                end_times.append(t_max[0])
        try:
            # Check that end_idx has been updated if the event is at the end of the sector/quarter
            if end_idx < start_idx:
                end_idx = len(time) - 1
            depths.append(1 - np.median(flux[start_idx:end_idx]))
        except UnboundLocalError:
            depths.append(1 - flux[t_mid_idx])
    return depths, widths, start_times, end_times


def monofind(time, flux, mad=3., var_mad=None):
    """Find individual TCEs in timeseries data based on the MAD of the light curve.

    Args:
        time (array): Time values of the light curve
        flux (array): Flux values of the light curve
        mad (float, optional): MAD multiplier for event detection. Defaults to 3.
        var_mad (array, optional): MAD values calculated using a rolling window. Defaults to None.

    Returns:
        tuple: Arrays of detected event indices and metadata dictionary
    """
    # Initialise metadata dictionary
    meta = {'threshold': [], 'depths': [], 'widths': [], 'start_times': [], 'end_times': []}
    # Calculate detection threshold
    threshold, meta_threshold = calculate_threshold(flux, mad, var_mad)
    meta['threshold'] = meta_threshold
    # Detect events
    events = detect_events(flux, threshold)
    # Return if no events found
    if len(events) == 0:
        return [], meta
    # Group events
    peaks = group_events(events)
    # Extract event metadata
    depths, widths, start_times, end_times = extract_event_metadata(time, flux, peaks)

    meta['depths'] = depths
    meta['widths'] = widths
    meta['start_times'] = start_times
    meta['end_times'] = end_times

    return peaks, meta
