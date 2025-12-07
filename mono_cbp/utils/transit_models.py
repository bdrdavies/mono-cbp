"""Transit model generation utilities for injection-retrieval testing."""

import numpy as np
import pymc as pm
import exoplanet as xo


def create_transit_models(
    depth_range=(1e-3, 1e-2),
    duration_range=(0.1, 1.0),
    num_depths=7,
    num_durations=7,
    time_range=(-1, 1),
    cadence_minutes=30,
    impact_parameter=0.0,
    period=10.0,
    limb_dark_coeffs=(0.3, 0.2)
):
    """Create a grid of synthetic transit light curves for injection-retrieval testing.

    This function generates transit models across a grid of transit depths and durations,
    which can be used to characterize pipeline completeness and detection efficiency.

    Args:
        depth_range (tuple, optional): (min, max) transit depth in fractional flux units.
                                       Defaults to (1e-3, 1e-2) = 0.1% to 1%.
        duration_range (tuple, optional): (min, max) transit duration in days.
                                         Defaults to (0.1, 1.0).
        num_depths (int, optional): Number of depth values to sample. Defaults to 7.
        num_durations (int, optional): Number of duration values to sample. Defaults to 7.
        time_range (tuple, optional): (start, end) time range in days centered on transit.
                                     Defaults to (-1, 1).
        cadence_minutes (float, optional): Observation cadence in minutes. Defaults to 30.
        impact_parameter (float, optional): Impact parameter (0 = center of limb). Defaults to 0.0.
        period (float, optional): Orbital period in days (arbitrary, >2×duration). Defaults to 10.0.
        limb_dark_coeffs (tuple, optional): Quadratic limb darkening coefficients (u1, u2).
                                           Defaults to (0.3, 0.2).

    Returns:
        dict: Dictionary containing:
            - 'time': Time array (same for all models)
            - 'models': List of transit model dictionaries, each containing:
                - 'flux': Normalized flux array
                - 'depth': Transit depth
                - 'duration': Transit duration in days
                - 'impact_parameter': Impact parameter
                - 'ror': Radius ratio (planet radius / star radius)
            - 'num_depths': Number of depth values
            - 'num_durations': Number of duration values
            - 'depth_range': Depth range tuple
            - 'duration_range': Duration range tuple
            - 'cadence_minutes': Cadence in minutes
    """
    # Define depth grid (log scale)
    depth_min, depth_max = depth_range
    depths = np.logspace(np.log10(depth_min), np.log10(depth_max), num=num_depths)

    # Define duration grid (linear scale)
    duration_min, duration_max = duration_range
    durations = np.linspace(duration_min, duration_max, num=num_durations)

    # Create time array
    time_start, time_end = time_range
    cadence = cadence_minutes / (60 * 24)  # Convert minutes to days
    time = np.arange(time_start, time_end + cadence, cadence)

    # Store transit models
    models = []

    # Loop over grid
    for depth in depths:
        for duration in durations:
            with pm.Model():
                # Approximate the radius ratio based on the transit depth
                ror = xo.LimbDarkLightCurve(limb_dark_coeffs).get_ror_from_approx_transit_depth(
                    depth, impact_parameter
                )

                # Obtain simple transit orbit
                orbit = xo.orbits.SimpleTransitOrbit(
                    period=period,
                    duration=duration,
                    t0=0,
                    ror=ror,
                    b=impact_parameter
                )

                # Generate the transit light curve
                light_curve = xo.LimbDarkLightCurve(limb_dark_coeffs).get_light_curve(
                    orbit=orbit, r=ror, t=time
                )

                # Compute the model flux
                flux = pm.math.sum(light_curve, axis=-1).eval()

            # Store model
            # Note: ror is a TensorVariable, so we need to evaluate it
            if hasattr(ror, 'eval'):
                ror_value = float(ror.eval())
            else:
                ror_value = float(ror)
            model_dict = {
                'flux': flux,
                'depth': depth,
                'duration': duration,
                'impact_parameter': impact_parameter,
                'ror': ror_value
            }
            models.append(model_dict)

    # Return dictionary with time and all models
    result = {
        'time': time,
        'models': models,
        'num_depths': num_depths,
        'num_durations': num_durations,
        'depth_range': depth_range,
        'duration_range': duration_range,
        'cadence_minutes': cadence_minutes
    }

    return result


def save_transit_models(models_dict, filepath):
    """Save transit models to an .npz file.

    Args:
        models_dict (dict): Dictionary from create_transit_models()
        filepath (str): Path to save .npz file
    """
    # Prepare data for saving
    save_dict = {
        'time': models_dict['time'],
        'num_depths': models_dict['num_depths'],
        'num_durations': models_dict['num_durations'],
        'depth_range': models_dict['depth_range'],
        'duration_range': models_dict['duration_range'],
        'cadence_minutes': models_dict['cadence_minutes']
    }

    # Add each model as a separate array
    for idx, model in enumerate(models_dict['models']):
        save_dict[f'model_{idx}_flux'] = model['flux']
        save_dict[f'model_{idx}_depth'] = model['depth']
        save_dict[f'model_{idx}_duration'] = model['duration']
        save_dict[f'model_{idx}_impact_parameter'] = model['impact_parameter']
        save_dict[f'model_{idx}_ror'] = model['ror']

    np.savez(filepath, **save_dict)


def load_transit_models(filepath):
    """Load transit models from an .npz file.

    Args:
        filepath (str): Path to .npz file created by save_transit_models()

    Returns:
        dict: Dictionary with same structure as create_transit_models() output
    """
    data = np.load(filepath, allow_pickle=True)

    # Extract metadata
    result = {
        'time': data['time'],
        'num_depths': int(data['num_depths']),
        'num_durations': int(data['num_durations']),
        'depth_range': tuple(data['depth_range']),
        'duration_range': tuple(data['duration_range']),
        'cadence_minutes': float(data['cadence_minutes']),
        'models': []
    }

    # Count models
    num_models = result['num_depths'] * result['num_durations']

    # Load each model
    for idx in range(num_models):
        model_dict = {
            'flux': data[f'model_{idx}_flux'],
            'depth': float(data[f'model_{idx}_depth']),
            'duration': float(data[f'model_{idx}_duration']),
            'impact_parameter': float(data[f'model_{idx}_impact_parameter']),
            'ror': float(data[f'model_{idx}_ror'])
        }
        result['models'].append(model_dict)

    return result
