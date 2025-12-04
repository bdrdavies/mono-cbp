"""Default configuration settings for mono-cbp pipeline."""

import copy

DEFAULT_CONFIG = {
    'transit_finding': {
        'edge_cutoff': 0.0,
        'mad_threshold': 3.0,
        'detrending_method': 'cb',
        'generate_vetting_plots': False,
        'generate_skye_plots': False,
        'generate_event_snippets': True,
        'save_event_snippets': True,
        'cadence_minutes': 30,
        'cosine': {
            'win_len_max': 12,
            'win_len_min': 1,
            'fap_threshold': 1e-2,
            'poly_order': 2,
        },
        'biweight': {
            'win_len_max': 3,
            'win_len_min': 1,
        },
        'pspline': {
            'max_splines': 25,
        },
        'filters': {
            'min_snr': 5,
            'max_duration_days': 1,
            'det_dependence_threshold': 18,
        },
    },
    'model_comparison': {
        'tune': 1000,
        'draws': 1000,
        'chains': 4,
        'cores': 4,
        'target_accept': 0.99,
        'sigma_threshold': 3,
        'aic_threshold': 2,
        'rmse_threshold': 1.2,
        'save_plots': False,
        'plot_dir': None,
    },
    'injection_retrieval': {
        'n_injections': 1000,
    },
}


def get_default_config():
    """Get a copy of the default configuration.

    Returns:
        dict: A deep copy of the default configuration
    """
    return copy.deepcopy(DEFAULT_CONFIG)


def merge_config(user_config, default_config=None):
    """Merge user configuration with default configuration.

    Args:
        user_config (dict): User-provided configuration
        default_config (dict, optional): Base configuration. Uses DEFAULT_CONFIG if None.

    Returns:
        dict: Merged configuration
    """
    if default_config is None:
        default_config = get_default_config()

    merged = copy.deepcopy(default_config)

    for key, value in user_config.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = merge_config(value, merged[key])
        else:
            merged[key] = value

    return merged
