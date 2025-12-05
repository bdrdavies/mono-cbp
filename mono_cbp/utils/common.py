"""Common utility functions used across the mono-cbp package."""

import logging
import numpy as np


def setup_logging(level=logging.INFO, log_file=None):
    """Set up logging for mono-cbp.

    Args:
        level (int, optional): Logging level. Defaults to logging.INFO (20).
        log_file (str, optional): Path to log file. If None, logs to console only.

    Returns:
        logging.Logger: Configured logger instance
    """
    logger = logging.getLogger('mono_cbp')
    logger.setLevel(level)

    # Clear existing handlers
    logger.handlers = []

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def get_snr(depth, error, duration, cadence=30):
    """Calculate signal-to-noise ratio for a threshold crossing event (TCE).

    Args:
        depth (float): TCE depth
        error (float): Combined measurement uncertainty (see documentation for details)
        duration (float): TCE duration in days
        cadence (int, optional): Cadence in minutes. Defaults to 30.

    Returns:
        float: Signal-to-noise ratio
    """
    cadence_days = cadence / 1440  # Convert minutes to days
    return (depth / error) * (np.sqrt(duration) / np.sqrt(cadence_days))
