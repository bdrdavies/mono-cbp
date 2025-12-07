# mono-cbp Documentation

Welcome to the documentation for **mono-cbp**, a pipeline for detecting and vetting circumbinary planets in TESS eclipsing binary light curves through single transit event identification.

## Table of Contents

1. [Installation](installation.md)
2. [Quick Start](quickstart.md)
3. [User Guide](user_guide.md)
4. [API Reference](api_reference.md)
5. [Configuration](configuration.md)

## Overview

The mono-cbp pipeline provides a complete workflow for:

- **Eclipse Masking**: Automated masking of primary and secondary eclipses in eclipsing binary systems using an input catalogue of the binary ephemerides and eclipse parameters (widths/positions)
- **Transit Detection**: Light curve detrending and single transit event detection using the `monofind` algorithm
- **Model Comparison**: Bayesian model comparison to classify detected events
- **Injection-Retrieval**: Completeness testing through synthetic transit injection

## Key Features

- Use individual components independently or as a complete pipeline
- Easily customize behavior through configuration files or dictionaries
- In-built filtering and vetting procedures
- Comprehensive error handling and logging

## Pipeline Workflow

```
Light Curve Data
      ↓
Eclipse Masking
      ↓
Transit Finding
      ↓
Model Comparison
      ↓
Candidate Transit Events
```

## Getting Help

- Review the Jupyter notebooks in examples/ for usage patterns
- Open an issue on [GitHub](https://github.com/bdrdavies/mono-cbp/issues)

## Citation

If you use mono-cbp in your research, please cite:

```bibtex
@software{mono_cbp,
  author = {Davies, Benjamin D R},
  title = {mono-cbp: Circumbinary Planet Detection Pipeline},
  year = {2025},
  url = {https://github.com/bdrdavies/mono-cbp}
}
```