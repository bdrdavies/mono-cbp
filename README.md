# mono-cbp: Search for Monotransits of Circumbinary Planets

A Python package for detecting circumbinary planets in TESS eclipsing binary light curves through the identification of single transit events ("monotransits").

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

**mono-cbp** is a pipeline designed to systematically search for circumbinary planets by detecting individual transit signatures in TESS eclipsing binary systems. The pipeline automates the complete workflow from masking stellar eclipses, threshold crossing event (TCE) detection, Bayesian vetting, and completeness analysis, making it easy to process large catalogues of eclipsing binaries.

## Key Features

- **Eclipse Masking**: Automatically mask primary and secondary eclipses in eclipsing binary light curves using eclipse positions and widths and binary ephemeris provided by an input catalogue
- **Transit Detection**: Removes unwanted trends from the input light curves and performs single-event detection using the by identifying Threshold Crossing Events (see [Hawthorn et al. 2024](https://academic.oup.com/mnras/article/528/2/1841/7589620?login=false))
- **Bayesian Model Comparison**: Event classification to discern transit-like events and systematics/detrending artefacts
- **Injection-Retrieval Testing**: Completeness analysis through synthetic transit injection and recovery statistics
- **Modular Architecture**: Use individual components independently or run the complete integrated pipeline
- **Configuration-Driven**: Easily customise parameters via Python dictionaries or JSON files without modifying code
- **Command-Line Interface**: Shell scripts and CLI subcommands for batch processing and reproducibility

## Installation

### Requirements

- Python 3.8 or higher (tested most rigorously with Python 3.9)

### From PyPI (Recommended)

The easiest way to install `mono-cbp` is from PyPI:

```bash
pip install mono-cbp
```

It is advisable to install `mono-cbp` into a Python environment using your favourite package manager, e.g. for `conda`:

```bash
conda create --name mono-cbp python=3.9
conda activate mono-cbp
pip install mono-cbp
```

This installs the package and creates the `mono-cbp` command-line tool.

### From Source

For development or to use the latest unreleased features:

```bash
git clone https://github.com/bdrdavies/mono-cbp.git
cd mono-cbp
pip install -e .
```

This installs the package in editable mode.

### Verify Installation

To check that the installation has been successful:

```bash
python -c "import mono_cbp; print(mono_cbp.__version__)"
mono-cbp --help
```

### Dependencies

All dependencies are automatically installed when you install `mono-cbp`.

See [pyproject.toml](pyproject.toml) or [requirements.txt](requirements.txt) for the complete dependency list and version constraints.

### Troubleshooting Installation

If you encounter issues:

- **Python version**: The package has been tested most thoroughly with Python 3.9.
- **Dependency conflicts**: The package pins specific versions of PyMC, ArviZ, and Bokeh for compatibility. If you have conflicts, create a fresh environment
- **Import errors**: If you see errors related to `bokeh` or `arviz`, ensure you have the correct versions installed (see [pyproject.toml](pyproject.toml))

## Examples & Tutorials

There are a series of Jupyter notebooks in the `examples/` directory to demonstrate how to use the package in your own code:

1. **[00_download_light_curves.ipynb](examples/00_download_light_curves.ipynb)** - Download TESS light curves in the `mono-cbp` format using [lightkurve](https://lightkurve.github.io/lightkurve/)
2. **[01_complete_pipeline.ipynb](examples/01_complete_pipeline.ipynb)** - End-to-end workflow on sample data
3. **[02_eclipse_masking.ipynb](examples/02_eclipse_masking.ipynb)** - Eclipse masking demo
4. **[03_transit_finding.ipynb](examples/03_transit_finding.ipynb)** - TCE detection example
5. **[04_model_comparison.ipynb](examples/04_model_comparison.ipynb)** - Bayesian model comparison example
6. **[05_injection_retrieval.ipynb](examples/05_injection_retrieval.ipynb)** - Completeness testing

## Documentation

Documentation is available in the `docs/` directory:

- **[docs/quickstart.md](docs/quickstart.md)** - Quickstart guide
- **[docs/data_formats.md](docs/data_formats.md)** - Input and output data format specifications
- **[docs/configuration.md](docs/configuration.md)** - Configuration system reference
- **[docs/api_reference.md](docs/api_reference.md)** - API documentation

## Support & Contact

For questions, issues, or feature requests:
- **Issues:** Open an issue on [GitHub Issues](https://github.com/bdrdavies/mono-cbp/issues)
- **Documentation:** Review the [full documentation](docs/)
- **Email:** ben.d.r.davies@warwick.ac.uk